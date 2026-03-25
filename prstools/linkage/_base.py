# __all__ = ['LinkageData'] # ok this does not seem to lead to the desired stuff
#| export
# %%writefile ../prstools/loaders.py
#!/usr/bin/env python  

"""
LinkageData
durr tst
"""

import os, re, sys
import scipy as sp
import numpy as np
import pandas as pd
from scipy import linalg
from sys import getsizeof
import warnings, importlib, json, os, glob, copy, uuid
from collections import OrderedDict, deque, defaultdict
try:
    from pysnptools.standardizer import Unit, UnitTrained
    import pysnptools as pst
except:
    Unit = defaultdict
    UnitTrained = defaultdict
try:
    from tqdm.cli import tqdm
except:
    tqdm = lambda x:x
# import prstools # this is needed below in regdef
import prstools as prst
from prstools.utils import suppress_warnings
try:
    import h5py
except:
    True
try:
    class SqrtNinv(Unit):
        def __init__(self):
            super(SqrtNinv, self).__init__()     
except:
    True 
try: import IPython as ip
except: pass

class _RemovalStagingLinkageData():
    


    @classmethod
    def from_prscsvars(cls, prscs_dt):
        _prscs2prstvarname = {'CHR':'chrom','SNP':'snp',None:'cm','BP':'pos','A1':'A1','A2':'A2','BETA':'beta_mrg','MAF':'maf','FLP':'flp'}
        # self._prscs2prstvarname = {'CHR':'chrom', 'SNP':'snp','BP':'pos','A1':'A1','A2':'A2','MAF':'maf','BETA':'beta_mrg'}
        def proc_sst(slice_df, *, i, **kwargs):
            map_dt = _prscs2prstvarname
            sst_df = slice_df.rename(columns=map_dt)
            if 'cm' not in sst_df.columns: sst_df['cm'] = 0.
            for key,item in kwargs.items():
                sst_df[key] = item
            maf = sst_df['maf']
            sst_df['std_ref'] = np.sqrt(2.0*maf*(1.0-maf))
            sst_df['i']=i
            return sst_df[list(map_dt.values())+[col for col in sst_df.columns if not col in map_dt.values()]]
        cnt = 0; i=0
        reg_dt = dict()
        for chrom, chr_dt in prscs_dt.items():
            ld_blk = chr_dt['ld_blk']
            sst_df = pd.DataFrame.from_dict(chr_dt['sst_dict']) 
            n_gwas = chr_dt['n_gwas']
            for D in ld_blk:
                if len(D) == 0:
                    continue
                start_j=cnt; stop_j=cnt+len(D)
                slice_df = sst_df.iloc[start_j:stop_j]
                slice_df['CHR'].unique()

                geno_dt = dict(regid=i,
                           chrom  =slice_df['CHR'].iloc[0],
                           start  =slice_df['BP'].iloc[0],
                           stop   =slice_df['BP'].iloc[-1]+1,
                           start_j=start_j,
                           stop_j =stop_j,
                           sst_df =proc_sst(slice_df, i=i, n_eff=n_gwas),
                           D      =D,
                              )
                reg_dt[i] = geno_dt
                i+=1; cnt+=len(D)
        clsattr_dt = dict(reg_dt=reg_dt, _n_snps_total=cnt)
        linkdata=cls(clsattr_dt=clsattr_dt)
        return linkdata

    @classmethod
    def from_cli_params(cls, *, ref, target, sst, n_gwas, chrom=range(1,23), verbose=True, **kwg):
        if type(chrom) is int: chrom = [chrom]
        if type(chrom) is str: chrom = chrom.split(',')
        from prstools import parse_genet
        default_param_dict = {'ref_dir': None, 'bim_prefix': None, 'sst_file': None, 'a': 1, 'b': 0.5, 'phi': None, 'n_gwas': None,
          'n_iter': 1000, 'n_burnin': 500, 'thin': 5, 'out_dir': None, 'chrom': range(1,23), 'beta_std': 'False', 'seed': None}
        def prscs_setup(ref_dir    = '../lnk/prscs/ldblk_1kg_eur__', # Adding the / here causes issues..
            sst_file   ='../prstools/PRScs/test_data/sumstats.txt__',
            bim_prefix ='../prstools/PRScs/test_data/test__',
            out_dir    ='./tst_out/weights_NOTREALLYNEEDDHERE', # This needs the exit slash "weights_" is a file prefix
            chrom      = 22,
            n_gwas     = 5000,
            n_iter     = 100,
            n_burnin   = 50,
            create_dfs = False,
            return_ld  = True):

            _excl_lst = ['self', 'param_dict', '_excl_lst','return_ld']
            param_dict = {key: item for key, item in locals().items() if not (key in _excl_lst)}
            param_dict['ref_dir'] = os.path.normpath(param_dict['ref_dir'])
            basename = os.path.basename(param_dict['ref_dir'])
            info = 'ukbb' if 'ukbb' in basename else ('1kg' if '1kg' in basename else None)
            if not info: raise Exception('No \'1kg\' or \'ukbb\' detected in dirname, which is required.'
                                         ' If it is present. Perhaps remove trailing /')
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + f'/snpinfo_{info}_hm3', int(chrom))
            vld_dict = parse_genet.parse_bim(param_dict['bim_prefix'], int(chrom))
            sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, param_dict['sst_file'], param_dict['n_gwas'])
            # import IPython as ip;print('woot'); ip.embed()

            if len(sst_dict['CHR']) == 0:
                print('Continuing to next chromosome, because there were zeros matching SNPs.')
                return False

            if create_dfs:
                ref_df = pd.DataFrame.from_dict(ref_dict)
                vld_df = pd.DataFrame.from_dict(vld_dict)
                sst_df = pd.DataFrame.from_dict(sst_dict)

            if return_ld:
                ld_blk, blk_size = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict, int(chrom))

            xtra_dt = {key:item for key,item in default_param_dict.items() if key in set(default_param_dict.keys()) - set(param_dict.keys())}
            param_dict.update(xtra_dt)

            return locals()


        # This takes a bit of time:
        prscs_dt = dict()
        for c in chrom: # the 21 here should be 1, its a POC.
            if verbose: print('##### process chromosome %d #####' % int(c))
            result = prscs_setup(chrom=c,ref_dir=ref, bim_prefix=target, sst_file=sst, n_gwas=n_gwas)
            if result: prscs_dt[c] = result
        linkdata = cls.from_prscsvars(prscs_dt)
        # import IPython as ip; print('here'); ip.embed()

        return linkdata
    
class GenotypeLinkageData():
    
    ###########################
    # Checking, Validation & Assertion
    
    def _check_xrd(self):

        if self.srd is not None:
            assert pst.snpreader.SnpReader in self.srd.__class__.__mro__

        if self.prd is not None:
            n_start = len(self.prd.iid)
            n_pheno = self.prd.shape[1]
            if n_pheno > 1: raise NotImplementedError(f'only one pheno in prd allowed {n_pheno} detected ({self.prd.col}), for now. remove other phenos')
            self.srd, self.prd = pst.util.intersect_apply([self.srd, self.prd])
            if len(self.prd.iid) != n_start:
                warnings.warn('Number of samples do not match up after internal intersection, samples were lost:' 
                              f'{n_start - len(self.prd.iid)}, start = {n_start}, after_intersection = {len(self.prd.iid)}')
                if not self.intersect_apply: raise Exception('Intersection was required, but may not performed. Hence raising this error.')

        if self.grd is not None:
            # Check alignment for now, auto alignment needs work cause iid stuffs:
            if self.srd is not None:
                if not np.all(self.grd.sid == self.srd.sid):
                    raise Exception('snps of grd and srd not matching up, align first,'
                                    ' auto align will be implemented later')
            else:
                raise NotImplementedError('Not sure what to do with grd if no srd is present. not implemented.')
                
    ###########################
    # Regions Administration:

    def init_regions(self):
        if not ('beta_mrg' in self.sst_df.columns):
            warnings.warn('No \'beta_mrg\' column detected in sst_df! This means that no summary stats were detected.')
        else:
            assert 'n_eff' in self.sst_df.columns
            if not ('std_sst' in sst_df.columns):
                warnings.warn('No \'std_sst\' column detected in sst_df! This will make it impossible to scale the computed PRS weights back to the original allele scale')

        cur_chrom = None
        i = 0; n_snps_cumsum = 0
        sst_df_lst = []
        for reg_cnt, (_, row) in enumerate(self.regdef_df.iterrows()):
            # Move region into specialized dictionary
            regid = row['regid']; chrom = row['chrom']
            start = row['start']; stop  = row['stop']

            # Map Variants to region
            ind = self.sst_df.chrom == chrom
            ind = (self.sst_df['pos'] >= start) & ind
            ind = (self.sst_df['pos'] < stop) & ind
            sid = self.sst_df['snp'][ind].values
            indices = self.srd.sid_to_index(sid)  # if sid not strickly present this will give an error!
            n_snps_reg = len(indices)
            if n_snps_reg == 0:
                continue
            else:
                geno_dt = dict(regid=regid,
                               chrom=chrom,
                               start=start,
                               stop=stop,
                               start_j=n_snps_cumsum)
                n_snps_cumsum += n_snps_reg
                geno_dt['stop_j'] = n_snps_cumsum
                sst_df = self.sst_df[ind].copy(); sst_df['i'] = i
                geno_dt['sst_df'] = sst_df
                assert geno_dt['start_j'] == sst_df.index[0]; sst_df_lst.append(sst_df)
                assert geno_dt['stop_j']  == sst_df.index[-1] + 1
                if self.srd is not None:
                    geno_dt['srd'] = self.srd[:, indices]
                    geno_dt['stansda'] = self.sda_standardizer() if self.sda_standardizer is not None else None
                else:
                    raise NotImplementedError()
                if self.grd is not None:
                    geno_dt['grd'] = self.grd[:, indices]
                    geno_dt['stangda'] = self.gda_standardizer() if self.gda_standardizer is not None else None
                # Count up if things are actually stored in reg_dt
                self.reg_dt[i] = geno_dt
                i += 1
        self._n_snps_total = n_snps_cumsum
        sst_df = pd.concat(sst_df_lst, axis=0)
        self.sst_df = sst_df
        
    ############################
    ## Compute : ###############

    # Local Linkage Stuff: ####
    if True:
    
        def compute_linkage_sameregion(self, *, i):
            return self.compute_linkage_shiftregion(i=i, shift=0)

        def regions_compatible(self, *, i, j):
            try:
                if self.reg_dt[i]['chrom'] == self.reg_dt[j]['chrom']:
                    res = True
                elif self._cross_chrom_ld:
                    res = True
                else:
                    res = False
            except Exception as e:
                if (not (i in self.reg_dt.keys())) or (not (j in self.reg_dt.keys())):
                    res = False
                else:
                    raise e
            return res

        def compute_linkage_shiftregion(self, *, i, shift):
            j = i + shift
            if self.regions_compatible(i=i, j=j):
                self_sda = self.get_sda(i=i)
                dist_sda = self.get_sda(i=j)
                n = len(self_sda.iid)
                S_shift = self_sda.val.T.dot(dist_sda.val) / (n) # old: - self.ddof)
                return S_shift
            else:
                self_sda = self.get_sda(i=i)
                return np.zeros((self_sda.val.shape[1], 0))

        def compute_linkage_cmfromregion(self, *, i, cm):
            geno_dt = self.reg_dt[i]; lst = []
            if cm < 0: # Doing left:
                stop_j   = geno_dt['start_j']
                cm_left  = geno_dt['sst_df'].loc[stop_j]['cm'] 
                slc_df = self.sst_df.loc[:stop_j-1]
                slc_df = slc_df[slc_df.chrom==geno_dt['chrom']]
                slc_df = slc_df[slc_df.cm > (cm_left + cm)]
                start_i = slc_df['i'].min()
                start_i = -7 if np.isnan(start_i) else start_i
                for cur_i in range(start_i, i):
                    lst.append(self.compute_linkage_shiftregion(i=i, shift=cur_i-i))
                    if start_i == -7: break
                L = np.concatenate(lst, axis=1)[:,-slc_df.shape[0]:] # concat & clip
                if self._setzero:
                    cms_reg    = geno_dt['sst_df']['cm'].values
                    cms_distal = slc_df['cm'].values
                    cms_L      =  cms_distal[np.newaxis,:] - cms_reg[:,np.newaxis]
                    setzero_L  = cms_L < cm
                    L[setzero_L] = 0
                    assert L.shape == setzero_L.shape
                return L
            else:
                start_j   = geno_dt['stop_j']
                cm_right  = geno_dt['sst_df'].loc[start_j-1]['cm']
                slc_df = self.sst_df.loc[start_j:]
                slc_df = slc_df[slc_df.chrom==geno_dt['chrom']]
                slc_df = slc_df[slc_df.cm < (cm_right + cm)]
                stop_i = slc_df['i'].max()
                stop_i = i+2 if np.isnan(stop_i) else stop_i + 1
                for cur_i in range(i+1, stop_i):
                    lst.append(self.compute_linkage_shiftregion(i=i, shift=cur_i-i))
                R = np.concatenate(lst, axis=1)[:,:slc_df.shape[0]] # concat & clip
                if self._setzero:
                    cms_reg    = geno_dt['sst_df']['cm'].values
                    cms_distal = slc_df['cm'].values
                    cms_R     =  cms_distal[np.newaxis,:] - cms_reg[:,np.newaxis]
                    setzero_R = cms_R > cm
                    R[setzero_R] = 0
                    assert R.shape == setzero_R.shape
                return R
        
        def compute_maf_region(self, *, i, store=True):
            raise Exception('This method needs to be looked at.')
            sst_df = self.reg_dt[i]['sst_df']
            def quadratic_formula(a, b, c):
                discriminant = b**2 - 4*a*c
                sqrt_discriminant = np.sqrt(np.maximum(0, discriminant))  # Taking the square root, ensuring non-negative values
                x1 = (-b + sqrt_discriminant) / (2*a)
                x2 = (-b - sqrt_discriminant) / (2*a)
                return x1, x2
            a=1;b=-1
            c = (self.get_allele_standev()**2)/2
            _, maf = quadratic_formula(a,b,c)
            sst_df['maf']=maf.flatten()
            return maf
    
    # Marginal Stuff: #############
    if True: # sumstats, allelefreq, ldscores
    
        def compute_sumstats_region(self, *, i, return_n_eff=False, return_allele_std=False):
            geno_dt = self.reg_dt[i]
            sda = self.get_sda(i=i)
            X = sda.val
            y = self.get_pda().val
            n_eff = len(y)
            beta_mrg = X.T.dot(y) / (n_eff) #old - self.ddof)
            output = [beta_mrg, n_eff] if return_n_eff else beta_mrg
            if return_allele_std: 
                if not self._save_s2sst: raise Exception('_save_s2sst has to be set to true for current setup to work')
                output.append(geno_dt['sst_df']['s'])
            return tuple(output)

        def compute_allelefreq_region(self, *, i):
            # Speed might be improved by using dot prod here, instead of sums
            # np.unique was way slower (5x)
            geno_dt = self.reg_dt[i]
            n, p_blk = sda.val.shape
            sst_df = geno_dt['sst_df'].copy()
            cnt0   = np.sum(sda.val==0, axis=0)
            cnt1   = np.sum(sda.val==1, axis=0)
            cnt2   = np.sum(sda.val==2, axis=0)
            cntnan = np.sum(np.isnan(sda.val), axis=0)
            assert np.allclose(cnt0 + cnt1 + cnt2 + cntnan, n)
            sst_df['altcnt=0']   = cnt0
            sst_df['altcnt=1']   = cnt1
            sst_df['altcnt=2']   = cnt2
            sst_df['altcnt=nan'] = cntnan
            sst_df['altfreq']    = (cnt1 + cnt2)/(n - cntnan)
            sst_df['missfreq']   = 1 - cntnan/n
            return sst_df

        def compute_ldscores_region(self, *, i):
            sst_df = self.reg_dt[i]['sst_df'].copy()
            L = self.get_left_linkage_region(i=i)
            D = self.get_auto_linkage_region(i=i)
            R = self.get_right_linkage_region(i=i)
            for k, j in enumerate(sst_df.index):
                slds = np.sum(L[k]**2) + np.sum(D[k]**2) + np.sum(R[k]**2)
                sst_df.loc[j, 'lds'] = np.sqrt(slds)
            return sst_df
        
    ############################
    ## init save load: ############### 
    
    def _load_all_snpdata(self):
        # load all regions
        for i, geno_dt in self.reg_dt.items():
            sda = geno_dt['srd'].read(dtype=self.dtype)
            stansda = sda.train_standardizer(apply_in_place=True,
                                             standardizer=geno_dt['stansda'])
            geno_dt['sda'] = sda
            geno_dt['stansda'] = stansda
        
        
    ############################
    ## Retrieve: ###############
    
    def retrieve_linkage_region(self, *, i):
        Yoootestreadcomment() # wait do we need to validate covariances here?
        
        shift = self.shift; cm = self.cm
        compute_sumstats = self.compute_sumstats

        if 'L' in geno_dt.keys():
            if 'D' in geno_dt.keys():
                if 'R' in geno_dt.keys():
                    return None  # everything is done now.
        
        
        # def something somsthing compute
        # This bit needs a bit of work to work well withitn the new class structure.
        if self.verbose: print(f'Computing LD for region #{i} on chr{geno_dt["chrom"]}', end='\r')
        # Refactor: if linkage is only in blocks this code will lead to recomputation...
        if (shift > 0):
            L_lst = []
            R_lst = []
            for cur_shift in range(1, shift + 1):
                L_lst.append(self.compute_linkage_shiftregion(i=i, shift=-cur_shift))
                R_lst.append(self.compute_linkage_shiftregion(i=i, shift=cur_shift))

            # Store Linkage in geno_dt
            geno_dt['L'] = np.concatenate(L_lst[::-1], axis=1)  # L stands for left
            geno_dt['D'] = self.compute_linkage_sameregion(i=i)  # Linkage within region, D is convention from LDpred 1
            geno_dt['R'] = np.concatenate(R_lst, axis=1)  # R stands for right

            # Indices needed for slicing and dicing matched variables (e.g. beta weights):
            geno_dt['start_j_L'] = geno_dt['start_j'] - geno_dt['L'].shape[1]
            geno_dt['stop_j_L'] = geno_dt['start_j']
            geno_dt['start_j_R'] = geno_dt['stop_j']
            geno_dt['stop_j_R'] = geno_dt['stop_j'] + geno_dt['R'].shape[1]

        elif (shift==0) and (cm is None):  # Only same region has to be done.
            geno_dt['D'] = self.compute_linkage_sameregion(i=i)

        elif (shift==0) and cm > 0:
            geno_dt['L'] = self.compute_linkage_cmfromregion(i=i, cm=-cm)
            geno_dt['D'] = self.compute_linkage_sameregion(i=i)
            geno_dt['R'] = self.compute_linkage_cmfromregion(i=i, cm=cm)

            # Indices needed for slicing and dicing matched variables (e.g. beta weights):
            geno_dt['start_j_L'] = geno_dt['start_j'] - geno_dt['L'].shape[1]
            geno_dt['stop_j_L'] = geno_dt['start_j']
            geno_dt['start_j_R'] = geno_dt['stop_j']
            geno_dt['stop_j_R'] = geno_dt['stop_j'] + geno_dt['R'].shape[1]

        if compute_sumstats:
            self.retrieve_sumstats_region(i=i)
    
    
    ############################
    ## Get:      ###############
    
    
    ###########################
    # Get methods genotype
        
    def get_sda(self, *, i):
        geno_dt = self.reg_dt[i]
        if 'sda' in geno_dt.keys():
            return geno_dt['sda']
        else:
            if 'srd' in geno_dt.keys():
                sda = geno_dt['srd'].read(dtype=self.dtype)
                sda, stansda = sda.standardize(standardizer=geno_dt['stansda'], return_trained=True)
                geno_dt['sda'] = sda
                geno_dt['stansda'] = stansda
                if self._save_s2sst:
                    geno_dt['sst_df']['s'] = stansda.stats[:,[1]]

                if 'loaded_sda' in geno_dt.keys():
                    self.reloaded_xda_cnt += 1
                    if self.reloaded_xda_cnt in [5, 20, 100, 400]:
                        warnings.warn(
                            f'Reloaded sda for the {self.reloaded_xda_cnt}\'th time. This causes memory swapping,'
                            ' that might make the computation of linkage quite slow.'
                            'Probably because memory limits and/or linkage size.')
                # Size determination and accounting:
                geno_dt['loaded_sda']=True
                self.cur_total_size_in_gb += getsizeof(sda.val) / 1024 ** 3
                self.xda_q.append((i,'sda'))  # put respective i in queue.
                while self.cur_total_size_in_gb > self.gb_size_limit:  # Keep removing till size is ok
                    i_2_rm, key = self.xda_q.popleft()
                    if i_2_rm == -1:
                        continue  # Continue to next iter if encountering a padding -1
                    rmgeno_dt = self.reg_dt[i_2_rm]
                    self.cur_total_size_in_gb -= getsizeof(rmgeno_dt[key].val) / 1024 ** 3
                    rmgeno_dt.pop(key)
                    if len(self.xda_q) <= 4:
                        raise Exception('The memory footprint of current settings is too high, '
                                        'reduce blocksize and/or correction windows or increase memory limits.')
                return sda
            else:
                raise Exception(f'No srd or sda found in region i={i}, this is not supposed to happen.')

    def get_pda(self):
        if not hasattr(self, 'pda'):
            pda = self.prd.read(dtype=self.dtype)
            pda, self.stanpda = pda.standardize(return_trained=True,
                            standardizer=self.pda_standardizer())
            self.pda = pda
        return self.pda
    
    
class _DiagnosticsPlusPlotting4LinkageData():
    
    def plot_manhattan(self, *args,**kwargs):
        sst_df = self.get_sumstats_cur()
        prst.utils.plot_manhattan(sst_df, *args, **kwargs)
    
class BaseLinkageData():
    
    _onthefly_retrieval=True # These underscore options are the advanced developer options 
    _save_vars = ['L','D','R','sst_df']
    _clear_vars = ['L','D','R','Ls','Ds','Rs','Z','Zs','Di','Dis','P','Ps'] # not sure Ps actually exists in linkdata.
    _uncache_vars = ['L','D','R','Z','Di','P'] # vars that should be uncached
    uncache = False
    _validate_linkage = True
    _side2varname = {'left':'Ls','center':'Ds','right':'Rs'}
    _cross_chrom_ld = False
    _save_s2sst = True
    _npd_cnt = 0

    def __init__(self, *, sst_df=None, regdef_df=None, clsattr_dt=None, #There should be sst_df or clsattr_dt
                 
                 srd=None, sda_standardizer=Unit,
                 prd=None, pda_standardizer=Unit,
                 lrd=None, lda_standardizer=None,
                 grd=None, gda_standardizer=False,
                 
                 shift=0, cm=None, _setzero=True,
                 
                 clear_xda=True,
                 clear_linkage=False,
                 compute_sumstats=False,
                 calc_allelefreq=False,
                 intersect_apply=True,
                 
                 gb_size_limit=10., 
                 dtype='float64', 
                 check=True, 
                 verbose=False):
        
        # bim and fam df have to be supplied because pysnptools halvely
        # implemented these portions of the genetic data into their object
        # meaning that srd cannot be relied uppon
        excl_lst = ['self','kwg_dt','excl_lst']
        kwg_dt = {key: item for key, item in locals().items() if not (key in excl_lst)}
        for key, item in locals().items():
            if not (key in excl_lst): 
                self.__setattr__(key, item)
        self._kwg_dt = copy.deepcopy(kwg_dt)
        # New rule: blx have to be created from the inside
        # Perhaps later it can be made into a special load instead of a compute

        # first-checks & inits:
        if cm is not None: assert cm > 0
        if lrd is not None: raise NotImplementedError('lrd not possible atm.')
        if grd is not None:
            assert gda_standardizer or (gda_standardizer is None)
        assert type(compute_sumstats) is bool
        self.reg_dt = dict()
        self.cur_total_size_in_gb = 0.0
        self.xda_q = deque()
        [self.xda_q.append((-1,'')) for _ in range(5)]  # put 5x -1 in queue
        self.reloaded_xda_cnt = 0
        self._fn_lst = []

        # Checks            
        if srd is not None:
            assert type(sst_df) is pd.DataFrame
            self._check_xrd()
            assert isinstance(sst_df, pd.DataFrame)
            if not isinstance(regdef_df, pd.DataFrame):
                self.regdef_df = load_regdef()
            self.init_regions()
        elif clsattr_dt is not None:
            # Fill attributes in case clsattr_dt is present:
            for key, item in clsattr_dt.items():
                setattr(self, key, item)
            reg_dt=dict()
            for pre_i, geno_dt in self.reg_dt.items(): reg_dt[int(pre_i)] = geno_dt
            self.reg_dt = reg_dt # An ugly type conversion hack cause json does not allow i to be integer, but forces it to be a string.
        elif not check:
            True # disabling checks allows for blind class generation.
        else:
            raise Exception('Essentials not present')
            
    

    ###########################
    # Checking, Validation & Assertion & Admin
    
    def clone(self):
        obj = self.__class__(**self.get_params())
        obj = self._add_attr_to(obj)
        return obj 
    
    def get_params(self):
        out=copy.deepcopy(self._kwg_dt); out.pop('_excl_lst',None)
        return out
    
    def _validate_sst_df(self, sst_df, extracols=list()):
        if not isinstance(sst_df, pd.DataFrame):
            raise TypeError("sst_df, which is the summary stats dataframe, must be a pandas DataFrame")
        required_columns = ['snp', 'A1', 'A2'] + extracols
        missing_columns = [col for col in required_columns if col not in sst_df.columns]
        if missing_columns:
            raise ValueError(f"Input sumstats dataframe is missing required columns: {', '.join(missing_columns)}")
        return sst_df.reset_index(drop=True)
    
    def summary(self, return_summary=False, show=True): #name inspired by keras
        if show: printfun = print
        else: 
            def printfun(*args,**kwargs): return None
        sum_dt = {}
        printfun('Summary placeholder string (given by linkdata.summary()), this should print'
                 ' thing like genomic inflation factor, ld-score reg result (if quick enough), ld-sst match score, etc')
        if return_summary:
            return sum_dt
        
    def get_i_list(self):
        return list(self.reg_dt.keys())
    
    def _add_attr_to(self, new_obj):
        for k, v in self.__dict__.items():
            if hasattr(self, '_copy_attrs'):
                if k in self._copy_attrs:
                    setattr(new_obj, k, v)
            elif not k == '_kwg_dt':
                setattr(new_obj, k, v)
        return new_obj

    ###########################
    ## Init,Save,Load: ########       
        
    def load_linkage_allregions(self):
        for i, geno_dt in self.reg_dt.items():
            self.load_linkage_region(i=i)
        if self.verbose: print('\nDone')

    def load_linkage_region(self, *, i):
        geno_dt = self.reg_dt[i]
        store_dt = geno_dt['store_dt']

        for varname, file_dt in store_dt.items():
            storetype = file_dt['storetype'] if 'storetype' in file_dt else 'pandascast'
            if storetype == 'pandacast':
                module = importlib.import_module('.'.join(file_dt['typestr'].split('.')[:-1]))
                cname  = file_dt['typestr'].split('.')[-1]
                CurClass = getattr(module, cname) # Retrieves module.submodule.submodule.. etc
                curfullfn = os.path.join(self.curdn, file_dt['fn'])
                geno_dt[varname] = CurClass(pd.read_hdf(curfullfn, key=file_dt['key']))
                if self.verbose: print(f'loading: fn={curfullfn} key={file_dt["key"]}'+' '*50, end='\r')
            elif storetype == 'prscs':
                if varname == 'D':
                    with h5py.File(file_dt['fn'], 'r') as f:
                        geno_dt[varname] = np.array(f[file_dt['key']]).copy()
                else: raise NotImplementedError()
            else: raise ValueError(f"Storage type \'{storetype}\' not recognized.")
                
    def save(self, fn, keyfmt='ld/chrom{chrom}/i{i}/{varname}', fmt='hdf5', mkdir=False, dn=None):
        self.curdn = os.path.dirname(fn) if (dn is None) else dn
        if mkdir: os.makedirs(self.curdn, exist_ok=True)
        if (fmt == 'hdf5'):
            True
        elif (fmt == 'prscs'): 
            self.save_prscsfmt(dn=fn)
            return None
        else:
            raise Exception(f'Only hdf5 or prscs file format supported atm, not {fmt}') 
        for i, geno_dt in self.reg_dt.items():
            self.save_linkage_region(i=i, fn=fn)

        # Saving of 'logistical' data for the object
        clsattr_lst = [ 'shift', 'cm', '_setzero',
         'clear_xda', 'clear_linkage', 'compute_sumstats', 'calc_allelefreq', 
         '_onthefly_retrieval', '_save_vars', '_clear_vars', 
         'gb_size_limit', 'dtype', 'verbose', 'n_snps_total']
        geno_lst = ['regid','chrom','start','stop','start_j','stop_j',
                    'start_j_L', 'stop_j_L', 'start_j_R', 'stop_j_R','store_dt']

        def caster(arg, types):
            if np.isscalar(arg):
                if isinstance(arg, np.integer): arg = int(arg)
            if type(arg) is int: return int(arg)
            assert type(arg) in types
            return arg

        #if hasattr(self,'s'): assert (self.s.shape == (self.n_snps_total,1))
        clsattr_dt = dict(); maxlen = 20
        for key in clsattr_lst:
            var = getattr(self, key)
            if type(var) is list:
                for item in var:
                    assert type(item) in (bool, str, float, int, type(None))
                    if type(item) is str: assert (len(item) < maxlen)
            elif type(var) is str:
                    assert len(var) < maxlen
            clsattr_dt[key] = caster(var, (list, bool, float, int, str, type(None)))

        reg_dt = dict()
        for i, geno_dt in self.reg_dt.items():
            newgeno_dt = dict()
            for key in geno_lst:
                if not (key in geno_dt.keys()): continue
                newgeno_dt[key] = caster(geno_dt[key], (str, int, dict))
            reg_dt[i] = newgeno_dt
        clsattr_dt['reg_dt'] = reg_dt     
        self._fn_lst = list(np.unique(self._fn_lst))
        for curfn in self._fn_lst:
            pd.DataFrame([json.dumps(clsattr_dt)]).to_hdf(os.path.join(self.curdn, curfn), key='clsattr_dt')

        if self.verbose: print('\nDone')

    def save_linkage_region(self, *, i, fn, keyfmt='ld/chrom{chrom}/i{i}/{varname}'): 
        # using 'store' instead of 'save' to indicate a connected relationship with 
        # the files used for this storage.
        geno_dt = self.reg_dt[i]
        chrom = geno_dt['chrom']
        curdn = self.curdn
        store_dt = dict() #geno_dt['store_dt']
        for varname, var in geno_dt.items():
            if varname in self._save_vars:
                curfn  = fn.format(**locals())
                key    = keyfmt.format(**locals())
                var    = geno_dt[varname]
                vartype = type(var)
                if vartype is np.ndarray: vartype = var.dtype.type
                curfullfn = os.path.join(curdn,curfn)
                try:
                    pd.DataFrame(var).to_hdf(curfullfn, key=key)
                except TypeError as e:
                    warnings.warn( #
                        'There was a TypeError during the saving of LinkageData data. This probably happened'
                        'because there were pyarrow types present. Hence now aiming to recast to numpy types.')
                    pd.DataFrame(var).apply(lambda col: col.to_numpy()).to_hdf(curfullfn, key=key)
                file_dt = dict(fn=curfn, key=key, 
                               typestr=vartype.__module__+'.'+vartype.__name__)
                store_dt[varname] = file_dt
                self._fn_lst.append(curfn)
                if self.verbose: print(f'saving: fn={curfullfn} key={key}'+' '*50,end='\r')
        geno_dt['store_dt'] = store_dt
        
    def save_prscsfmt(self, *, out_dn=None, base_dn=None, cohort='ref', pop='adj', overwrite=True, allow_multi=False, verbose=None):
        if verbose is None: verbose = self.verbose
        print(out_dn)
        assert (out_dn is not None) ^ (base_dn is not None), 'Need to supply out_dn OR base_dn'
        if out_dn is not None: out_dn = os.path.normpath(out_dn)
        base_dn = os.path.dirname(out_dn) if base_dn is None else base_dn
        if base_dn == '': base_dn = './'
        assert os.path.exists(base_dn), f'Dir {base_dn} does not exist. You should create it.'
        assert overwrite, 'only allowed with overwrite=True'
        #print('inputs:',cohort,pop)
        join = os.path.join; split=os.path.split
        ref_bn = 'ldblk_{cohort}_{pop}' if out_dn is None else os.path.basename(out_dn)
        fnfmt0 = join(base_dn,ref_bn+'/ldblk_{cohort}_chr{chrom}.hdf5')
        fnfmt1 = join(base_dn,ref_bn+'/snpinfo_{cohort}_hm3')
        fnfmt2 = join(base_dn,ref_bn+'/snpextinfo_{cohort}_{pop}_hm3.tsv')
        # fnfmt = './ldblk_1kg_chr{chrom}.hdf5'
        ## A check check:
        if not allow_multi:
            ind = self.get_sumstats_cur()['snp'].duplicated()
            assert ind.sum() == 0, 'snp ids are not unique, this is not allowed atm. contact dev.'
        self.file_dt=dict()
        import gc; gc.collect()
        if overwrite:
            prefn = fnfmt1.format_map(locals())
            pat = join(split(prefn)[0],'ldblk_[a-z]*')
            print(pat); fn_lst = glob.glob(pat)
            if len(fn_lst) > 30: Exception(f'Using {pat} trying to remove mor than 30 files, probably too much to remove, please clear that target directory manually')
            else: print(f'Removing/Overwriting {len(fn_lst)} files.')
            for oldfn in fn_lst: os.remove(oldfn) # os.remove only does files here (fails on dir input)
        def fun(chrom, **kwg):
            chrom=int(chrom)
            kwg['chrom'] = chrom
            fn = fnfmt0.format_map(kwg)
            os.makedirs(os.path.dirname(fn), exist_ok=True)
            self.file_dt[chrom] = h5py.File(fn, 'w')
            return self.file_dt[chrom]
        blkcnt_dt = defaultdict(lambda:1)
        for i, geno_dt in prst.utils.get_pbar(list(self.reg_dt.items()), colour='yellow'):
            if not 'chrom' in geno_dt:
                clst = geno_dt['sst_df']['chrom'].unique()
                assert len(clst) == 1
                geno_dt['chrom'] = clst[0]
            chrom = int(geno_dt['chrom'])
            #if i ==1: jkergjk
            f = self.file_dt[chrom] if chrom in self.file_dt.keys() else fun(**locals())
            group = f.create_group(f'blk_{blkcnt_dt[chrom]}')
            blkcnt_dt[chrom] += 1
            D = self.get_linkage_region(i=i)
            #group.create_dataset('ldblk', data=D) # old
            #group.create_dataset('ldblk', data=D, compression="gzip", compression_opts=7, shuffle=True, chunks=True) ## <<<---------- keep looking here
            group.create_dataset('ldblk', data=D, compression='lzf', shuffle=True, chunks=True) ## <<<---------- keep looking here
            arr = np.array(geno_dt['sst_df']['snp'].to_numpy(),dtype='S')
            group.create_dataset('snplist', data=arr) #, dtype='S')
        for f in self.file_dt.values(): f.close()
        sst_df = self.get_sumstats_cur()
        cols = ['chrom','snp','pos','A1','A2','maf_ref']
        save_df = sst_df[cols]
        save_df.columns = ['CHR','SNP','BP','A1','A2','MAF']
        save_df = pd.DataFrame(save_df)
        save_df.to_csv(fnfmt1.format_map(locals()), sep='\t', index=False)
        if 'af_A1_ref' in sst_df.columns:
            sst_df[cols+['af_A1_ref']].to_csv(fnfmt2.format_map(locals()), sep='\t', index=False)
            

    ########################### 
    ## Compute: ############### 
    

    ############################
    ## Retrieve: ###############
    
    # Local Linkage: ############
    if True:
    
        def retrieve_linkage_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.retrieve_linkage_region(i=i)
            if self.verbose:   print('\nDone')
            if self.clear_xda: self.clear_all_xda()

        def retrieve_linkage_region(self, *, i, force=True):
            geno_dt = self.reg_dt[i] 
            if 'store_dt' in geno_dt.keys():
                self.load_linkage_region(i=i)
                return None
            
        def retrieve_slicedlinkage_region(self, *, i, varname='Ds'):
            geno_dt = self.reg_dt[i]
            sst_df = geno_dt['sst_df']
            if 'Ds' == varname:
                D=self.get_specified_data_region(i=i, varname=varname.rstrip('s'), checkdims=False)
                bidx = sst_df['bidx'] if 'bidx' in sst_df.columns else np.arange(len(D)) # This will fail if D and sst_df dont match in dims
                preDs = D[bidx][:,bidx]; info=False
                Ds, info = prst.io.validate_linkage(preDs, return_info=True) if self._validate_linkage else preDs
                if info: self._npd_cnt = self._npd_cnt + 1 # Like this to prevent class attribute to get used (not +=)
                #msg = 'Non posidefinite LD matrices detected, applying correction.'
                #if info: warnings.warn(msg)
                geno_dt['Ds']=Ds
            else: raise NotImplementedError()

        def retrieve_precision_region(self, *, i, store_cnum=True, perc=1e-6, maxcond=1e6):
            D = self.get_linkage_region(i=i)
            U,s,Vt = linalg.svd((D+np.eye(len(D))*perc)/(1+perc)) 
            cond = s.max()/s.min()
            if store_cnum: self.reg_dt[i]['cond'] = cond
            if not cond <= maxcond:
                self.reg_dt[i]['Di'] = linalg.pinv(D)
            else:
                self.reg_dt[i]['Di'] = (U*(1/s))@Vt
            
        def retrieve_halfmatrix_region(self, *, i, varname):
            if varname == 'Dihalf':
                Di = self.get_precision_region(i=i)
                Dh = self.get_specified_data_region(i=i, varname='Dhalf') 
                #U,S,Vt=np.linalg.svd(Di,full_matrices=False, hermitian=True); 
                #Dihalf=(U*np.sqrt(S))
                self.reg_dt[i]['Dihalf'] = Di@Dh
            elif varname == 'Dhalf':
                D = self.get_linkage_region(i=i)
                U,S,Vt=np.linalg.svd(D,full_matrices=False, hermitian=True); 
                Dhalf=(U*np.sqrt(S))
                self.reg_dt[i]['Dhalf'] = Dhalf
            else:
                raise Exception(f'Option \'{varname}\' not recognized')
                
        def retrieve_startstopjs_region(self, *, i, varname):
            geno_dt = self.reg_dt[i]
            sst_df = geno_dt['sst_df']
            if varname in ['start_j','stop_j']:
                if 'idx' in sst_df:
                    idx = sst_df['idx'].to_numpy()
                    assert (np.arange(idx[0], idx[0]+len(idx)) == idx).all(), 'Something wrong with reference.'
                    geno_dt['start_j'] = idx[0]; geno_dt['stop_j'] = idx[-1]+1
                    
                          
    # SumStat: ##############
    if True:

        def retrieve_sumstats_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.retrieve_sumstats_region(i=i)

        def retrieve_sumstats_region(self, *, i):
            geno_dt = self.reg_dt[i] 
            sst_df  = geno_dt['sst_df']
            if not 'std_ref' in sst_df:
                maf = sst_df['maf_ref']
                sst_df['std_ref'] = np.sqrt(2.0*maf*(1.0-maf))
            if 'beta_mrg' in geno_dt.keys():
                return None # Sumstat present so no need to compute anything.
            elif 'beta_mrg' in sst_df.columns:
                geno_dt['beta_mrg'] = sst_df[['beta_mrg']].to_numpy()
                return None
            sst_df['beta_mrg'], sst_df['n_eff'], sst_df['std_sst'] = self.compute_sumstats_region(i=i, 
                                                        return_n_eff=True, return_allele_std=True)

        def retrieve_ldscores_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.retrieve_ldscores_region(i=i)

        def retrieve_ldscores_region(self, *, i):
            geno_dt = self.reg_dt[i]
            sst_df = geno_dt['sst_df']
            if not 'lds' in sst_df.columns:
                newsst_df = self.compute_ldscores_region(i=i)
                geno_dt['sst_df'] = newsst_df
            if self.clear_linkage:
                self.clear_linkage_region(i=i)

    # Clearing/uncache Functions: #####
    if True:

        def clear_all_xda(self):
            while len(self.xda_q) != 0:
                i_2_rm, key = self.xda_q.popleft()
                if i_2_rm == -1:
                    continue  # Continue to next iter if encountering a padding -1
                rmgeno_dt = self.reg_dt[i_2_rm]
                self.cur_total_size_in_gb -= getsizeof(rmgeno_dt[key].val) / 1024 ** 3
                rmgeno_dt.pop(key)
            [self.xda_q.append((-1,'')) for _ in range(5)]  # put 5x -1 in queue
            
        def clear_linkage_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.clear_linkage_region(i=i)
            prst.utils.clear_memory()
            if self.verbose: print('\nDone') 

        def clear_linkage_region(self, *, i):
            geno_dt = self.reg_dt[i]
            key_lst = list(geno_dt.keys())
            for key in key_lst:
                if key in self._clear_vars:
                    geno_dt.pop(key, None)
                    
        def uncache_region(self,*,i):
            geno_dt = self.reg_dt[i]
            for key in list(geno_dt.keys()): # list() bit is needed to stop size-changed-during-iteration error
                if key in self._uncache_vars:
                    geno_dt.pop(key, None)
        
        def uncache_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.uncache_region(i=i)
            prst.utils.clear_memory()
            if self.verbose: print('\nDone') 
            
            
    ############################ 
    ## Get: #################### 
    
    # Local Linkage: ###########
    if True:

        def get_specified_data_region(self, *, i, varname, checkdims=True):
            try:
                return self.reg_dt[i][varname]
            except KeyError as e:
                if '_glocal' in varname:
                    self.retrieve_linkage_region_glocalshiftwindow(i=i)
                elif varname in 'LDR':
                    self.retrieve_linkage_region(i=i)
                elif varname in 'Ls-Ds-Rs': # comparable or faster than list() formulaton.
                    self.retrieve_slicedlinkage_region(i=i,varname=varname)
                elif varname == 'Z':
                    self.retrieve_linkage_region_global(i=i)
                elif varname == 'Di':
                    self.retrieve_precision_region(i=i)
                elif 'half' in varname:
                    self.retrieve_halfmatrix_region(i=i, varname=varname)
                elif '_j' in varname:
                    self.retrieve_startstopjs_region(i=i, varname=varname)
                else:
                    raise Exception(f'varname={varname}, on-the-fly not possilbe for this variable.')
                try:
                    # Being here means a retrieval was nessicary: 
                    var = self.reg_dt[i][varname]
                    msg = (f'Retrieval of {varname} was nessicary, but the result did not have '
                    'the same dimension as the sumstat for the region, something went wrong, '
                    'probably linkage i.e. LD was retrieved and later linkage data was sliced/merged. '
                    '(can make that work. just did not implement it yet.)')
                    if not any(elem.isupper() for elem in varname): checkdims=False
                    if checkdims: assert var.shape[0] == self.reg_dt[i]['sst_df'].shape[0], msg
                    if self.uncache: self.uncache_region(i=i)
                    return var
                except KeyError as e:
                    print('Failed, eventough on-the-fly retrieval was attempted')
                    raise e
            else:
                    raise Exception('on-the-fly retrieval blocked, set _onthefly_retrieval=True if desired')
        
        def get_range_region(self,*, i, side='center'):
            suffix = '_'+self._side2varname[side] if side != 'center' else ''
            k0='start_j'+suffix; k1='stop_j'+suffix
            try:
                return self.reg_dt[i][k0], self.reg_dt[i][k1]
            except:
                a=self.get_specified_data_region(i=i, varname=k0)
                b=self.get_specified_data_region(i=i, varname=k1)
                return a,b
        
        def get_linkage_region(self, *, i, side='center'):
            #_side2varname = {'left':'Ls','center':'Ds','right':'Rs'} # als geheugen steuntje
            return self.get_specified_data_region(i=i, varname=self._side2varname[side])
        
        #@something # This does nothing, just some small decorator xps
        def get_precision_region(self, *, i, side='center'):
            if not side == 'center': raise Exception(f'Option not valid: {side}')
            return self.get_specified_data_region(i=i, varname='Di')
            
    # Sumstats: #################
    if True:
        
        def get_s(self):
            sst_df = self.get_sumstats_cur()
            try: 
                s = sst_df[['s']].values
                assert np.isnan(s).sum() == 0
                return s
            except:
                stansda = self.get_stansda()
                s = self.get_stansda().stats[:,[1]]
                self.s = s
                return s
            
        def get_sumstats_cur(self, maf=False):
            if maf: self.retrieve_maf_allregions()
            sst_df_lst = []
            for i, geno_dt in self.reg_dt.items():
                sst_df = geno_dt['sst_df']
                sst_df_lst.append(sst_df)
            if len(sst_df_lst) == 0:
                msg = ('The (matched) reference is completely empty (i.e. not a single snp). This could mean an empty '
                       'reference or no overlaping variants between inputs (e.g. target-bim, sumstat & LD reference) or '
                      'forinstance that --chrom 3 was selected, but that chromosome 3 is not in the reference.')
                raise RuntimeError(msg)
            sst_df = pd.concat(sst_df_lst, axis=0)
            return sst_df

        def get_stansda(self, standardizer='unit'):
            if not standardizer=='unit': raise NotImplementedError('contact dev')
            
            if hasattr(self, 'stansda'):
                if type(self.stansda) is UnitTrained:
                    return self.stansda
                else:
                    raise NotImplementedError('contact dev')
                    
            standardizer_lst = []
            for i, geno_dt in self.reg_dt.items():
                #(not 'stansda' in geno_dt.keys())
                if (not type(geno_dt['stansda']) is UnitTrained) & self._onthefly_retrieval:
                    self.retrieve_linkage_region(i=i)
                if type(geno_dt['stansda']) is UnitTrained:
                    standardizer_lst.append(geno_dt['stansda'])
                else:
                    raise Exception('No standardizer detected. Compute this first. Contact dev if issue persists.')

            assert np.all([type(stan) is UnitTrained for stan in standardizer_lst])            
            sid = np.concatenate([stan.sid for stan in standardizer_lst])
            assert np.unique(sid).shape[0] == sid.shape[0]

            stats = np.concatenate([stan.stats for stan in standardizer_lst])
            combined_unit_standardizer = UnitTrained(sid, stats)
            self.stansda = combined_unit_standardizer
            return combined_unit_standardizer
        
        def get_allele_standev(self, source='ref'):
            self.retrieve_sumstats_allregions()
            sst_df = self.get_sumstats_cur()
            if source == 'ref':
                return sst_df['std_ref'].to_numpy()[:,np.newaxis]
            elif source == 'sst':
                return sst_df['std_sst'].to_numpy()[:,np.newaxis]
            elif source == 'afsst':
                af = sst_df['af_A1_sst']
                maf = np.minimum(af, 1.0 - af)
                sst_df['std_afsst'] = np.sqrt(2.0*maf*(1.0-maf))
                return sst_df['std_afsst'].to_numpy()[:, None]
            else:
                raise Exception(f'Allele-standev data source {ref} not found.')
        
        def get_beta_marginal(self):
            beta_mrg_lst = []
            for i, geno_dt in self.reg_dt.items():
                beta_mrg = self.get_beta_marginal_region(i=i)
                beta_mrg_lst.append(beta_mrg)
            beta_mrg_full = np.concatenate(beta_mrg_lst)
            return beta_mrg_full
        
        def get_beta_marginal_region(self, *, i):
            if not 'beta_mrg' in self.reg_dt[i]:
                self.retrieve_sumstats_region(i=i)
            return self.reg_dt[i]['beta_mrg']
        
        @property
        def shape(self):
            return (self.n_snps_total, len(self.reg_dt.keys()))
        
        @property
        def n_snps_total(self):
            if not hasattr(self, '_n_snps_total'): 
                self._n_snps_total = self.get_sumstats_cur().shape[0]
            return self._n_snps_total
    
    ############################
    ## Set: ####################
    
    def set_sumstats(self, sst_df, merge=False, check=True, extracols=['i']):
        if not 'SparseLinkageData' in str(self.__class__): raise NotImplementedError('this is not implemented yet; non sparse linkagedata, has pre-sliced LD, which wont match up.')
        sst_df = self._validate_sst_df(sst_df, extracols=extracols)
        if merge:
            cur_sst_df = self.get_sumstats_cur()
            if len(set(cur_sst_df) - set(sst_df)) != 0 or check:
                raise NotImplementedError('set_sumstats method is in an alpha state, please contact dev.') 
        else:
            if check:
                cur_sst_df = self.get_sumstats_cur()
                assert len(sst_df) == len(cur_sst_df), "Lengths of sst_df (sumstat) and current sst_df do not match"
                assert all(sst_df['snp'] == cur_sst_df['snp']), "snp column values do not match, also A1 and A2 columns must match"
                assert all(sst_df['A1']  == cur_sst_df['A1']),  "A1 column values do not match, also A2 columns must match"
                assert all(sst_df['A2']  == cur_sst_df['A2']),  "A2 column values do not match"

        for i, df in sst_df.groupby('i'):
            self.reg_dt[i]['sst_df'] = df

        return self
    
    def set_population(self, pop):
        assert type(pop) is str
        self.pop = pop.upper()
        return self
    

class RefLinkageData(BaseLinkageData, _DiagnosticsPlusPlotting4LinkageData):
    
    uncache=True
    _extradropdupcols = ['i','bidx', 'blkid'] 
    def get_extradropdupcols(self):
        return self._extradropdupcols
    
    @classmethod
    def _save_snpregister(cls, *, ref, reg_bn, chrom='all', verbose=False):

        self = cls(check=False)
        assert chrom == 'all'
        if chrom == 'all': chrom = '*'
        
        pat0 = os.path.join(ref,'snpextinfo*')
        pat1 = os.path.join(ref,'snpinfo*')
        lst = glob.glob(pat0)+glob.glob(pat1)
        if len(lst) == 0:
            raise FileNotFoundError(f'No snpinfo file(s) found using the prefix {pat1}, are you sure the --ref directory is a proper prstools reference?')
        ref_fn = lst[0]
        ref_df = prst.io.load_ref(ref_fn, verbose=False)
        out_fn = os.path.join(ref,reg_bn)
        if verbose: print(f'Creating snp register for efficient operations (only once)[{os.path.basename(ref_fn)}] @ {out_fn} ', end='', flush=True)

        # Load all snps-ids if no blk # if not 'blk' in ref_df.columns and not blk:
        chr_dt = {}; reg_dt = {}; tlst = []; i = 0
        for cur_chrom in list(range(1,24)): # SOMETHING NEEDS TO GET FIXED HERE!!!!!!!!!
            h5lst = glob.glob(os.path.join(ref,f'*chr{cur_chrom}.hdf5'))
            if len(h5lst) == 0: continue
            assert len(h5lst) == 1, 'Issue with reference'
            hdf_chr = h5py.File(h5lst[0], 'r')
            chr_dt[cur_chrom] = hdf_chr
            for num in range(1,len(hdf_chr)+1):
                blkid = f'blk_{num}'
                snps = np.array(hdf_chr[blkid]['snplist'][:])
                #snps = np.array(hdf_chr[blkid]['snplist'], dtype=str) # Old version no 3.6 compat.
                dt = {0:snps,'i':i,'blkid':blkid,'bidx':np.arange(len(snps))}
                sst_df = pd.DataFrame(dt)
                h5 = hdf_chr[blkid]['ldblk']
                geno_dt = dict(blkid=blkid, sst_df=sst_df, store_dt=dict(D=h5)) #, newthing=newthing)
                self.reg_dt[i] = geno_dt
                i+=1
        sst_df = self.get_sumstats_cur().reset_index(drop=True) 
        assert (ref_df['snp'].to_numpy().astype(str) == sst_df[0].to_numpy().astype(str)).all()
        ref_df[['i','bidx','blkid']]=sst_df[['i','bidx','blkid']]
        if not 'std_ref' in ref_df.columns:
            maf = ref_df['maf_ref']
            ref_df['std_ref'] = np.sqrt(2.0*maf*(1.0-maf))

        # Put the staged data back into reg_dt
        new_reg_dt={}
        for i, (old_i, cur_df) in enumerate(ref_df.groupby('i', sort=True)):
            geno_dt = self.reg_dt[old_i]
            cur_df['i'] = i
            geno_dt['sst_df'] = cur_df
            new_reg_dt[i] = geno_dt
        self.reg_dt = new_reg_dt

        ref_df=self.get_sumstats_cur()
        #ref_df.to_csv(out_fn, sep='\t', index=False)
        prst.io._pd_to_atomizer(to_file=ref_df.to_csv, fn=out_fn, sep='\t', index=False)
        if verbose: print('-> Done')

    @classmethod
    def from_ref(cls, ref, chrom='all', return_locals=False, reg_bn='snpregister.tsv', storetype='prscs', verbose=False, **kwg):
        ref = prst.utils.validate_path(ref=ref, must_exist=True, handle_prstdatadir='allow', verbose=False)
        reg_fn = os.path.join(ref, reg_bn)
        if not os.path.isfile(reg_fn): cls._save_snpregister(ref=ref, reg_bn=reg_bn, verbose=verbose)
        lst=glob.glob(os.path.join(ref,'snpextinfo*'))
        if len(lst)>0:
            assert len(lst) == 1, f"Only put 1 snpextinfo file in {ref}. Found {len(lst)} : {lst}"
            ## Not actually doing this bit just yet:
#             ext_df = prst.io.load_sst(lst[0], n_gwas=None, calc_beta_mrg=False, nrows=5, check=False, verbose=False)
#             tst_df = prst.io.load_sst(reg_fn, n_gwas=None, calc_beta_mrg=False, nrows=5, check=False, verbose=False)
#             if not all(col in tst_df for col in ext_df.columns):
#                 cls._save_snpregister(ref=ref, reg_bn=reg_bn, verbose=verbose)

        ref_df = prst.io.load_ref(reg_fn, chrom=chrom, verbose=verbose) ## <--- here the magic for chrom slicing takes place..
        if not 'check' in kwg: kwg['check']=False
        self = cls(**kwg)
        chrom_fn_dt={}
        for chrom in ref_df['chrom'].unique():
            lst = glob.glob(os.path.join(ref,f'*chr{chrom}.hdf5'))
            assert len(lst) == 1, 'Something wrong with reference, consider redownload.'
            chrom_fn_dt[chrom] = lst[0]
            
        # make new reg_dt and put back: ################################################################ COMBINE
        new_reg_dt = {}
        for i, (old_i, cur_df) in enumerate(ref_df.groupby('i', sort=True)):
            #row = cur_df.iloc[0] # slow line! 
            chrom = cur_df['chrom'].iloc[0] # better than slicing the first row
            blkid = cur_df['blkid'].iloc[0]
            cur_df['i'] = i # if stuff was dropped already from the reference i needs to be updated.
            file_dt = dict(fn=chrom_fn_dt[chrom], key=f'{blkid}/ldblk', storetype=storetype)
            store_dt=dict(D=file_dt)
            geno_dt = dict(sst_df=cur_df,store_dt=store_dt)
            new_reg_dt[i] = geno_dt
            i += 1
        self.reg_dt=new_reg_dt
        return self
    
    @classmethod
    def from_cli_params(cls, *, ref, target, sst, n_gwas=None, chrom='*', verbose=False, colmap=None, pop=None, cli=True, 
                        sstrename_dt=dict(maf='maf_sst',af_A1='af_A1_sst'), **kwg): 
        # Basic checks:
        if target is not None: tsttarget = '.'.join(target.split('.')[:-1])+'.bim' if (target.split('.')[-1] in ('bim','fam','bed')) else target+'.bim'
        else: tsttarget=None
        ref, sst, tsttarget = prst.utils.validate_path(ref=ref, sst=sst, tsttarget=tsttarget, must_exist=True, handle_prstdatadir='allow', verbose=verbose)
        msg=f'Population argument specified (pop={pop}), but for this approach this information is currently not used.'
        if pop is not None and pop != 'pop': warnings.warn(msg)
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        
        # Loading:
        orisst_df = prst.load_sst(sst, calc_beta_mrg=True, n_gwas=n_gwas, colmap=colmap, verbose=verbose, cli=cli)
        target_df, _ = prst.load_bimfam(target, fam=False, chrom=chrom, start_string = 'Loading target file.    ', verbose=verbose) if target else (None,None)
        linkdata = cls.from_ref(ref, chrom=chrom, verbose=verbose, sst_df=orisst_df, **kwg)
        ref_df   = linkdata.get_sumstats_cur()
        msg = (f'\033[1;31mWARNING: The size of the reference (={ref_df.shape[0]} snps) is much smaller than the sumstat (={orisst_df.shape[0]} snps). '
               'Are you sure you are using the right reference and not the reference example?\033[0m')
        if ref_df.shape[0] < 1e4 and orisst_df.shape[0] > 1e5: warnings.warn(msg) 
        orisst_df.rename(columns=sstrename_dt, inplace=True)

        # Matching:
        if verbose: print('Matching sumstat & reference ', end='', flush=True)
        ddups = linkdata.get_extradropdupcols()
        sst_df = prst.merge_snps(ref_df, orisst_df, flipcols=['beta_mrg','beta'], handle_missing='filter', extradropdupcols=ddups)
        if verbose and target: print('& target. ', end='', flush=True)
        sst_df = prst.merge_snps(sst_df, target_df, flipcols=[], handle_missing='filter', extradropdupcols=ddups, warndupcol=True) if target else sst_df
        n_match = sst_df.shape[0]
        msg = (f'-> {n_match:,} common variants after matching ' +
                          f'reference ({(n_match/max(ref_df.shape[0],1))*100:.1f}% incl.), ' +
                          (f'target ({(n_match/max(target_df.shape[0],1))*100:.1f}% incl.) and ' if target else 'and ') +
                           f'sumstat ({(n_match/max(orisst_df.shape[0],1))*100:.1f}% incl.).')
        if verbose: print(msg)
         # generate spacing between loading and fit()
        if verbose and hasattr(orisst_df, 'msg') and cli and type(orisst_df.msg) is str: print(orisst_df.msg, '\n')
        elif verbose: print('\n')
        linkdata = linkdata.merge(sst_df, warndupcol=False, inplace=True)
        return linkdata
        
    def merge(self, sst_df, inplace=False, flipcols='auto', drop=True, aligned=False, check=True, handle_missing='filter', warndupcol=True, 
              extradropdupcols='auto',dropalldupcols=False):
        from prstools import merge_snps
        assert extradropdupcols == 'auto', 'only avail option atm is \'auto\'' 
        assert drop is True
        assert check is True
        ddups = self.get_extradropdupcols()
        if flipcols == 'auto':
            flipcols = [col for col in sst_df.columns if col in ['beta','beta_mrg', 'allele_weight']]
        cur_sst_df = self.get_sumstats_cur().drop(flipcols+['n_eff'], errors='ignore', axis=1)
        if not aligned:
            new_sst_df = merge_snps(cur_sst_df, sst_df, flipcols=flipcols, handle_missing=handle_missing, 
                extradropdupcols=ddups, warndupcol=warndupcol, dropalldupcols=dropalldupcols)
        else:
            raise NotImplementedError
            check(); new_sst_df = pd.concat() # maybe also works with filter.
            
        # make new reg_dt and put back: ################################################################ COMBINE
        new_sst_df = prst.io.validate_dataframe_index(new_sst_df)
        new_sst_df['idx'] = new_sst_df.index # not sure what i put this in again, explain please
        nreg_dt = {}
        for i_new, (i_old, df) in enumerate(new_sst_df.groupby('i', sort=True)):
            geno_dt = self.reg_dt[i_old]
            if not inplace: 
                df = df.copy()
                geno_dt = copy.deepcopy(geno_dt)
            df['i'] = i_new
            geno_dt['sst_df'] = df
            geno_dt.pop('beta_mrg', None)
            nreg_dt[i_new] = geno_dt
        flinkdata = self if inplace else self.clone()
        flinkdata.reg_dt = nreg_dt
        return flinkdata
    
    def xs(self, keys, on='i', sort=True, makecopy=True):
        assert on=='i', "For now only i is allowed for linkdata slicing/ xs\'ing"
        assert sort is True, 'atm input keys and all stuff needs to be sorted'
        assert makecopy is True, 'assuming copy only for now'
        keys = np.sort(keys)
        new_sst_df = pd.concat([self.reg_dt[key]['sst_df'] for key in keys])
        new_sst_df = new_sst_df.reset_index(drop=True)
        new_sst_df = prst.io.validate_dataframe_index(new_sst_df)
        new_sst_df['idx'] = new_sst_df.index # not sure what i put this in again, explain please
        nreg_dt = {}
        for i_new, (i_old, df) in enumerate(new_sst_df.groupby('i', sort=True)):
            geno_dt = self.reg_dt[i_old]
            if makecopy:
                df = df.copy()
                geno_dt = copy.deepcopy(geno_dt)
            df['i'] = i_new
            geno_dt['sst_df'] = df
            geno_dt.pop('beta_mrg', None)
            nreg_dt[i_new] = geno_dt
        newlinkdata = self.clone()
        newlinkdata.reg_dt = nreg_dt
        return newlinkdata
        
    def groupby(self, by=None, sort=True, warndupcol=True, skipempty=True, needmerge=None):
        assert skipempty, 'Only option is to skip the empty groupbys for now.'
        import time, itertools
        sst_df = self.get_sumstats_cur()
        groupings = list(sst_df.groupby(by, sort=sort))
        sst_df=[]
        sets = [set(cdf['i'].unique()) for grp, cdf in groupings]
        if needmerge is None: needmerge = any(s1 & s2 for (i, s1), (j, s2) in itertools.combinations(enumerate(sets), 2))
        for grp, cdf in groupings:
            if needmerge:
                nlink = self.merge(cdf.reset_index(), warndupcol=warndupcol, dropalldupcols=True, inplace=False)
            else:
                keys=np.sort(cdf['i'].unique())
                nlink = self.xs(keys, on='i')
            if len(nlink.get_i_list()) > 0: yield grp, nlink
            else: warnings.warn(f'Grouping by {by} specifically for {by}={grp} led to an empty LD + sumstat (i.e. not data), so skipping {by}={grp}')

    
class SparseLinkageData(BaseLinkageData):
    
    @classmethod
    def from_cli_params(cls, *, ref, target, sst, n_gwas, chrom='*', pop=None, verbose=True, return_locals=False, pyarrow=True, colmap=None, **kwg):
        if pop is None: raise Exception('Population not specified. please specify population.')
        bim, _ = prst.io.load_bimfam(target,fam=False)
        pop = pop.upper()
        reg_dt, sst_df, _extra = prst.io._load_sparse_data(chrom=chrom, ref=ref, sst=sst, pop=pop, n_gwas=n_gwas, target=target, pyarrow=pyarrow,
                                           return_locals=return_locals, colmap=colmap, verbose=verbose)
        linkdata = cls(check=False, verbose=verbose)
        linkdata.reg_dt = reg_dt
        linkdata._extra = _extra
        linkdata.set_population(pop.split('-')[-1])
        assert np.all(linkdata.get_allele_standev('ref') != 0)
        return linkdata

    def retrieve_linkage_region(self, *, i, regu=0):
        # Your additional logic before calling the parent method
        geno_dt = self.reg_dt[i]
        if 'P' in geno_dt and (not 'D' in geno_dt):

            # Load required indices & data:
            rsst_df = self.get_specified_data_region(i=i, varname='sst_df')
            idx_reg = rsst_df.index.to_numpy()
            idx_inc = rsst_df['pindex'].to_numpy()
            P = self.get_specified_data_region(i=i, varname='P')
            # Compute LD matrix & store:
            nonzero_ind = P.diagonal() != 0
            Ps = P[nonzero_ind][:,nonzero_ind] #.toarray()
            from sksparse.cholmod import cholesky
            I = sp.sparse.diags(np.ones(Ps.shape[0])).tocsc()
            solver = cholesky(Ps)
            Dfull = solver(I).toarray()
            #Dfull = linalg.pinv(Ps+np.eye(len(Ps))*regu) 
            inc_ind = np.zeros(P.shape[0], dtype='bool')
            inc_ind[idx_inc] = True
            ind = inc_ind[nonzero_ind]
            D = Dfull[ind][:,ind]
            geno_dt['D'] = D
            
    def retrieve_precision_region(self, *, i, store_cnum=True, perc=1e-6, maxcond=1e10):
        D = self.get_linkage_region(i=i)
        U,s,Vt = linalg.svd((D+np.eye(len(D))*perc)/(1+perc))
        cond = s.max()/s.min()
        if store_cnum: self.reg_dt[i]['cond'] = cond
        assert cond <= maxcond
        self.reg_dt[i]['Di'] = (U*(1/s))@Vt
        
    def retrieve_sumstats_region(self, *, i):
        pop = self.pop
        geno_dt = self.reg_dt[i] 
        sst_df  = geno_dt['sst_df']
        def fun(arg):
            if arg > 0.5: return 1-arg
            else: return arg
        maf = sst_df[pop].apply(fun)
        sst_df['maf_ref'] = maf 
        sst_df['std_ref'] = np.sqrt(2.0*maf*(1.0-maf))
        return super().retrieve_sumstats_region(i=i)
        
        
if not '__file__' in locals():
    import sys
    if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+['ipykernel_launcher.py' in sys.argv[0]]):
        with open('../prstools/linkage/_base.py', 'w') as loadrf: loadrf.write(In[-1])
        print('Written to:', loadrf.name)