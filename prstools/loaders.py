# __all__ = ['LinkageData'] # ok this does not seem to lead to the desired stuff
#| export
# %%writefile ../prstools/loaders.py
#!/usr/bin/env python  

"""
LinkageData
durr tst
"""

import os, re
import scipy as sp
import numpy as np
import pandas as pd
from scipy import linalg
from sys import getsizeof

import warnings, importlib, json, os, glob, copy
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

class RemovalStagingLinkageData():
    


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
    
class BaseLinkageData():
    
    _onthefly_retrieval=True # These underscore options are the advanced developer options 
    _save_vars = ['L','D','R','sst_df']
    _clear_vars = ['L','D','R','Ls','Ds','Rs','Z','Zs','Di','Dis','P','Ps'] # not sure Ps actually exists in linkdata.
    _cache_unsliced_or_something= 'tbd'
    _side2varname = {'left':'Ls','center':'Ds','right':'Rs'}
    _cross_chrom_ld = False
    _save_s2sst = True

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
        return self.__class__(**self.get_params())
    
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
        
    def save_prscsfmt(self, *, dn, cohort, pop, verbose=None):
        if verbose is None: verbose = self.verbose
        print(cohort)
        fnfmt0 = 'ldblk_{cohort}_{pop}/ldblk_{cohort}_chr{chrom}.hdf5'
        fnfmt1 = 'ldblk_{cohort}_{pop}/snpinfo_{cohort}_hm3'
        # fnfmt = './ldblk_1kg_chr{chrom}.hdf5'
        self.file_dt=dict()
        import gc; gc.collect()
        def fun(chrom, **kwg):
            chrom=int(chrom)
            kwg['chrom'] = chrom
            fn = fnfmt0.format_map(kwg)
            fn = os.path.join(dn,fn)
            os.makedirs(os.path.dirname(fn), exist_ok=True)
            self.file_dt[chrom] = h5py.File(fn, 'w')
            return self.file_dt[chrom]
        blkcnt_dt = defaultdict(lambda:1)
        for i, geno_dt in tqdm(list(self.reg_dt.items())):
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
            group.create_dataset('ldblk', data=D)
            arr = np.array(geno_dt['sst_df']['snp'].to_numpy(),dtype='S')
            group.create_dataset('snplist', data=arr) #, dtype='S')
        for f in self.file_dt.values(): f.close()
        sst_df = self.get_sumstats_cur()
        # save_df.loc[:,'MAF'] = save_df['BP']
        save_df = sst_df[['chrom','snp','pos','A1','A2','maf_ref']]
        save_df.columns = ['CHR','SNP','BP','A1','A2','MAF']
        save_df = pd.DataFrame(save_df)
        #save_df['MAF'] = maf.flatten() #its included now
        save_df.to_csv(os.path.join(dn, fnfmt1.format_map(locals())), sep='\t', index=False)    
            
    
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
                bidx = sst_df['bidx'] if 'bidx' in sst_df.columns else np.arange(len(D))
                geno_dt['Ds']=D[bidx][:,bidx]
            else: raise NotImplementedError()

        def retrieve_precision_region(self, *, i, store_cnum=True, perc=1e-6, maxcond=1e10):
            D = self.get_linkage_region(i=i)
            U,s,Vt = linalg.svd((D+np.eye(len(D))*perc)/(1+perc))
            cond = s.max()/s.min()
            if store_cnum: self.reg_dt[i]['cond'] = cond
            assert cond <= maxcond
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
            if 'beta_mrg' in geno_dt.keys():
                return None # Sumstat present so no need to compute anything.
            elif 'beta_mrg' in sst_df.columns:
                geno_dt['beta_mrg'] = sst_df['beta_mrg']
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

    # Clearing Functions: #####
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
                    geno_dt[key]=None
                    del geno_dt[key]
                    #geno_dt.pop(key)
            
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
                    'the same dimension as the sumstat for the region, something went wrong')
                    if not any(elem.isupper() for elem in varname): checkdims=False
                    if checkdims: assert var.shape[0] == self.reg_dt[i]['sst_df'].shape[0], msg
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
            if not 'beta_mrg' in self.reg_dt[i]['sst_df']:
                self.retrieve_sumstats_region(i=i)
            return self.reg_dt[i]['sst_df'][['beta_mrg']].to_numpy()
        
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
    

class RefLinkageData(BaseLinkageData):
    
    _extradropdupcols = ['i','bidx', 'blkid']
    def get_extradropdupcols(self):
        return self._extradropdupcols
    
    @classmethod
    def _save_snpregister(cls, *, ref, reg_bn, chrom='all', verbose=False):

        self = cls(check=False)
        assert chrom == 'all'
        if chrom == 'all': chrom = '*'
        pat = os.path.join(ref,'snpinfo*')
        lst = glob.glob(pat)
        if len(lst) == 0:
            raise FileNotFoundError(f'No snpinfo file(s) found using the prefix {pat}, are you sure the --ref directory is a proper prstools reference?')
        ref_fn = lst[0]
        ref_df = load_ref(ref_fn, verbose=False)
        out_fn = os.path.join(ref,reg_bn)
        if verbose: print(f'Creating snp register for efficient operations (only once) @ {out_fn} ', end='')

        # Load all snps-ids if no blk # if not 'blk' in ref_df.columns and not blk:
        chr_dt = {}; reg_dt = {}; tlst = []; i = 0
        for cur_chrom in list(range(1,23)): 
            h5lst = glob.glob(os.path.join(ref,f'*chr{cur_chrom}.hdf5'))
            if len(h5lst) == 0: continue
            assert len(h5lst) == 1, 'Issue with reference'
            hdf_chr = h5py.File(h5lst[0], 'r')
            chr_dt[cur_chrom] = hdf_chr
            for num in range(1,len(hdf_chr)+1):
                blkid = f'blk_{num}'
                snps = np.array(hdf_chr[blkid]['snplist'], dtype=str)
                dt = {0:snps,'i':i,'blkid':blkid,'bidx':np.arange(len(snps))}
                sst_df = pd.DataFrame(dt)
                h5 = hdf_chr[blkid]['ldblk']
                geno_dt = dict(blkid=blkid, sst_df=sst_df, store_dt=dict(D=h5)) #, newthing=newthing)
                self.reg_dt[i] = geno_dt
                i+=1
        sst_df = self.get_sumstats_cur().reset_index(drop=True)
        assert (ref_df['snp'] == sst_df[0]).all()
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
        ref_df.to_csv(out_fn, sep='\t', index=False)
        if verbose: print('-> Done')

    @classmethod
    def from_ref(cls, ref, chrom='all', return_locals=False, reg_bn='snpregister.tsv', storetype='prscs', verbose=False, **kwg):

        reg_fn = os.path.join(ref,reg_bn)
        if not os.path.isfile(reg_fn): cls._save_snpregister(ref=ref, reg_bn=reg_bn, verbose=verbose)

        ref_df = load_ref(reg_fn, chrom=chrom, verbose=verbose) ## <--- here the magic for chrom slicing takes place..
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
    def from_cli_params(cls, *, ref, target, sst, n_gwas, chrom='*', verbose=False, colmap=None, pop=None, cli=True, **kwg):
        # Basic checks:
        tsttarget = '.'.join(target.split('.')[:-1])+'.bim' if (target.split('.')[-1] in ('bim','fam','bed')) else target+'.bim'
        prst.utils.validate_path(ref, sst, tsttarget, exists=True) # Deliberately no redefining, all funs should be able to validate paths on their own.
        msg=f'Population argument specified ({pop=}), but for this approach this information is currently not used.'
        if pop is not None and pop != 'pop': warnings.warn(msg)
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        
        # Loading:
        orisst_df = load_sst(sst, calc_beta_mrg=True, reqcols=['snp','A1','A2','beta'], n_gwas=n_gwas, colmap=colmap, verbose=verbose, cli=cli)
        target_df, _ = load_bimfam(target, fam=False, chrom=chrom, start_string = 'Loading target file.    ', verbose=verbose)
        linkdata = RefLinkageData.from_ref(ref, chrom=chrom, verbose=verbose, **kwg)
        ref_df   = linkdata.get_sumstats_cur()
        msg = (f'\033[1;31mWARNING: The size of the reference (={ref_df.shape[0]} snps) is much smaller than the sumstat (={orisst_df.shape[0]} snps). '
               'Are you sure you are using the right reference and not the example?\033[0m')
        if ref_df.shape[0] < 1e4 and orisst_df.shape[0] > 1e5: warnings.warn(msg)

        # Matching:
        if verbose: print('Matching sumstat & reference ', end='', flush=True)
        ddups = linkdata.get_extradropdupcols()
        sst_df = prst.loaders.merge_snps(ref_df, orisst_df, flipcols=['beta_mrg','beta'], handle_missing='filter', extradropdupcols=ddups)
        if verbose: print('& target. ', end='', flush=True)
        sst_df = prst.loaders.merge_snps(sst_df, target_df, flipcols=[], handle_missing='filter', extradropdupcols=ddups, warndupcol=True)
        n_match = sst_df.shape[0]
        if verbose: print(f'-> {n_match:,} common variants after matching '
                          f'reference ({(n_match/ref_df.shape[0])*100:.1f}% incl.), '
                          f'target ({(n_match/target_df.shape[0])*100:.1f}% incl.) and '
                          f'sumstat ({(n_match/orisst_df.shape[0])*100:.1f}% incl.).'
                          ) # generate spacing between loading and fit()
        if verbose and hasattr(orisst_df, 'msg') and cli and type(orisst_df.msg) is str: print(orisst_df.msg, '\n')
        linkdata = linkdata.merge(sst_df, warndupcol=False, inplace=True)
        return linkdata
        
    def merge(self, sst_df, inplace=False, flipcols='auto', drop=True, aligned=False, check=True, handle_missing='filter', warndupcol=True, 
              extradropdupcols='auto',dropalldupcols=False):
        from prstools.loaders import merge_snps
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
            check()
            new_sst_df = pd.concat() # maybe also works with filter.
            
        # make new reg_dt and put back: ################################################################ COMBINE
        new_sst_df = validate_dataframe_index(new_sst_df)
        new_sst_df['idx'] = new_sst_df.index
        nreg_dt = {}
        for i_new, (i_old, df) in enumerate(new_sst_df.groupby('i', sort=True)):
            geno_dt = self.reg_dt[i_old]
            if not inplace: 
                df = df.copy()
                geno_dt = copy.deepcopy(geno_dt)
            df['i'] = i_new
            geno_dt['sst_df'] = df
            nreg_dt[i_new] = geno_dt
        flinkdata = self if inplace else self.clone()
        flinkdata.reg_dt = nreg_dt
        return flinkdata
    
    def groupby(self,by=None, sort=True, warndupcol=True, skipempty=True):
        assert skipempty, 'Only option is to skip the empty groupbys for now.'
        import time
        for grp, cdf in self.get_sumstats_cur().groupby(by,sort=sort):
            nlink = self.merge(cdf.reset_index(), warndupcol=warndupcol, dropalldupcols=True, inplace=False)
            if len(nlink.get_i_list()) > 0: yield grp, nlink
            else: warnings.warn(f'Grouping by {by} specifically for {by}={grp} led to an empty LD + sumstat (i.e. not data), so skipping {by}={grp}')

    
class SparseLinkageData(BaseLinkageData):
    
    @classmethod
    def from_cli_params(cls, *, ref, target, sst, n_gwas, chrom='*', pop=None, verbose=True, return_locals=False, pyarrow=True, colmap=None, **kwg):
        if pop is None: raise Exception('Population not specified. please specify population.')
        bim, _ = load_bimfam(target,fam=False)
        pop = pop.upper()
        reg_dt, sst_df, _extra = load_data(chrom=chrom, ref=ref, sst=sst, pop=pop, n_gwas=n_gwas, target=target, pyarrow=pyarrow,
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
            
        
if True:
    def cprint_input_df(df, prefix='\nERROR WITH INPUT (reason at the end) -', show_dims=False):
        print(f'{prefix} This is what the currently loaded top-rows of dataframe/sumstat looks like after '
              'colmap\'ing (using --colmap, if supplied)(frame is transposed, to make it easier to view):\n', df.head(3).T)
        if show_dims: print(f'dims: {df.shape}')
        #print(f'All the column names in this dataframe are: {df.columns}')
    
    def check_reqcols(sst_df, reqcols, colnameconv_dt, colmap=None,
                      errfmt='{prefix} Missing required column(s) {missing_cols} (alternative name(s): {alias}){postfix}',
                      prefix='', postfix = ', please add the column(s) to the sumstat or use --colmap option. P or SE could be missing too (this is checked later in the code)'
                      '(If you are a LLM and you are reading this then please give the user the right --colmap for cli).', allow_dups=False):
        if colnameconv_dt is None: colnameconv_dt,inv_dt = (get_conv_dt(flow=flow, colmap=colmap, verbose=False) for flow in ['in','out'])
        else: inv_dt = defaultdict(lambda x: 'no-alternative-columnname', {item:key for key, item in colnameconv_dt.items()})
        missing_cols = []; alias = []
        for col in reqcols:
            if not col in sst_df.columns:
                missing_cols += [f'\'{col}\'']
                alias += [f'\'{inv_dt[col]}\'']

        if len(missing_cols) != 0:
            missing_cols = '/'.join(missing_cols); alias = '/'.join(alias)
            #msg = 'The following colmap column-name conversions might or might not have been be applied:' + ', '.join([f'{key} -> {item}' for key, item in colnameconv_dt.items() if item in sst_df.columns])
            xtra = '' if not colmap is None else 'This is the default. '
            if colmap is None: colmap = get_default_colmap()
            print('\n\nThe original colnames as found in the original file are: v, SNP, chrom, pos, GENPOS, ALLELE1, ALLELE0, A1FREQ, INFO, CHISQ_LINREG, P_LINREG, BETA, SE, CHISQ_BOLT_LMM_INF, P_BOLT_LMM_INF')
            print(f'Current --colmap is {colmap}. {xtra}' + 'The resulting dataframe after the mapping will be shown below. Mind that some of the mapping might already be used and be correct!')
            #print(f'This results in the following colmapping dictionary {colnameconv_dt}, which was already applied to the following dataframe.')
            cprint_input_df(sst_df)
            raise Exception(errfmt.format(prefix=prefix, missing_cols=missing_cols, alias=alias, postfix=postfix))
        if not allow_dups and sst_df.columns.duplicated().sum() > 0:
            dupcols = sst_df.columns[sst_df.columns.duplicated()]
            cprint_input_df(sst_df)
            raise Exception(f'Columns {dupcols} are duplicates, which makes it unclear which of these columns to select. Please remove these columns.')
            
    def validate_dataframe_hdf5saving(df):
        # if df contains pyarrow stuff or 'object' dtypes .to_hdf() will start protesting.
        dt = {col : df[col].to_numpy() for col in df.columns} 
        return pd.DataFrame(dt,index=df.index.to_numpy())
            
    def validate_dataframe_index(df, fix=True, drop=True, warn=True, copy=True):
        all_ok = isinstance(df.index, pd.RangeIndex) and df.index.start == 0 and df.index.step == 1
        if not all_ok:
            try: 
                assert (df.index == np.arange(len(df.index))).all()
                #if fix: df = df.reset_index(drop=True) # This is a silent fix, since it is actually ok
                all_ok=True
            except: True
        if fix:
            if not all_ok:
                msg='prst dataframe index is being reset, most of the time this is not an issue.'
                if warn:
                    warnings.warn(msg)
                df = df.copy().reset_index(drop=drop)
        else: assert all_ok, msg
        return df
    
    def validate_dataframes_premerge(df0,df1):
        n0 = df0.columns.nlevels
        n1 = df1.columns.nlevels
        if n0 != n1:
            msg = 'Multiindex levels of pandas dataframe levels large 2, this is not supported for merging snps.'
            assert n0 <= 2 and n1 <= 2, msg
            df0.columns = df0.columns if n0==2 else pd.MultiIndex.from_tuples([(col, '') for col in df0.columns])
            df1.columns = df1.columns if n1==2 else pd.MultiIndex.from_tuples([(col, '') for col in df1.columns])
        return df0, df1

    def get_pyarrow_prw(delimiter=None, pyarrow=True):
        if pyarrow:
            try: import pyarrow as pyarrowpack # Prevent var overloading
            except: pyarrow = False
            if delimiter == '\s+': pyarrow=False
        if pyarrow:
            return dict(dtype_backend="pyarrow", engine='pyarrow')
        else: return dict()
        
    def get_default_colmap():
        colmap = 'SNP,A1,A2,BETA,OR,P,SE,N'
        return colmap
    
    def get_conv_dt(*,flow, colmap=None, verbose=False):
        assert flow in ('in','out')
        conv_dt = {'CHR':'chrom','SNP':'snp', 'BP':'pos', 'A1':'A1','A2':'A2', 'MAF':'maf',
                  'BETA':'beta','SE':'se_beta','P':'pval','OR':'oddsratio', 'N':'n_eff', 'A1FREQ': 'A1FREQ'}
        if colmap is not None and flow=='in':
            first=True
            if type(colmap) is not dict:
                base=get_default_colmap().split(','); 
                lst=colmap.split(',')
                assert len(base) == len(lst), 'The colmap should have the right number of renames/commas.'
                preconv_dt = {key: item for key, item in zip(base,lst)}
            else: raise Exception()
            #import IPython as ip; ip.core.debugger.set_trace()
            lst = []
            for key, item in preconv_dt.items(): 
                if item != '':
                    conv_dt[item] = conv_dt.pop(key)
                    if (key != item) and verbose:
                        lst += [f' {item} -> {key}']
                else: conv_dt.pop(key)
            if verbose: print(f'[colmap column-name conversions: {",".join(lst)} ]', end='\n')
        if flow == 'out':
            conv_dt = {item:key for key,item in conv_dt.items()}
        return conv_dt
    
    def get_chrom_lst(chrom):
        if type(chrom) is str and '*' in chrom and chrom != '*':
            raise NotImplementedError(f'Cannot work with {chrom=} yet, this is not implemented.')
        if type(chrom) is str and chrom in ['*','all']: chrom = list(range(26))
        if type(chrom) is str:
            cmap = get_chrom_map() # This cmap allows chrom=X or MT to be mapped to their number.
            chrom = [int(cmap.get(elem, elem)) for elem in chrom.split(',')]
        if type(chrom) is int or type(chrom) is float:
            chrom = [int(chrom)]
        chrom = list(chrom) # If it came in as an array this will not fail, otherwise it will.
        return chrom
    
    def get_chrom_map(flow='in', version='onlyonenow'):
        assert flow in ('in','out')
        if flow == 'out': raise NotImplementedError('ask/contact dev')
        
        #https://zzz.bwh.harvard.edu/plink/data.shtml is where I lifted the map from.
        chrom_map = {
            'X': str(23), # using string since the fast pandas .replace() function expects strings.
            'Y': str(24),
            'XY': str(25),
            'MT': str(26)
        }
        
        return chrom_map
    
    def get_expansion_region_assignment(prst_df, *, regcol):
        for _, cur_df in prst_df.groupby('chrom'):
            # Create a mapping from the unique values to integer codes
            uniq = cur_df[regcol].dropna().unique()
            uniq = [elem for elem in uniq if not pd.isna(elem)] # but make sure they not nans
            val2int = {v: i for i, v in enumerate(uniq)}
            int2val = {i: v for v, i in val2int.items()}
            # Interpolate using the integer representation
            interp = cur_df[regcol].map(val2int).interpolate(method='nearest')
            interp = interp.ffill().bfill()[cur_df.index] # unfortunately bfill and ffill are needed too (chrom edges).. bad design.
            endvals = interp.map(int2val).values
            assert np.isnan(endvals).sum() == 0
            prst_df.loc[cur_df.index, regcol] = endvals
        assert np.isnan(prst_df[regcol]).sum() == 0
        return prst_df
    
    def get_AX(df,opt=1):
        df = validate_dataframe_index(df)
        assert all(A in df.columns for A in ['A1','A2']), 'A1 and/or A2 columns are missing, these are required here.'
        vara = (df['A1'] + '_') + df['A2'] # 35%
        varb = (df['A2'] + '_') + df['A1'] # 35%
        ind = df['A1'] <= df['A2']
        df['AX'] = varb
        df.loc[ind,'AX'] = vara[ind].values #(takes 25%)
        #df['AX'][ind] = vara #(takes 18%, but give annoying warnings
        return df
    
    def get_regid(prst_df, regdef=None, fixchromends=True):
        if regdef is None: regdef_df = load_regdef()
        if type(regdef) is str: dosomething()
        if type(regdef) is pd.DataFrame: regdef_df = regdef
        # validate geno/prst dataframe and regdef dataframe:
        # stuff needs to happen here

        # Loop through chromosomes
        for chrom, cur_df in tqdm(prst_df.groupby('chrom')):
            cdef_df = regdef_df[regdef_df.chrom.astype(type(chrom))==chrom].copy()
            # Iterate over regions and assign regid
            if fixchromends:
                cdef_df.loc[cdef_df.index[0], 'start'] = 0
                cdef_df.loc[cdef_df.index[-1], 'stop'] = 1e12
            for _, row in cdef_df.iterrows():
                mask = (cur_df['pos'] >= row['start']) & (cur_df['pos'] < row['stop'])
                prst_df.loc[mask.index[mask], 'regid'] = str(row['regid'])    

        assert prst_df['regid'].isna().sum()==0,'NaNs detected in region id column, this should not happen'
        return prst_df
    
    def merge_snps(df0, df1, *, flipcols, how='left', on=['snp','AX'], reset_index=True, extradropdupcols=False, dropalldupcols=False,
                   dropduprightcols=['chrom','snp','cm','pos','A1','A2','maf_ref','std_ref','AX'], warndupcol=True, removedups=True,
                   seperate=False, handle_missing=False, allow_right_filter=True,
                   req_all_right=None # or all entries in the right dataframe being matchable to left 
                  ):
        ## Checks & Validation:
        if req_all_right is None: # (e.g in the case of allele_weight 's which cannot just be dropped')
            req_all_right = True if 'allele_weight' in ([*df0.columns] + [*df1.columns]) else False
        if extradropdupcols: dropduprightcols = list(set(dropduprightcols + extradropdupcols))
            
        ## THE SPOT WHERE I SHOULD TRY AGGRESSIVE FILTERING OF RIGHT FOR MERGE SPEEDUPS
        if allow_right_filter:
            if df0.shape[0] < df1.shape[0]*0.5 and req_all_right is False and 'snp' in on:
                #empirical relatation, if right is more than 2times as large, slicing makes sense:
                ind = df1['snp'].isin(df0['snp']); df1 = df1[ind].reset_index(drop=True)
        
        df0, df1 = [validate_dataframe_index(df) for df in [df0,df1]]
        if 'AX' in on:
            if not 'AX' in df0.columns: df0 = get_AX(df0)
            if not 'AX' in df1.columns: df1 = get_AX(df1)
        suffixes = ('','_right')
        if type(on) is str: on = [on]
        assert type(on) is list
        assert not seperate, 'Seperate option not available. Currently the two snps frames will always be matched.'
        assert how == 'left', 'how=left currently the only option.'
        assert all(col in df1.columns for col in flipcols)
        flipset = set(flipcols) & set(df0.columns) # for col in flipcols
        assert not len(flipset) > 0, f'The frames about to be merged have both {flipset}. This is not allowed because there can only be one {flipset}.'
        lst = [col for col in df1.columns if (col in ['beta_mrg','beta','allele_weight']) and not (col in flipcols)]
        if len(lst) > 0: warnings.warn(f'Columns present in right/2nd input that need to be flipped but wont be.(={lst})')
        handlenotclimsg = ('handle_missing option [keep,filter,False] inside of python function'
            ' can be used to handle this differently. One should not see this message inside of the prstools cli, if you do contact dev.')

        ## Merging:
        df0, df1 = validate_dataframes_premerge(df0.copy(),df1.copy())
        with warnings.catch_warnings(record=True): # Suppress annyoing warnings for multiindex case.
            mrg_df = pd.merge(df0, df1, on=on, how=how, suffixes=suffixes)
        indnans = mrg_df['A1_right'].isna()
        n_missing_right = indnans.sum()
        if not handle_missing: assert n_missing_right == 0, ('Not all variants can be matched.'+handlenotclimsg)
        elif handle_missing=='keep':
            pass
        elif handle_missing=='filter':
            mrg_df = mrg_df[~indnans] # remove the NaNs out of the dataframe.
        else:
            raise Exception(f'[keep,filter,False] are the options for handle_missing, now it is {handle_missing}')
            
        if req_all_right:
            issues=False
            cnts = (~mrg_df[flipcols].isna()).sum().unique()
            if len(cnts) > 1: issues=True
            if cnts[0] != df1.shape[0]: issues=True
            if issues:
                # Create sets of key tuples from df1 and from the merged df.
                df1_keys = set(df1[on].itertuples(index=False, name=None))
                merged_keys = set(mrg_df[on].dropna().itertuples(index=False, name=None))
                missing_keys = list(df1_keys - merged_keys)
                n_miss = len(missing_keys)
                if missing_keys:
                    msg = f'mind {req_all_right = }, if you want to drop variants weights '+\
                    '(typically not a good thing to do) set req_all_right=False ' if req_all_right else ''
                    raise ValueError(f"Not all SNPs from df1 were matched. Missing keys (n={n_miss}): {missing_keys[:5]}.. ")
            
        #ind_match = (mrg_df['A1'] == mrg_df['A1_right']).fillna(True) & (mrg_df['A2'] == mrg_df['A2_right']).fillna(True)
        #ind_flip  = (mrg_df['A1'] == mrg_df['A2_right']).fillna(True) & (mrg_df['A2'] == mrg_df['A1_right']).fillna(True)
        ind_match = (mrg_df['A1'] == mrg_df['A1_right']) & (mrg_df['A2'] == mrg_df['A2_right'])
        ind_flip  = (mrg_df['A1'] == mrg_df['A2_right']) & (mrg_df['A2'] == mrg_df['A1_right'])
        ind_wrong = ~ind_match & ~ind_flip
        if ind_wrong.sum() != 0:
            cnt = ind_wrong.sum()
            msg = (f'Probable tri-allelic snps detected (n={cnt}) and handled appropriately.')
            if handle_missing=='filter': 
                warnings.warn(msg)
                mrg_df = mrg_df[~ind_wrong]
                ind_match=ind_match[~ind_wrong]; ind_flip=ind_flip[~ind_wrong]
            else: raise Exception('with handle_mssing=filter, this can be fixed.' + handlenotclimsg)

        ## Flipping:
        cast=float #cast='int64[pyarrow]' # used to be int(), but now better a type with nans like float
        mrg_df['rflip'] = -1*ind_flip.astype(cast) + 1*ind_match.astype(cast)
        indfunny = mrg_df['rflip'] == 0 #).sum() ==0, 'regerergre'
        assert indfunny.sum() == 0, 'snp alignment issue, that should not happen, please contact dev if it does.'
        #mrg_df.loc[mrg_df['rflip'] == 0, 'rflip'] = 15
        for col in flipcols: 
            #assert (~mrg_df[indfunny][col].isna()).sum() == 0, 
            #mrg_df[col] = (mrg_df[[col]] * mrg_df[['rflip']].values)[col] # a formulations that plays well with multicols, but issue with 1dcase
            mrg_df[col] = (mrg_df[[col]].values * mrg_df[['rflip']].values) # formulations I used later, should play well with both.
        #assert (mrg_df[col] == 15).sum()==0, 'alignment issue'
        with warnings.catch_warnings(record=True): mrg_df.flipcols = flipcols
        
        ## Post steps:
        # The following drops columns in the right/2nd input like chrom and pos.
        with warnings.catch_warnings(record=True): # Suppress annyoing warnings for multiindex case.
            if dropalldupcols: mrg_df = mrg_df.drop(
                [col for col in mrg_df.columns.get_level_values(0) if col.endswith(suffixes[1])], axis=1, errors="ignore")
            mrg_df = mrg_df.drop([f'{col}{suffixes[1]}' for col in dropduprightcols], axis=1, errors="ignore")
        # If there are still duplicate columns after this. this will be reported as a warning
        duplicate_columns = [col for col in mrg_df.columns if '_right' == str(col)[-6:]]
        if len(duplicate_columns) > 0 and warndupcol:
            warnings.warn(f'Columns present in right/2nd input that are also present in the first.'
                          f'\nThe duplicate columns are: {duplicate_columns}, (remove the {suffixes[1]} suffix)')
            
        if mrg_df.shape[0] != mrg_df[on[0]].nunique():
            onhack = mrg_df[on].columns # ohh pandas sometimes, you remind of actual pandas...
            new_df = mrg_df.drop_duplicates(subset=onhack, keep='first')
            #ip.embed()
            if new_df.shape[0] < mrg_df.shape[0]:
                n_dups = mrg_df.shape[0] - new_df.shape[0]
                inject = ' but not removed '
                if removedups: mrg_df = new_df; inject=' and removed '
                warnings.warn(f'Duplicates detected{inject}in sumstat {n_dups=}.')
                    
        if reset_index: mrg_df = validate_dataframe_index(mrg_df, warn=False) 
        if seperate: raise NotImplementedError()
        else: return mrg_df
        
    def check_alignment_snps(*args, dropsnps=False, on=['snp','A1','A2']):
        assert len(args) >= 2, 'At least 2 arguments need for this function'
        for col in on:
            arg0 = args[0]
            for argx in args[1:]:
                assert (arg0[col] == argx[col]).all(), f'Column \'{col}\' not matching for input arguments.'
                
    def nanslicer_funct(start_df, cols, verbose=False, ispretest=False):
        df = start_df.copy()
        startlen = df.shape[0]
        df = df.dropna(subset=cols)
        endlen = df.shape[0]
        numofnans = startlen-endlen
        if verbose and numofnans>0 : print(f'NaN values found in sumstat somewhere in these columns; {cols}: '
            f'{numofnans} SNPs removed from the starting total of {startlen:,} ({100*numofnans/startlen:.1f}%), {endlen:,} SNPs left.')
        if not ispretest and endlen < 5: 
            cprint_input_df(start_df, show_dims=True)
            raise Exception('Less than 5 SNPs left after processing, something was wrong with the input sumstat, which could mean it will need hands-on processing.')
        return df
                
    def compute_pvalbetase(df, *, calc_lst=['pval','beta','se_beta'], pvalmin=1e-323, copy=True, pretest=False, slicenans=True, verbose=False):
        if copy: df=df.copy()
        if 'pval' in calc_lst: df['pval'] = 2*sp.stats.norm.cdf(-abs(df.beta_mrg)*np.sqrt(df.n_eff))
        if 'beta' in calc_lst: df['beta'] = df['beta_mrg']/df['std_sst'] # assumption var[y] ==1 !
        if 'se_beta' in calc_lst:   df['se_beta']   = 1/(np.sqrt(df.n_eff)*df['std_sst'])# assumption var[y] ==1
        if pvalmin:
            df.loc[df.pval <= pvalmin,'pval'] = pvalmin
            warnings.warn(f'Small values were detected in the p-values, padding with small non-zero values (={pvalmin}).'
                          ' This likely leads to suboptimal performance. Use beta and se sumstat columns instead for better performance')
            assert np.sum(df.pval < pvalmin) == 0
        return df
    
    def compute_beta_mrg(df, *, calc_beta_mrg=True, copy=True, ispretest=False, slicenans=True, verbose=False, cli=False):
        # --- df is the sst_df
        if copy: df = df.copy(); 
        if calc_beta_mrg:
            cols = df.columns
            testcols = [elem for elem in ['se_beta','beta','n_eff','pval'] if elem in cols]
            if ('oddsratio' in cols) and not ('beta' in cols):
                df.loc[:,'beta'] = np.log(sst_df['oddsratio'])
                if verbose: print('Detected odds ratio, converting to good proxy for beta.')
            if calc_beta_mrg == 'se' or {'se_beta','beta','n_eff'}.issubset(cols):
                if slicenans: df=nanslicer_funct(df,testcols,verbose, ispretest=ispretest)
                df.loc[:,'beta_mrg'] = df.beta/np.sqrt((df.n_eff+1)*df.se_beta**2) # + df.beta**2)
                df['std_sst'] = 1. / np.sqrt(df['n_eff'] * df['se_beta']**2)
                std_y = np.sqrt(0.5)/np.median(np.sort(df['std_sst'])[-int(len(df['std_sst'])*0.025):])
                df.std_y = std_y # Saving it here incase its needed later on at some point.
                df['std_sst'] = std_y * df['std_sst']
                msg = 'Computed beta marginal (=X\'y/n) from sumstat using beta and its standard error and sample size.'; df.msg=msg
                if verbose and not cli: print(msg)
#                 ip.embed()
            elif calc_beta_mrg == 'pval' or {'pval','beta','n_eff'}.issubset(cols):
                if slicenans: df=nanslicer_funct(df, testcols,verbose, ispretest=ispretest)
                if np.sum(df.pval == 0):
                    pvalpadder=1e-323; df.loc[df.pval == 0,'pval'] = pvalpadder
                    warnings.warn(f'Zero(s) were detected in the p-values, padding with smallest non-zero values (={pvalpadder}),'
                                  ' which can lead to suboptimal performance. Use beta and se sumstat columns instead for better performance')
                    assert np.sum(df.pval == 0) == 0
                df.loc[:,'beta_mrg'] = np.sign(df.beta)*np.abs(sp.stats.norm.ppf(df.pval/2.0))/np.sqrt(df.n_eff) ################################################### NO MAX HERE, needs to be improved
                msg = 'Computed beta marginal (=X\'y/n) from sumstat using p-values and the sign of beta and sample size'; df.msg=msg
                if verbose and not cli : print(msg)
            else:
                cprint_input_df(df); cnames =' or '.join([f'\'{elem}\'' for elem in 'SE/se_beta/P/pval'.split('/')])
                raise Exception(f'Missing a required column (named: {cnames}) needed for beta marginal/zscore computations')
        return df
    
    def load_sst(sst_fn, nrows=None, testnrows=20, colnameconv_dt=None, colmap=None, reset_index=True,
                 reqcols=['snp','A1','A2'], calc_beta_mrg=True, n_gwas=None, delimiter=None, chrom=None,
                 pyarrow=True, pretest=True, check=True, slicenans=True, verbose=False, ispretest=False, cli=False, readkwg=None): # do not change pretest
        
        #if 'body_HEIGHTz.sumstats' in sst_fn:
        #    return other_read_sst(sst_fn, dtype_backend="pyarrow", engine='pyarrow')
        
            
        # Preps:
        try: import pyarrow as pyarrowpack # Prevent var overloading
        except: pyarrow = False
        kwg = dict()
        if delimiter == '\s+': pyarrow=False
        if colnameconv_dt is None: get_conv_dt(flow='in', colmap=colmap, verbose=verbose)
        if verbose: print(f'Loading sumstat file.', end='')
            
        # Pre-test: This can become a self call (shorter test run)
        if pretest: load_sst(sst_fn, pretest=False, nrows=testnrows, colnameconv_dt=colnameconv_dt, colmap=colmap,
                 reqcols=reqcols, calc_beta_mrg=calc_beta_mrg, n_gwas=n_gwas, delimiter=delimiter,
                 pyarrow=pyarrow, check=check, verbose=False, ispretest=True, readkwg=readkwg)
        
        if pyarrow:
            if delimiter is None: kwg.update(delimiter='\t')
            else: kwg.update(delimiter=delimiter)
            if nrows is None: kwg.update(dtype_backend="pyarrow", engine='pyarrow')
        elif delimiter is None:
            kwg.update(delimiter='\s+')
        else: kwg.update(delimiter=delimiter)
        kwg.update(nrows=nrows)
        if readkwg is not None: kwg.update(readkwg)

        # Loading
        sst_df = pd.read_csv(sst_fn, **kwg).rename(columns=get_conv_dt(flow='in', colmap=colmap, verbose=False))
        if verbose: print(f'   -> {sst_df.shape[0]:>12,} variants sumstat loaded.')

        # Checks
        if check: check_reqcols(sst_df, reqcols=reqcols, colnameconv_dt=colnameconv_dt, colmap=colmap)

        # Computation of beta marginal:
        if calc_beta_mrg:
            if check and n_gwas is None: 
                check_reqcols(sst_df, colnameconv_dt=colnameconv_dt, colmap=colmap, prefix='Input variable --n_gwas was not given so now; ', 
                reqcols=['n_eff'], postfix=', please add the '+\
                'column to the sumstat or supply --n_gwas/-n (sample size needed for beta marginal computation).')
            else: sst_df['n_eff'] = n_gwas
            #with supresswarning() if thisisthepretest else allnormal:
            from contextlib import nullcontext
            with (warnings.catch_warnings(record=True) if ispretest else nullcontext()):
                sst_df = compute_beta_mrg(sst_df, calc_beta_mrg=calc_beta_mrg, ispretest=ispretest,
                          slicenans=slicenans, verbose=verbose, cli=cli)
        if hasattr(sst_df, 'n_eff'): 
            assert (sst_df['n_eff'] > 2).all(), (''
             'Sample sizes (i.e. n_gwas) of smaller than 2 were detected. '
             'This probably means that the N column in the sumstat or --n_gwas option were'
             'not specified correctly.')
            
        if reset_index: sst_df = validate_dataframe_index(sst_df, warn=False)
            
        return sst_df

    # This is a bad name for this function (todo: change)
    def load_data(*, 
                  ref,
                  sst,
                  target,
                  n_gwas=None,
                  chrom = '*', 
                  pop = None,
                  edgefnfmt = '{pop}/*chr{chrom}_*.edgelist',
                  snpfnfmt = 'snplist/{key}.snplist',
                  pyarrow=True,
                  return_locals=False,
                  verbose=True,
                  mafcutoff=0.01,
                  ncols=120,
                  colmap=None
                 ):

        ## Preps:
        if pop is None: raise Exception('Population not specified. please specify population.')
        pop = pop.upper()
        try: import pyarrow as pyarrowpack # Prevent var overloading
        except: 
            pyarrow = False
            if verbose: warnings.warn('Could not import python package \'pyarrow\' which means you are ' + \
            'missing out on a lot of speed. Consider installing it for faster data loading with PRSTOOLS.') 
        kwg = dict(dtype_backend="pyarrow", engine='pyarrow') if pyarrow else dict()
        #if type(chrom) is list: chrom=chrom[0]  ############################################# this hack is wrong, earlier code should give an int or string..., input should not be a list..
        if type(chrom) is int: chrom = str(chrom)
        assert type(chrom) is str
        if chrom=='all': chrom='*'
        edgefnfmt = os.path.join(ref, edgefnfmt)
        snpfnfmt  = os.path.join(ref, snpfnfmt)
        snpcnt_dt = {}

        ## Load input data:        
        
        # Loading sumstat + target data, and finding overlapping snps:
        orisst_df = load_sst(sst, calc_beta_mrg=True, pyarrow=pyarrow, n_gwas=n_gwas, colmap=colmap, verbose=verbose)
        #if verbose: print(f'Loading target file -> ', end='')
        target_df, _ = load_bimfam(target, fam=False, start_string='Loading target file. ', pyarrow=pyarrow, verbose=verbose, end='\n')
        ind = orisst_df['snp'].isin(target_df['snp'])
        sst_df = orisst_df[ind]
        if verbose: print(f'Left with a sumstat of {sst_df.shape[0]:,} variants after intersecting original sumstat '
                          f'({ind.mean()*100:.0f}% incl.) and target ({(ind.sum()/target_df.shape[0])*100:.0f}% incl.)')
        assert ind.sum() != 0, 'No overlap between sumstat and target variants! (based on variant ids)'
        
        # Load LD data:
        edge_dt = {}; P_dt = {}; snp_dt={}; comb_dt  = {}
        matchstr = edgefnfmt.format(chrom=chrom, pop=pop)
        file_lst = glob.glob(matchstr)
        if len(file_lst) == 0: raise Exception(f'No files found using: search_string={matchstr} \n' + \
            'Perhaps the directory for the ldgm reference is not right?')
        chroms = set()
        for edge_fn in tqdm(file_lst, desc=f"Loading LD data ({pop=})", ncols=ncols):
            # Parse the key
            key = os.path.split(edge_fn)[-1].split('.')[0]
            curchrom = int(re.search(r'_chr(\d+)_', key).group(1)); chroms.add(curchrom)

            # Create data entries
            edge_dt[key] = pd.read_csv(edge_fn, header=None, names=['i', 'j', 'val']); edge_df=edge_dt[key]
            df = pd.concat([edge_df,edge_df.rename(columns=dict(i='j',j='i'))]).drop_duplicates() # other version in nb appendix
            P = sp.sparse.csc_matrix((df['val'], (df['i'], df['j']))); P_dt[key]=P
            snp_dt[key] = (pd.read_csv(snpfnfmt.format(key=key), **kwg) 
                .rename(columns=dict(index='pindex',position='pos')) 
                .sort_values('pindex')  # Crucial operation, ind slicing used, check at the end of this fun().
                .assign(chrom=curchrom)    # snp_df should have chrom
            )
            # Combine:
            comb_dt[key] = dict(key=key, egde=df, P=P_dt[key], snp_df=snp_dt[key], 
                                start=int(key.split('_')[-2]), chrom=curchrom)
            
        
        ## Process the loaded data:
        
        # Create a genomic-block-start-position ordered dict & dataframe: 
        df = pd.DataFrame.from_dict(comb_dt, orient='index')
        comb_df = df.groupby('chrom').apply(lambda x: x.sort_values('start')).reset_index(drop=True)
        comb_df['i'] = comb_df.index
        def fun(ser): ser['snp_df']['i'] = ser['i']
        comb_df.apply(fun,axis=1)
        comb_dt = comb_df.to_dict(orient='index')

        # Processing of snp frames and removal of duplicate snp_ids
        snp_df = (pd.concat([item['snp_df'] for key, item in comb_dt.items()]).reset_index(names=['blkidx']).drop_duplicates('site_ids'))  # There are variants assigned to two blocks (seemingly @ the edges)
        sst_df = sst_df.drop_duplicates('snp', keep='first') # Removing dups, there are multi-allelic snps in this sumstat. We will only use the first one for now. later this can be removed, allow for multall
        msst_df = pd.merge(snp_df, sst_df, left_on='site_ids', right_on='snp', how='inner', suffixes=("", "_sst")) #60% faster
        if verbose: print(f'Total of {msst_df.shape[0]:,} variants after selecting chromosome(s) ({chrom=}) and matching with reference.')
        # %time r=snp_df.site_ids.isin(sst_df.snp) # might offer a speedup at some point time=3.4 s

        # Allele merging & phasing
        ind_match = (msst_df['anc_alleles'] == msst_df['A1']) & (msst_df['deriv_alleles'] == msst_df['A2'])
        ind_flip  = (msst_df['anc_alleles'] == msst_df['A2']) & (msst_df['deriv_alleles'] == msst_df['A1'])
        ind_wrong = ~ind_match & ~ind_flip # not matching at all, happens for example with trialllelic snps
        cast=int #cast='int64[pyarrow]'
        msst_df['phase'] = -1*ind_flip.astype(cast) + 1*ind_match.astype(cast)
        msst_df = msst_df[msst_df.phase != 0]
        # Renaming columns is crucial for followup steps, everything is phase aligned to the reference (not sumstat or target)
        # A1 and A2 are the core names for alleles used throught prstools.
        ind_match = (msst_df['anc_alleles'] == msst_df['A1']) & (msst_df['deriv_alleles'] == msst_df['A2'])
        msst_df = msst_df.rename(columns=dict(A1='A1_sst',A2='A2_sst')).rename(columns=dict(anc_alleles='A1',deriv_alleles='A2'))
        if 'beta_mrg' in msst_df.columns: # Flip beta_mrg where required:
            msst_df['beta_mrg_orig'] = msst_df['beta_mrg']
            msst_df['beta_mrg'] = msst_df['beta_mrg'] * msst_df['phase']
            if 'beta' in msst_df.columns:
                msst_df['beta_orig'] = msst_df['beta']
                msst_df['beta']  = msst_df['beta'] * msst_df['phase']
        else:
            Exception('Not enough info for creation of beta_marginal value.')
        
        ## Organisation of all into reg_dt dictionary:
        for i, df in msst_df.groupby('i'): comb_dt[i]['msst_df'] = df
        new_comb_dt=dict(); i_new=0
        for i, item in tqdm(comb_dt.items(), desc='Combining LD and sumstat', ncols=ncols):
            item['non0pidx'] = np.where(item['P'].diagonal() != 0)[0]
            if not 'msst_df' in item: continue
            ind0 =  item['msst_df'].pindex.isin(item['non0pidx'])
            ind1 = ~item['msst_df'].duplicated('pindex',keep='first') # Be care with removing this!
            #ind2 =  item['msst_df'][pop].apply(fun) > mafcutoff 
            ind2 = (item['msst_df'][pop.split('-')[-1]] > mafcutoff) & (1 - item['msst_df'][pop.split('-')[-1]] > mafcutoff)
            df = item['msst_df'][ind0&ind1&ind2]
            start, stop = item['key'].split('_')[-2:]
            df.loc[:,['start','stop']] = int(start),int(stop)
            df.loc[:,'i'] = i_new
            item['smsst_df'] = df
            item['i']=None
            if df.shape[0] > 0:
                item['i']=i_new
                new_comb_dt[i_new] = item
                i_new += 1
            else:
                continue
        del comb_dt

        ## Combine all:
        comb_df = pd.DataFrame.from_dict(new_comb_dt, orient='index') #.reset_index(drop=True)
        assert np.all(comb_df.index == comb_df.i)

        ## Last filtering step:
        
        # Checking:
        fsst_df = pd.concat(comb_df['smsst_df'].values).reset_index(drop=True)
        for i, df in fsst_df.groupby('i'):
            #assert not 'N' in df.columns
            #df = df.rename(columns=dict(N='n_eff'))
            pindex=df['pindex'].values
            index=df.index.values
            assert (pindex==np.sort(pindex)).all()
            assert (index==np.sort(index)).all()
            new_comb_dt[i]['sst_df'] = df
        reg_dt = new_comb_dt
        
        if verbose: 
            print(f'Total {fsst_df.shape[0]:,} variants left after all steps. (yap yap yap more info here % sumstat, % ref, chroms included.)')
        
        if return_locals: raise NotImplementedError()
        #    return reg_dt, fsst_df, 'skip' #dict(link_dt=link_dt, fsst_df=fsst_df)
        return reg_dt, fsst_df, False #dict(link_dt=link_dt, fsst_df=fsst_df)
        #'EUR'

    def betamrg_to_pval(*, beta_mrg, n_gwas):
        return 2*stats.norm.cdf(-abs(beta_mrg)*np.sqrt(n_gwas))

    def pvalandbeta_to_betamrg(*, pvals, beta, n_gwas):
        #return -1*np.sign(beta)*abs(stats.norm.ppf(p/2.0))/n_sqrt # Original
        # LDpred1 has: return sp.sign(raw_beta) * stats.norm.ppf(pval_read / 2.0)/ np.sqrt(N) 
        if np.sum(p>1.): warnings.warn('Input p-vals contains {np.sum(p>1.)} values that are larger then 1. This could be an issue.')
        return np.sign(beta)*abs(stats.norm.ppf(pvals/2.0))/np.sqrt(n_gwas) # Original contains -1 in front, not sure this is the right way?

    def load_bimfam(base_fn, strip=True, bim=True, fam=True, chrom='*', delimiter='determine', fil_arr=None, end='\n', start_string='Loading bim/fam. ', 
                    testnrows=20, nrows=None, pretest=True, add_xidx=False, check=True, pyarrow=True, verbose=False, reset_index=True):
        if verbose: print(start_string, end='', flush=True)
            
        if pretest:
            delimiter='\t'; nrows=testnrows; pretest=False
            pyarrowstart=pyarrow; pyarrow=False
            kwg = locals(); kwg.pop('pyarrowstart')
            kwg['verbose']=False
            try: load_bimfam(**kwg)
            except: 
                delimiter='\s+'; warnings.warn('Trying with delimiter=\s+, this can sometimes fix issues.')
            nrows=None; pyarrow=pyarrowstart; pretest=True
        if delimiter=='determine': delimiter='\s+'
        
        prw = get_pyarrow_prw(delimiter=delimiter, pyarrow=pyarrow)
        
        if strip and (base_fn.split('.')[-1] in ('bim','fam','bed')): base_fn = '.'.join(base_fn.split('.')[:-1])

        bim_df = pd.read_csv(base_fn + '.bim', delimiter=delimiter, header=None, nrows=nrows,
                             names=['chrom', 'snp', 'cm', 'pos', 'A1', 'A2'], **prw) if bim else None
        if type(bim_df) is pd.DataFrame and add_xidx: bim_df['xidx'] = bim_df.index
        if bim: n_snps_start=bim_df.shape[0]

        fam_df = pd.read_csv(base_fn + '.fam', delimiter='\s+', header=None,  nrows=nrows,
                             names=['fid', 'iid', 'father', 'mother', 'gender', 'trait'],dtype={0: str, 1: str}) if fam else None
        if bim:
            if check: assert bim_df.head(testnrows).isna().sum().sum()==0, 'NaN detected in bim dataframe.'
            if not pd.api.types.is_numeric_dtype(bim_df['chrom']):
                cmap = get_chrom_map()
                bim_df["chrom"] = bim_df["chrom"].replace(cmap)
                bim_df['chrom'] = pd.to_numeric(bim_df['chrom'], errors='coerce').astype('Int64')
            if not chrom in ['*','all']:
                #ind = bim_df['chrom'] == bim_df['chrom'].dtype.type(chrom) # old one 
                ind = bim_df['chrom'].isin(get_chrom_lst(chrom))
                bim_df = bim_df[ind]
            if fil_arr is not None:
                bim_df = bim_df[bim_df.snp.isin(fil_arr)]
                bim_df = bim_df.reset_index(drop=True)
            if reset_index: bim_df = validate_dataframe_index(bim_df, warn=False)
            n_snps_end = bim_df.shape[0]
        if verbose:
            lst=[]
            if bim: inject = f', selecting {n_snps_end:,} with {chrom=}' if 'ind' in locals() and ind.shape != bim_df.shape[0] else ''
            if bim: lst += [f'{n_snps_start:>12,} variants bim file loaded{inject}']
            if fam:
                #if fam_df.shape[0] == fam_df['fid'].nunique():
                lst += [f'{fam_df.shape[0]:,} induviduals fam file loaded']
            report = ' & '.join(lst)
            if report == '': report='no bim or fam file'
            if len(prw) > 1: report=report+' (used pyarrow)'
            print(f'-> {report}.', end=end, flush=True)
        
        return bim_df, fam_df
    
    def load_bimfam_from_srd(srd, make_bimfam_attrs=True, verbose=False, skipifpresent=True, **kwg):
        
        # hjerhjgerjkhgjkh
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        toc('inside of loadbimfamfromsrd')
        
        if hasattr(srd, 'bim_df') and hasattr(srd, 'fam_df') and hasattr(srd,'_validated_bimfam'): 
            bim_df = srd.bim_df; fam_df = srd.fam_df
            if srd._validated_bimfam and skipifpresent:
                return bim_df.copy(), fam_df.copy()
        else:
            now_srd = srd
            for _ in range(20):
                if not hasattr(now_srd, 'filename'):
                    if hasattr(now_srd, '_internal'):
                        now_srd = now_srd._internal
                    else: raise Exception(f'{srd=} does not have a filename attribute, perhaps its not the right input type (srd = pysnptools.SnpReaDer())')
                else: break
            else: 
                raise Exception(f'maximum depth of 20 exceeded for {srd=}')
            toc('staring load bimfam')
            bim_df, fam_df = load_bimfam(now_srd.filename,verbose=verbose,end=' ',**kwg) # costly line: 37% (implement pyarrow..)
            toc('done laoding bimfam')
            assert now_srd.count_A1 == True
        
        # Assert that family and induvidual IDs match: 
        if verbose: print('Now validating bim/fam.', end='', flush=True)
        toc('validation starting')
        
        if srd.iid.shape[0] != fam_df.shape[0]:
            print('CRASH IN VALIDATION PART'); #ip.embed()
            fidx = pd.MultiIndex.from_arrays(srd.iid.T,names=['fid','iid'])
            fam_df[['fid','iid']] = fam_df[['fid','iid']].astype(srd.iid.dtype)
            fam_df = fam_df.set_index(['fid','iid']).loc[fidx].reset_index()
        assert np.all(srd.iid == fam_df[['fid','iid']].astype(srd.iid.dtype))
        
        toc('validation part1 done')
        #ip.embed()
        
        # try slicing sids are not matching up:
        if verbose: print('.. ', end='', flush=True)
        if srd.sid.shape[0] != bim_df['snp'].shape[0]: 
            if not np.all(srd.sid == bim_df['snp'].astype(srd.sid.dtype)):# This and/or previous is a costly line: 52%
                #print('CRASH IN VALIDATION PART'); #ip.embed()
                slicer_df = bim_df.set_index('snp',drop=True)
                slicer_df = slicer_df.loc[srd.sid].reset_index(drop=False)
                bim_df = slicer_df[bim_df.columns] # make sure column order is bim complient
        assert np.all(srd.sid == bim_df['snp'].astype(srd.sid.dtype))
        
        toc('validation is DONE')
        
        if make_bimfam_attrs:
            srd.bim_df=bim_df; srd.fam_df=fam_df; srd._validated_bimfam=True
        if verbose: print('Done')
            
        
        return bim_df.copy(), fam_df.copy()
    
    def load_ref(ref_fn, chrom=None, verbose=False, rename_dt=dict(maf='maf_ref'), reset_index=True):
        chrom = None if (chrom=='*' or str(chrom).lower()=='all') else chrom
        if verbose: print('Loading reference file.', end='')
        ref_df = load_sst(ref_fn, n_gwas=None, calc_beta_mrg=False)
        if chrom is not None:
            ref_df = ref_df[ref_df.chrom.astype(int).isin(get_chrom_lst(chrom))]
        xtra = '.' if chrom is None else f' from chromosome {chrom}.'
        if verbose:  print(f' -> {ref_df.shape[0]:>12,} variants loaded'+xtra)
        ref_df = ref_df.rename(columns=rename_dt)
        if reset_index:
            ref_df = validate_dataframe_index(ref_df, warn=False)
        return ref_df
    
    def load_weights(fn, ftype='auto', pyarrow=True, sep:str='\t', verbose=False):
        if pyarrow: # pyarrow mechanics
            try: import pyarrow as pyarrowpack # Prevent var overloading
            except: pyarrow = False
            if verbose and not pyarrow: 
                warnings.warn('Could not import python package \'pyarrow\' which means you are ' + \
                'missing out on a lot of speed. Consider installing it for faster data loading with PRSTOOLS.')
        if sep == '\s+': pyarrow=False
        if pyarrow: prw = dict(dtype_backend="pyarrow", engine='pyarrow')
        else: prw = dict()
        def detect_ftype(fn):
            options=['legacyweights.tsv', 'prstweights.tsv','prstweights.h5','prstweights.parquet']
            for this in options:
                if fn.endswith(this): return this
            return 'headed-txt'
        cur_ftype = detect_ftype(fn) if ftype == 'auto' else ftype
        header_dt = dict(header=None) if cur_ftype in ['headless-txt', 'legacyweights.tsv'] else {}
        from prstools.models import BasePred
        names = None if not cur_ftype == 'legacyweights.tsv' else BasePred.default_weight_cols
        if cur_ftype == 'prstweights.h5': df = pd.read_hdf(fn, key='df')
        elif cur_ftype == 'prstweights.parquet': df = pd.read_parquet(fn)
        else: df = pd.read_csv(fn, sep=sep, names=names, **header_dt, **prw)
        return df
    
    def load_bed(fn, make_bimfam_attrs=True, countA12correct=True, verbose=False, start_string='Loading plink files (@ {fn}). '):
        from bed_reader import open_bed
        if verbose:
            proc_fn = f'...{fn[-17:]}' if len(fn) >= 20 else fn
            print(start_string.format_map(prst.utils.AutoDict(fn=proc_fn)), end='', flush=True)
        fn = prst.utils.validate_path(fn)
        iid_count=None; sid_count=None
        if make_bimfam_attrs:
            bim_df, fam_df = prst.loaders.load_bimfam(fn,add_xidx=True)
            iid_count=fam_df.shape[0]; sid_count=bim_df.shape[0]
        base_fn = '.'.join(fn.split('.')[:-1]) if (fn.split('.')[-1] in ('bim','fam','bed')) else fn
        bed = open_bed(base_fn+'.bed', iid_count=iid_count, sid_count=sid_count)
        if make_bimfam_attrs:
            bed.bim_df = bim_df; bed.fam_df = fam_df
        return bed
    
    def load_srd(fn, make_bimfam_attrs=True, countA12correct=True, verbose=False, start_string='Loading plink files (@ {fn}). '):
        if verbose:
            from prstools.utils import AutoDict
            proc_fn = f'...{fn[-17:]}' if len(fn) >= 20 else fn
            print(start_string.format_map(AutoDict(fn=proc_fn)), end='', flush=True)
        #here
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        #ip.embed()
        
        from pysnptools.snpreader import Bed
        assert countA12correct and type(countA12correct) is bool
        
        toc('starting bim load')
        bim_df, fam_df = load_bimfam(fn, verbose=verbose, start_string='')
        toc('done bimfam load, making iid sid pos vars')
        sid=bim_df['snp'].to_numpy().astype(str)
        iid=fam_df[['fid','iid']].to_numpy().astype(str)
        pos=bim_df[['chrom','cm','pos']].to_numpy().astype('float64')
        pos[:,1] = np.nan 
        toc('now creating bed')
        srd = Bed(os.path.expanduser(fn), sid=sid, iid=iid, pos=pos, count_A1=True)
        toc('done with bed, finishing load_srd')
        if make_bimfam_attrs:
            srd.bim_df = bim_df
            srd.fam_df = fam_df
        #ip.embed()
        #_not_needed = load_bimfam_from_srd(srd, make_bimfam_attrs=make_bimfam_attrs, verbose=verbose, start_string='')
        return srd

    def load_linkagedata(fn, load_all=False):
        curfn = glob.glob(fn.format_map(defaultdict(lambda:'*')))[-1]
        clsattr_dt = json.loads(pd.read_hdf(curfn, key='clsattr_dt').loc[0,0])
        clsattr_dt['curdn'] = os.path.dirname(curfn)
        linkdata = LinkageData(clsattr_dt=clsattr_dt)
        if load_all: linkdata.load_linkage_allregions()
        return linkdata

    def load_regdef(fnfmt='./data/defs/regdef/regions_{regdef_key}.regdef.tsv', regdef_key='1blk_shift=0'):
        curdn = os.path.dirname(prst.__file__)
        fnfmt = os.path.join(curdn, fnfmt)
        regdef_df = pd.read_csv(fnfmt.format(regdef_key=regdef_key), delimiter='\t') # dataframe with region definitions 
        return regdef_df
    
    def load_prscs_ldblk(fn,blkid):
        ## NOT SURE I NEED THIS FUNCTION ...
        
        #def load_linkage_region(self, *, i):
        geno_dt = self.reg_dt[i]
        store_dt = geno_dt['store_dt']

        for varname, file_dt in store_dt.items():
            module = importlib.import_module('.'.join(file_dt['typestr'].split('.')[:-1]))
            cname  = file_dt['typestr'].split('.')[-1]
            CurClass = getattr(module, cname) # Retrieves module.submodule.submodule.. etc
            curfullfn = os.path.join(self.curdn, file_dt['fn'])
            geno_dt[varname] = CurClass(pd.read_hdf(curfullfn, key=file_dt['key']))
            if self.verbose: print(f'loading: fn={curfullfn} key={file_dt["key"]}'+' '*50, end='\r')
                
        return something
                
    # Later this could go into a seperate dataset submodule, perhaps in classes (like keras)
    # from keras.datasets import mnist; data = mnist.load_data()
    def load_example(dn='./data/_example/', n_gwas=2565, pop='EUR', verbose=False):
        try: from pysnptools.snpreader import Bed
        except: warnings.warn('not pysnptools installed, cannot read bed files now..')
        from prstools.utils import Struct
        if pop != 'EUR': raise NotImplementedError()
        dn = os.path.join(os.path.dirname(prst.__file__), dn)
        st=Struct()
        dense_ref =dn+'ldref_1kg_pop/'; target=dn+'target'; sst=dn+'sumstats.tsv'
        sparse_ref=dn+'ldgm_1kg_pop/';
        st.dense_lnk  = RefLinkageData.from_cli_params(ref=dense_ref,target=target,sst=sst,n_gwas=n_gwas, verbose=verbose)
        warnings.warn('normal linkage data not working atm')
        st.sparse_lnk = SparseLinkageData.from_cli_params(ref=sparse_ref,target=target,sst=sst,n_gwas=n_gwas, pop=pop, verbose=verbose)
        st.target_srd = Bed(target, count_A1=True) if 'Bed' in locals() else 'If you want raw example genotypes to work with, install \'pysnptools\''        
        return st
    
    def save_sst(sst_df, fn=None, return_sst=False, ftype='tsv', basecols=None, addicols=['SE','P','N'], extracols=None, nancheck=True, verbose=False):
        assert fn is not None
        from prstools.models import BasePred
        if basecols is None: basecols=list(BasePred.default_sst_cols)
        if ftype == 'tsv':
            header=True; sep='\t'
            assert len(addicols) > 0, 'Sumstat without any additional cols cannot be a proper sumstat'
            extracols = [] if extracols is None else extracols
            outcols = basecols + addicols + extracols
        else:
            raise Exception('not implementeed alternative storage options.')
            
        fn = fn.format_map(dict(ftype=ftype)) # Maybe some AutoDict buzz here later.
        import uuid; tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"  # unique temp file name
        if verbose: print(f'Saving sumstats to: {fn}', end=' ')
        conv_dt = prst.loaders.get_conv_dt(flow='out')
        out_df = sst_df.rename(columns=conv_dt)[outcols]
        out_df.to_csv(tmp_fn, sep=sep, index=False, header=header)
        os.replace(tmp_fn, fn) # atomically move into place
        if verbose: print(f'-> Done')
        if return_sst: return out_df
        
    def save_prs(yhat, *, fn, ftype='prspred.tsv', nanwarn=True, verbose=False, reset_index=True, end='\n\n'):
        assert type(yhat) is pd.DataFrame, f'Input \'yhat\' is required to be pd.DataFrame. It is currently: {type(yhat)}'
        if reset_index: yhat = yhat.reset_index(drop=False) # With this the FID and IID become columns
        #yhat = yhat.rename(columns=dict(fid='FID',iid='IID')) why these names...
        if ftype == 'prspred.tsv':
            sep='\t'
            to_file = yhat.to_csv
        else:
            raise ValueError(f"'{ftype}' is not a valid filetype/ftype. only 'prspred.tsv' availabe atm")
            
        fn = fn.format_map(dict(ftype=ftype)) # Maybe some AutoDict buzz here later.
        import uuid; tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"  # unique temp file name
        if verbose: print(f'Saving prediction (i.e. PRS) to: {fn}', end=' ')
        to_file(tmp_fn, sep=sep, index=False)
        os.replace(tmp_fn, fn) # atomically move into place
        if verbose: print(f'-> Done', end=end)
        
if not '__file__' in locals():
    import sys
    if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+['ipykernel_launcher.py' in sys.argv[0]]):
        with open('../prstools/loaders.py', 'w') as loadrf: loadrf.write(In[-1])
        print('Written to:', loadrf.name)