from abc import ABC, abstractmethod
import copy, time, warnings, math, traceback, sys, os, glob
import scipy as sp
import numpy as np
import pandas as pd
from scipy import linalg, stats
import prstools as prst
from prstools.models._compute import dpsi, gigrnd, g

# __all__ = [k for k, v in globals().items() if not k.startswith('_') and v.__module__ == __name__]

# __all__ = [
#     name for name, val in globals().items()
#     if not name.startswith('_')
#     and callable(val)
#     and getattr(val, '__module__', None) == __name__
# ]

__all__ = ['BasePred','MultiPred','GroupByModel','PredPRS','PRSCS2']


class BasePred(ABC):

    weight_filetypes = ['extprst.tsv','prscs.tsv']
    default_weightfile_type = 'prstweights.tsv'
    default_weight_cols = ['chrom','snp','pos','A1','A2','allele_weight']
    extra_weight_cols   = False
    default_sst_cols = ['SNP','A1','A2','BETA']
    _nancheck = True
    scaling = 'ref'
    shuffle = False
    _clear_cache=True
    _close_pbar = True
    _default_n_jobs=4
    
    def remove_cache(self):
        if hasattr(self,'cache_dt'):
            del self.cache_dt
    
    def _checktype(self, obj, classname): # This methods needs some work
        #if not (type(obj).__name__ in list(classnames)): raise TypeError(f'{type(obj)} not allowed as linkdata input. Must be {classnames}')        
        if not classname in type(obj).__name__: raise TypeError(f'{type(obj)} not a {classname} type (or subtypes thereof).')
        
    def _order(self, input_lst):
        assert type(input_lst) is list
        if self.shuffle:
            return input_lst[np.random.permutation(len(input_lst))]
        else:
            return input_lst
        
    def _return_stddev(self, _stansda):
        _stddev = _stansda.stats[:,1][:,np.newaxis]
        if np.any(np.isinf(_stddev)): warnings.warn('inf value detected in standardizer. implies maf=0')
        _stddev[np.isinf(_stddev)] = 1
        return _stddev
    
    @classmethod
    def _get_cli_epilog(cls, commentccode='32'):
        #commentccode='32;2' 
        from textwrap import dedent
        def format_color(text, color_code): return f"\033[{color_code}m{text}\033[0m"
        string = format_color('test',5)
        insert = len(cls.__name__.lower())*' '
        cmdname=cls.__name__.lower() #,string=string, insert=insert
        ldrefname='ldgm_1kg_pop' if 'sparse' in cls.__doc__.lower() else 'ldref_1kg_pop'
        chromopt='--chrom \'*\' ' if 'sparse' in cls.__doc__.lower() else ''
        epilog=f'''\
        Examples --> can be directly copy-pasted (:
         prst downloadutil --pattern example --destdir ./; cd example  {insert}                                        # Makes \'example\' dir in current path.
         prstools {cmdname} --ref {ldrefname} --target target --sst sumstats.tsv {chromopt}--n_gwas 2565 --out ./result-{cmdname} # Run the model with example data.
         prst {cmdname} -r {ldrefname} -t target -s sumstats.tsv -n 2565 {chromopt}-o ./result-{cmdname} --pred                  # A shorter version of previous that also does the predictions.
        '''
        #plink --bfile target --out prspred --keep-allele-order --score ./result-{cmdname}_* 2 4 6 # Make predictions from weights (plink must be installed).
        newepi = []
        for elem in epilog.split('\n'):
            elems = elem.split('#')
            elems[-1] = format_color('#'+elems[-1],commentccode)
            if len(elems)>1:
                newepi.append(''.join(elems))
            else:
                newepi.append(elem)
        epilog = '\n'.join(newepi)
        
        return dedent(epilog)
    
    @classmethod
    def _get_cli_spkwg(cls, basic_pkwargs='figurethisout'):
        
        placeholdersection = thisstuffisnotactiveyet
        if basic_pkwargs and type(basic_pkwargs) is bool:
            basic_pkwargs = dict(basics=dict(
                    args=['-h', '--help'], 
                    kwargs=dict(action='help', help='Show this help message and exit.')))
        from textwrap import dedent
        doc = dedent(cls.__doc__)
        epilog = getattr(cls, '_get_cli_epilog', lambda: None)()
        spkwg = dict(
            cmdname    = cls.__name__.lower(), #command name
            clsname    = cls.__name__, #class name
            description= doc,
            help       = doc.split('\n')[0],
            epilog     = epilog,
            modulename = cls.__module__,
            groups     = dict(
                general=dict(
                    grpheader='Options',
                    pkwargs={**basic_pkwargs, **retrieve_pkwargs(cls.from_cli_params_and_run)})
            ),
            subtype    = 'PRSTCLI'
        )
        
        ##### THIS STUFF IS ACTUALLY NOT BEING USED, BUT THIS IS HOW IT SHOULD BE USED.
        verifythisisnotused
        
        # Add data related kwargs:
        data_group = model_parser.add_argument_group('Data Arguments (first 5 required)')
        data_group.add_argument("--ref","--ref_dir","-r", 
                **prc(dict(required=True, metavar='<dir/refcode>', 
                    help="Path to the directory that contains the LD reference panel. You can download this reference data "
                           f"manually with \'{basecmd} downloadref\'. Soon it will be possible to automatically download it on the fly.")))
        data_group.add_argument("--target", "--bim_prefix", "-t", 
                **prc(dict(required=True, metavar='<bim-prefix>', 
                    help="Specify the directory and prefix of the bim file for the target dataset.")))
        data_group.add_argument("--sst","--sst_file","-s", **prc(dict(required=True, metavar='<file>', 
                    help="The summary statistics file from which the model will be created. The file should contain columns SNP, A1, A2, P and BETA or OR (in this order)."
                         f" SNP should contain rsid\'s. See {format_color('https://tinyurl.com/sstxampl','34')} for an example.")))
        data_group.add_argument("--n_gwas","-n", 
                **prc(dict(required=True, type=int, metavar='<num>', 
                    help="Sample size of the GWAS")))
        data_group.add_argument("--out","--out_dir","-o", 
                **prc(dict(required=True, metavar='<dir>', 
                    help="Output prefix for the results (variant weights). This should be a combination of the desired output dir and file prefix.")))
        data_group.add_argument("--chrom", 
                    type=lambda x: x.split(','), metavar='<chroms>', default=range(1, 23), 
                    help="Optional: Chromosomes to include. You can specify multiple by using a separting comma e.g. \"--chrom 1,2,3\".")

        return spkwg
    
    @classmethod
    def from_params(cls, groupby=False, **kwg):
        #import IPython as ip; ip.embed() 
        if str(groupby) == '-1': groupby=False
        if groupby:
            model = GroupByModel(cls(**kwg), groupby=groupby, verbose=kwg.get('verbose',False),
                    **{key:item for key,item in kwg.items() if not key in ['verbose','groupby']})
        else:
            model = cls(**kwg, groupby=groupby)
        return model
    
    @classmethod
    def from_cli_params_and_run(cls, *, ref, target, sst, n_gwas=None, chrom='all', fnfmt='_.{ftype}', ftype='prstweights.tsv', groupbydefault=False,
                                verbose=True, pkwargs=None, out=None, return_models=True, fit=True, pop=None, colmap=None, pred=False, command=None, **kwargs):
        from prstools.loaders import RefLinkageData
        
        # Initialize model object(s) (multiple since hyperparam ranges, and maybe chroms):
        def testkey(key): return (key in pkwargs) if pkwargs else True # The verbose in the next line overwrites the verbose in the 'kwargs' dict
        groupby = kwargs.get('groupby', pkwargs.get('groupby',{}).get('kwargs',{}).get('default', groupbydefault))
        model = cls.from_params(**dict({key: item for key, item in kwargs.items() if testkey(key)}, verbose=verbose, groupby=groupby))
        
#         # Loop through different models, likely a parameter grid:
#         for model in models: # not sure about this atm
        
        # Gen output file name format and do quick check if output file can be saved before a lot of work is done:
        if out: out_fnfmt = model.create_output_fnfmt(**locals())

        # Initialize data objects, fit the model & predict:
        linkdata = RefLinkageData.from_cli_params(ref=ref, target=target, sst=sst, 
                        n_gwas=n_gwas, chrom=chrom, colmap=colmap, pop=pop, verbose=verbose)
        model.set_linkdata(linkdata)
        if fit: model.fit()
        if out:  model.save_weights(out_fnfmt, ftype=ftype) # Store fitting result
            
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        
        if pred: # Prediction
            #import IPython as ip; ip.embed() # save these lines for later..
            #model.remove_linkdata(); linkdata.clear_linkage_allregions # seems to do pretty much nothing.. anyway xp was 5% mem, which jumped to 20 and 60 later
            bed = prst.loaders.load_bed(target, verbose=verbose)
            yhat = model.predict(bed); 
            prst.loaders.save_prs(yhat, fn=out_fnfmt, verbose=verbose) # Store prediction result
        if return_models: 
            return model
    
    @staticmethod
    def basenaming(item):
        if type(item) is str:
            newitem = os.path.basename(item)
            if newitem == '': newitem = os.path.basename(item.rstrip('/\\'))
            if newitem == '': newitem = item 
        else: newitem = item
        return newitem
    
    @staticmethod
    def create_output_fnfmt(*, cls, out, fnfmt, prstlogs=True, testsave=True, ftype=None, **kwg):
        assert testsave and prstlogs, 'testsave must be enable at this point'
        from prstools.utils import AutoDict
        mname = cls.__name__.lower()
        out_fnfmt = out + fnfmt
        kwgkwg = {} if not 'kwargs' in kwg else kwg['kwargs']
        format_dt = AutoDict({key: cls.basenaming(item) for key, item in {**locals(), **kwg, **kwgkwg}.items()})
        if 'ftype' in format_dt: format_dt.pop('ftype')
        out_fnfmt = out_fnfmt.format_map(format_dt)
        if testsave: # saving quick check, before lots of work is done
            #out_fn = out_fnfmt.format(ext='tmp') +'.tmp'
            out_fn = out_fnfmt.format_map(dict(ftype='tmp')) #+ f'{np.random.randint(0,10**6):07}' + '.tmp'
            pd.DataFrame(['Currently being computed']) \
            .to_csv(out_fn, index=False, header=False);
            os.remove(out_fn) # briefly uncommented this to see doulbe slurm submission issue on mgh cluster.
        mainout_fn = out_fnfmt.format_map(dict(ftype=ftype))
        if os.path.isfile(mainout_fn):
#             msg = f"\033[1;31mWARNING:\033[0m The file {mainout_fn} already exists! If this code finishes, it will be overwritten."
            msg = f"\033[1;31mWARNING: {mainout_fn} already exists! If you let this code finish, it will be overwritten.\033[0m"
            #msg = f'WARNING: The file {mainout_fn} already exists! If this code finishes it will be overwritten.'
            warnings.warn(msg)
        if prstlogs:
            prstlogs_fn = out_fnfmt.format(ftype='json'); dn, fn = os.path.split(prstlogs_fn)
            prstlogs_fn = os.path.join(dn, '.prstoolslogs', fn)
            prstlogs = prst.utils.get_prstlogs()
            prstlogs.set_prstlogs_fn(prstlogs_fn, save=True)
        # out_fnfmt can be a completed file name or a string that f'{still}{has}{things}{that_have_to_be_filled_in}'
        # However, {ftype} (==filetype) will never be filled in, so you can have file.log and file.results.
        return out_fnfmt
        
    def save_weights(self, fn, return_weights=False, ftype='auto', extra_weight_cols=None, nancheck=None, end='\n\n'):
        options=['legacyweights.tsv','prstweights.tsv','prstweights.h5','prstweights.parquet']# give one of these extensions for auto
        nancheck = self._nancheck if nancheck is None else nancheck
        if extra_weight_cols is None: extra_weight_cols = self.extra_weight_cols
        cols=list(self.default_weight_cols)
        if ftype=='auto':
            for opt in options: 
                if fn.endswith(opt): ftype=opt
            if ftype=='auto': ftype=self.default_weightfile_type
        if ftype == 'legacyweights.tsv': 
            ext = 'legacyweights.tsv'; header=False; selcols = list(cols)
        elif ftype.split('.')[0] == 'prstweights':
            header=True
            ewc = extra_weight_cols
            if not ewc: extcols = []
            else: extcols = [col for col in self.get_weights().columns if col not in cols] if type(ewc) is bool and ewc else ewc
            selcols = cols + extcols
        else:
            raise Exception('Model weight file saving format could not be properly determined.')
            
        fn = fn.format_map(dict(ftype=ftype)) # Maybe some AutoDict buzz here later.
        import uuid; tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"  # unique temp file name
        if self.verbose: print(f'Saving model weights (filetype={ftype}) to: {fn}', end=' ')
        if nancheck: assert np.sum(self.get_weights()['allele_weight'].isna().sum()) == 0
        fin_df = self.get_weights()[selcols]
        if ftype.endswith('prstweights.h5'): pd.DataFrame(fin_df.to_numpy(), index=fin_df.index, columns=fin_df.columns).to_hdf(tmp_fn, key='df')
        elif ftype.endswith('prstweights.parquet'): fin_df.to_parquet(tmp_fn)
        else: fin_df.to_csv(tmp_fn, sep='\t', index=False, header=header)
        os.replace(tmp_fn, fn) # atomically move into place
        if self.verbose: print(f'-> Done', end=end)
        if return_weights: return df
        
    def save_sst(self, fn, return_sst=False, ftype='tsv', basecols=None, addicols=None, nancheck=False):
        kwg = {}
        if addicols is not None: kwg['addicols'] = addicols
        if basecols is not None: kwg['basecols'] = basecols
        out_df = prst.loaders.save_sst(sst_df=self.sst_df, fn=fn, return_sst=return_sst, ftype='tsv', 
                    nancheck=nancheck, verbose=self.verbose, **kwg)
        if return_sst: return out_df
    
    def get_pbar(self, iterator, *, make_range_var=True, **kwg):
        # Maybe import funny wrapper class
        msg = f"Object implement iterator has no length. This is required. More info: {type(iterator)=} , {iterator=}"
        assert hasattr(iterator, '__len__'), msg
        if make_range_var: # This is something we really want because we dont want a part of,
            # the linkdata to remain stuck in a pbar (and hence stuck in memory...)
            iterator = range(len(iterator))
        pbar = prst.utils.get_pbar(iterator, **kwg)
        return pbar
    
    def get_iterator(self, init_iterator, pbar=None, **kwg):
        if not pbar:
            for elem in init_iterator:
                yield elem
        else:
            if type(pbar) is bool:
                self.pbar = self.get_pbar(init_iterator, **kwg)
            assert hasattr(self.pbar,'update'), 'pbar not an tqdm-like pbar object, while it should be.' 
            
            try:
                for elem in init_iterator:
                    yield elem
                    self.pbar.update()
            except Exception as e:
                # Optional: log or print here
                raise  # re-raise the original exception
            finally:
                if self._close_pbar:
                    self.pbar.close(); self.pbar=True
    
    def get_params(self):
        out=copy.deepcopy(self._kwg_dt); out.pop('_excl_lst',None)
        return out
    
    def get_linkdata(self):
        if hasattr(self,'_linkdata'): return self._linkdata
        else: raise Exception(f'The model being run ({self}) is requesting linkdata, but it is not present. use model.fit(linkdata), or use model.set_linkdata(linkdata)')
    
    @property
    def linkdata(self):
        return self.get_linkdata()
        
    def get_weights(self, return_frame=True):
        if not return_frame: raise NotImplementedError('contact dev') 
        if hasattr(self, 'weights_df'): return self.weights_df
        else: raise Exception('No weights present for this model, which probably means it was not run.')
            
    def set_linkdata(self, linkdata, requires_attrs=None, warn=True, ignore_none=True):
        if linkdata is None and ignore_none: return None
        self._checktype(linkdata, 'LinkageData')
        if requires_attrs is not None:
            for elem in requires_attrs:
                assert hasattr(linkdata, elem) and getattr(linkdata, elem) is not None
        if getattr(self, '_linkdata', None) is not None and warn:
            warnings.warn(f'There is already linkdata present for {self}, but is now overwitten with new linkdata input.')
        self._linkdata = linkdata
        return self
            
    def _set_weights(self, weights_df, sort=True, reset_index=True, silentsort=False):
        if not isinstance(weights_df, pd.DataFrame): # can this be done with decorator?
            raise TypeError("Input must be a DataFrame.")
        if not weights_df.shape[0] == len(weights_df['snp'].unique()): # 200 ms (== ok, longest step of this method)
            raise Exception('Duplicate snp-ids present. This is not allowed at the '
                            'moment (meaning multiallelic snp are not possible).')
        assert all(col in weights_df.columns for col in self.default_weight_cols)
        nan_ser = weights_df[self.default_weight_cols].isna().sum()
        assert nan_ser.sum() == 0, (f'There are nan values in the weights, this is not right! All NaNs '
                                    f'need to be removed. Printing NaN counts for respective columns: \n {nan_ser}')
        # If there are errors here it could be because there is no 
        # implemented X/Y/MT chrom functionality (mapping needed). Contact dev.
        possiblysortedweights_df = weights_df.sort_values(['chrom','pos']) if sort else weights_df
        if not (possiblysortedweights_df.index == weights_df.index).all():
            if not silentsort:
                warnings.warn('Input weights were sorted on chromosome (chrom) and position (pos), since inputs weren\'t. '
                              'This is usually only important if using this code inside of python and not if using the prstools commandline.')
        weights_df = possiblysortedweights_df
        self.weights_df = weights_df.reset_index(drop=True) if reset_index else weights_df
        
    def remove_linkdata(self):
        self.linkdata.clear_linkage_allregions()
        del self._linkdata
        
    def clone(self):
        return self.__class__(**self.get_params())
    
    def _compute_sst_inside_pred(**kwg):
        raise NotImplementedError('whoops not implemented this yet')
        cols = ['chrom', 'snp', 'pos', 'A1', 'A2']
        chunk_sst_df = wchunk_df[cols].copy()
        chunk_sst_df.columns = chunk_sst_df.columns.get_level_values(0)
        chunk_sst_df['std_psd'] = s
        X = chunk_sda.val ### OKE, the standard deviation is not always gonna be 1 here, it should be
        #aa=X.mean(axis=0)
        #bb=X.std(axis=0)
        #Xs = X-aa
        #Xs = Xs/bb
        Xs = (X-X.mean(axis=0))/X.std(axis=0)
        if not 'nsamps' in locals() or not 'ready_df' in locals():   
            nsamps = trait_df.shape[0] - trait_df.isna().sum()
            ready_df = trait_df.fillna(0).astype(float)
            # ready_df = (ready_df-ready_df.mean())/ready_df.std()# hey! this should NOT be done here
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = wchunk_df['flipper'].values[:,None]*(Xs.T@ready_df)/nsamps
            r.index = wchunk_df.index
        chunk_sst_df = pd.concat([chunk_sst_df,r], axis=1)
        sst_lst += [chunk_sst_df]
        return stuff
    
    def predict(self, bed, *, n_inchunk=1000, groupby=None, validate=True,
                localdump=False, weight_type='allele', trait_df=None, colour='#7f00ff'): # <-- The more esotheric stuff on this line
        
        if 'pysnptools' in str(type(bed)):
            srd = bed; del bed
            return self.srdpredict(**locals())
        
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        if self.verbose: print(f'Predicting phenotypes i.e. generating PRS (in chucks of {n_inchunk} snps).', flush=True)
        assert weight_type in ('allele','standardized'); 
        if trait_df is not None: assert localdump
        weights_df = self.get_weights()
        
        if validate:
            bim_df = bed.bim_df.copy() # For the next line did it the other way around for mem-footprint.
            weights_df = prst.loaders.merge_snps(weights_df, bim_df, req_all_right=False, handle_missing='filter', flipcols=[])
            weights_df['allele_weight']=weights_df['allele_weight']*weights_df['rflip'] # This is a bit of a hack
            weights_df = weights_df.sort_values('xidx')
        
        # Loop through Genome:
        yhat_dt = dict(); sst_dt = {}
        if len(weights_df['allele_weight'].shape) == 1: n_traits = 1
        else: n_traits = weights_df['allele_weight'].shape[1]
        #for itr in self.get_iterator(range(n_iter), pbar=self.pbar)
        for grp, wgrp_df in self.get_iterator(weights_df.groupby(groupby), pbar=self.pbar, colour=colour) if groupby is not None else [(None, weights_df)]:
            yhat = np.zeros((bed.iid_count, n_traits)); sst_lst = []
            inner_pbar = self.pbar if grp is None else None
            for start in self.get_iterator(range(0, wgrp_df.shape[0], n_inchunk), pbar=inner_pbar, colour=colour):
                wchunk_df = wgrp_df.iloc[start:start+n_inchunk]
                cxidx = wchunk_df['xidx']
                X = bed.read(index=np.s_[:,cxidx])
                m = np.nanmean(X, axis=0)
                idx = np.where(np.isnan(X))
                s = np.nanstd(X, axis=0) if weight_type == 'standardized' else None
                X[idx] = np.take(m, idx[1])
                w = wchunk_df['allele_weight']; w=w if type(w) is pd.DataFrame else w.to_frame(name='prs')
                if weight_type == 'standardized': w = s*w
                yhat += X@w.values.astype(X.dtype) #chunk_df['allele_weight']
                if trait_df is not None: # Compute beta marginal too if required
                    self._compute_sst_inside_pred(**locals())

            columns = w.columns # considering doing something special with ('prs',f'{colname}') here.. \newline
            # , but multiindex will give funny/bad-4-users prs pred files downstream so..
            yhat = pd.DataFrame(yhat, index=pd.MultiIndex.from_arrays(bed.fam_df[['fid','iid']].values.T, names=["fid", "iid"]), columns=columns)
            if trait_df is not None: sst_dt[grp] = pd.concat(sst_lst, axis=0)
            yhat_dt[grp] = yhat

        output = yhat if groupby is None else yhat_dt
        if localdump: output=locals()
        return output
            
    def srdpredict(self, srd, *, n_inchunk=1000, groupby=None, check='depreciated-arg', validate=True, 
                localdump=False, weight_type='allele', trait_df=None, colour=None): # <-- The more esotheric stuff on this line
        
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        
        string = '(matching inputs now)' if validate else ''
        if self.verbose: print(f'Predicting phenotypes for given snp inputs i.e. generating PRS (done in chucks of {n_inchunk} snps){string}.', flush=True)
        toc('here now')
        print('egrreg',flush=True)
        assert weight_type in ('allele','standardized'); 
        if trait_df is not None: assert localdump
        weights_df = self.get_weights()
        print('after getweights',flush=True)
        if validate:
            toc('really staring validation')
            #import IPython as ip; ip.embed()
            subsrd = srd
            for _ in range(5): subsrd = subsrd if hasattr(subsrd, 'count_A1') else subsrd._internal
            assert subsrd.count_A1 is True, 'Have to set count_A1=True, for this function to work, reinitialize snpreader with count_A1=True'
            ind = weights_df.snp.isin(srd.sid)
            if not ind.all(): raise Exception('Weight snps not in target. This is required, consider imputing target.')
            from prstools.loaders import load_bimfam_from_srd
            bim_df, fam_df = load_bimfam_from_srd(srd, skipifpresent=True)
            # First we validate SNP id alignment:
            
            if weights_df.snp.shape[0] != srd.sid.shape[0] or not np.all(weights_df.snp == bim_df.snp): # Do things to make this true
                toc('statring sid_to_index')
                weights_df['idx_srd'] = srd.sid_to_index(weights_df.snp) # ask a small # of snps to sid_to_index asks about 10G for 22M 1kg snpset..
                # corresponding bim is about 700MB and bed 2.8G (but should not be loaded..), not great...
            toc('continue after sid2index')
                #srd = srd[:,srd.sid_to_index(weights_df.snp)]
            #assert np.all(srd.sid == weights_df.snp)
            xbim_df = bim_df.set_index('snp',drop=True).loc[weights_df.snp].reset_index(drop=False)
            # Second-ly we validate Allele alignment:
            ind_match = (xbim_df['A1'] == weights_df['A1']) & (xbim_df['A1'] == weights_df['A1'])
            ind_flip  = (xbim_df['A1'] == weights_df['A2']) & (xbim_df['A1'] == weights_df['A2'])
            ind_wrong = ~ind_match & ~ind_flip # not matching at all, happens for example with trialllelic snps
            assert ind_wrong.sum()==0, 'There are SNPs that cannot be matched. This is probably because trialllelic snps'
            cast=int
            weights_df['flipper'] = -1*ind_flip.astype(cast) + 1*ind_match.astype(cast)
            if trait_df is not None:
                assert (trait_df.index.to_frame(index=False).to_numpy().astype(srd.iid.dtype) == srd.iid).all()
            toc('done with validation')
        
        # Loop through Genome:
        yhat_dt = dict(); sst_dt = {}
        if len(weights_df['allele_weight'].shape) == 1: n_traits = 1
        else: n_traits = weights_df['allele_weight'].shape[1]
        toc('really starting pRS prediction looop')
        for grp, wgrp_df in self.pbar(weights_df.groupby(groupby)) if groupby is not None else [(None, weights_df)]:
            yhat = np.zeros((srd.shape[0], n_traits)); sst_lst = []
            inner_pbar = self.pbar if grp is None else lambda x: x
            for start in inner_pbar(range(0, wgrp_df.shape[0], n_inchunk)):
                wchunk_df = wgrp_df.iloc[start:start+n_inchunk]
                chunk_srd = srd[:,wchunk_df['idx_srd']]
                chunk_sda = chunk_srd.read()
                chunk_sda, chunk_stansda = chunk_sda.standardize(return_trained=True) 
                # Reason for standardisation here is to allow for straightforward dealing with NaN values, 
                # since they are automatically mean imputed by pysnptools.
                s=chunk_stansda.stats[:,[1]]; m=chunk_stansda.stats[:,[0]]; f=wchunk_df[['flipper']].values
                wstan = np.ones(s.shape) if weight_type == 'standardized' else s
                w = wchunk_df['allele_weight']; w=w if type(w) is pd.DataFrame else w.to_frame(name='prs')
                cur_wtilde = f*wstan*w
                yhat += chunk_sda.val@cur_wtilde.values.astype(chunk_sda.val.dtype)
                
                if trait_df is not None: # Compute beta marginal too if required
                    cols = ['chrom', 'snp', 'pos', 'A1', 'A2']
                    chunk_sst_df = wchunk_df[cols].copy()
                    chunk_sst_df.columns = chunk_sst_df.columns.get_level_values(0)
                    chunk_sst_df['std_psd'] = s
                    X = chunk_sda.val ### OKE, the standard deviation is not always gonna be 1 here, it should be
#                     aa=X.mean(axis=0)
#                     bb=X.std(axis=0)
#                     Xs = X-aa
#                     Xs = Xs/bb
                    Xs = (X-X.mean(axis=0))/X.std(axis=0)
                    if not 'nsamps' in locals() or not 'ready_df' in locals():   
                        nsamps = trait_df.shape[0] - trait_df.isna().sum()
                        ready_df = trait_df.fillna(0).astype(float)
                        # ready_df = (ready_df-ready_df.mean())/ready_df.std()# hey! this should NOT be done here
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        r = wchunk_df['flipper'].values[:,None]*(Xs.T@ready_df)/nsamps
                        r.index = wchunk_df.index
                    chunk_sst_df = pd.concat([chunk_sst_df,r], axis=1)
                    sst_lst += [chunk_sst_df]
                 

            columns = cur_wtilde.columns # considering doing something special with ('prs',f'{colname}') here.. \newline
            # , but multiindex will give funny/bad-4-users prs pred files downstream so..
            yhat = pd.DataFrame(yhat, index=pd.MultiIndex.from_arrays(srd.iid.T, names=["fid", "iid"]), columns=columns)
            if trait_df is not None: sst_dt[grp] = pd.concat(sst_lst, axis=0)
            yhat_dt[grp] = yhat

        output = yhat if groupby is None else yhat_dt
        if localdump: output=locals()
        return output
    
class MultiPred():
    
    @classmethod
    def from_weights(cls, weights, sort=True, **kwg):
        if not isinstance(weights, pd.DataFrame): # can this be done with decorator?
            raise TypeError("Input must be a DataFrame.")
        if type(weights) is dict: raise NotImplementedError()
        # A check for the required columns is also needed, somewhere.
        model = cls(**kwg)
        model._set_weights(weights, sort=sort)
        return model
    
    @classmethod
    def from_dict(cls, weights_dt, ref_df=None, verbose=False, greedy=False, on=None, remove_allnan=False, **kwg):
        assert not greedy, 'Greedy options not implemented yet.'
        assert len(weights_dt) > 0, 'weights_dt is empty'
        if ref_df is not None: allweights_df=ref_df.copy()
        elif len(weights_dt) == 1: return cls.from_weights(list(weights_dt.values())[0].copy(), verbose=verbose, **kwg)
        else: raise NotImplementedError(f'... contact dev if this option (ref_df != None) is desired.')
        on_dt = dict() if on is None else dict(on=on)
        if len(on_dt) == 0: on=['snp','A1','A2']
        assert allweights_df[on].isna().sum().sum() == 0, f'cannot have Nans in starter columns on={on}'
        msg = 'Duplicated snp present, this is probably cause of multi-allelic snps, this needs to be implemented contact dev'
        assert allweights_df['snp'].duplicated().sum() == 0, msg
        
        # merging mechanics:
        for wname, curweights_df in tqdm(weights_dt.items()):
            curweights_df = curweights_df.copy() # This line is crucial for the PandasMimic type used to load from disk.
            cmissing = set(cls.default_weight_cols) - set(curweights_df.columns)
            allweights_df = prst.loaders.merge_snps(allweights_df, curweights_df, flipcols=['allele_weight'], how='left', handle_missing='keep', **on_dt)
            allweights_df = allweights_df.rename(columns=dict(allele_weight=f"allele_weight_{wname}")).drop('rflip',axis=1)
#             missing = set(cls.default_weight_cols) - set(allweights_df.columns)
#             if missing: raise ValueError(f"Missing columns: {missing}")

        #Create multi-index columns:
        lst = [(col, '') if not col.startswith('allele_weight') else ('allele_weight', col.split('allele_weight_')[-1]) for col in allweights_df.columns]
        newcols = pd.MultiIndex.from_tuples(lst)
        allweights_df.columns = newcols
        assert (allweights_df['allele_weight'] == 0).sum().sum() == 0
        if remove_allnan:
            ind = allweights_df['allele_weight'].isna().sum(axis=1) < allweights_df['allele_weight'].shape[1]
            allweights_df = allweights_df[ind]
        allweights_df = allweights_df.fillna(0)
        
        model = cls.from_weights(allweights_df, verbose=verbose, **kwg)
        return model
            
    @classmethod
    def from_path(cls, path_or_list, 
                  ref_df=None,
                  ftype='auto', pyarrow=True, remove_allnan=False,
                  sep:str='\t', # Seperator for the inputs files, default is \t (tab).
                  verbose=False, on=None, **kwg):
            
        if type(path_or_list) is str: fn_lst=glob.glob(path_or_list)
        else: fn_lst=path_or_list; assert type(path_or_list) is list, 'Input to this function should be path or list.'
        msg = f'No file(s) {path_or_list} found (yes this code matches multiple files e.g. ./thesis_v*.tex)'
        assert len(fn_lst) > 0, msg
        if type(ref_df) is str: ref_df, _ = prst.loaders.load_bimfam(ref_df, fam=False, verbose=verbose)
        class PandasMimic():
            def __init__(self,fn, ftype):
                # potentially we can do a rapid check here for the weight loading.
                assert type(fn) is str and type(ftype) is str
                self.fn = fn; self.ftype=ftype; self.pyarrow=pyarrow
                self.sep=sep
            def copy(self):
                return prst.loaders.load_weights(**vars(self))
        weights_dt=dict()
        for fn in fn_lst:
            weights_dt[fn]=PandasMimic(fn=fn, ftype=ftype)
            
        model = cls.from_dict(weights_dt, ref_df=ref_df, on=on, remove_allnan=remove_allnan, verbose=verbose, **kwg)

        return model
    
from joblib import Parallel, delayed 
class GroupByModel(MultiPred, BasePred):
    
    def __init__(self, _model, *, groupby, n_jobs=BasePred._default_n_jobs, pbar:bool=True, verbose=False, **xtras):
        
        # Stuff all the args into fields.
        _excl_lst = ['self', 'kwg_dt']
        kwg_dt = {key: item for key, item in locals().items() if not (key in _excl_lst)}
        for key, item in locals().items():
            if not (key in _excl_lst): 
                self.__setattr__(key, item)
        self._kwg_dt = copy.deepcopy(kwg_dt) 
        
    def clone(self):
        raise Exception(f'{self.__class__} cannot be cloned.')
        
    def get_model_clone(self):
        return self._model.clone()
        
    def fitold(self):

        linkdata = self.get_linkdata()
        self.model_dt = dict()
        
        def passthrough(arg):
            return arg
        
        if self.verbose: print('Starting iterations of model(s):')
        assert type(self.groupby) is str, 'groupby must be string, if you want to use multiple columns combined then contact dev.'
        
        
        
        nuniq = linkdata.get_sumstats_cur()[self.groupby].nunique()
        tot_iters = getattr(self._model,'n_iter',1)*nuniq
        pbar = self.get_pbar(iterator=range(tot_iters))
        contents = {key: item for key, item in linkdata.groupby(self.groupby, sort=True, skipempty=True)}
        del linkdata
        self.remove_linkdata()
        
        #  tqdm(linkdata.groupby(self.groupby, sort=True, skipempty=True), total=22)
        # this loop in a multi processed way?
        #prst.utils.get_ip().embed()
        print(sys.argv)
        if 'testmulti' in ' '.join(sys.argv) or '/opt/conda/lib/python3.11/site-packages/ipykernel_launcher.py' in ' '.join(sys.argv):
            from prstools.utils import save_to_interactive; save_to_interactive(dict(loc_dt=locals()))
            crash()
         
            
        for grp, cur_linkdata in linkdata.groupby(self.groupby, sort=True, skipempty=True):
            model = self.get_model_clone()
            model.verbose = False; model.pbar = pbar
            model._close_pbar = False
            model.fit(cur_linkdata)
            self.model_dt[grp] = model
        pbar.close(); self.pbar=True
        self.combine_set_weights()
        return self
    
    
    def fit(self, linkdata=None):
        
        ## Prep portion
        def worker(model, cur_linkdata, pbar, grp):
            model.verbose = False
            model.pbar = pbar
            model.fit(cur_linkdata)
            return grp, model
        self.set_linkdata(linkdata, ignore_none=True)
        linkdata = self.get_linkdata()
        self.model_dt = dict()
        if self.verbose: print('Starting iterations of model(s):')
        assert type(self.groupby) is str, 'groupby must be string, if you want to use multiple columns combined then contact dev.'
        nuniq = linkdata.get_sumstats_cur()[self.groupby].nunique()
        tot_iters = getattr(self._model,'n_iter',1)*nuniq
        contents = {key: item for key, item in linkdata.groupby(self.groupby, sort=True, skipempty=True)}
        del linkdata; self.remove_linkdata()
        
        # MultiProcessing portion:
        #with Manager() as manager:
        fakebar = prst.utils.FakeMultiprocPbar(manager=None) if self.pbar else False
        real_pbar = self.get_pbar(iterator=range(tot_iters), fakebar=fakebar, deamon=True) if self.pbar else False
        if self.pbar: # Do fakebar hacks to make everything work
            mgr = fakebar._mgr; fakebar._mgr=None
        results = Parallel(n_jobs=self.n_jobs, max_nbytes=None)(
            delayed(worker)(self.get_model_clone(), cur_linkdata, fakebar, grp)
            for grp, cur_linkdata in contents.items())
        if self.pbar: 
            real_pbar.close(); real_pbar=None
            prst.utils.clear_memory(); # Crucial line because gc.collect() inside, else thing go wrong later.
            fakebar.close(); mgr.shutdown(); mgr=None
        for grp, model in results: self.model_dt[grp] = model
        self.combine_set_weights()
        return self

    def combine_set_weights(self):
        assert hasattr(self,'model_dt'), f'No models present, so cannot create a working weights set for {self}.'
        weights_df = pd.concat([model.get_weights() for grp, model in self.model_dt.items()], axis=0) #for grp, model in self.model_dt.items():
        self._set_weights(weights_df, silentsort=True)
        
class PredPRS(MultiPred, BasePred):
    
    "(not implemented yet) PredPRS: It predict polygenic risk scores if you give it weights (aka. give it the files produced by the methods in PRSTOOLS)(internal: PredPRS.from_weights())."
        
    def __init__(self, *,
         pbar:bool=True,               # [not functional yet] Display a progress bar during optimization. 
         verbose:bool=False):           # Verbose mode when fitting the model.
        
        # Stuff all the args into fields.
        _excl_lst = ['self', 'kwg_dt']
        kwg_dt = {key: item for key, item in locals().items() if not (key in _excl_lst)}
        for key, item in locals().items():
            if not (key in _excl_lst): 
                self.__setattr__(key, item)
        self._kwg_dt = copy.deepcopy(kwg_dt)
        if not self.pbar or not verbose: self.pbar = lambda x: x  
        else: self.pbar = tqdm
            
    @classmethod
    def from_cli_params_and_run(cls, *,  
        ref, # this is ref bkabkjfjkegkhjr kjhergjkerjhk
        chrom=[1,2,3], # blablablbal lablablab moar blablbal
        fnfmt='_.{ftype}', ftype='prstweights.tsv',                
        verbose=True,  # jkerkjgerkjghjk jkhergkjh
        pkwargs=None, 
        out=None, 
        return_model=True, 
        fit=True ,
        colmap=None
        ):
        
        from prstools.loaders import LinkageData
        
        # Gen output file name format and do quick check if output file can be saved before a lot of work is done:
        if out: out_fnfmt = cls.create_output_fnfmt(**locals())
        
        # Initialize data & model objects and fit the model:
        assert colmap is None, 'colmap not implemented for these non sparse prs cli tools yet, complain to the developer!'
        linkdata = LinkageData.from_cli_params(ref=ref, target=target, sst=sst, n_gwas=n_gwas, chrom=chrom)
        def testkey(key): return (key in pkwargs) if pkwargs else True
        model = cls(**dict({key: item for key, item in kwargs.items() if testkey(key)}, verbose=verbose))
        model.set_linkdata(linkdata)
        if fit: model.fit()
        
        # Store the results:
        if out: model.save_weights(out_fnfmt, ftype=ftype)

        return model
       
try:
    profile
except NameError:
    def profile(func):
        return func

class PRSCS2(BasePred):
    
    "PRS-CS v2: A polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors."
    _gig = None
    
    def __init__(self, *,
         n_iter=1000,              # Total number of MCMC iterations.
         n_burnin=0.5,             # Number of burn-in iterations if larger than 1 or fraction of n_iter if smaller then 1.
         n_slice=1,                # Thinning of the Markov chain.
         shuffle=False,
         seed=-1,                  # Random seed for reproducibility.
         a=1.0,                    # Parameter a in the gamma-gamma prior.
         b=0.5,                    # Parameter b in the gamma-gamma prior. 
         phi=-1.,                  # Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a Bayesian approach
         clip=1.,                  # Clip parameter. The default works best in pretty much all cases.
         sampler='Rue',            # Sampler algorithm. Rue sampling is the original sampler, which gives good results.
         groupby:str='chrom',      
         local_rm:bool=False,    
         compute_score:bool=False,   
         clear_linkdata:bool=True,
         pop='pop',
         n_jobs=BasePred._default_n_jobs, # This sets the number of jobs for parallel processing.
         pbar:bool=True,              
         verbose:bool=False): 
        
        # Stuff all the args into fields.
        _excl_lst = ['self', 'kwg_dt']
        kwg_dt = {key: item for key, item in locals().items() if not (key in _excl_lst)}
        for key, item in locals().items():
            if not (key in _excl_lst): 
                self.__setattr__(key, item)
        self._kwg_dt = copy.deepcopy(kwg_dt)
        
        if self.seed == -1: self.seed = None
        #if not self.pbar: self.pbar = lambda x: x
        #else: self.pbar = tqdm if pbar is None or type(pbar) is bool else pbar
        if self.phi == -1: self.phi=None
        self.do_phi_updt, self.phi = (True, 1.0) if self.phi is None else (False, self.phi)
        if self.phi is not None: assert self.phi > 0
        n_burnin = int(n_burnin*n_iter) if n_burnin < 1 else int(n_burnin)
        self.n_burnin = n_burnin
        assert (n_iter-n_slice) > n_burnin
        self.sampler=str(sampler).lower()
        self.pop = self.pop.upper()
        assert self.sampler in ['rue','bhat','sld']
    
    
#     def _gig(self, p,a,b, psi=None):
#         x = np.zeros(b.shape) if psi is None else psi
#         for j in range(b.shape[0]): # This loop gets everything back in shape. 
#             x[j] = gigrnd(p, a[j], b[j])
#             #psi[j] = gigrnd(a-0.5, 2.0*delta[j], n_eff*beta[j]**2/sigma)#, seed=seed)
#         # else: raise ValueError(f"Option not recognized: {self.gigsampler}")
#         return x
    
    def _compute_beta_tilde(self, *, beta, i_reg, linkdata):
        beta_tilde = linkdata.get_beta_marginal_region(i=i_reg)
        if self.local_rm: # RM
            raise NotImplementedError()
        return beta_tilde
    
    @profile
    def fit(self, linkdata=None):
        
        # Loading variables:
        self.set_linkdata(linkdata, ignore_none=True)
        s=self; linkdata=s.linkdata; 
        n_burnin=s.n_burnin; n_slice=s.n_slice; 
        n_iter=s.n_iter; n_pst=(n_iter-n_burnin)/n_slice
        a=s.a; b=s.b; phi=s.phi
        verbose=s.verbose; do_phi_updt=self.do_phi_updt
        beta_mrg = linkdata.get_beta_marginal()
        p        = len(beta_mrg)
        n_eff    = linkdata.get_sumstats_cur()['n_eff'].median()
        i_lst    = linkdata.get_i_list()
            
        # Initalisations:
        if self.seed != None: np.random.seed(self.seed)
        beta=np.zeros((p,1)); beta_est=np.zeros((p,1)); beta_ml=np.zeros((p,1))
        psi=np.ones((p,1)); psi_est=np.zeros((p,1)); self.scores=[]
        sigma=1.; sigma_est=0.; phi_est=0.;
        #if self.pbar and type(self.pbar)is bool self.pbar = tqdm
        
        # Sampling Loops:
        if verbose: print('Starting iterations of Sampler:')
        for itr in self.get_iterator(range(n_iter), pbar=self.pbar):
            quad = 0; i_reg=None
            if not self.pbar:
                do_show = ((itr % 10 == 0) | (itr<3)) & verbose
                if do_show: print(f'-> itr={itr}, i_reg={i_reg} <-  ', end='\r')
            for i_reg in self._order(i_lst):

                # Compute beta_tilde, a corrected GWAS sumstat zscore [=RM]:
                beta_tilde = self._compute_beta_tilde(beta=beta, i_reg=i_reg, linkdata=linkdata)

                # Sample beta from MVN:
                s2 = sigma; s=np.sqrt(s2)
                idx_reg = range(*linkdata.get_range_region(i=i_reg));
                if self.sampler == 'rue':
                    D = linkdata.get_linkage_region(i=i_reg)
                    dinvt = D + np.diag(1.0/psi[idx_reg].T[0])
                    test = dinvt@beta_tilde
                    dinvt_chol = linalg.cholesky(dinvt)
                    beta_tmp = (linalg.solve_triangular(dinvt_chol, beta_tilde, trans='T') +
                                np.sqrt(sigma/n_eff)*np.random.randn(len(D), 1))
                    beta[idx_reg] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                    quad += np.dot(np.dot(beta[idx_reg].T, dinvt), beta[idx_reg])              
                elif self.sampler == 'bhat':
                    D     = linkdata.get_linkage_region(i=i_reg)
                    Dh    = linkdata.get_specified_data_region(i=i_reg, varname='Dhalf') 
                    T     = np.diag(psi[idx_reg].T[0])*s2/n_eff
                    u     = np.random.randn(len(D),1)*np.sqrt(T.diagonal()[:,np.newaxis])
                    seedk = np.random.randn(len(D),1)
                    beta_tilde_samp = beta_tilde/s2 - D@u/s2 - (Dh@seedk)/(s*np.sqrt(n_eff)) #(D@Dih@seedk)
                    Q = D@T/s2 + np.eye(len(D))/n_eff
                    beta[idx_reg] = u + T@linalg.solve(Q, beta_tilde_samp)
                    quad += np.dot(np.dot(beta[idx_reg].T, D+np.diag(1.0/psi[idx_reg].T[0])), beta[idx_reg])
                elif self.sampler == 'sld':
                    D      = linkdata.get_linkage_region(i=i_reg) # This is not ok, mod later.
                    Di     = linkdata.get_precision_region(i=i_reg) #auto_linkage_region(i=i_reg)
                    Dih    = linkdata.get_specified_data_region(i=i_reg, varname='Dihalf') 
                    #if not np.any(beta_ml[idx_reg]): beta_ml[idx_reg] = Di@beta_tilde
                    T     = np.diag(psi[idx_reg].T[0])*s2/n_eff
                    u     = np.random.randn(len(Di),1)*np.sqrt(T.diagonal()[:,np.newaxis])
                    seedk = np.random.randn(len(Di),1)
                    k     = Dih@seedk
                    alpha = Di@beta_tilde-u-(s/np.sqrt(n_eff))*k
                    W     = (T + (s2/n_eff)*Di)
                    beta[idx_reg] = u + T@linalg.solve(W, alpha)
                    quad += np.dot(np.dot(beta[idx_reg].T, D+np.diag(1.0/psi[idx_reg].T[0])), beta[idx_reg])
                else:
                    raise Exception('Sampler not recognized:', self.sampler)
                
            if self.compute_score:
                if callable(self.compute_score): score = self.compute_score(**locals())
                else: score = n_eff/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad)
                self.scores.append(score)
                
            # Stuffs: (more tweaking prob needed)
            err = max(n_eff/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n_eff/2.0*sum(beta**2/psi))
            sigma = 1.0/np.random.gamma((n_eff+p)/2.0, 1.0/err)
            delta = np.random.gamma(a+b, 1.0/(psi+phi))

            # Sample Variance of the Weight prior:
            if self._gig: psi = self._gig(a-0.5, 2.0*delta, n_eff*beta**2/sigma, psi=psi)
            else:
                for j in range(p): psi[j] = gigrnd(a-0.5, 2.0*delta[j], n_eff*beta[j]**2/sigma)
            if self.clip: psi[psi>self.clip] = self.clip #Clipping.

            # Sample Phi or continue with set value:
            if self.do_phi_updt == True: # Could be tweaked with range_p_filter for speed.
                w = np.random.gamma(1.0, 1.0/(phi+1.0))
                phi = np.random.gamma(p*b+0.5, 1.0/(sum(delta)+w))

            # Posterior:
            if (itr>n_burnin) and ((itr%n_slice)==0):
                beta_est = beta_est + beta/n_pst
                psi_est = psi_est + psi/n_pst
                sigma_est = sigma_est + sigma/n_pst
                phi_est = phi_est + phi/n_pst
                
        #for me not run should
        #Post proc & storage:
        weights_df = linkdata.get_sumstats_cur().copy()
        weights_df['raw_weight'] = beta_est
        weights_df['allele_weight'] = beta_est/linkdata.get_allele_standev(source=self.scaling)
        self.weights_df = weights_df
        if self.clear_linkdata: self.remove_linkdata()
        if callable(self.compute_score): itr=-1; self.compute_score(**locals())
        if verbose: print('----- Done with Sampling -----')
        return self


if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+
          ['ipykernel_launcher.py' in sys.argv[0]] + 
          [not '__file__' in locals()]):

    if 'In' in locals() and _isdevenv_prstools:
        code = In[-1] 
        with open('../prstools/models/_base.py', 'w') as f: f.write(code)
        print('Written to:', f.name); time.sleep(0.03)
        print('starting here in models:')
        get_ipython().system('time python ../prstools/_cmd.py --dev')
        #!prst --dev | head -3
    print('Done')