from abc import ABC, abstractmethod
import copy, time, warnings, math, traceback, sys, os, glob
import scipy as sp
import numpy as np
import pandas as pd
from scipy import linalg, stats
import prstools as prst
from prstools.models._compute import dpsi, gigrnd, g
from prstools.utils import PRSTCLI
try:
    from fastcore.script import call_parse, Param
except:
    def Param(*args, **kwg):
        return None
# PRSTCLI = None  
# AutoPRSTCLI 

# __all__ = [k for k, v in globals().items() if not k.startswith('_') and v.__module__ == __name__]

# __all__ = [
#     name for name, val in globals().items()
#     if not name.startswith('_')
#     and callable(val) 
#     and getattr(val, '__module__', None) == __name__
# ]

__all__ = ['BasePred','BaseMulti','GroupByModel','MultiPRS','PRSCS2']
# __all__ = ['BasePred','MultiPred','GroupByModel','PredPRS','PRSCS2']


class BasePred(ABC):

    weight_filetypes = ['extprst.tsv','prscs.tsv']
    default_weightfile_type = 'prstweights.tsv'
    default_weight_cols = ['chrom','snp','pos','A1','A2','allele_weight']
    extra_weight_cols   = False
    default_sst_cols = ['SNP','A1','A2','BETA']
    dtype_pred = 'float32'
    _nancheck = True
    scaling = 'ref'
    shuffle = False
    _clear_cache=True
    _close_pbar = True
    _default_n_jobs=4
    _allow_missing=True
    algo_pred = 'i8fast'
    _display_info = True
    
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
         prst downloadutil --pattern example --destdir ./; cd example  {insert}                                       # Makes \'example\' dir in current path.
         prstools {cmdname} --ref {ldrefname} --target target --sst sumstats.tsv {chromopt}--n_gwas 2565 --out ./result-{cmdname} # Run the model with example data.
         prst {cmdname} -r {ldrefname} -t target -s sumstats.tsv -n 2565 {chromopt}-o ./result-{cmdname}                          # A shorter version of previous.
        '''
        # prst {cmdname} -r {ldrefname} -t target -s sumstats.tsv -n 2565 {chromopt}-o ./result-{cmdname} --pred # A shorter version of previous that also does the predictions.
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
    def _get_cli_spkwg(cls, basic_pkwargs=True):
        from prstools.utils import retrieve_pkwargs
        from textwrap import dedent
        
        if basic_pkwargs and type(basic_pkwargs) is bool:
            basic_pkwargs = dict(
                basics=dict(args=['-h','--help'], kwargs=dict(action='help', help='Show this help message and exit.')),
                cpus=dict(args=['--cpus','-c'], kwargs=dict(metavar='<number-of-cpus>', default=prst._cmd.get_default_cpus(), type=int, 
                                                help='The number of CPUs to use (1–5 is generally most efficient). It is generally best to first maxout --n_jobs before increasing the number of '
                                                            'CPUs above 1. To disable set to -1.')))

#         intcaster = lambda x: int(float(x))
#         def str_or_none(x): return None if x.lower() == "none" else x
#         intcaster = int
#         str_or_none = str
        def str_or_none(arg): return str(arg)
        def intcaster(x): return int(float(x))
#          \033[34mhttps://tinyurl.com/sstxampl\033[0m 
        data_pkwargs = dict(
        ref=dict(args=['--ref','-r'], kwargs=dict(required=True, metavar='<dir/refcode>', 
                help="Path to the directory that contains the LD reference panel. You can download this reference data manually too "
                "with '{basecmd} downloadutil'. Current config relevant to this functionality (prstdatadir: {prstdatadir}, auto_download: {auto_download}).")),
        target=dict(args=['--target','-t'], kwargs=dict(required=True, type=str_or_none, metavar='<bim-prefix>', 
                help="Specify the bim file or its prefix of desired target dataset. The bim file should be in plink format. You can also set "
                "the target to 'none' in which case the variant set will not be filtered for variants present in the target set.")),
        sst=dict(args=['--sst','-s'], kwargs=dict(required=True, metavar='<file>', 
                help="The summary statistics file from which the model will be created. The file should contain columns: SNP, A1, A2, BETA or OR, P or SE information. "
                "At the moment, the file is assumed to be tab-seperated, if you like other formats please let devs know. Alternative column names can be specified with "
                "--colmap (more info below). SNP column should contain rsid's, but now these can be filled from See \033[34mhttps://tinyurl.com/sstxampl\033[0m for a sumstat example.")),
        out=dict(args=['--out','-o'], kwargs=dict(required=True, metavar='<dir+prefix>', 
                help="Output prefix for the results (variant weights). This should be a combination of the desired output dir + file prefix.")),
        n_gwas=dict(args=['--n_gwas','-n'], kwargs=dict(required=False, type=intcaster, metavar='<num>', default=None, 
                help="Sample size of the GWAS. Not required if sumstat has a 'N' column and overrules column data if specified.")),
        chrom=dict(args=['--chrom'], kwargs=dict(required=False, type=str, metavar='<chroms>', default='all', 
                help="Optional: Select specific chromosome to work with. You can specify a specific chromosome as e.g. \"--chrom 3\". All chromosomes are used by default.")),
        colmap=dict(args=['--colmap'], kwargs=dict(type=str, metavar='<colnames>', default='{default_colmap}', 
                help="Optional: Allows one to specify an alterative column name for the internally used columns snp,A1,A2,beta,or,pval,se_beta,n_eff,af_A1,, "
                "(in that order). Forinstance \"--colmap rsid,a1,a2,beta_gwas,,pvalue,beta_standard_error,,,,\" (OR, N, FRQA1, are excluded in this example). "
                "When the command is run a quick this_column -> that_column conversion table will be shown. Additionaly prstools has many internal checks to make "
                "sure a good PRS will be generated! The original default colmap works with the PRS-CS standard sumstat formatting.")),
        pred=dict(args=['--pred','-p'], kwargs=dict(required=False, metavar='<yes/no>', type=str, default='auto', 
                help="Optional: Add this argument to set behavior for PRS generation for the induviduals in the target dataset. "
                "With the 'auto' option (which is the default) the tool tries to generate a prediction unless the --chrom option is set. Available options: (yes/no/auto)."))
        )

        from textwrap import dedent
        doc = dedent(cls.__doc__)
        epilog = getattr(cls, '_get_cli_epilog', lambda: None)()
        display_info = getattr(cls, '_display_info', False)
        groups = dict(
            general=dict(grpheader='General Options', pkwargs={**basic_pkwargs}),
            data   =dict(grpheader='Data Arguments',  pkwargs=data_pkwargs),
            model  =dict(grpheader='Model Arguments (all optional)', pkwargs=retrieve_pkwargs(cls))
        )
        spkwg = dict(
            cmdname     = cls.__name__.lower(), #command name
            clsname     = cls.__name__, #class name
            description = doc, # order of desc and help good here, dont check again.
            display_info = display_info,
            help        = doc.split('\n')[0],
            epilog      = epilog,
            modulename  = cls.__module__,
            groups      = groups,
            subtype     = 'PRSTCLI'
        )
        
#         @classmethod, ignore comment
#         def _get_cli_spkwg(cls, basic_pkwargs=True): ## This badboi wraps the super method to enhance it.
#             nargskeys = ['input','selectcols','sortcols','assertunique','antiglobs']
#             reqkeys = ['input','out']
#             spkwg = super()._get_cli_spkwg(basic_pkwargs=basic_pkwargs) 
#             for key in nargskeys: spkwg['groups']['general']['pkwargs'][key]['kwargs'].update(nargs='+')
#             for key in reqkeys: spkwg['groups']['general']['pkwargs'][key]['kwargs'].update(required=True)
#             return spkwg
 
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
    
    @staticmethod
    def _get_pkwargs_for_class(cls): # fyi, making a classmethod here, gave issue since. 
        # it would get BasePred as class, did not find other fix
        # pkwargs are the kwargs and the defaults of the prstools cli.
        from prstools._parser_vars import get_subparserkwg_lst
        for elem in get_subparserkwg_lst():
            if elem['clsname'] == cls.__name__:
                pkwargs = elem['pkwargs']; break
        return pkwargs
    
    @classmethod
    def from_cli_params_and_run(cls, *, ref, target, sst, n_gwas=None, chrom='all', fnfmt='_.{ftype}', ftype='prstweights.tsv', groupbydefault=False,
                                verbose=True, pkwargs=None, out=None, return_models=True, fit=True, pop=None, colmap=None, pred='auto', regdef=None, 
                                command=None, **kwargs):
        try: from prstools.linkage import AutoLinkageData
        except: from prstools.linkage import RefLinkageData as AutoLinkageData
        
        # Initialize model object(s) (multiple since hyperparam ranges, and maybe chroms):
        if pkwargs is None: pkwargs = cls._get_pkwargs_for_class(cls)
        def testkey(key): return (key in pkwargs) if pkwargs else True # The verbose in the next line overwrites the verbose in the 'kwargs' dict 
        groupby = kwargs.get('groupby', pkwargs.get('groupby',{}).get('kwargs',{}).get('default', groupbydefault))
        model = cls.from_params(**dict({key: item for key, item in kwargs.items() if testkey(key)}, verbose=verbose, groupby=groupby))
        
        ## Loop through different models, likely a parameter grid:
        #for model in models: # not sure about this atm
        
        # Gen output file name format and do quick check if output file can be saved before a lot of work is done:
        if out: out_fnfmt = model.create_output_fnfmt(**locals()); prstlogs=prst.utils.get_prstlogs()

        # Initialize data objects, fit the model & predict:
        linkdata = AutoLinkageData.from_cli_params(ref=ref, target=target, sst=sst, 
                        n_gwas=n_gwas, chrom=chrom, colmap=colmap, pop=pop, verbose=verbose, regdef=regdef, out_fnfmt=out_fnfmt, **kwargs)
        model.set_linkdata(linkdata)
        if fit: model.fit()
        if out: model._save_results(out_fnfmt, out=out, ftype=ftype) # Store fitting result, most often this will be the weights.
        prstlogs['times']['methodstop'] = pd.Timestamp.now() # Save model endtime
        
        if pred == 'auto' and not hasattr(model, 'weights_df'): pred = 'no'
        if pred == 'auto' and chrom != 'all': pred='no'
        if pred and pred != 'no': # Prediction
            #model.remove_linkdata(); linkdata.clear_linkage_allregions # seems to do pretty much nothing.. anyway xp was 5% mem, which jumped to 20 and 60 later
            try: 
                bed = prst.io.load_bed(target, verbose=verbose)
                yhat = model.predict(bed); 
                prst.io.save_prs(yhat, fn=out_fnfmt, verbose=verbose) # Store prediction result
            except Exception as e:
                msg = (f"Could not generate prediction (e.g. plink file missing)" 
                       f" so since --pred='auto' the prediction step will be skipped (target={target})")
                if pred == 'auto': print(msg)
                else: raise e

        if return_models: 
            return model
    
    def _save_results(self, fn, *, out, ftype):
        res = self.save_weights(fn, ftype=ftype)
        return res
    
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
        #import uuid; tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"  # unique temp file name
        if self.verbose: print(f'Saving model weights (filetype={ftype}) to: {fn}', end=' ', flush=True)
        if nancheck: assert np.sum(self.get_weights()['allele_weight'].isna().sum()) == 0
        fin_df = self.get_weights()[selcols]
        #if ftype.endswith('prstweights.h5'): to_file=pd.DataFrame(fin_df.to_numpy(), index=fin_df.index, columns=fin_df.columns).to_hdf; kwg=dict(key='df')
        if ftype.endswith('prstweights.h5'): to_file=fin_df.astype(object).infer_objects().to_hdf; kwg=dict(key='df')
        elif ftype.endswith('prstweights.parquet'): to_file=fin_df.to_parquet; kwg={}
        else: to_file=fin_df.to_csv; kwg=dict(sep='\t', index=False, header=header)
        prst.io._pd_to_atomizer(fn=fn, to_file=to_file, **kwg)
        #os.replace(tmp_fn, fn) # atomically move into place
        if self.verbose: print(f'-> Done', end=end, flush=True)
        if return_weights: return df
        
    def save_sst(self, fn, return_sst=False, ftype='tsv', basecols=None, addicols=None, nancheck=False):
        kwg = {}
        if addicols is not None: kwg['addicols'] = addicols
        if basecols is not None: kwg['basecols'] = basecols
        out_df = prst.io.save_sst(sst_df=self.sst_df, fn=fn, return_sst=return_sst, ftype='tsv', 
                    nancheck=nancheck, verbose=self.verbose, **kwg)
        if return_sst: return out_df
    
    def get_pbar(self, iterator, *, make_range_var=True, **kwg):
        # Maybe import funny wrapper class
        msg = f"Object implement iterator has no length. This is required. More info: type(iterator)={type(iterator)} , iterator={iterator}"
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
        
    def remove_cache(self):
        if hasattr(self,'cache_dt'):
            del self.cache_dt
        
    def remove_linkdata(self, clean_mem=True):
        self.linkdata.clear_linkage_allregions()
        del self._linkdata
        for attr in ['_linkdata', 'cache_dt']:
            try: delattr(self, attr)
            except Exception: pass
        if clean_mem: prst.utils.clear_memory(postsleep=1.)
        
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
    
    def predict(self, bed, *, n_inchunk=1000, groupby=None, validate=True, dtype=None, algo=None,
                localdump=False, weight_type='allele', trait_df=None, colour='#7f00ff'): # <-- The more esotheric stuff on this line
        
        if 'pysnptools' in str(type(bed)):
            srd = bed; del bed
            return self.srdpredict(**locals())
        if type(bed) is str:
            if self.verbose: print('Input to predict() is a string and perhaps a filepath. Trying to load..')
            bed = prst.io.load_bed(bed,verbose=self.verbose)
        
        prstlogs = prst.utils.get_prstlogs()
        tic, toc = prstlogs.get_tictoc()
        assert weight_type in ('allele','standardized');
        if algo != 'ori': assert weight_type == 'allele'
        if trait_df is not None: assert localdump
        weights_df = self.get_weights()
        if len(weights_df['allele_weight'].shape) == 1: n_traits = 1
        else: n_traits = weights_df['allele_weight'].shape[1]  
        if self.verbose: print(f'Predicting {n_traits} phenotype(s) i.e. generating PRS, in chucks of {n_inchunk} snps. ', flush=True, end='')
        dtype = self.dtype_pred if dtype is None else dtype
        algo  = self.algo_pred if algo is None else algo
        msg = ''
        
        if validate:
            bim_df = bed.bim_df.copy() # For the next line did it the other way around for mem-footprint.
            weights_df = weights_df.drop(columns='xidx', errors='ignore') # Make sure it does not have xidx by chance.
            p_pre = weights_df.shape[0]
            weights_df = prst.merge_snps(weights_df, bim_df, req_all_right=False, handle_missing='filter', flipcols=[]) ## HEY WAIT... WHAT!!!
            #prst.utils.get_ip().embed()
            #weights_df['allele_weight']=weights_df['allele_weight']*weights_df['rflip'] # 20TB crash..
            #weights_df['allele_weight']=weights_df[['allele_weight']]*weights_df[['rflip']] # nans
            weights_df['allele_weight']=weights_df['allele_weight'].mul(weights_df['rflip'], axis=0) # This is a bit of a hack
            weights_df['xidx'] = weights_df['xidx'].astype('int64', errors='ignore')
            assert weights_df['xidx'].dtype == 'int64', 'xidx contained a nan, this should not happen for regular users, contact dev'
            weights_df = weights_df.sort_values('xidx')
            p_post = weights_df.shape[0]; n_missing = p_pre - p_post; perc = (n_missing/p_pre) * 100.
            inject = ', which is acceptable depending on use-case (<10%)' if perc < 10 else ''
            if n_missing > 0: msg += f'\nMissing {n_missing:,} variants ({perc:.0f}%) in the target that are in the weights{inject}.'
            if n_missing > 0 and not self._allow_missing: raise RuntimeError(msg)
        if self.verbose: print(msg)
        
        # Loop through Genome:
        yhat_dt = dict(); sst_dt = {}; X=None
        if len(weights_df['allele_weight'].shape) == 1: n_traits = 1
        else: n_traits = weights_df['allele_weight'].shape[1]  
        #for itr in self.get_iterator(range(n_iter), pbar=self.pbar)
        for grp, wgrp_df in self.get_iterator(weights_df.groupby(groupby), pbar=self.pbar, colour=colour) if groupby is not None else [(None, weights_df)]:
            yhat = np.zeros((bed.iid_count, n_traits)); sst_lst = []
            inner_pbar = self.pbar if grp is None else None
            for start in self.get_iterator(range(0, wgrp_df.shape[0], n_inchunk), pbar=inner_pbar, colour=colour):
                wchunk_df = wgrp_df.iloc[start:start+n_inchunk]
                cxidx = wchunk_df['xidx']
                if algo == 'ori':
                    X = bed.read(index=np.s_[:,cxidx], dtype=dtype)
                    m = np.nanmean(X, axis=0)
                    idx = np.where(np.isnan(X))
                    s = np.nanstd(X, axis=0) if weight_type == 'standardized' else None
                    X[idx] = np.take(m, idx[1])
                elif algo == 'i8fast':
                    if X is None: X = bed.read(index=np.s_[:,cxidx], dtype=dtype)
                    nsamp,_= X.shape
                    X8 = bed.read(index=np.s_[:,cxidx], dtype='int8')
                    mask = (X8 == -127)
                    m = X8.sum(axis=0).astype(dtype)
                    msksum = mask.sum(axis=0).astype(dtype)
                    m += msksum*127
                    m /= (nsamp - msksum)
                    try:
                        assert X.shape == X8.shape
                        np.copyto(X, X8)
                    except: X = X8.astype(dtype)
                    for j in range(X.shape[1]): # <--- This one is faster!
                        X[mask[:, j], j] = m[j]
                    # for i in range(X.shape[0]):
                    #     X[i, mask[i,:]] = m[mask[i,:]]
                else: raise ValueError(f"Unknown algorithm: {algo}")
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
                localdump=False, weight_type='allele', trait_df=None, colour=None, dtype=None): # <-- The more esotheric stuff on this line
        
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
            from prstools.io import _load_bimfam_from_srd
            bim_df, fam_df = _load_bimfam_from_srd(srd, skipifpresent=True)
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
    
class BaseMulti(): ## This is a base class so should Not generate objects i.e. instances.
    
    @classmethod
    def from_weights(cls, weights, sort=True, verbose=False, **kwg):
        if verbose: print('Initializing multi-weight model using input weights. ',flush=True, end='')
        if not isinstance(weights, pd.DataFrame): # can this be done with decorator?
            raise TypeError("Input must be a DataFrame.")
        if type(weights) is dict: raise NotImplementedError()
        # A check for the required columns is also needed, somewhere.
        model = cls(**kwg, verbose=verbose)
        model._set_weights(weights, sort=sort)
        if verbose: print('-> Done')
        return model
    
    @classmethod
    def from_dict(cls, weights_dt, ref_df=None, verbose=False, greedy=False, on=None, remove_allnan=False, **kwg):
        
        assert not greedy, 'Greedy options not implemented yet.'
        assert len(weights_dt) > 0, 'weights_dt is empty'
        if ref_df is not None: allweights_df=ref_df.copy()
        elif len(weights_dt) == 1: return cls.from_weights(list(weights_dt.values())[0].copy(), verbose=verbose, **kwg)
        else: raise NotImplementedError(f'... contact dev if this option (ref_df is None) is desired.')
        on_dt = dict() if on is None else dict(on=on)
        if len(on_dt) == 0: on=['snp','A1','A2']
        assert allweights_df[on].isna().sum().sum() == 0, f'cannot have Nans in starter columns on={on}'
        msg = 'Duplicated snp present, this is probably cause of multi-allelic snps, this needs to be implemented contact dev'
        cnt = allweights_df['snp'].duplicated().sum(); msg+=f'(dupcount={cnt:,})'
        if cnt != 0: warnings.warn(msg)
        if len(weights_dt)==1:
            if verbose: print(f'{"Combining":<12}Since number of input weights is 1, we do not have to combine weight files.')
            return cls.from_weights(list(weights_dt.values())[0], verbose=verbose, **kwg)
        
        # merging mechanics:
        pre_self = cls(**kwg, verbose=verbose) ## <---- self is init twice, see cls.from*() line below
        for wname, curweights_df in pre_self.get_iterator(weights_dt.items(), pbar=pre_self.pbar):
            pre_self.pbar.set_description(f"{'Combining':<12}")
            curweights_df = curweights_df.copy() # This line is crucial for the PandasMimic type used to load from disk.
            cmissing = set(cls.default_weight_cols) - set(curweights_df.columns)
            allweights_df = prst.merge_snps(allweights_df, curweights_df, flipcols=['allele_weight'], how='left', handle_missing='keep', **on_dt)
            allweights_df = allweights_df.rename(columns=dict(allele_weight=f"allele_weight_{wname}")).drop('rflip',axis=1)

        #Create multi-index columns:
        lst = [(col, '') if not col.startswith('allele_weight') else ('allele_weight', col.split('allele_weight_')[-1]) for col in allweights_df.columns]
        newcols = pd.MultiIndex.from_tuples(lst)
        allweights_df.columns = newcols
        assert (allweights_df['allele_weight'] == 0).sum().sum() == 0
        if remove_allnan:
            ind = allweights_df['allele_weight'].isna().sum(axis=1) < allweights_df['allele_weight'].shape[1]
            allweights_df = allweights_df[ind]
            
        # Some postprocessing related to nans (dropping nanonly rows)
        ind = ~allweights_df['allele_weight'].isna().all(axis=1)
        allweights_df = allweights_df[ind].reset_index(drop=True)
        allweights_df = allweights_df.fillna(0)
        
        model = cls.from_weights(allweights_df, verbose=verbose, **kwg)
        if verbose: ''
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
        if type(ref_df) is str: ref_df, _ = prst.load_bimfam(ref_df, fam=False, verbose=verbose)
        class PandasMimic():
            def __init__(self,fn, ftype):
                # potentially we can do a rapid check here for the weight loading.
                assert type(fn) is str and type(ftype) is str
                self.fn = fn; self.ftype=ftype; self.pyarrow=pyarrow
                self.sep=sep
            def copy(self):
                return prst.io.load_weights(**vars(self))
        weights_dt=dict()
        for fn in fn_lst:
            weights_dt[fn]=PandasMimic(fn=fn, ftype=ftype)
            
        model = cls.from_dict(weights_dt, ref_df=ref_df, on=on, remove_allnan=remove_allnan, verbose=verbose, **kwg)

        return model

    
class GroupByModel(BaseMulti, BasePred):
    
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
        from joblib import Parallel, delayed
        
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
        #print(1); prst.utils.get_memory_usage() # remove me later
        contents = {key: item for key, item in linkdata.groupby(self.groupby, sort=True, skipempty=True)}
        del linkdata; self.remove_linkdata()
        #print(2); prst.utils.get_memory_usage() # remove me later
        # MultiProcessing portion:
        #with Manager() as manager:
        fakebar = prst.utils.FakeMultiprocPbar(manager=None) if self.pbar else False
        real_pbar = self.get_pbar(iterator=range(tot_iters), fakebar=fakebar, deamon=True) if self.pbar else False
        if self.pbar: # Do fakebar hacks to make everything work
            mgr = fakebar._mgr; fakebar._mgr=None
        with warnings.catch_warnings(): # Bit annoying this catch warning is needed... but it is, else freaky warnings for my users
            warnings.filterwarnings("ignore", category=UserWarning,
                message=r".*worker stopped while some jobs were given to the executor.*")
            prst.utils.clear_memory()
            #print(5,'before doing parallel -- clear'); prst.utils.get_memory_usage() 
            results = Parallel(n_jobs=self.n_jobs, max_nbytes=None, mmap_mode=None)(delayed(worker)(self.get_model_clone(), cur_linkdata, fakebar, grp) for grp, cur_linkdata in contents.items())
        if self.pbar: 
            real_pbar.close(); real_pbar=None
            prst.utils.clear_memory(); # Crucial line because gc.collect() inside, else things go wrong later.
            fakebar.close(); mgr.shutdown(); mgr=None
        for grp, model in results: self.model_dt[grp] = model
        self.combine_set_weights()
        return self

    def combine_set_weights(self):
        assert hasattr(self,'model_dt'), f'No models present, so cannot create a working weights set for {self}.'
        weights_df = pd.concat([model.get_weights() for grp, model in self.model_dt.items()], axis=0) #for grp, model in self.model_dt.items():
        self._set_weights(weights_df, silentsort=True)
    
class MultiPRS(BaseMulti, BasePred, PRSTCLI):
    """\
    MultiPRS: It generates polygenic risk scores if you give it weights (
    Note: currently one needs to run "prst config" first.
    """
    # PredPRS Does not exist anymore
    
    _multiprs_mode_helpmsg =  ("Different modes for file handling. If 'strict' is present an error in forinstance "
        "the loading of a weight will lead to a stop. ")
    
    def __init__(self, *,
        mode:Param(_multiprs_mode_helpmsg, str)='strict',     
        ref='1kg_hm3',
        # A reference to which all will be joined. <-DeveloperNotes for later (next line must be empty to prevent this note to end up in cli)
        # Determines how missing variants are handled.
                 
        missing='flex',
        val=None,
        groupby:str=None,
        pbar:bool=True,
        verbose:bool=False
        ):
        # Input sumstats. Can be specified with * operator. For example "--sst CoolCohort_trait=*_hapmap3.tsv"
        # Stuff all the args into fields.
        _excl_lst = ['self', 'kwg_dt']
        kwg_dt = {key: item for key, item in locals().items() if not (key in _excl_lst)}
        for key, item in locals().items():
            if not (key in _excl_lst): 
                self.__setattr__(key, item)
        self._kwg_dt = copy.deepcopy(kwg_dt)
    
    @classmethod
    def _get_cli_epilog(cls, commentccode='32'):
        # Put an example here when u ready.
        return None

    @classmethod
    def _get_cli_spkwg(cls, basic_pkwargs=True): ## This badboi wraps the super method to enhance it.
        spkwg = super()._get_cli_spkwg(basic_pkwargs=basic_pkwargs)
        weights = {'args': ['--weights', '-w'],
             'kwargs': dict(help='The weight file(s) for the PRS. Can be specified with * operator. For example "-w greatprsweights_trait=*.tsv".  ergerg',
              required=True, metavar='<file>', type=str, nargs='+')} 
        drops = ['ref','sst','n_gwas','chrom','colmap']
        for drop in drops:
            if drop in spkwg['groups']['data']['pkwargs']:
                spkwg['groups']['data']['pkwargs'].pop(drop)
            if drop in spkwg['groups']['model']['pkwargs']:
                spkwg['groups']['model']['pkwargs'].pop(drop)
        pkwargs = spkwg['groups']['data']['pkwargs']
        pkwargs = dict(weights=weights, **pkwargs)
        pkwargs['out']['kwargs']['help'] += 'Suggestion: "--out {trimweights}_{target}".' # <- Additions
        spkwg['groups']['data']['pkwargs'] = pkwargs #put back
        return spkwg
    
    @classmethod
    def from_cli_params_and_run(cls,
            weights=None,
            target=None, # mind this stuff here now does not populate the cli arguments
            out:str=None, # Output file name. A good suggestion: "--out {trimsst}.prstsst.tsv" 
            pred='auto',
            mode=None,     
            ref='1kg_hm3',
            fnfmt='.{ftype}',
            prs_ftype='prstprs.tsv',
            weights_ftype='prstweights.parquet',
            return_models=True,
            pkwargs=None,
            verbose=True,            
            **kwg # this kwg catches command and func for a smooth run
            ):
        
        # Initialize model object(s) (multiple since hyperparam ranges, and maybe chroms):
        assert type(ref) is str, f'variable "ref" should be string! it is {type(ref)}'
        if pkwargs is None: pkwargs = cls._get_pkwargs_for_class(cls)
        def testkey(key): return (key in pkwargs) if pkwargs else True # The verbose in the next line overwrites the verbose in the 'kwargs' dict 
        #groupby = kwargs.get('groupby', pkwargs.get('groupby',{}).get('kwargs',{}).get('default', groupbydefault))
        # mind that in current version,
        groupby = False ## going for this untill its more clear
        modelkwg = dict({key: item for key, item in kwg.items() if testkey(key)}, verbose=verbose, groupby=groupby)
        self = cls.from_params(**modelkwg)        
        ref = prst.utils.validate_path(ref=ref,  handle_prstdatadir='allow')
        if os.path.isdir(ref):
            ref_lst = glob.glob(os.path.join(ref, '*.bim'))
            msg = f'Found multiple or no .bim files in target {ref}. There should be one. perhaps --ref input was not proper?'
            assert len(ref_lst) == 1, msg
            ref = ref_lst[0]
        ref_df, _ = prst.load_bimfam(ref, fam=False)
        mode = self.mode
        
        # Prepping & validating input & output filenames & load target bim if needed:
        if type(weights) is str: fn_lst = [weights]
        elif type(weights) is list: fn_lst = weights
        else: ValueError(f'Input weights i.e. --weights are of the type {type(weights)}. This is unexpected and the following code is not designed for it.')   
        new_lst=[]
        for fn in fn_lst:
            if any(elem in fn for elem in ['*','[']): new_lst += glob.glob(fn)
            else: new_lst += [fn]
        fn_lst = new_lst
        uniqfns  = np.unique(fn_lst)
        nuniqfns = len(uniqfns)
        fn_lst = uniqfns
        out_lst = []
        for fn in fn_lst: # Quickly check files exist:
            if not os.path.isfile(fn): open(fn)
        for fn in fn_lst:
            weights = os.path.basename(fn); ftype=prs_ftype
            trimweights = prst.io.get_fn_trimmed(weights, make_bn=True)
            if out: out_fnfmt = cls.create_output_fnfmt(**locals())
            out_lst += [out_fnfmt]
        nuniqouts = len(np.unique([os.path.basename(out_fn) for out_fn in out_lst]))
        msg = f'\nGot {nuniqfns} unique input filenames and generated {nuniqouts} unique output filename(s).'
        if (nuniqouts) != 1: 
            msg += ('\nBecause there is not 1 unique output we cannot proceed, '
                'consider changing/removing bracket ({}), if you have them in --out.'); raise ValueError(msg)
        #else:  msg += '\nBecause there is not 1 output for every 1 input we cannot proceed! Set e.g. "--out {trimweights}.prs.tsv" to fix.' 
        #msg += 'Will be combining the weights and dont worry.. before the prediction is starting a combined version will be stored which can be reloaded quickly.'
        if verbose: print(msg+'\n')
        target_df, _ = prst.load_bimfam(target, fam=False, start_string = 'Loading target file.', verbose=False) if target else (None,None) # Load the target to make sure it work
            
        # Loop through files:
        if verbose: print('Loading & Combining weights:')
        weights_dt = {}
        pbar = prst.utils.get_pbar(list(zip(fn_lst,out_lst))) 
        pbar.set_description(f"{'Loading':<12}")
        lasterr=False; first = True
        for fn, out_fn in pbar:
            bn = os.path.basename(fn)
            pbar.set_postfix(bn=bn)
            try:
                #print(f'\nIt seems {os.path.basename(out_fn)} does not exist (or a redo was requested) so we are making it from {fn}')
                weights_df = prst.load_weights(fn)
                weights_dt[bn] = weights_df
            except Exception as e:
                lasterr = e
                if 'strict' in self.mode:
                    raise e from Exception('Set --mode flex to have the procedure skip input files that give issues.')
        if lasterr != False:
            print('This is the last error from the bunch:')
            if not 'flex' in mode:
                raise e
            else: print(lasterr)
        
        # Combine the weights into one frame:    
        model = cls.from_dict(weights_dt, **modelkwg, ref_df=ref_df)
        
        if len(weights_dt) > 1: 
            if prst.io.get_pyarrowinstalled_bool():
                model._save_results(out_fnfmt, out=out, ftype=weights_ftype)
            else: 
                if verbose: print('Skipping multi-weights saving since pyarrow is not installed\n')
        elif verbose: print("Not re-saving weights since there was only 1 input weight.\n")
        
        if pred and pred != 'no': # Prediction
            try: 
                bed = prst.io.load_bed(target, verbose=verbose)
                yhat = model.predict(bed); 
                prst.io.save_prs(yhat, fn=out_fnfmt, verbose=verbose, ftype=prs_ftype) # Store prediction result
            except Exception as e:
                msg = (f"Could not generate prediction (e.g. plink file missing)" 
                       f" so since --pred='auto' the prediction step will be skipped (target={target})")
                if pred == 'auto': print(msg)
                else: raise e
                    
        if return_models: 
            return model
    
try:
    profile
except NameError:
    def profile(func):
        return func

class PRSCS2(BasePred):
    
    "PRS-CS v2: A polygenic prediction method that infers posterior SNP effect sizes under continuous shrinkage (CS) priors."
    _gig = None
    _default_sampler='rue'
    
    def __init__(self, *,
         n_iter=1000,              # Total number of MCMC iterations.
         n_burnin=0.5,             # Number of burn-in iterations if larger than 1 or fraction of n_iter if smaller than 1.
         n_slice=1,                # Thinning of the Markov chain.
         shuffle=False,
         seed=-1,                  # Random seed for reproducibility.
         a=1.0,                    # Parameter a in the gamma-gamma prior.
         b=0.5,                    # Parameter b in the gamma-gamma prior. 
         phi=-1.,                  # Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a Bayesian approach.
         clip=1.,                  # Clip parameter. The default works best in pretty much all cases.
         sampler='default',        # Sampler algorithm. The default is Rue sampling, which is the original sampler and gives good results.
         groupby:str='chrom',
         local_rm:bool=False,    
         compute_score:bool=False,
         clear_linkdata:bool=True,
         scaling='ref',
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
        if self.sampler == 'default': self.sampler = self._default_sampler
        #assert self.sampler in ['rue','bhat','sld']
        self.pop = self.pop.upper()

    
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
    
#     @profile 
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
        get_ipython().system('python ../prstools/_cmd.py --dev-secret')
        #!prst --dev | head -3
    print('Done')