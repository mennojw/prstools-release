#| export
import os, json, time
from prstools._ext_utils import *
import warnings, threading
from abc import ABC, abstractmethod
try:
    from fastcore.script import call_parse, Param
except:
    def Param(**kwg):
        return None
try: import IPython as ip
except: pass
#     Param = fun
def optional_import(path, name=None, default=None):
    try:
        mod = __import__(path, fromlist=[''])
        return getattr(mod, name or path.split('.')[-1]) if name else mod
    except ImportError:
        return default

def plot_manhattan(data_df, x=None, y='-logp', regcol='chrom', pvalmin=1e-323, palette='bright', aspect=4, s=6., **snskwg):
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from prstools.loaders import compute_pvalbetase
    data_df = data_df.reset_index(drop=True).reset_index()
    if not '-logp' in data_df.columns:
        if not 'pval' in data_df.columns:
            data_df = compute_pvalbetase(data_df, calc_lst=['pval'], pvalmin=pvalmin)
        data_df['-logp'] = -np.log10(data_df['pval'])
    if x is None: x = 'index'
    plot = sns.relplot(data=data_df, x=x, y=y, aspect=aspect, hue=regcol, palette=palette, legend=None, s=s, **snskwg)
    chrom_df=data_df.groupby(regcol)[x].median()
    plot.ax.set_xlabel(regcol); plot.ax.set_xticks(chrom_df); plot.ax.set_xticklabels(chrom_df.index)
    plot.fig.suptitle('Manhattan plot')
    plt.show()
    
manhattan_plot = plot_manhattan
    
def validate_path(*args, exists=False):
    newargs=[]
    for arg in args:
        arg = os.path.expanduser(arg)
        if exists and not os.path.exists(arg):
            open(
                arg)
            #import errno
            #raise FileNotFoundError( errno.ENOENT, os.strerror(errno.ENOENT), arg)
            #raise FileNotFoundError(f'No such file or directory: \'{arg}\'')
        newargs += [arg]
    return newargs[0] if len(newargs) == 1 else tuple(newargs)

def get_tqdm():
    try:
        from tqdm.cli import tqdm
    except:
        tqdm = lambda x:x
    return tqdm

def get_timestring_from_td(td):
    s = td.total_seconds()
    s = int(s)
    h,s=divmod(s,3600)
    m,s=divmod(s,60)
    return f'Completed in {h}h, {m}m and {s}s'

# class FakeMultiprocPbar:
#     """
#     Manager-backed, multiprocess-safe lightweight progress counter.
#     """
#     def __init__(self):
#         # Create the manager immediately so 'with' is optional
#         from multiprocessing import Manager
#         #self._manager = Manager()
#         manager = Manager()
#         self._val = manager.Value('i', 0)
#         self._lock = manager.Lock()
# #         self._val = self._manager.Value('i', 0)
# #         self._lock = self._manager.Lock()

#     # Context-manager hooks
#     def __enter__(self):
#         return self

#     def __exit__(self, exc_type, exc, tb):
#         self.close()           # ensure manager shuts down

#     # Same API as before
#     def update(self, n=1):
#         with self._lock:
#             self._val.value += n

#     @property
#     def n(self):
#         with self._lock:
#             return int(self._val.value)

#     def close(self):
#         """Release the Manager process (safe to call more than once)."""
#         1
# #         if self._manager is not None:
# #             self._manager.shutdown()
#             self._manager = self._val = self._lock = None

class FakeMultiprocPbar:
    def __init__(self, manager=None):
        from multiprocessing import Manager
        if manager is None:
            mgr = Manager()
            self._mgr = mgr # Internal management, which does not seem to work..
        else: 
            mgr = manager
            self._mgr = None
        self._val  = mgr.Value('i', 0)
        self._lock = mgr.Lock()
#         self._shutdown = mgr.shutdown      # plain function, not the Manager itself

    def update(self, n=1):
        with self._lock:
            self._val.value += n
            
    def __exit__(self, exc_type, exc, tb):
        self.close()           # ensure manager shuts down

    @property
    def n(self):
        with self._lock:
            return self._val.value

    def close(self):
        if self._mgr:
            self.mgr.shutdown()
            self._mgr = None


def get_pbar(iterator, ncols=200, colour='green', bar_format='{l_bar}{bar:70}{r_bar}',
             mininterval=0.5, desc=None, fakebar=None, deamon=False, **kwg):
    assert not 'total' in kwg, 'Need to use iterator, cannot use total argument.'
    create_pbar = get_tqdm()
    pbar = create_pbar(iterator, ncols=ncols, colour=colour,
          bar_format=bar_format,
          mininterval=mininterval, desc=desc, **kwg) 
    
    if deamon:
        try: sleeptime = float(deamon)
        except: sleeptime = 1.0
        tot_iters = len(iterator)
        stop_event = threading.Event()
        def do_update():
            if fakebar:
                delta = (fakebar.n - pbar.n)
            else: delta = -1
            if delta > 0: pbar.update(delta)
            else: pbar.refresh()
        def monitor():
            while not stop_event.is_set() and pbar.n < tot_iters:
                do_update()
                time.sleep(sleeptime)
        t = threading.Thread(target=monitor, daemon=True)
        t.start()
        pbar._reall_close = pbar.close
        def altclose(*args,**kwargs):
            stop_event.set()
            time.sleep(sleeptime)
            do_update()
            t.join(0.1)
            return pbar._reall_close(*args,**kwargs)
        pbar.close = altclose
    return pbar

def save_to_interactive(dct=None, maxframes=20):
    import sys
    if dct is None: raise NotImplementedError('automatic retrieval of vars not ready yet, you have to give an argument like'
            ' save_to_interactive(dict(loc_dt=locals()))')
    # Walk up the stack looking for '__name__'
    # with a value of '__main__' in frame globals
    for n in range(maxframes):
        cur_frame = sys._getframe(n)
        name = cur_frame.f_globals.get('__name__')
        if name == '__main__':
            # Yay - we're in the stack frame of the interactive interpreter!
            # So we update its frame globals with the dict containing our data
            cur_frame.f_globals.update(dct)
            break
            
## CLI mechanics functionality:
def process_subparserkwgs(subparserkwg_lst):
    from textwrap import dedent
    for i, spkwg in enumerate(subparserkwg_lst):
        if type(spkwg) is type({}): continue
        elif 'PRSTCLI' in str(getattr(spkwg,'__bases__','')):
            new_spkwg = spkwg._get_cli_spkwg()
            subparserkwg_lst[i] = new_spkwg
        elif 'BasePred' in str(getattr(spkwg,'__bases__','')):
            from .utils import retrieve_pkwargs
            doc = dedent(spkwg.__doc__)
            new_spkwg = dict(
                cmdname    = spkwg.__name__.lower(), #command name
                clsname    = spkwg.__name__, #class name
                description= doc,
                help       = doc.split('\n')[0],
                epilog     = spkwg._get_cli_epilog(),
                module     = spkwg.__module__,
                pkwargs    = retrieve_pkwargs(spkwg),
                subtype    = 'BasePred'
            )
            subparserkwg_lst[i] = new_spkwg
        else: 
            print(str(getattr(spkwg,'__bases__','')))
            raise NotImplementedError('Contact dev.') 
    return subparserkwg_lst
            
def store_argparse_dicts(subparserkwg_lst, show=False, sort_dicts=False, lst=None, store=True):
    from ._cmd import parse_args
    if not lst: lst = []
    import json
    from pprint import PrettyPrinter
    NoneType = type(None)
    subparserkwg_lst = parse_args(argv=[], subparserkwg_lst=subparserkwg_lst, return_spkwg=True)
    lst.append(subparserkwg_lst)
    proc_dt = dict()
    def fun(obj):
        if type(obj) is type:
            proc_dt[repr(obj)] = obj.__name__
            return repr(obj)
        else:raise Exception('This is mjw prst code, it seems a particular argument cannot be transformed into a argparse dict')
    if store:
        string = json.dumps(subparserkwg_lst, default=fun, indent=2)
        string = PrettyPrinter(indent=1,width=200,sort_dicts=sort_dicts).pformat(subparserkwg_lst)
        stringproc = string + ''
        for key, item in proc_dt.items():
            stringproc = stringproc.replace(key,item)
        test_lst = eval(stringproc)
        assert test_lst == subparserkwg_lst
        finstring = '# Warning: Dont edit here, This file was generated automatically, using a developer tool.'
        finstring += '\nNoneType = type(None)\nsubparserkwg_lst = '+stringproc
        with open('../prstools/_parser_vars.py', 'w') as f: 
            f.write(finstring)
        if show: print(finstring)
    return subparserkwg_lst


try:
    from fastcore.script import Param
except:
    class Param:
        def __init__(self, *args, **kwg):
            self.args =args
            self.all = kwg

        def __call__(self,*args, **kwg):
            return self
    
# from fastcore.docments import docments
# from fastcore.script import Param
def retrieve_pkwargs(cls):
    from fastcore.docments import docments
    from fastcore.script import Param
    import argparse
    
    def anno_parser2(func,  # Function to get arguments from
                    prog:str=None,  # The name of the program
                    return_pkwargs=False):
        #assert 'docments' in locals().keys()
        "Look at params (annotated with `Param`) in func and return an `ArgumentParser`"
        p = argparse.ArgumentParser(description=func.__doc__, prog=prog) #, formatter_class=_HelpFormatter) 
        pkwargs = {}
        for k,v in docments(func, full=True, returns=False, eval_str=True).items():
#             if k == 'verbose' and "XPRS" in str(cls): ergergger
            param = v.anno
            if not isinstance(param,Param): 
                param = Param(v.docment, v.anno, default=v.default)
                param.default = v['default']
            else:
                param = Param(param.help,param.type) #, default=v.default)
                param.default = v['default']
#                 print(k,v,v['default'])
#                 param.set_default(v.default)
#             param.default=None
#             param.set_default(v.default)
            
#             print(param.help)
            args = [f"{param.pre}{k}"]
            kwargs = param.kwargs
            pkwargs[k] = dict(args=args, kwargs=kwargs, v=v, param=param)
            p.add_argument(*args, **kwargs)
        p.add_argument(f"--pdb", help=argparse.SUPPRESS, action='store_true')
        p.add_argument(f"--xtra", help=argparse.SUPPRESS, type=str)
        if return_pkwargs:
            return p, pkwargs
        else:
            return p

    pkwargs = {}
    for argname, item in anno_parser2(cls, return_pkwargs=True)[1].items():
        if argname == 'kwg': continue
        kwargs = {}; v = item['v']; param=item['param']
        hvar = param.help
        if hvar == '': hvar=None
        kwargs['help'] = hvar
        if v.anno is type(False):
            ctype=v.anno # Fastcore's param does weird stuff with bool typing Param(type=bool).type-> None !
        else: ctype=param.type
        kwargs['type'] = ctype
        kwargs['default'] = param.default
        
        if kwargs['type'] in [str,int] and kwargs['default'] is None:
            kwargs['default'] = 'SUPPRESS'
        
#         akwargs = {};
#         akwargs['help'] = v['docment']  # old way of doing it
#         maybetype = v['anno']
#         akwargs['type'] = maybetype.type if not type(maybetype) is type else maybetype
#         akwargs['default'] = v['default']
        
        
#         if not (akwargs == kwargs):
#             print(kwargs)
#             print(akwargs)
#         if akwargs['default'] != kwargs['default']: jkergkj
#         if argname == 'verbose': 
#             if len(str(hvar)) > 6:
#                 print(cls)
#                 ergergerg

#         import IPython as ip
#         if 'sep' in argname: ip.embed(); ergegr
            
#         if argname == 'noheader': jkhergkjherk

#         if argname == 'verbose' and "XPRS" in str(cls):
#             print(kwargs)
#             print(akwargs)
#             ergergerg
            
        pkwargs[argname] = dict(args=['--'+argname], kwargs=kwargs)
        
    return pkwargs
            
class PRSTCLI(ABC):
    
    @classmethod
    @abstractmethod
    def _get_cli_spkwg(cls):
        pass
        """Generate CLI subparser kwargs."""
        
    @classmethod
    @abstractmethod
    def from_cli_params_and_run(self):
        """Method that runs the class from cli, can be access from the cli and from inside of python env."""
        pass
            
class AutoPRSTCLI(PRSTCLI):
    
    @classmethod
    def _get_cli_spkwg(cls, basic_pkwargs=True):
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
        return spkwg
        #print(cls.from_cli_params_and_run)
        
def extract_fn(url):
    import urllib
    return urllib.parse.urlparse(url).path.split('/')[-1]

def download_tar(url, dn='', desc_prefix=''):
    import urllib.request
    try:
        from tqdm.auto import tqdm
    except:
        tqdm = lambda x:x
    fn = extract_fn(url)
    fn = os.path.join(os.path.expanduser(dn), fn)
    with tqdm(total=0, unit='B', unit_scale=True, ncols=120, colour='green',
              bar_format='{l_bar}{bar:35}{r_bar}',
              mininterval=1, desc='{:<22}'.format(desc_prefix+os.path.basename(fn))) as pbar:
        def download_progress(count, block_size, total_size):
            if count == 0: pbar.total = total_size + 2;
            pbar.update(min(pbar.total-pbar.n,block_size))
        if not ((os.path.isfile(fn) or os.path.isdir(fn.replace('.tar.gz','')))):
            filename, headers = urllib.request.urlretrieve(url, fn+'.tmp', reporthook=download_progress)
            if os.path.isfile(fn+'.tmp'): os.rename(fn+'.tmp', fn)
        else:
            pbar.total = 1; pbar.update(1)
        
        
class DownloadUtil(AutoPRSTCLI): #, AutoPRSTSubparser):
    
    '''\
    Download and unpack LD reference panels and other data.
    The files that can be downloaded and unpacked with thirs command includes the standard reference files for PRS-CS and PRS-CSx.
    Additional information can be found at https://github.com/getian107/PRScsx
    '''
    
    @classmethod
    def from_cli_params_and_run(cls,
    destdir:str=None, # Directory in which all the data will be downloaded. Option required if you want the download to start.
    pattern:str='ALL',# A string pattern that retrieves every file that it matches. Matches everything by default. Without --destdir option (required to start downloading) one can see which files get matched.
    list=False, # Show which files can optionally be downloaded and exit.
    mkdir=False,
    command=None,
    func=None,
    keeptar=False, # Keep the tar.gz files in the destdir. If this option is not given they will be deleted automatically to save space.
    **kwg
                   ):

        # Defs & Inits:
        import tarfile, contextlib, pandas as pd
        try:
            from tqdm.auto import tqdm
        except:
            tqdm = lambda x:x
        data = [
            ["snpinfo_mult_1kg_hm3",  "https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=1", "1000G multi-ancestry SNP info (for PRS-CSx) (~106M)"],
            ["snpinfo_mult_ukbb_hm3", "https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=1", "UKBB multi-ancestry SNP info (for PRS-CSx) (~108M)"],
            ["ldblk_1kg_afr.tar.gz",  "https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=1", "1000G AFR Population LD panel (~4.44G)"],
            ["ldblk_1kg_amr.tar.gz",  "https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=1", "1000G AMR Population LD panel (~3.84G)"],
            ["ldblk_1kg_eas.tar.gz",  "https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=1", "1000G EAS Population LD panel (~4.33G)"],
            ["ldblk_1kg_eur.tar.gz",  "https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=1", "1000G EUR Population LD panel (~4.56G)"],
            ["ldblk_1kg_sas.tar.gz",  "https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=1", "1000G SAS Population LD panel (~5.60G)"],
            ["ldblk_ukbb_afr.tar.gz", "https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=1", "UKBB AFR Population LD panel (~4.93G)"],
            ["ldblk_ukbb_amr.tar.gz", "https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz?dl=1", "UKBB AMR Population LD panel (~4.10G)"],
            ["ldblk_ukbb_eas.tar.gz", "https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=1", "UKBB EAS Population LD panel (~5.80G)"],
            ["ldblk_ukbb_eur.tar.gz", "https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=1", "UKBB EUR Population LD panel (~6.25G)"],
            ["ldblk_ukbb_sas.tar.gz", "https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=1", "UKBB SAS Population LD panel (~7.37G)"],
            ["example.tar.gz", "https://www.dropbox.com/scl/fi/yi6lpbp0uhqiepayixvtj/example.tar.gz?rlkey=kvd7r17wuory9ucqdk4rh55jw&dl=1", "PRSTOOLS Example data (3.8M)"],
            ["g1000.tar.gz",'https://www.dropbox.com/scl/fi/97lsbtoomhti3q6x2wttf/g1000.tar.gz?rlkey=9hd85oytgnpv6wvbapvu2rk2m&st=3k4fq9ub&dl=1', "European 1kg plink dataset for hapmap3 (~64M)"]
            #["example.tar.gz","https://www.dropbox.com/scl/fi/7fg6c9e5dnmb0n4cdfquz/example.tar.gz?rlkey=31u2948paz539uw61jq37oe8s&dl=1", "PRSTOOLS Example data (70mb)"] 
        ]
        columns = ["filename", "url", "description"]
        links_df = pd.DataFrame(data, columns=columns)[['filename','description','url']]

        # Preprocessing:
        if list: # Overloading list is bad practice, I know.
            print('\nFiles available for downloading & unpacking (if following is difficult to read consider making terminal temporarily wider):\n')
            print(links_df.to_string(index=False, justify="left"))
            return True
        if not pattern == 'ALL':
            print('\nA pattern was used, which matched the following files:')
            ind = links_df['filename'].str.contains(pattern)
            links_df = links_df[ind]
            print(', '.join(links_df['filename'].to_list()))
        if not destdir: print('\n--destdir was not specified so download will not start.\n'); return True
        if not os.path.isdir(os.path.expanduser(destdir)):
            if mkdir: raise Exception(f'\nIt appears the supplied --destdir does not exist, please create: {destdir}')

        # Downloading:
        print(f'\nDownloading data, which might take some time. Data will be stored in: {destdir}')
        lst = []
        for idx, row in links_df.iterrows():
            fn = os.path.join(os.path.expanduser(destdir), row['filename']); dn=fn.replace('.tar.gz','')
            if os.path.isdir(dn): print(f'For {fn} the associated directory {dn} already exists,'
                                        ' therefore download & unpack is skipped. Remove directory for a redownload.'); lst+=[False]; continue
            download_tar(row['url'], dn=destdir); lst+=[True]
        links_df = links_df[lst]

        # Untarring:
        print('\nFinished downloading all LD data. Now we need to unpack all the tar.gz files (takes some time to start):') 
        def untar_file(archive_path, destination):
            with tarfile.open(archive_path, 'r:gz') as tar:
                file_names = tar.getnames()
                progress_bar = tqdm(total=len(file_names), 
                    ncols=120, colour='green', desc='Extracting')
                for file in tar:
                    tar.extract(file, destination)
                    progress_bar.update(1)
                    progress_bar.set_postfix(file=file.name)
                progress_bar.close()

        for idx, row in links_df.iterrows():
            curfn = row['filename']
            if 'tar.gz' in curfn:
                untar_file(os.path.join(destdir,curfn), destdir) 

        if not keeptar: print('Deleting the following files: ', end='')
        for idx, row in links_df.iterrows():
            curfn = row['filename']
            if not keeptar and 'tar.gz' in curfn:
                with contextlib.suppress(BaseException):
                    os.remove(os.path.join(destdir, curfn))
                    print(curfn, end=', ')
                    
        # Done: 
        print('\nCompletely done with downloading & unpacking\n')
        return None

class Combine(AutoPRSTCLI): #, AutoPRSTSubparser):
# class CombinerUtil(DownloadUtil): #, AutoPRSTSubparser):
#     pass
    
    '''\
    A tool to combine genetics-related text files.
    '''
    
    
    @classmethod
    def _get_cli_spkwg(cls, basic_pkwargs=True):
        nargskeys = ['input','selectcols','sortcols','assertunique','antiglobs']
        reqkeys = ['input','out']
        spkwg = super()._get_cli_spkwg(basic_pkwargs=basic_pkwargs)
        for key in nargskeys: spkwg['groups']['general']['pkwargs'][key]['kwargs'].update(nargs='+')
        for key in reqkeys: spkwg['groups']['general']['pkwargs'][key]['kwargs'].update(required=True)
        return spkwg
    
    insert = 'a long text text'*20
    
    @classmethod
    def from_cli_params_and_run(cls,
            #testarg:Param(f"A notebook name or glob {insert} to convert", str, required=True)='defaultstr', # something
            input:str=None, # Input files to be read and combined. Inputs assumed to be delimited text files,, for now.
            out:str=None, # Output file name.
            reqlen:int=None, # require the number of input files to be a specific number, else procedures stops. (e.g. 22 so it combine 22 chroms)
            sep:str='\t', # Seperator for the inputs files, default is \t (tab).
            antiglobs:list=['*.log','*.tmp'], # Globs/patterns that should be removed from the input files. e.g. ['*.log'] removes files ending with .log
            intype:Param('Type of input file. The \'auto\' default option detects this'
                         ' automatically based on extensions. If extension is not recognized it assumes'
                         ' \'headed-txt\'. Other options are: \'headless-txt\', \'legacyweights.tsv\', \'prstweights.tsv\'', str)='auto', # kjherkjgekrj
            assertunique:str=["SNP,A1,A2"], # Check if the row is unique for these columns. Ignored if columns are not present. Empty list (=[]) is you want this check to not happen.
            sortcols:str=None, # The columns that should be used for sorting.
            selectcols:str=None, # In case not specified the columns of the input files are used.
            noheader:bool=None, # Dont use a header for the output file. For plink and other tools this is sometimes needed.
            outtype:str='tsv', # Specify the filetype to save. For now the option is tab separated file (tsv).
            pyarrow=True,
            verbose=True,
            #     # Following two are needed to not have the args.func(var(args)) not crash: 
            #     command=None,
            #     func=None, 
            **kwg # this kwg catches command and func for a smooth run
            ):
        
        # Preprocessing & Checks:
        assert outtype in ['tsv','prstweights.tsv','legacyweights.tsv'], f"outtype given ({outtype}) is not a valid option."
        import pandas as pd; from tqdm.cli import tqdm; import fnmatch; from prstools.models import BasePred
        assert input is not None and out is not None, '--input/out are required arguments. Supply these arguments.'
        def splitter(lst):
            if lst: lst=[chunk for elem in lst for chunk in elem.split(',')]
            return lst
        antiglobs=splitter(antiglobs); assertunique=splitter(assertunique);
        sortcols=splitter(sortcols); selectcols=splitter(selectcols);
        fn_lst = input; oldlen = len(fn_lst); assert type(fn_lst) is list
        assert not any('*' in fn for fn in fn_lst), 'Found \'*\' in input string, this is not supported yet, you probably typed --input=*, --input * does work.'
        fn_lst = [fn for fn in fn_lst if not any(fnmatch.fnmatch(fn, antiglob) for antiglob in antiglobs)]
        #print(fn_lst)
        extra = f'Removed {oldlen-len(fn_lst)} with antiglob (e.g. with *.log extension)' if oldlen-len(fn_lst) > 0 else ''
        #print(fn_lst)
        if verbose: print(f'The number of input is {len(fn_lst)}, now processing them. {extra}')
        if reqlen: assert len(fn_lst) == reqlen, f'The required number of files (reqlen={reqlen}) is not present, so this code is exiting.'
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
            options=['legacyweights.tsv', 'prstweights.tsv']
            for this in options:
                if fn.endswith(this): return this
            return 'headed-txt'
        # Loading of files from disk: 
        pbar = tqdm(fn_lst, miniters=1); df_dt = {} #, refresh=True)
        for fn in pbar:
            pbar.set_postfix_str(f'loading: {fn}')
            cur_intype = detect_ftype(fn) if intype == 'auto' else intype
            header_dt = dict(header=None) if cur_intype in ['headless-txt', 'legacyweights.tsv'] else {}
            names = None if not cur_intype == 'legacyweights.tsv' else BasePred.default_weight_cols
            df = pd.read_csv(fn, sep=sep, names=names, **header_dt, **prw)
            df_dt[fn] = df
            #print(df.head())
        totrows = sum(df.shape[0] for df in df_dt.values())
        
        # Checks:
        ncols = np.unique(list(df.shape[1] for df in df_dt.values()))
        if len(ncols) > 1: raise Exception(f'The input files have different numbers of columns ({ncols}).'
                                           'Make sure all the input files have the same number of columns')
        if verbose: print(f'Read {len(fn_lst)} files with a total of {totrows} rows and {ncols} columns.')
        
        # Combining dataframes:
        if verbose: print(f'Concatenating the {len(fn_lst)} frames',end=' ')
        ordered_keys = df_dt.keys()
        df = pd.concat([df_dt[key] for key in ordered_keys], axis=0) # axis=0 is default.
        if verbose: print('-> concatenation done. A view of the first and last rows:')
        if verbose: print(df)
        
        # More checks:
        assert df.shape[1] == ncols[0], (f'Oke so the df.shape[1]={df.shape[1]}, which means the input '
            'files must have different column names. modify the column names s.t. they are the same.')
        if assertunique:
            if not (set(assertunique) - set(df.columns)):
                assert not df[assertunique].duplicated().any(), (f'Rows for columns {assertunique} the combined frame are not unique. '
                'There are duplicates. Change assertunique options or remove duplicates.')
            else: warnings.warn(f'{set(assertunique) - set(df.columns)}, not present in frame columns. cannot perform assertion of uniques for columns {assertunique}')
        
        # Sorting & Selecting
        if sortcols: df=df.sort_values(sortcols)
        if selectcols: df=df[selectcols]

        # Save the combined dataframe:
        if verbose: print(f'Saving the combined frame to: {out}')
        header = True if not noheader else False
        if not selectcols: selectcols = df.columns
        dcols = BasePred.default_weight_cols
        if outtype == 'legacyweights.tsv': selectcols = dcols; header=False
        if outtype == 'prstweights.tsv': selectcols = dcols + [col for col in df.columns if col not in dcols]; header=True # other columns in the back
        df[selectcols].to_csv(out, sep='\t', index=False, header=header)
        if verbose: print(f'----- All done with combining and saving data -----') 
            
        return df
    
# class PRSTLogs(dict):
#     _instance = None
#     def __new__(cls):
#         if cls._instance is None:
#             cls._instance = super().__new__(cls)
#             cls._instance._prstlogs_fn = None
#         return cls._instance

#     def set_prstlogs_fn(self, filename, save=True):
#         """Set log filename and ensure the directory (one level deep) exists."""
#         dir_name = os.path.dirname(filename)
#         if dir_name:
#             os.makedirs(dir_name, exist_ok=True)
#         self._prstlogs_fn = filename
#         if save:
#             self.save()

#     def save(self):
#         """Save the dictionary to the log file if set."""
#         def _safe_json(obj):
#             """Ensure JSON serialization does not fail by converting non-serializable objects to strings."""
#             try:
#                 json.dumps(obj)  # Try serializing directly
#                 return obj
#             except (TypeError, OverflowError):
#                 return str(obj)  # Fallback to string representation
        
#         if self._prstlogs_fn:
#             with open(self._prstlogs_fn, "w") as f:
#                 json.dump(self, f, indent=2, default=_safe_json)
                
#     def finish(self):
#         if hasattr(self,'_prstlogs_fn'):
#             self.save()

            
# import json
# import os

class CycleDict(dict):

    def __getitem__(self, key):
        """Override to return a new PRSTLogs instance for missing keys."""
        if key not in self:
            self[key] = CycleDict()  # Automatically create a new PRSTLogs for missing keys
        return super().__getitem__(key)

class PRSTLogs(CycleDict):
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            tic, toc = cls._instance.get_tictoc(); tic('')
            cls._instance._prstlogs_fn = None
        return cls._instance

    def set_prstlogs_fn(self, filename, save=True):
        """Set log filename and ensure the directory (one level deep) exists."""
        dir_name = os.path.dirname(filename)
        if dir_name:
            os.makedirs(dir_name, exist_ok=True)
        self._prstlogs_fn = filename
        if save:
            self.save()

    def save(self):
        """Save the dictionary to the log file if set."""
        def _safe_json(obj):
            """Ensure JSON serialization does not fail by converting non-serializable objects to strings."""
            try:
                json.dumps(obj)  # Try serializing directly
                return obj
            except (TypeError, OverflowError):
                return str(obj)  # Fallback to string representation
        
        if self._prstlogs_fn:
            with open(self._prstlogs_fn, "w") as f:
                json.dump(self, f, indent=2, default=_safe_json)
                
    def get_tictoc(self):
        if not hasattr(self,'timer'):
            self.timer = Timer(get_memory_usage)
        return self.timer.tic, self.timer.toc

    def finish(self):
        if hasattr(self, '_prstlogs_fn'):
            self.save()

# Global access function
def get_prstlogs(): return PRSTLogs()

def get_argnames(fun):
    import inspect
    return inspect.signature(fun).parameters.keys()


if not '__file__' in locals():
    import sys
    if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+['ipykernel_launcher.py' in sys.argv[0]]):
        with open('../prstools/utils.py', 'w') as loadrf: loadrf.write(In[-1])
        print('Written to:', loadrf.name)
        if 'In' in locals() and _isdevenv_prstools:
            print('starting here in models:') 
            get_ipython().system('prst --dev | head -3')
            