#| export
import os, json, time, shutil
from prstools._ext_utils import *
import prstools as prst
import warnings, threading
from abc import ABC, abstractmethod
try:
    from fastcore.script import call_parse, Param
except:
    def Param(*args, **kwg):
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
    
# from prstools._cmd import load_config, save_config, remove_config
_devonlymsg =  "If you see this message something went wrong in an unexpected way. Please contact the developer."

def load_config(*a, **kw):
    from prstools._cmd import load_config as _f; return _f(*a, **kw)
def save_config(*a, **kw):
    from prstools._cmd import save_config as _f; return _f(*a, **kw)
def remove_config(*a, **kw):
    from prstools._cmd import remove_config as _f; return _f(*a, **kw)

def plot_manhattan(data_df, x=None, y='-logp', regcol='chrom', pvalmin=1e-323, palette='bright', aspect=4, s=6.,title='Manhattan plot', **snskwg):
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from prstools.io import compute_pvalbetase
    data_df = data_df.reset_index(drop=True).reset_index()
    if not '-logp' in data_df.columns:
        if not 'pval' in data_df.columns:
            data_df = compute_pvalbetase(data_df, calc_lst=['pval'], pvalmin=pvalmin)
        data_df['-logp'] = -np.log10(data_df['pval'])
    if x is None: x = 'index'
    plot = sns.relplot(data=data_df, x=x, y=y, aspect=aspect, hue=regcol, palette=palette, legend=None, s=s, **snskwg)
    chrom_df=data_df.groupby(regcol)[x].median()
    plot.ax.set_xlabel(regcol); plot.ax.set_xticks(chrom_df); plot.ax.set_xticklabels(chrom_df.index)
    plot.fig.suptitle(title)
    plt.show()
    
manhattan_plot = plot_manhattan


def validate_path(*, must_exist=True, handle_prstdatadir=False, verbose=False, **kwg):
    if not type(handle_prstdatadir) is list: handle_prstdatadir = [handle_prstdatadir for _ in range(len(kwg))]
    else: len(handle_prstdatadir) == len(kwg), 'input paths for validation must match handle_prstdatadir list in length'
    for chandle_prstdatadir in handle_prstdatadir: assert chandle_prstdatadir in ('only','allow',False) 
    
    newargs=[]
    prstcfg = load_config()
    links_df = _get_linksprst()
    prstdatadir = prstcfg.get('prstdatadir', None)
    auto_download = prstcfg['auto_download']
    #if prstdatadir is not None: assert type(prstdatadir) is str, f'prstdatadir must be of type string! now it is: {type(prstdatadir)}'
    
    for chandle_prstdatadir, (key, arg) in zip(handle_prstdatadir, kwg.items()):
        
        if not arg: newargs += [arg]; continue # Skipp and continue if None or False
        assert type(arg) is str, 'Inputs must be strings to be valid and be validated as paths.'
        arg = os.path.expanduser(arg)
        argexists = os.path.exists(arg)
        do_prstdd = must_exist and (not argexists or (chandle_prstdatadir=='only')) and chandle_prstdatadir in ('allow','only')
        
        # prstdatadir part: see if it is in the pdd
        if do_prstdd: # note: a lot of this determination of status_prstdatadir is not subsequently used.
            msg = f'prstdatadir is not set, but required to retrieve "{arg}" from prstdatadir, for option {key}. Please set prstdatadir with "prst config".'
            if not prstdatadir and (chandle_prstdatadir=='only'): raise ValueError(msg)

            # Determine status_prstdatadir and potentially correct arg
            is_prstdatadir_compliant = len(arg)>0 and (not arg[0] in '/.') and (arg.count('/') <= 2)
            if prstdatadir and not os.path.exists(prstdatadir):
                msg = (f'\033[1;31mWARNING: It seems prstools data storage location '
                           f'(i.e. prstdatadir) does not exist (anylonger) @ {prstdatadir}\033[0m')
                if not os.path.exists(prstdatadir): warnings.warn(msg)
                status_prstdatadir = 'dirdisappeared'
            elif prstdatadir and is_prstdatadir_compliant:
                tstarg = os.path.join(prstdatadir, arg)
                argexists = os.path.exists(tstarg)
                if argexists: arg=tstarg
                if key=='ref' and not argexists:
                    tstarg = os.path.join(prstdatadir, 'ldblk_'+arg);
                    if os.path.exists(tstarg): arg=tstarg
                argexists = os.path.exists(arg)
                status_prstdatadir= 'argexists' if argexists else 'argabsent'
            elif prstdatadir: status_prstdatadir = 'argnoncompliant'
            else: status_prstdatadir='notset'
        
        # links_df part: see if it can be downloaded, and download if possible
        do_prstlink = do_prstdd and (not argexists and is_prstdatadir_compliant)
        ind = links_df['filename'].str.replace('.tar.gz','') == arg
        ind = ind | (links_df['filename'].str.replace('.tar.gz','').str.replace('ldblk_','') == arg)
        if do_prstlink and ind.sum() == 1:
            df = links_df[ind]
            print(f'The data for "{arg}" was not yet available in the prstools data directory, but the data is downloadable.')
            msg = (f'Requested {arg} is available for download, but your prstdatadir is not set or auto_download disabled so cannot download {arg}. ' 
                    'You can manually download and set prstdatadir (using "prst config").')
            if auto_download and prstdatadir: DownloadUtil.from_cli_params_and_run(destdir=prstcfg['prstdatadir'], pattern=arg)
            else: raise RuntimeError(msg)
            bn = df['filename'].str.replace('.tar.gz','').iloc[0]
            tstarg = os.path.join(prstdatadir, bn)
            assert os.path.exists(tstarg), (f'Although based on input "{arg}" a download was performed. '
                f'This did not result in a file or directory at {tstarg}. This should not happen and is a bug. please contact developer')
            arg = tstarg
        elif do_prstlink and ind.sum() > 1:
            links_df['matching'] = ind
            msg = ''
            msg += f'Issue with input argument "--{key}"\n'
            msg += 'It was not present in prstools data dir, but perhaps it can be downloaded.\n'
            msg += 'However it did not match uniquely with the available data downloads:\n'
            msg += str(links_df[['filename','matching','description']]) + '\n'
            msg += 'It should match uniquely!\n\n'
            raise ValueError(msg+f'"--{key}" should match uniquely for data download.')
        elif do_prstlink:
            msg = f'Could not find "{arg}" no such file or directory. '
            msg += f'Argument "{arg}" for option "--{key}" was not found in prstools data dir, current dir or available in downloads.\n'
            msg += f'Have your argument for "--{key}" match a fitting option from the table above, or specify a path that exists.'
            disp_df = links_df[['filename','description']]
            disp_df.loc[:,'filename'] = disp_df['filename'].str.replace('.tar.gz','').str.replace('ldblk_','')
            disp_df = disp_df.rename(columns=dict(filename=f'possible arguments for --{key}'))
            msg = '\n' + str(disp_df) + '\nFileNotFoundError: ' + msg
            raise FileNotFoundError(msg)
            
        if not os.path.exists(arg) and must_exist: # <- Ok so after handling it still does not exist, we will be throwing a custom error
            msg = f'Could not find "{arg}" no such file or directory.' # you could make the messages a bit better
            xtra_maybelater = '\nNote: the argument does start with . or / or has 2+ of slashes. ' if do_prstdd and not is_prstdatadir_compliant else ''
            xtra=''; msg = msg + xtra
            raise FileNotFoundError(msg)
        newargs += [arg]

    return newargs[0] if len(newargs) == 1 else tuple(newargs)
    
# def validate_path(*, must_exist=True, handle_prstdatadir=False, verbose=False, **kwg):
#     assert handle_prstdatadir in ('only','allow',False)
#     newargs=[]
#     prstcfg = load_config()
#     links_df = _get_linksprst()
#     prstdatadir = prstcfg.get('prstdatadir', None)
    
#     for key, arg in kwg.items():
        
#         if not arg: newargs += [arg]; continue # Skipp and continue if None or False
#         assert type(arg) is str, 'Inputs must be strings to be valid and be validated as paths.'
#         arg = os.path.expanduser(arg)
        
#         if must_exist and (not os.path.exists(arg) or (handle_prstdatadir=='only')):
            
#             checked_prstdatadir=False
            
#             msg = f'prstdatadir is not set, but required to retrieve "{arg}" from prstdatadir. Please set prstdatadir with prst config.'
#             if not prstdatadir and (handle_prstdatadir=='only'): raise ValueError(msg)
                
#             if handle_prstdatadir and type(prstdatadir) is str and len(arg)>0 and (not arg[0] in '/.') and (arg.count('/') <= 2):
                
#                 msg = f'\033[1;31mWARNING: It seems prstools data storage location (i.e. prstdatadir) does not exist (anylonger) @ {prstdatadir}\033[0m'
#                 if not os.path.exists(prstdatadir):
#                     warnings.warn(msg)
#                 else:
#                     tstarg = os.path.join(prstdatadir, arg)
#                     prstargexists = os.path.exists(tstarg)
#                     if os.path.exists(tstarg): arg=tstarg
#                     if key=='ref':
#                         tstarg = os.path.join(prstdatadir, 'ldblk_'+arg)
#                         if os.path.exists(tstarg): arg=tstarg
#                     checked_prstdatadir=True
                    
#                 if (not os.path.exists(arg) or (handle_prstdatadir in ('only','allow') and not prstargexists)):
#                     #links_df = prst.utils._get_linksprst() 
                    
#                     # See if present in links_df
#                     if key == 'ref':
#                         ind = links_df['filename'].str.startswith('ldblk_')
#                         links_df = links_df[ind].reset_index(drop=True)
#                     ind = links_df['filename'].str.replace('.tar.gz','') == arg
#                     ind = ind | (links_df['filename'].str.replace('.tar.gz','').str.replace('ldblk_','') == arg)
                    
#                     if ind.sum() == 1:
#                         df = links_df[ind]
#                         print(f'The data for "{arg}" was not yet available in the prstools data directory, so it is downloaded now.')
#                         if prstcfg['auto_download']: DownloadUtil.from_cli_params_and_run(destdir=prstcfg['prstdatadir'], pattern=arg)
#                         else: raise RuntimeError(f'auto_download disabled so cannot download {arg}')
#                         bn = df['filename'].str.replace('.tar.gz','').iloc[0]
#                         tstarg = os.path.join(prstdatadir, bn)
#                         assert os.path.exists(tstarg), (f'Although based on input "{arg}" a download was performed. '
#                             f'This did not result in a file or directory at {tstarg}. This should not happen and is a bug. please contact developer')
#                         arg = tstarg
#                     elif ind.sum() > 1:
#                         links_df['matching'] = ind
#                         msg = ''
#                         msg += f'Issue with input argument "--{key}"\n'
#                         msg += 'It was not present in prstools data dir, but perhaps it can be downloaded.\n'
#                         msg += 'However it did not match uniquely with the available data downloads:\n'
#                         msg += str(links_df[['filename','matching','description']]) + '\n'
#                         msg += 'It should match uniquely!\n\n'
#                         raise Exception(msg+f'"--{key}" should match uniquely for data download.')
#                     else:
#                         msg = ''
#                         msg += f'Argument "{arg}" for option "--{key}" was not found in prstools data dir, current dir or available in downloads.\n'
#                         msg += f'Have your argument for "--{key}" match one of these options:\n'
#                         disp_df = links_df[['filename','description']]
#                         disp_df.loc[:,'filename'] = disp_df['filename'].str.replace('.tar.gz','').str.replace('ldblk_','')
#                         msg += str(disp_df) + '\n'
#                         raise Exception(msg + f"-> So in the end could not find '{arg}'")
                    
#             if not os.path.exists(arg): # <- Ok so after handling it still does not exist, we will be throwing a custom error
#                 if checked_prstdatadir: print(f'\n--> Also looked for \'{arg}\' in prstdatadir: {prstdatadir} <--\n')
#                 notsetbuthere_prstdatadir = True 
#                 if (handle_prstdatadir in ('only', 'allow') and prstdatadir is None) else False:
                    
#                 try:
#                     open(arg) # this is to make it throw an error, because the file/dir does not exist. 
#                 except Exception as e:
#                     msg = f'Could not find "{arg}" no such file or directory.'
#                     xtra = ', '.join([elem for elem in os.listdir(prstdatadir) if elem[:1] != '.']) if checked_prstdatadir else ''
#                     xtra = f'Looked in prstdatadir too.\nFor this the options are: {xtra}' if checked_prstdatadir else ''
#                     xtra += '' if denied_prstdatadir else ''
#                     if 
#                     msg = msg + xtra
#                     raise FileNotFoundError(msg) from e
#             #import errno
#             #raise FileNotFoundError( errno.ENOENT, os.strerror(errno.ENOENT), arg)
#             #raise FileNotFoundError(f'No such file or directory: \'{arg}\'')
#         newargs += [arg]
#     return newargs[0] if len(newargs) == 1 else tuple(newargs)

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
        #self._shutdown = mgr.shutdown      # plain function, not the Manager itself

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


def get_pbar(iterator, maxncols=200, colour='green', bar_format='{l_bar}{barstr}{r_bar}',
            mininterval=0.5, desc=None, fakebar=None, barlen='auto', deamon=False, **kwg):
    # Prep:
    assert not 'total' in kwg, 'Need to use iterator, cannot use total argument.'
    create_pbar = get_tqdm()
    detcols = shutil.get_terminal_size((999, None)).columns
    ncols = min(maxncols, detcols)
    barlen = max(min(70, detcols-65), 10)
    bar_format = bar_format.format_map(prst.utils.AutoDict(barstr=f'{{bar:{barlen}}}'))
    
    # Pbar creation:
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
            gbs = prst.utils.get_memory_usage(show=False)
            if gbs: 
                peak = max(gbs, getattr(fakebar,'_peak', 0.0))
                pbar.set_postfix({'RAM': f"{gbs:.2f}G",'Peak': f"{peak:.2f}G"}, refresh=False)
                fakebar._peak = peak
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

def clear_pbars(verbose=True):
    try:
        from tqdm import tqdm
        for inst in list(tqdm._instances): inst.close()
        if verbose: print('Succesfully cleared all pbars.')
    except:
        if verbose: print('Issue with clearing pbars.')

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
    from textwrap import dedent, indent
    from pprint import PrettyPrinter
#     NoneType = type(None)
    subparserkwg_lst = parse_args(argv=[], subparserkwg_lst=subparserkwg_lst, return_spkwg=True)
    lst.append(subparserkwg_lst)
    proc_dt = dict()
    def fun(obj):
        special_notactuallytypes = ['intcaster', 'str_or_none']
        if type(obj) is type:
            proc_dt[repr(obj)] = obj.__name__
            return repr(obj)
        elif obj.__name__ in special_notactuallytypes:
#             proc_dt[obj.__name__] = obj.__name__
            proc_dt[repr(obj)] = obj.__name__
            return obj.__name__ #repr(obj)
        else:raise Exception('This is mjw prst code, it seems a particular argument cannot be transformed into a argparse dict')
    if store:
        string = json.dumps(subparserkwg_lst, default=fun, indent=2)
        string = PrettyPrinter(indent=1, width=200, sort_dicts=sort_dicts).pformat(subparserkwg_lst)
        stringproc = string + ''
        stringproc = '\n'.join([' '*4 + elem for elem in stringproc.split('\n')]) # do pyhton indent
        for key, item in proc_dt.items():
            stringproc = stringproc.replace(key,item)
        specstring = ''
        specstring += 'NoneType = type(None)\n'
        specstring += 'def str_or_none(x): return None if x.lower() == "none" else x\n'
        specstring += 'def intcaster(x): return int(float(x))\n'
        exec(specstring)
        test_lst = eval(stringproc)
#         try:
#             assert test_lst == subparserkwg_lst
#         except: get_ip().embed()
        finstring = ''
        finstring += '# Warning: Dont edit here, This file was generated automatically, using a developer tool.\n'
        finstring += "def get_subparserkwg_lst():\n"
        finstring += '\n'.join([4*' '+elem for elem in specstring.split('\n')]) + '\n'
        finstring += '    subparserkwg_lst = '+stringproc + '\n'
        finstring += "    return subparserkwg_lst"
        ns = {}; exec(finstring, ns)
        lst = ns["get_subparserkwg_lst"]() # small code exe test
        fn = '../prstools/_parser_vars.py'
        with open(fn, 'w') as f: f.write(finstring)
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
            param = v.anno
            if not isinstance(param,Param): 
                param = Param(v.docment, v.anno, default=v.default)
                param.default = v['default']
            else:
                param = Param(param.help,param.type) #, default=v.default)
                param.default = v['default']
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
        display_info = getattr(cls, '_display_info', False)
        spkwg = dict(
            cmdname    = cls.__name__.lower(), #command name
            clsname    = cls.__name__, #class name
            description= doc, # Help and discription are in the right spot, naming is counter intuitive, save urself time and dont check 
            help       = doc.split('\n')[0], # and dont check again.
            epilog     = epilog,
            display_info = display_info,
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
        
class Config(AutoPRSTCLI): #, AutoPRSTSubparser):
    '''\
    Interactively set configuration for prstools (run this first if you are new to prstools).
    This command always runs interactively and prompts the user to configure
    the prstools data directory and related settings.
    '''
    @classmethod
    def from_cli_params_and_run(cls,
        setup:bool=False,    # Run the interactive setup. Currently the config is -> auto_download: {auto_download}, prstdatadir: {prstdatadir}
        noprompt:bool=False, 
        # Under development: Currently it is only possible to run the prstools config command with the interactive setup or use --yes2all
        yes2all:bool=False,  # Run setup and say yes to all the questions being asked
        prstdatadir='make in current directory',
        auto_download=True,
        **kwg
        ):
        
        assert prstdatadir == 'make in current directory', 'This is currently the only valid option!'
        print('You have activated the prstools configuration setup!\n')
        print("Press ctrl+c to cancel any time.\n")
        prstcfg = prst.utils.load_config()
        oldprstdatadir = prstcfg.get('prstdatadir',None)
        if oldprstdatadir: print(f"The current prstdatadir is: {oldprstdatadir})\n")
        if yes2all: print('--yes2all option activated!\n')
        # prompt prstdatadir current? 
        # prompt autodownload 
        
        def ask(prompt, default=True):
            if yes2all: return True
            yn = "[Y/n]" if default else "[y/N]"
            while True:
                reply = input(f"{prompt} {yn}: ").strip().lower()
                if reply == "" and default is not None:
                    return default
                if reply in ("y", "yes"):
                    return True
                if reply in ("n", "no"):
                    return False
                print("Please answer y or n.")

        # --- prompts ---
        prstdatadir = os.path.join(os.getcwd(), "prstdata")
        answer = ask("Create the prstools data directory (prstdatadir) in the current path?\n"
                    "(Used forinstance to store LD reference data. Run this command in a different dir if you'd like another location)\n"
                    f"@ {prstdatadir}")
        if not answer: prstdatadir=None
        if answer: os.makedirs(prstdatadir, exist_ok=True)
        answer = ask("Enable automatic download of reference datasets?", default=auto_download)
        auto_download = answer

        updater=dict(prstdatadir=prstdatadir,
                            auto_download=auto_download)
        print('New configuration:')
        for key, item in updater.items():
            print(f'  {key:15} = {item}')
        print()
        
        if not all(updater.values()): 
            print('\033[1;38;2;179;125;19mWarning: If prstdatadir and auto_download are not set you will miss out on several prstools features.\033[0m\n')
        
        # populate config dictionary
        prstcfg = load_config()
        prstcfg.update(updater)
        print('Saving config to ~/.prstools/config.json ', end='')
        save_config(prstcfg)
        print('-> Done')
        return None
        
def extract_fn(url):
    import urllib
    return urllib.parse.urlparse(url).path.split('/')[-1]

def download_tar(url, dn='', desc_prefix=''):
    import urllib.request
    try: from tqdm.cli import tqdm
    except: tqdm = lambda x:x
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
            
def _get_linksprst():
    import pandas as pd
    data = [
    ["snpdb_mini.tsv.gz", 'https://www.dropbox.com/scl/fi/q54612a8zs1guoo39xj2z/snpdb_mini.tsv.gz?rlkey=ujc60us9j4j0bzuctrh0fwgw3&st=d67dx3eo&dl=1', 
     'tiny snpdb_full used for build detection (~0.2M)'],
    ["snpdb_full.tsv.gz", 'https://www.dropbox.com/scl/fi/wujxo9rmwcfltoljx9hlc/snpdb_full.tsv.gz?rlkey=83mmcxwylff17h9k67xdz2zec&st=1dgiy1te&dl=1', 'SNP db with rsids & hg19+38 positions in 1kg (~0.93G)'],
    ["snpinfo_mult_1kg_hm3",  "https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=1", "1000G multi-ancestry SNP info (for PRS-CSx) (~106M)"],
    ["snpinfo_mult_ukbb_hm3", "https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=1", "UKBB multi-ancestry SNP info (for PRS-CSx) (~108M)"],
    ["ldblk_1kg_afr.tar.gz",  "https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=1", "1000G AFR Population LD panel (~4.44G)"],
    ["ldblk_1kg_amr.tar.gz",  "https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=1", "1000G AMR Population LD panel (~3.84G)"],
    ["ldblk_1kg_eas.tar.gz",  "https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=1", "1000G EAS Population LD panel (~4.33G)"],
    ["ldblk_1kg_eur.tar.gz",  "https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=1", "1000G EUR Population LD panel (~4.56G)"],
    ["ldblk_1kg_sas.tar.gz",  "https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=1", "1000G SAS Population LD panel (~5.60G)"],
    ["1kg_hm3.tar.gz",        "https://www.dropbox.com/scl/fi/0nn9za9wbg6n0ki3e6371/1kg_hm3.tar.gz?rlkey=7v1hqr4jnacvvfqv6jj13lcd2&st=yv4fcssd&dl=1",
     "A 1kg plink dataset for hapmap3 snps with all pops (2.5K induv) (~348M)."],
    ["ldblk_ukbb_afr.tar.gz", "https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=1", "UKBB AFR Population LD panel (~4.93G)"],
    ["ldblk_ukbb_amr.tar.gz", "https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz?dl=1", "UKBB AMR Population LD panel (~4.10G)"],
    ["ldblk_ukbb_eas.tar.gz", "https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=1", "UKBB EAS Population LD panel (~5.80G)"],
    ["ldblk_ukbb_eur.tar.gz", "https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=1", "UKBB EUR Population LD panel (~6.25G)"],
    ["ldblk_ukbb_sas.tar.gz", "https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=1", "UKBB SAS Population LD panel (~7.37G)"],
    ["example.tar.gz", "https://www.dropbox.com/scl/fi/yi6lpbp0uhqiepayixvtj/example.tar.gz?rlkey=kvd7r17wuory9ucqdk4rh55jw&dl=1", "PRSTOOLS Example data (3.8M)"],
#     ["g1000.tar.gz",'https://www.dropbox.com/scl/fi/97lsbtoomhti3q6x2wttf/g1000.tar.gz?rlkey=9hd85oytgnpv6wvbapvu2rk2m&st=3k4fq9ub&dl=1', "European 1kg plink dataset for hapmap3 (~64M)"]
        #["example.tar.gz","https://www.dropbox.com/scl/fi/7fg6c9e5dnmb0n4cdfquz/example.tar.gz?rlkey=31u2948paz539uw61jq37oe8s&dl=1", "PRSTOOLS Example data (70mb)"] 
    ]
    columns = ["filename", "url", "description"]
    links_df = pd.DataFrame(data, columns=columns)[['filename','description','url']]
    return links_df
        
class DownloadUtil(AutoPRSTCLI): #, AutoPRSTSubparser):
    
    '''\
    Download and unpack LD reference panels and other data.
    The files that can be downloaded and unpacked with this command includes the standard reference files for PRS-CS and PRS-CSx.
    Additional information can be found at https://github.com/getian107/PRScsx
    '''
    
    @classmethod
    def from_cli_params_and_run(cls,
    destdir:str=None, # Directory in which all the data will be downloaded. Option required if you want the download to start. If you specify 'prstdatadir' it will try to save to the prstdatadir (if configured).
    pattern:str='ALL',# A string pattern that retrieves every file that it matches. Matches everything by default. Without --destdir option (required to start downloading) one can see which files get matched.
    list=False, # Show which files can optionally be downloaded and exit.
    listfull=False, # Show complete filelist with long urls and exit.
    mkdir=False,
    command=None,
    func=None,
    keeptar=False, # Keep the tar.gz files in the destdir. If this option is not given they will be deleted automatically to save space.
    **kwg
                   ):

        # Defs & Inits:
        prstcfg = load_config()
        import tarfile, contextlib, pandas as pd
        try: from tqdm.cli import tqdm
        except: tqdm = lambda x:x
        links_df = _get_linksprst()
        def truncate_middle(s, maxlen=60):
            if s is None or len(s) <= maxlen:
                return s
            half = (maxlen - 3) // 2
            return s[:half] + "..." + s[-half:]

        # Preprocessing:
        if list or listfull: # Overloading list is bad practice, I know.
            print('\nFiles available for downloading & unpacking (if difficult to read consider making terminal wider):\n')
            if not listfull: links_df["url"] = links_df["url"].map(lambda x: truncate_middle(x, 40))
            print(links_df.to_string(index=False, justify="left"))
            return True
        if not pattern == 'ALL':
            print('\nA pattern was used, which matched the following files:') 
            ind = links_df['filename'].str.contains(pattern)
            links_df = links_df[ind]
            print(', '.join(links_df['filename'].to_list()))
        else: print('Downloading all files (no filter --pattern was used):\n',', '.join(links_df['filename'].to_list()))
        if destdir == 'prstdatadir':
            prstdatadir = prstcfg['prstdatadir']
            print(f'Set destination is the prstdatadir! @ {prstdatadir}')
            if prstdatadir is None: raise Exception('prstdatadir was not set in config! run prst config to set it')
            destdir = prstdatadir
        if not destdir: print('\n--destdir was not specified so download will not start.\n'); return True
        if not os.path.isdir(os.path.expanduser(destdir)):
            if mkdir: raise Exception(f'\nIt appears the supplied --destdir does not exist, please create: {destdir}')

        # Downloading:
        print(f'\nDownloading data, which might take some time. Data will be stored in: {destdir} \n')
        lst = []
        for idx, row in links_df.iterrows():
            fn = os.path.join(os.path.expanduser(destdir), row['filename']); dn=fn.replace('.tar.gz','')
            if os.path.isdir(dn): print(f'For {fn} the associated directory {dn} already exists,'
                                        ' therefore download & unpack is skipped. Remove directory for a redownload.'); lst+=[False]; continue
            download_tar(row['url'], dn=destdir); lst+=[True]
        links_df = links_df[lst]

        # Untarring:
        print('\nFinished downloading all data. Now we need to unpack all the tar.gz files (takes some time to start):') 
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
    
    '''\
    A tool to combine genetics-related text files. If all is well you won't need this.
    '''
    
    
    @classmethod
    def _get_cli_spkwg(cls, basic_pkwargs=True): ## This badboi wraps the super method to enhance it.
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
    

class Transform(AutoPRSTCLI): #, AutoPRSTSubparser):
    
    '''\
    Transform messy sumstats and other inputs into neat standardized formats
    It can process whole directories at a time, using stars (--sst sst-prefix_*.tsv)
    '''
    
    _this_msg = 'yoyoyoy lololo ergergerg'
    
    _transform_mode_helpmsg =  ("yooo Different modes for file handling. If 'strict' is present an error in forinstance "
        "the loading of a sumstat will lead to a stop. 'redo' will force redo if output exists (else it will skip). with --mode strictredo you will combine options.")
    
    @classmethod
    def _get_cli_spkwg(cls, basic_pkwargs=True): ## This badboi wraps the super method to enhance it.
        nargskeys = ['sst']
        reqkeys = ['sst','out']
        spkwg = super()._get_cli_spkwg(basic_pkwargs=basic_pkwargs)
        for key in nargskeys: spkwg['groups']['general']['pkwargs'][key]['kwargs'].update(nargs='+')
        for key in reqkeys: spkwg['groups']['general']['pkwargs'][key]['kwargs'].update(required=True)
        return spkwg
    
    @classmethod
    def from_cli_params_and_run(cls,
            sst:str=None, # Input sumstats. Can be specified with * operator. For example "--sst CoolCohort_trait=*_hapmap3.tsv"
            colmap:str='{default_colmap}', # something here try default_colmap here too: {default_colmap}
            n_gwas:float=None, # Manually specify gwas sample size of the input sumstats. Most often it can be extracted from the sumstat itself so you do not need to specify this.
            out:str=None, # Output file name. A good suggestion: "--out {trimsst}.prstsst.tsv" 
            mode:Param(_transform_mode_helpmsg, str)='strict',
            # jergjkerkj jkergjkhegrjk
            sep:str='\t', # Seperator for the inputs files, default is \t (tab).
            n_eff_handling='raw',
            pyarrow=True,
            verbose=True,
            **kwg # this kwg catches command and func for a smooth run
            ):
        
        #dry:bool=False, # Specify to do a dry run to get an idea of what the command will do
        if type(sst) is str: fn_lst = [sst]
        elif type(sst) is list: fn_lst = sst
        else: ValueError(f'Input --sst is unrecognized type: {type(sst)}')
        new_lst=[]
        for fn in fn_lst:
            if any(elem in fn for elem in ['*','[']): new_lst += glob.glob(fn)
            else: new_lst += [fn]
        fn_lst = new_lst
        if colmap == '{default_colmap}': colmap = prst.io._get_default_colmap()
        out_lst = []
        uniqfns  = np.unique(fn_lst)
        nuniqfns = len(uniqfns)
        fn_lst = uniqfns
        for fn in fn_lst:
            sst = os.path.basename(fn)
            trimsst = prst.io.get_fn_trimmed(sst)
            out_fn = out.format(**locals())
            out_lst += [out_fn]
        nuniqouts = len(np.unique(out_lst))
        msg = f'\nGot {nuniqfns} unique input filenames and generated {nuniqouts} unique output filenames from them.'
        if (nuniqfns - nuniqouts) == 0: msg += '\nBecause there is 1 output for every 1 output we can proceed.'
        else:  msg += '\nBecause there is not 1 output for every 1 input we cannot proceed! Set e.g. "--out {trimsst}.tsv" to fix.'; raise ValueError(msg)
        print(msg)
        print(f'\nPath of first file that we will try to generate: {out_lst[0]}')
        print(f'It will be generated from: {fn_lst[0]}\n')
        
        pbar = prst.utils.get_pbar(list(zip(fn_lst,out_lst)))
        e=False; first = True; firstwarn=True
        for fn, out_fn in pbar:
            if verbose and first: print('\n'); first=False
            if not os.path.exists(out_fn) or 'redo' in mode:
                inject = '(i.e. already ready for subsequent steps) Note: this warning is only shown for the first occurance.'
                msg = f"\033[1;31m  '.prst' string detected in input filename {fn}. Are you sure you are not accidentally using prstools-ready homogenized sumstats as inputs? {inject}\033[0m"
                if firstwarn and '.prst' in fn: warnings.warn(msg); firstwarn=False
                try:
                    print(f'\nIt seems {os.path.basename(out_fn)} does not exist (or a redo was requested) so we are making it from {fn}')
                    sst_df = prst.load_sst(fn, colmap=colmap, verbose=verbose, delimiter=sep, n_gwas=n_gwas, n_eff_handling=n_eff_handling)
                    sst_df['snp'] = sst_df['snp'].fillna('NA')
                    prst.io.save_sst(sst_df, out_fn, verbose=verbose)
                except Exception as e:
                    if 'strict' in mode:
                        raise e from Exception('Set --mode flex to have the procedure skip input files that give issues.')
            else: 
                if verbose: print(f'Skipping input {fn} because output {out_fn} already exists.', flush=True)
        if e != False:
            print('This is the last error from the bunch:')
            raise e
            
    
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
            def fun():
                gbs = get_memory_usage(show=False)
                return f'{gbs:.3}G'
            self.timer = Timer(fun)
        return self.timer.tic, self.timer.toc
    
    @property
    def tic(self):
         return self.get_tictoc()[0]
    
    @property
    def toc(self):
        return self.get_tictoc()[1]

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
        print('Written to:', loadrf.name, flush=True)
        if 'In' in locals() and _isdevenv_prstools:
            print('starting here in models:') 
            get_ipython().system('prst --dev-secret | head -3')
            