import os, time, sys, argparse, json
# import gc,tracemalloc,objgraph,psutil

def get_default_cpus():
    cpus = 1 # This is the default number of cpus
    return cpus

def set_cpu_envvars(cpus, output_string=True):
    if cpus != -1: # disabling.
        #MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
        #print(f'Setting environmental variables to control number of cpus used (={cpus}).')
        n_cores = cpus # Correction, this does do something, although its quite opaque
        os.environ['OPENBLAS_NUM_THREADS'] = str(n_cores)
        os.environ['MKL_NUM_THREADS'] = str(n_cores)  # For Intel MKL
        os.environ['OMP_NUM_THREADS'] = str(n_cores)
        os.environ['NUMEXPR_NUM_THREADS'] = str(n_cores)  # For NumExpr
        #os.environ['OMP_DYNAMIC'] = 'FALSE' # After some really odd behavior seen on elyn where sum(proc-util) = 100% (set cpu count)
        #os.environ['OMP_THREAD_LIMIT'] = str(n_cores) # I decided to not include these last 4...
        #os.environ['OMP_NESTED'] = 'FALSE'
        #os.environ['OMP_PROC_BIND'] = 'TRUE'    
    else:
        True
    cpu_display_string = cpus if cpus != -1 else 'Not specified'
    return str(cpu_display_string)
        
# This functioning is needed here, a lazy retrieval s.t. imports only happen when stuff is actually run:
def retrieve_classmethod(*, clsname, methodname, modulename='prstools.models'):
    def pipefun(*args,**kwargs):
        import importlib
        cls = getattr(importlib.import_module(modulename), clsname)
        return getattr(cls,methodname)(*args,**kwargs)
    return pipefun

class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"

def process_argkwargs(kwargs, setverbosetrue=False):
    prstcfg = SafeDict(get_config())
    kwargs = kwargs.copy() 
    kwargs['help'] = argparse.SUPPRESS if kwargs.get('help', None) is None else kwargs['help']
    if 'default' in kwargs and type(kwargs['default']) is str and '{' in kwargs['default']:
        kwargs['default'] = kwargs['default'].format_map(prstcfg)
    if type(kwargs['help']) is str:
        kwargs['help'] = kwargs['help'].format_map(prstcfg)
    if kwargs.get('type',None) is bool:
        if not kwargs['default']:
            kwargs.pop('type')
            kwargs['action']='store_true'
            kwargs['default'] = argparse.SUPPRESS
        if setverbosetrue:
            raise NotImplementedError()
    if 'default' in kwargs and kwargs['default'] == 'SUPPRESS':
        kwargs['default'] = argparse.SUPPRESS
    if (kwargs['help'] != argparse.SUPPRESS):
        if not (kwargs.get('default', 'idontexist') in ['idontexist', argparse.SUPPRESS]):
            kwargs['help'] += f" (default: {kwargs.pop('default')})"
        if 'required' in kwargs:
            if kwargs['required'] is True:
                kwargs['help'] = kwargs['help'] + ' (required)'
    return kwargs
def format_color(text, color_code): return f"\033[{color_code}m{text}\033[0m"
def format_bold(text): return f"\033[1m{text}\033[0m"

## Formatting functionality:
class CustomArgumentParser(argparse.ArgumentParser):
    
        def format_usage(self):
            usage = super().format_usage()
            return self.usage_optionals_formattingx(usage)
        def format_help(self):
            string = super().format_help()
            string = self.usage_optionals_formattingx(string)
            return string.replace('usage:','\nUsage:\n')
        
        def usage_optionals_formattingx(self, string):
            string = string.replace('] [','  ')
            return string
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,argparse.RawDescriptionHelpFormatter): #, argparse.RawDescriptionHelpFormatter):
    def __init__(self, prog, indent_increment=1, max_help_position=60, width=None, minwidth=0):
        if not width: 
            import shutil # Width processing
            width = shutil.get_terminal_size().columns-2
            width = max(width, minwidth)
        _excl_lst = ['self', 'kwg_dt','_excl_lst','CustomFormatter','__class__','minwidth', 'shutil']
        kwg_dt = {key: item for key, item in locals().items() if not (key in _excl_lst)}
        super(CustomFormatter, self).__init__(**kwg_dt) #prog, indent_increment=indent_increment, max_help_position=max_help_position, width=width)      

    def _get_default_metavar_for_optional(self, action):
        return '<'+action.dest+'>' #.upper()

    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = '<'+action.dest+'>'#.upper()
                args_string = self._format_args(action, default)
                parts.extend(action.option_strings)
                parts = parts[::-1]
                parts[-1] += ' %s' % args_string
            return ', '.join(parts)
        
def get_defaultprstcfg():
    from prstools import __version__, _date
    defaultprstcfg = {
        'auto_download': False, 'prstdatadir': None, 
        '__version__':__version__, '_date': _date, 'basecmd': 'prst',
        'default_colmap':'SNP,A1,A2,BETA,OR,P,SE,N,FRQA1,CHR,BP',
        'default_conv_dt': {'CHR':'chrom','SNP':'snp', 'BP':'pos', 'A1':'A1','A2':'A2', 'MAF':'maf',
              'BETA':'beta','SE':'se_beta','P':'pval','OR':'oddsratio', 'N':'n_eff', 'FRQA1': 'af_A1'}
    }
    return defaultprstcfg

def save_config(prstcfg, fn='default', failsilent=False):
    if fn == 'default':
        fn = os.path.expanduser("~/.prstools/config.json")
    try: # Make config file if it does not exist, but pass if fails.
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        msg = ("# This is a prstools config file, it can be safely deleted to "
            "reset the config, if possible use \'prstools config\' to change the config\n")
        save_dt = dict(_comment=msg)
        save_dt.update(prstcfg)
        with open(fn, "w") as f: 
            f.write(json.dumps(save_dt, indent=2))
    except Exception as e: 
        if failsilent: pass
        else: raise e
            
def remove_config():
    fn = os.path.expanduser("~/.prstools/config.json")
    if os.path.exists(fn):
        os.remove(fn)

def load_config():
    fn = os.path.expanduser("~/.prstools/config.json")
    if os.path.exists(fn):
        try:
            with open(fn, "r") as f:
                prstcfg = json.loads(f.read())
        except Exception as e:
            msg = f'Something went wrong with loading config file: {fn} Deleting the file/reseting the config can fix the issue.'
            raise Exception(msg) from e
        ## Validations and assertions:
        defaultprstcfg = get_defaultprstcfg(); update=False
        for key, item in prstcfg.items():
            if type(item) is str and item[:1] == '~':
                prstcfg[key]=os.path.expanduser(item)
        if not all(key in prstcfg for key in defaultprstcfg.keys()):
            for key in defaultprstcfg.keys():
                if not key in prstcfg:
                    prstcfg[key] = defaultprstcfg[key]
                    update=True
        if update: save_config(prstcfg, fn, failsilent=True)
    else:
        prstcfg = get_defaultprstcfg()
        save_config(prstcfg, fn, failsilent=True) 
    
    return prstcfg

_prstcfg = None
def get_config(reload=False):
    global _prstcfg
    if _prstcfg is None and not reload:
        _prstcfg = load_config() 
    return _prstcfg

def parse_args(argv=None, description="Convenient and powerfull Polygenic Risk Score creation. \n\'prst\' is a commandline shorthand for \'prstools\'",
               subparserkwg_lst=None, basecmd='prstools', return_spkwg=False, reload=False):
    
    # Prepare parsing params:
    if argv is None: argv=sys.argv[1:]; basecmd=sys.argv[0]
    if subparserkwg_lst is None:
        if reload: print('Used reload=True, reload is not needed anymore if all is well.')
        from prstools._parser_vars import get_subparserkwg_lst
        subparserkwg_lst = get_subparserkwg_lst()
    else:
        from prstools.utils import process_subparserkwgs
        subparserkwg_lst=process_subparserkwgs(subparserkwg_lst)
        # subparserkwg_lst contains the information to construct parsers
        # This is generated from model code by the developer and saved
    prstcfg = get_config()
        
    # Construct Parser:
    parser = CustomArgumentParser(
        description=description,
        argument_default=argparse.SUPPRESS,
        add_help=False,
        formatter_class=CustomFormatter
    )
    general_group = parser.add_argument_group('General Options')
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help=argparse.SUPPRESS) # , help='Show help')
    
    if return_spkwg: return subparserkwg_lst
    def str_or_none(x): return None if x.lower() == "none" else x
    def prscsx_linkfun(**kwg): from prstools.PRScsx.PRScsx import main; sys.argv = sys.argv[1:]; main()
    def prscs_linkfun(**kwg): from prstools.PRScs.PRScs import main; sys.argv = sys.argv[1:]; main()        
    xtr = dict(subtype='external')
    ext_lst = [dict(cmdname='prscsx',help="PRS-CSx (original): A cross-population polygenic prediction method with continuous shrinkage "
                    "(CS) priors trained with multiple GWAS summary statistics.",func=prscsx_linkfun, **xtr),
               dict(cmdname='prscs',help="PRS-CS (original): A polygenic prediction method with continuous shrinkage (CS) priors trained with GWAS summary statistics.", func=prscs_linkfun, **xtr)]
    ext_lst = []
    subparserkwg_lst += ext_lst
        
    # Add Subparser structure:
    subparser = parser.add_subparsers(title="Models & Utility Commands", dest="command", metavar='<command>'+' \b'*2, help="") #trick: +' \b'*0
    extcmds = []; prc = process_argkwargs
    for i, spkwg in enumerate(subparserkwg_lst):
        if spkwg['subtype'] == 'BasePred':
            # Create Model parser and add basic help:
            basemodels = ['prscs2','prscs2a']
            
            model_parser = subparser.add_parser(spkwg['cmdname'],
                                                help=spkwg['help'] if spkwg['cmdname'] in basemodels else argparse.SUPPRESS,
                                                description=spkwg['description'],
                                                epilog=spkwg['epilog'],
                                                formatter_class=CustomFormatter,
                                                argument_default=argparse.SUPPRESS,
                                                add_help=False)
            
            if not spkwg['cmdname'] in basemodels:
                subparser._choices_actions = [ca for ca in subparser._choices_actions if ca.dest != spkwg['cmdname']]
            modelgeneral_group = model_parser.add_argument_group('General Options')
            modelgeneral_group.add_argument('-h', '--help', action='help', help='Show this help message and exit.') #default=argparse.SUPPRESS
            # WARNING!: there is still something wrong with the default of this --cpus cli argument, it does not seem to get pushed into the setting of the number of cores function.
            modelgeneral_group.add_argument('--cpus', '-c', **prc(dict(metavar='<number-of-cpus>', default=get_default_cpus(), type=int, 
                        help='The number of CPUs to use (1–5 is generally most efficient). Set to -1 to disable.')))
            
            # Add data related kwargs:
            intcaster = lambda x: int(float(x))
            data_group = model_parser.add_argument_group('Data Arguments (first 5 required)')
            data_group.add_argument("--ref","-r", 
                    **prc(dict(required=True, metavar='<dir/refcode>', 
                        help="Path to the directory that contains the LD reference panel. You can download this reference data "
                               "manually too with \'{basecmd} downloadutil\'. "
                               "Current config relevant to this functionality (prstdatadir: {prstdatadir}, auto_download: {auto_download}).")))
            data_group.add_argument("--target", "-t", type=str_or_none,
                    **prc(dict(required=True, metavar='<bim-prefix>',
                        help="Specify the bim file or its prefix of desired target dataset. The bim file should be in plink format.")))
            data_group.add_argument("--sst","-s", **prc(dict(required=True, metavar='<file>', 
            help="The summary statistics file from which the model will be created. The file should contain columns: SNP, A1, A2, BETA or OR, P or SE information. "
                  "At the moment, the file is assumed to be tab-seperated, if you like other formats please let devs know. "
                 f"Alternative column names can be specified with --colmap (more info below). SNP column should contain rsid\'s, but now these can be filled from  "
                 f"See {format_color('https://tinyurl.com/sstxampl','34')} for a sumstat example.")))
            data_group.add_argument("--out","-o", 
                    **prc(dict(required=True, metavar='<dir+prefix>', 
                        help="Output prefix for the results (variant weights). This should be a combination of the desired output dir + file prefix.")))
            data_group.add_argument("--n_gwas","-n",
                    **prc(dict(required=False, type=intcaster, metavar='<num>', default=None,
                        help="Sample size of the GWAS. Not required if sumstat has a 'N' column and overrules column data if specified.")))
            data_group.add_argument("--chrom", #lambda x: x.split(',')
                    **prc(dict(required=False,type=str, metavar='<chroms>', default='all', 
                        help="Optional: Select specific chromosome to work with. You can specify a specific chromosome as e.g. \"--chrom 3\". All chromosomes are used by default.")))
            data_group.add_argument("--colmap", #lambda x: x.split(',')
                    **prc(dict(type=str, metavar='<colnames>', default='{default_colmap}',# Allows one to specify an alterative column name for columns SNP,A1,A2,BETA,OR,P,SE,N (in that order). "
                        help="Optional: Allows one to specify an alterative column name for the internally used columns snp,A1,A2,beta,or,pval,se_beta,n_eff,af_A1,, (in that order). "
                                    "Forinstance \"--colmap rsid,a1,a2,beta_gwas,,pvalue,beta_standard_error,,,,\" (OR, N, FRQA1, are excluded in this example). "
                                    "When the command is run a quick this_column -> that_column conversion table will be shown. "
                                    "Additionaly prstools has many internal checks to make sure a good PRS will be generated! The original default colmap works with "
                                    "the PRS-CS standard sumstat formatting.")))
            data_group.add_argument("--pred", "-p", 
                    **prc(dict(required=False,metavar='<yes/no>',type=str, default='auto',
                                    help="Optional: Add this argument to set behavior for PRS generation for the induviduals in the target dataset (yes/no/auto).")))

            # Add model-related arguments (hyper parameters and such):
            modelargs_group = model_parser.add_argument_group('Model Arguments (all optional)')
            for argname, item in spkwg['pkwargs'].items():
                #cargs, ckwargs = process_pkwargs(item) << -- underconstruction, need 2 add verbose=True
                modelargs_group.add_argument(*item['args'], **process_argkwargs(item['kwargs']))
            # Below 'func' contains a delayed import of a classmethod that can run the full method from arg
            # This means the import (which can be slow) will take place only when run, leading to big speedups 
            # and a snappy cmdline tool, which is nice for users. In case the model is PRSCS2 it in effect says:
            # ..  .set_defaults(func=PRSCS2.run_from_cli_params_thiswholenamecouldchange_imasocalledclassmethod, .. etc
            func = retrieve_classmethod(clsname=spkwg['clsname'], methodname='from_cli_params_and_run')
            model_parser.set_defaults(func=func, pkwargs=spkwg['pkwargs'], model_parser=model_parser)
        elif spkwg['subtype'] == 'external':
            linked_parser = subparser.add_parser(spkwg['cmdname'], help=spkwg['help'], description='unneeded', add_help=False) #, description='descption prscx')
            linked_parser.set_defaults(func=spkwg['func'])
            extcmds+=[spkwg['cmdname']] 
        elif spkwg['subtype'] == 'PRSTCLI':
            subcmd_parser = subparser.add_parser(spkwg['cmdname'], help=spkwg['help'],
                                    description=spkwg['description'],
                                    epilog=spkwg['epilog'],
                                    formatter_class=CustomFormatter,
                                    argument_default=argparse.SUPPRESS,
                                    add_help=False)
            
            # Add groups and respective arguments:
            for grpname, grpkwg in spkwg['groups'].items():
                cur_group = subcmd_parser.add_argument_group(grpkwg['grpheader'])
                for argname, item in grpkwg['pkwargs'].items():
                    cur_group.add_argument(*item['args'], **process_argkwargs(item['kwargs']))
            
            # Set func to be linked to all the args:
            func = retrieve_classmethod(modulename=spkwg['modulename'], clsname=spkwg['clsname'], methodname='from_cli_params_and_run')
            subcmd_parser.set_defaults(func=func) # Yes it needs to b
        else:
            raise Exception('Subparser subtype not recognized, Contact dev.')
    
    # Commence actual parsing:
    if len(argv)<2 and not 'config' in argv: argv+=['-h']
    knargs, _ = parser.parse_known_args(argv)
    if knargs.command in extcmds:
        return knargs
    else:
        args = parser.parse_args(argv)
        return args

def main(argv=None):
    
    if argv is None: argv=sys.argv[1:]
    args = parse_args(argv)
    display_info = True if 'pkwargs' in args else False
    from prstools import __version__, _date# as version, date
    timestampfmt = "%a, %d %b %Y %H:%M:%S %z"
    if display_info:
        param_dt = vars(args)
        topstr = '\n'.join([
            f'PRSTOOLS v{__version__} ({_date})',
            f'Running command: {args.command}',
            f'Options in effect:'
        ])
        print(topstr)
        for act in  args.model_parser._actions:
            key = act.dest
            if not key in param_dt: continue
            item = param_dt[key]
            fun = lambda: item is args.pkwargs[key]['kwargs'].get('default',None)
            cond0 = fun() if key in args.pkwargs else False
            if cond0 or key in ['func','pkwargs','command','model_parser']: continue
            pitem = item if type(item) is not bool else ''
            print('  --%s %s' % (key, pitem))
        print()
    else: print(f'PRSTOOLS v{__version__} ({_date})')

    args_dt = vars(args)
    # Need to happen here because it needs to happen before numpy/scipy is imported: 
    #**{key:item for key, item in args_dt.items() if key=='cpus'}if 'cpus' in args_dt: 
    cpus = args_dt.get('cpus', get_default_cpus())
    cpu_display_string = set_cpu_envvars(cpus, output_string=True)

    # Initialize logs and grab certain parts:
    from prstools.utils import get_prstlogs; 
    import pandas as pd; import socket
    start = pd.Timestamp.now(); hostname = socket.gethostname(); cwd=os.getcwd()
    if 'seed' in args_dt.get('pkwargs',''):
        if not 'seed' in args_dt:
            args_dt['seed']=int(time.time()) % (2**32)  
        seed = args_dt['seed']
        displayseed = f'Random seed:       {seed}'
    else: displayseed=None
#     print( args_dt.get('pkwargs',''))
    if 'n_jobs' in args_dt.get('pkwargs',''):
#         import prstools as prst
#         prst.utils.get_ip().embed()
        if not 'n_jobs' in args_dt:
            n_jobs = str(args_dt['pkwargs']['n_jobs']['kwargs']['default'])
        else: n_jobs = str(args_dt['n_jobs'])
        displaynjobs = f'Number of workers: {n_jobs} [--n_jobs]'
        cpu_display_string = str(cpu_display_string) + ' (per worker)'
    else: displaynjobs=None
    
    prstlogs = get_prstlogs()
    tic, toc = prstlogs.get_tictoc()
    prstlogs['__version__'] = __version__
    prstlogs['times']['start'] = start
    prstlogs['cwd']      = cwd
    prstlogs['hostname'] = hostname
    prstlogs['argv'] = argv
    #prstlogs['args_dt'] = args_dt 
    if display_info:
        envstr = '\n'.join(elem for elem in [
        f'Hostname:          {hostname}',
        f'Working directory: {cwd}',
        f'Start time:        {start.strftime(timestampfmt)}',
        displayseed,        # skipping this if no seed
        displaynjobs,
        f'Number of cpus:    {cpu_display_string}',
        f' '] if elem)
        print(envstr)
    ''
    # Run the actual task:
    args_dt['return_models'] = False
    result = args.func(**args_dt)
    import prstools as prst
    # Some post main-task things:
    stop = pd.Timestamp.now(); prstlogs['times']['stop'] = stop; prstlogs.finish()
    timestr = prst.utils.get_timestring_from_td(stop-start)
    if display_info: print(f'End time: {stop.strftime(timestampfmt)} -  {timestr}',)

if '_isdevenv_prstools' in locals() or '--dev-secret' in sys.argv:
    
    if 'In' in locals():
        with open('../prstools/_cmd.py', 'w') as f: f.write(In[-1])
        print('Written to:', f.name); 
        
    from prstools import models, utils; import importlib
    importlib.reload(models); importlib.reload(utils)
    try:
        from prstools.models import PRSCS2, MultiPRS
        from prstools.utils import DownloadUtil, store_argparse_dicts, Combine, Config, Transform
        try: from prstools.models._ext import _ext_cli_selection
        except: _ext_cli_selection = []
        extra = [getattr(models,elem) for elem in _ext_cli_selection]
        subparserkwg_lst = [Config, DownloadUtil, Transform, Combine, PRSCS2, MultiPRS] + extra
        store_argparse_dicts(subparserkwg_lst)
        print('Saved new argparse dict. (mind: dont forget the suppress mechanism, this is something in the argparse-dict processing)') 
    except Exception as e: 
        print(e, 'Import issue.'); raise e

    if 'In' in locals():
        get_ipython().system('time python ../prstools/_cmd.py | tail -1 #l2pred # afster a  normal pip install, !time prst seems to have parity with this faster approach.')
else:
    if __name__ == '__main__':
        main()