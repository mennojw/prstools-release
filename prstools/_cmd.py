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

##################################################### BEST IF THIS IS ALREADY DONE BY THE TIME THIS CODE RUNS, BELOW ############
def process_argkwargs(kwargs, setverbosetrue=False):
    kwargs = kwargs.copy() 
    kwargs['help'] = argparse.SUPPRESS if kwargs.get('help', None) is None else kwargs['help']
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
##################################################### BEST IF THIS IS ALREADY DONE BY THE TIME THIS CODE RUNS, ABOVE ############

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

def parse_args(argv=None, description="Convenient and powerfull Polygenic Risk Score creation. \n\'prst\' is a commandline shorthand for \'prstools\'",
               subparserkwg_lst=None, basecmd='prstools', return_spkwg=False, reload=False):
    
    # Prepare parsing params:
    if argv is None: argv=sys.argv[1:]; basecmd=sys.argv[0]
    #from prstools.models import L2Pred, PRSCS
    #subparserkwg_lst = [L2Pred, PRSCS]
    if subparserkwg_lst is None:
        if reload: import importlib; from prstools import _parser_vars; importlib.reload(_parser_vars)
        from prstools._parser_vars import subparserkwg_lst
    else:
        from prstools.utils import process_subparserkwgs
        subparserkwg_lst=process_subparserkwgs(subparserkwg_lst)
        # subparserkwg_lst contains the information to construct parsers
        # This is generated from model code by the developer and saved

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
            basemodels = ['prscs2']
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
                        help='The number of cpus to use. Generally most efficient if chosen to be between 1 and 5. \
                              Functionality can be turned-off completely by setting it to -1.')))
            
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
                        help="The summary statistics file from which the model will be created. The file should contain columns: SNP, A1, A2, BETA or OR, P or SE information. "
                              "At the moment, the file is assumed to be tab-seperated, if you like other formats please let devs know. "
                             f"Alternative column names can be specified with --colmap (more info below). SNP column should contain rsid\'s. "
                             f"See {format_color('https://tinyurl.com/sstxampl','34')} for an example.")))
            data_group.add_argument("--out","--out_dir","-o", 
                    **prc(dict(required=True, metavar='<dir+prefix>', 
                        help="Output prefix for the results (variant weights). This should be a combination of the desired output dir and file prefix.")))
            data_group.add_argument("--n_gwas","-n",
                    **prc(dict(required=False, type=int, metavar='<num>', default=None,
                        help="Sample size of the GWAS. Not required if sumstat has a 'N' column and overrules column data if specified.")))
            data_group.add_argument("--chrom", #lambda x: x.split(',')
                    **prc(dict(required=False,type=str, metavar='<chroms>', default='all', 
                        help="Optional: Select specific chromosome to work with. You can specify a specific chromosome as e.g. \"--chrom 3\". All chromosomes are used by default.")))
            data_group.add_argument("--colmap", #lambda x: x.split(',')
                        type=str, metavar='<colnames>', # Allows one to specify an alterative column name for columns SNP,A1,A2,BETA,OR,P,SE,N (in that order). "
                        help="Optional: Allows one to specify an alterative column name for the internally used columns snp,A1,A2,beta,or,pval,se_beta,n_eff (in that order). "
                                    "Forinstance \"--colmap rsid,a1,a2,beta_gwas,,pvalue,beta_standard_error,\" (OR & N are excluded in this example). "
                                    "When the command is run a quick this_column -> that_column conversion table will be shown. "
                                    "Additionaly prstools has many internal checks to make sure a good PRS will be generated! By default this is the PRS-CS standard. (default: SNP,A1,A2,BETA,OR,P,SE,N)")
            data_group.add_argument("--pred", "-p", **prc(dict(required=False,metavar='<yes/no>',type=str, default='yes',
                                    help="Optional: Add this argument to set behavior for PRS generation for the induviduals in the target dataset (yes/no).")))

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
        elif spkwg['subtype'] == 'function':
            funct_parser = subparser.add_parser(spkwg['cmdname'], help=spkwg['help'],
                                                description=spkwg['description'],
                                                formatter_class=CustomFormatter,
                                                argument_default=argparse.SUPPRESS,
                                                add_help=False)
            general_group = funct_parser.add_argument_group('Options')
            general_group.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
            for argname, item in spkwg['pkwargs'].items():
                general_group.add_argument(*item['args'], **process_argkwargs(item['kwargs']))
            func = globals()[spkwg['cmdname']]
            funct_parser.set_defaults(func=func)
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
    if len(argv)<2: argv+=['-h']
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

if '_isdevenv_prstools' in locals() or '--dev' in sys.argv:
    
    if 'In' in locals():
        with open('../prstools/_cmd.py', 'w') as f: f.write(In[-1])
        print('Written to:', f.name); 
        
    from prstools import models, utils; import importlib
    importlib.reload(models); importlib.reload(utils)
    try:
        from prstools.models import PRSCS2, PredPRS
        from prstools.utils import DownloadUtil, store_argparse_dicts , Combine
        try: from prstools.models._ext import _ext_cli_selection
        except: _ext_cli_selection = []
        extra = [getattr(models,elem) for elem in _ext_cli_selection]
        store_argparse_dicts([DownloadUtil, Combine,
                              PRSCS2 , PredPRS 
                              ] + extra)
        print('Saved new argparse dict. (mind: dont forget the suppress mechanism, this is something in the argparse-dict processing)')
    except Exception as e: 
        print(e, 'Import issue.'); raise e

    if 'In' in locals():
        get_ipython().system('time python ../prstools/_cmd.py | tail -1 #l2pred # afster a  normal pip install, !time prst seems to have parity with this faster approach.')
else:
    if __name__ == '__main__':
        main()