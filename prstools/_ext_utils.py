#| export
# %%writefile ../mjwt/utils.py

# Something python bash python something (test)dur ergerg -- new nb overwrite (oke now another overwrite in nb mjwt)
# Doing something to make all the import suggestions in the code etc work nicely will come later
# this requires writing fancy __init__.py files and having subdirs, _stuff.py files and more..
# I put this here, because I had a hunch it would generate desired python behavior.
# __all__ = ['Tree','fullvars','sizegb','beep','Timer','implot', 'redo_all_above', 'do_all_above', 'Struct'] 

import time, os, pickle, inspect, functools, contextlib, warnings
from sys import getsizeof
from collections import defaultdict

import numpy as np
import scipy as sp
import pandas as pd

from IPython.display import Audio, display, Javascript
from IPython import get_ipython

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
except: warnings.warn('Issue with importing matplotlib, now working around it.')

_rbg_cdict = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),
         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 0.0, 0.0))}
try:
    plt.colormaps.register('rbg', mcolors.LinearSegmentedColormap('rbg', _rbg_cdict, 100))
except:
    True
    
# Perfect for format_map()
class AutoDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'
    
def printfun(arg='ergergreg'):
    print(arg+'  - something')

def suppress_warnings(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = func(*args, **kwargs)
        return result
    return wrapper

# Tree:
def Tree():
    return defaultdict(Tree)

#fullvars:
def fullvars(obj):
    return {key: getattr(obj, key) for key in dir(obj)}

# Size computer:
def sizegb(var, verbose=False):
    size_in_gb = getsizeof(var)/1024**3
    if verbose:
        print('GB: ', size_in_gb)
    return size_in_gb

# Beeper:
def beep(n_beeps=1, seconds=1.,framerate = 4410,v=.8, rate=1, pitch=1., fn='completion_2.mp3'):
    from IPython.display import Audio, clear_output
    t = np.linspace(0, seconds, int(framerate*seconds))
    audio_data = np.sin(2*np.pi*300*t*pitch)+np.cos(2*np.pi*240*t*pitch)
    audio_data = v*audio_data/(audio_data.max()*1.01)
    fn = fn=os.path.join(os.path.dirname(__file__), fn)
    for bcnt in range(n_beeps):
        if os.path.exists(fn):
            a = Audio(filename=fn, rate=10, autoplay=True)
            a.rate=10.
            display(a)
        else:
            display(Audio(audio_data, rate=framerate, autoplay=True, normalize=False))
        time.sleep(seconds+0.25); 
        if bcnt < (n_beeps-2): clear_output()
    #return 

# Timer:
class Timer():
    def __init__(self, fun=None):
        self.time_dt = {}
        if fun is None: fun=lambda : ''
        self.fun = fun
        self.timeepoch_cnt = -1

    def tic(self, dashline=True):
        if dashline:
            print('-'*42)
        self.timeepoch_cnt += 1
        cnt = self.timeepoch_cnt
        self.time_dt[cnt] = [time.time()]
        self.time_lst = self.time_dt[cnt]

    def toc(self, rep=''):
        cnt = self.timeepoch_cnt
        time_lst = self.time_lst
        t = time.time()
        current_line = inspect.currentframe().f_back.f_lineno
        funres = self.fun()
        print(f'Step({len(time_lst)}): {t-time_lst[0]:.3f}, {t-time_lst[-1]:.3f} {funres} <-- line({current_line}) -- {rep}', flush=True)
        time_lst.append(t)
        

def find_varname(var):
    lcls = inspect.stack()[2][0].f_locals
    for name in lcls:
        if id(var) == id(lcls[name]):
            return name
    return None

        
def legimplot(M, show=True, cmap=None, interpolation='none', title=None, figsize=None, **kwargs):
    if title is None:
        #title = find_varname(M) # so so..
        frame = inspect.getouterframes(inspect.currentframe())[1]
        string = inspect.getframeinfo(frame[0]).code_context[0].strip()
        args = string[string.find('(') + 1:-1].split(','); names = []
        for i in args:
            if i.find('=') != -1:
                names.append(i.split('=')[1].strip())
            else:
                names.append(i)    
        title = names[0]
    
    # Determining plot aspect ratios and figsizes:
    aspect = 'auto' if (('aspect' not in kwargs.keys()) & (figsize is not None)) else defaultdict(lambda:None,kwargs)['aspect']
    if figsize is not None: plt.figure(figsize=figsize)
    if aspect == 'auto': title = title + ', aspect=\'auto\''
        
    # Do the plot:
    plt.imshow(M, cmap=cmap, aspect=aspect, interpolation=interpolation); plt.colorbar(); plt.title(title); 
    if show: plt.show() #optional .show supression
    

def implot(*args, cmap='seismic', midpoint=np.median, title=None, figsize=None, f=1.4, 
            ncols=4, nrows=1, vmin=None, vmax=None, show=True, interpolation='none', **kwargs):
    
    # snippets for later dev:
    #fig, axes = plt.subplots(3,3, figsize=[4,4])
    #r = plt.gcf()
    if len(args)==1 and type(args[0]) is dict:
        title=list(args[0].keys()); args=list(args[0].values()); 
    
    # Args call parsing:
    if title is None: # Process the function call for plot titles:
        frame = inspect.getouterframes(inspect.currentframe())[1]
        code = inspect.getframeinfo(frame[0]).code_context[0].strip()
        arr=np.array(list(code)); admin=np.zeros(len(arr))
        admin[arr=='[']= 1; admin[arr=='(']= 1
        admin[arr==']']=-1; admin[arr==')']=-1
        commas=arr==','; commas[admin.cumsum()>1]=0
        for num in np.where(commas)[0]: arr[num] = '$'
        title = map(str.strip,''.join(list(arr[admin.cumsum()>0][1:])).split('$'))      
    if type(title) is str: title = [title]
    titles_dt = defaultdict(lambda : '', {num : titlehere for num, titlehere in enumerate(title)})
    
    # Determining plot aspect ratios and figsizes:
    aspect = None if (('aspect' not in kwargs.keys()) & (figsize is not None)) else defaultdict(lambda:None,kwargs)['aspect']
    #if figsize is not None: plt.figure(figsize=figsize)
    suffix_title = ' aspect=\'auto\'' if aspect == 'auto' else ''
    height = plt.rcParams['figure.figsize'][1]
    if figsize is None: figsize = [f*height*ncols, height]
    
    # Plot all the heatmaps:
    cnt=0
    while cnt < len(args):            
        fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=figsize)
        if not hasattr(axes,'__len__'): axes = [axes]
        for nn, ax in enumerate([*axes]):
            if len(args) > cnt:
                if midpoint is not None:
                    try:
                        vals = args[cnt].flatten()
                    except:
                        vals = args[cnt].values.flatten()
                    vals = vals[~np.isnan(vals)]
                    mid = midpoint if not callable(midpoint) else midpoint(vals)
                    delta  = max(np.abs([min(vals)-mid, max(vals)-mid]))
                    vmin = mid-delta; vmax=mid+delta
                im = ax.imshow(args[cnt], cmap=cmap, aspect=aspect, interpolation=interpolation, vmin=vmin, vmax=vmax); 
                plt.colorbar(im, ax=ax); ax.title.set_text(titles_dt[cnt]+suffix_title)
            else: 
                if show: ax.axis('off')
            cnt += 1
        if show: plt.show()
    
# Redo all above:
def redo_all_above():
    display(Javascript('Jupyter.notebook.kernel.restart(); IPython.notebook.execute_cells_above();'))
    # Best 2 run in middle of browser page, wait and visual view will end up in roughly same place.
    
# Redo all above:
def do_all_above():
    display(Javascript('IPython.notebook.execute_cells_above();'))
    # Best 2 run in middle of browser page, wait and visual view will end up in roughly same place.

class Struct(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)
        
    def __setattr__(self, key, value):
        self[key] = value
        
    def __dir__(self):
        return self.keys()
    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)

# class Struct(dict):
#     def __init__(self, dt=None):
#         """Convert a dictionary to a class
#         @param :adict Dictionary
#         """
#         if dt is not None:
#             self.__dict__.update(dt) 
                
#     def update(self, dt):
#         """Convert a dictionary to a class
#         @param :adict Dictionary
#         """
#         assert type(dt) is dict
#         super().update(dt)
#         self.__dict__.update(dt)
        
def psrc(obj, return_source=False):
    """Print the code of a Python object
    """
    src = inspect.getsource(obj)
    print(src)
    if return_source:
        return src

def jobinfo(job, return_string=False):
    string = get_ipython().system('jobinfo {job.job_id}')
    print(string)
    if return_string:
        return string
    
    
def corr(X, Y=None, ddof=0): # so a correlation might not actually have a ddof.
    if Y is None: Y = X 
    assert type(X) is pd.DataFrame
    assert X.shape[0] == Y.shape[0]
    X = (X-X.mean())/X.std(ddof=ddof)
    Y = (Y-Y.mean())/Y.std(ddof=ddof)
    arr = X.values.T.dot(Y.values)/(X.shape[0]-ddof)
    C = pd.DataFrame(arr, index=X.columns, columns=Y.columns)
    return C

def pcorr(X, Y, ignore_nans=False):

    if np.isnan(X).any()  or np.isnan(Y).any():
        if not ignore_nans:
            raise Exception('Input contains NaNs, this is not allowed')
            
    n = X.shape[0]

    m_X = np.mean(X, axis=0)[np.newaxis,:]
    m_Y = np.mean(Y, axis=0)[np.newaxis,:]
    # sn_X = np.sum( ,axis=0)
    # sn_Y = np.sum( ,axis=0)
    s_X = np.std(X, axis=0)[np.newaxis,:]
    s_Y = np.std(Y, axis=0)[np.newaxis,:]
    s_XY = s_X.T.dot(s_Y)

    R = (X-m_X).T.dot(Y-m_Y)
    R = (R/s_XY)/n

    # Compute P-values matrix:
    dist = sp.stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
    P = 2*dist.cdf(-abs(R))
    
    return R, P

class ResultsStorage():

    def __init__(self, *, timer, suffix='xp', res_base_dn='./results/', verbose=True):

        self.timer=timer
        self.timer.tic(dashline=False) # Start timer.
        self.suffix = suffix
        self.res_base_dn = res_base_dn
        self.ts = pd.Timestamp.now()
        self.ts_str = self.ts.floor('s').isoformat().replace('T','_')
        self.res_dn = os.path.join(res_base_dn, f'{self.ts_str}_{suffix}/')
        self.res_dn = os.path.abspath(os.path.expanduser(self.res_dn))
        if verbose:
            print('ResultStorage dir : ', self.res_dn)
        assert not os.path.exists(self.res_dn)
        os.makedirs(self.res_dn)


    def write_pickle(self, sub_dn='', mode='wb', **kwargs):
        assert len(kwargs) == 1
        key = list(kwargs.keys())[0]
        fn = os.path.join(self.res_dn, sub_dn, key + '.pickle')
        fn = self._prepare_path(fn)
        with open(fn, mode) as f:
            pickle.dump(kwargs[key], f)

    def write_nbformat(self, nb, *, fn, mode='w'):
        fn = self._prepare_path(os.path.join(self.res_dn, fn))
        with open(os.path.join(self.res_dn, fn), mode) as f:
            nbformat.write(nb, f)

    def write_string(self, string, *, fn, mode='w'):
        fn = self._prepare_path(os.path.join(self.res_dn, fn))
        with open(fn, mode) as f:
            f.write(string)

    def _prepare_path(self, fn):
        full_dn = ntpath.dirname(fn)
        os.makedirs(full_dn, exist_ok=True)
        return fn

    def report_completion(self):

        t = self.timer.time_dt[0][0]
        dt = str(time.time()-t) + 's'
        report_string = f"""
        Completed!
        Total-time-taken: {str(pd.Timedelta(dt))}
        """
        self.write_string(report_string, fn='completed.txt')
        
def get_ip():
    import IPython as ip
    return ip

def get_memory_usage(show=True, prefix=''):
    # Checking for Open File Handles !!! 4 prst reflinkagedata issues
    #https://medium.com/brexeng/debugging-and-preventing-memory-errors-in-python-e00be55e7cf2
    import psutil, os #, gc, time
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    gbs=mem_info.rss / (1024 ** 3) # Convert to GBs
    if show:
        print(f"{prefix}Memory usage: {gbs:.3} GB",)
    return gbs

def clear_memory(malloc_trim=True):
    if malloc_trim:
        try:
            import ctypes
            libc = ctypes.CDLL("libc.so.6")
            libc.malloc_trim(0)  # aggressively free heap back to OS #<- this line appears slowest
        except:
            pass
    import gc
    gc.collect()
    
# True
if '_isdevenv_prstools' in locals():
    if _isdevenv_prstools:
        with open('../mjwt/utils.py', 'w') as f: f.write(In[-1])
        print('Written to:', f.name)