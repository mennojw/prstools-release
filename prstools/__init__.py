__version__ = "0.0.11"
_date = "18-02-2025"
# from . import *

# from . import models
# from . import loaders
# from . import utils
# from . import _cmd as cmd

# import .models
# from models import __
# from . import models.L2Pred
# 
import importlib as _importlib  # Import takes around 6 microseconds

# List of submodules to be included
_submodules = [
    '_cmd',
    'loaders',
    'models',
    'utils'
]

__all__ = _submodules + [
    # 'LowLevelCallable',
    # 'test',
    # 'show_config',
    '__version__',
]

def __dir__():
    return __all__

def __getattr__(name):
    if name in _submodules:
        return _importlib.import_module(f'prstools.{name}')
    try:
        return globals()[name]
    except KeyError:
        raise AttributeError(f"Module 'prstools' has no attribute '{name}'")
