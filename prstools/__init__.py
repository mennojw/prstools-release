__version__ = "0.0.68"
_date = "21-04-2026"

import importlib as _importlib

# Existing lazy submodules (UNCHANGED) 
_submodules = [
    "_cmd",
    "io",
    "linkage",
    "models",
    "utils",
    "errors",
]

# NEW: lazy top-level symbol map
# name -> module path
_lazy_exports = {
    # io / loaders-style API
    #"get_build": "prstools.io",
    #"get_rsids": "prstools.io",
    
    # load_*
    "load_bed": "prstools.io",
    "load_bimfam": "prstools.io",
    "load_example": "prstools.io",
    "load_linkagedata": "prstools.io",
    "load_prscs_ldblk": "prstools.io",
    "load_ref": "prstools.io",
    "load_regdef": "prstools.io",
    "load_snpdb": "prstools.io",
    "load_srd": "prstools.io",
    "load_sst": "prstools.io",
    "load_weights": "prstools.io",

    # save_*
    "save_bim": "prstools.io",
    "save_fam": "prstools.io",
    "save_prs": "prstools.io",
    "save_sst": "prstools.io",
    
    # An extra:
	"merge_snps": "prstools.io",
}

__all__ = (
    _submodules
    + list(_lazy_exports.keys())
    + ["__version__"]
)

def __dir__():
    return __all__

def __getattr__(name):
    # OLD BEHAVIOR: lazy submodules
    if name in _submodules:
        return _importlib.import_module(f"prstools.{name}")

    # NEW BEHAVIOR: lazy functions
    if name in _lazy_exports:
        mod = _importlib.import_module(_lazy_exports[name])
        try:
            return getattr(mod, name)
        except AttributeError:
            raise AttributeError(
                f"Module '{_lazy_exports[name]}' does not define '{name}'"
            )

    # Fallback (unchanged)
    try:
        return globals()[name]
    except KeyError:
        raise AttributeError(
            f"Module 'prstools' has no attribute '{name}'"
        )



# _ergegr_version__ = "0.0.58"
# _date = "22-10-2025"
# # from . import *
# 
# # from . import models
# # from . import loaders
# # from . import utils
# # from . import _cmd as cmd
# 
# # import .models
# # from models import __
# # from . import models.L2Pred
# # 
# import importlib as _importlib  # Import takes around 6 microseconds
# 
# # List of submodules to be included
# _submodules = [
#     '_cmd',
#     'loaders',
#     'models',
#     'utils'
# ]
# 
# __all__ = _submodules + [
#     # 'LowLevelCallable',
#     # 'test',
#     # 'show_config',
#     '__version__',
# ]
# 
# def __dir__():
#     return __all__
# 
# def __getattr__(name):
#     if name in _submodules:
#         return _importlib.import_module(f'prstools.{name}')
#     try:
#         return globals()[name]
#     except KeyError:
#         raise AttributeError(f"Module 'prstools' has no attribute '{name}'")
