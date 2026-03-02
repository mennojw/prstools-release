import sys # very basic import

## Base exception
class BaseError(Exception):
    """Base class for all prstools user-facing errors."""
    pass

## Specific semantic / user-facing errors
class SumstatSchemaError(BaseError):
    """Raised when required sumstat columns are missing or misnamed."""
    pass

class LoadError(RuntimeError):
    __slots__ = ('inner', 'stage')
    def __init__(self, inner, *, stage):
        super().__init__(str(inner))
        self.inner = inner
        self.stage = stage

if not '__file__' in locals():
    import sys
    import numpy as np
    if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+['ipykernel_launcher.py' in sys.argv[0]]):
        with open('../prstools/errors.py', 'w') as loadrf: loadrf.write(In[-1])
        print('Written to:', loadrf.name)
        if 'In' in locals() and _isdevenv_prstools:
            print('starting here in models:') 
            get_ipython().system('prst --dev | head -3')
            