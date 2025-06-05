
from prstools.models._base import *
try:
    from prstools.models._ext import *
except:
    pass

# if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+
#           ['ipykernel_launcher.py' in sys.argv[0]] + 
#           [not '__file__' in locals()]):

#     if 'In' in locals() and _isdevenv_prstools:
#         code = In[-1] 
#         with open('../prstools/models/__init__.py', 'w') as f: f.write(code)
#         print('Written to:', f.name); time.sleep(0.03)
#         !time python ../prstools/_cmd.py --dev  
#         #!prst --dev | head -3
#     print('Done')
