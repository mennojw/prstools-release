import os, re, sys, time
import scipy as sp
import numpy as np
import pandas as pd
from scipy import linalg
from sys import getsizeof
import warnings, importlib, json, os, glob, copy, uuid
from collections import OrderedDict, deque, defaultdict
try:
    from pysnptools.standardizer import Unit, UnitTrained
    import pysnptools as pst
except:
    Unit = defaultdict
    UnitTrained = defaultdict
try:
    from tqdm.cli import tqdm
except:
    tqdm = lambda x:x
# import prstools # this is needed below in regdef
import prstools as prst
from prstools.utils import suppress_warnings, _devonlymsg
try:
    import h5py
except:
    True
try:
    class SqrtNinv(Unit):
        def __init__(self):
            super(SqrtNinv, self).__init__()     
except:
    True 
try: import IPython as ip
except: pass

def validate_dataframe_hdf5saving(df):
    # if df contains pyarrow stuff or 'object' dtypes .to_hdf() will start protesting.
    dt = {col : df[col].to_numpy() for col in df.columns} 
    return pd.DataFrame(dt,index=df.index.to_numpy())

def validate_dataframe_index(df, fix=True, drop=True, warn=True, inplace=True, must_exist=None):
    all_ok = isinstance(df.index, pd.RangeIndex) and df.index.start == 0 and df.index.step == 1
    if not inplace: df = df.copy()
    if not all_ok:
        try: 
            assert (df.index == np.arange(len(df.index))).all()
            #if fix: df = df.reset_index(drop=True) # This is a silent fix, since it is actually ok
            all_ok=True
        except: True
    if fix:
        if not all_ok:
            msg='prst dataframe index is being reset, most of the time this is not an issue.'
            if warn:
                warnings.warn(msg)
            df.reset_index(drop=drop, inplace=True)
    else: assert all_ok, msg
    return df

def validate_dataframe_A1A2(df, warn=True, must_exist=False):
    assert not must_exist, 'Not implemented must_exist=True'
    if all(col in df.columns for col in ['A1','A2']):
        head = df.head(200)
        if any(head[col].str.contains('[atcg]').any() for col in ['A1','A2']):
            msg = 'Lower case a,t,c,g letters found in A1 & A2. Casting them to upper-case.'
            if warn: warnings.warn(msg)
        df["A1"] = df["A1"].str.upper()
        df["A2"] = df["A2"].str.upper()
    return df

def validate_dataframe_chrompos(df, warn=True, must_exist=False):
    assert not must_exist, 'Not implemented must_exist=True'
    if 'chrom' in df and not pd.api.types.is_numeric_dtype(df['chrom']):
        cmap = get_chrom_map()
        df["chrom"] = df["chrom"].replace(cmap) # the to_numeric() step can be slow, speedup is involved.
        df['chrom'] = pd.to_numeric(df['chrom'], errors='coerce').astype('Int64') # rrx= bim_df['chrom'].unique() appears fast, so perhaps fix, mod cmap
    if 'pos' in df and not pd.api.types.is_numeric_dtype(df['pos']):
        msg = 'Data in position/bp column is not a number, trying to cast it as a number. This could indicate issues with your inputs e.g. sumstat.'
        if warn: warnings.warn(msg)
        df['pos'] = pd.to_numeric(df['pos'], errors='coerce').astype('Int64') # rrx= bim_df['chrom'].unique() appears fast, so perhaps fix, mod cmap
    return df

def validate_dataframe_n_eff(df, warn=True, must_exist=False):
    if hasattr(df, 'n_eff'):
        #### This logic should not be here, but in diagnostics
        assert (df['n_eff'] > 2).all(), (''
         'Sample sizes (i.e. n_gwas) of smaller than 2 were detected. '
         'This probably means that the N column in the sumstat or --n_gwas option were'
         'not specified correctly.')
    return df
        
def validate_dataframe_select(df, select=None, **kwg):
    msg = 'The argument select was not properly set. This is required. options are ["index","A1A2","chrompos"]. Supply as list.'
    if select is None: raise ValueError(msg)
    badsel = set(select) - set(["index","A1A2","chrompos","n_eff"])
    if len(badsel) > 0: raise ValueError(f'bad selection: {badsel} <--'+msg)
    if 'index'    in select: df = validate_dataframe_index(df,**kwg)
    if 'A1A2'     in select: df = validate_dataframe_A1A2(df,**kwg)
    if 'chrompos' in select: df = validate_dataframe_chrompos(df,**kwg)
    if 'n_eff'    in select: df = validate_dataframe_n_eff(df,**kwg)
    return df

def validate_dataframes_premerge(df0,df1):
    n0 = df0.columns.nlevels
    n1 = df1.columns.nlevels
    if n0 != n1:
        msg = 'Multiindex levels of pandas dataframe levels large 2, this is not supported for merging snps.'
        assert n0 <= 2 and n1 <= 2, msg
        df0.columns = df0.columns if n0==2 else pd.MultiIndex.from_tuples([(col, '') for col in df0.columns])
        df1.columns = df1.columns if n1==2 else pd.MultiIndex.from_tuples([(col, '') for col in df1.columns])
    return df0, df1

def validate_linkage(D, frac=0.001, atol=1e-6, clip=0., return_info=False):
    assert 0 < frac < 1, 'frac needs to be between 0 and 1'
    diag = D.diagonal()
    ndiag = diag*(1+frac)
    np.fill_diagonal(D, ndiag)
    try: L = linalg.cholesky(D, lower=True, check_finite=False)
    except: L = None
    np.fill_diagonal(D, diag)
    if L is None: # cholesky decomp failed, need to fix
        assert np.allclose(D, D.T, atol=atol), 'Input matrix is not symmetric, or contains nans'
        Dm = (D + D.T)/2
        lg, Q = np.linalg.eigh(Dm)
        lg[lg<clip]=clip
        Dm = (Q*lg)@Q.T
    else: Dm = D    
    info = L is None
    out = Dm, info if return_info else Dm
    return out

def get_fn_trimmed(fn, exts=[".gz", ".tsv.gz", ".tsv", '.txt','.csv','csv.gz']):
    for ext in exts:
        if fn.endswith(ext):
            return fn[:-len(ext)]
    return fn

def get_diagnostics(df, verbose=True, allchecks=False):
    assert all(col in df.columns for col in ['A1','A2']), 'Input frame does have to have A1 and A2 columns!'
    if verbose: print(f'{verbose=} and verbose is True by default. And frame has A1 and A2')
    df = get_AX(df)
    issue = False; issue_results = []; n_mlt = -1
    if 'snp' in df.columns:
        if verbose: print('Input frame has "snp" column')
        on = ['snp','AX']
        ind_dups = df.duplicated(subset=on)
        n_dups = ind_dups.sum()
        if n_dups > 0:
            issue = True
            msg = (f'There are duplicates in the input frame using {on=}! {n_dups=}. It could be a simple flipped copy '
                    'in a plink file. returning the rows of the snps that had the offending rows.')
            if verbose: warnings.warn(msg)
            cdf = df[df['snp'].isin(df[ind_dups]['snp'])]
            issue_results += [cdf]
        on = ['snp']
        ind_mlt = df.duplicated(subset=on)
        n_mlt = ind_mlt.sum()
        msg = f'There are multi-allelic snps present,' if n_mlt > 0 else 'All snp-ids are unique,'
        if verbose: print(msg + f'{n_mlt=} {on=}')
        ind = df['snp'].str.startswith('rs')
        perc = ind.mean()*100
        msg = 'For the "snp" column'
        if verbose: print(msg)
    else: 
        if verbose: print('Frame does not have "snp" column. This is odd.')
    
    if ((not issue) or allchecks) and ('pos' in df.columns):
        assert 'chrom' in df.columns, 'Should have chrom if one has pos column!'
        on = ['chrom','pos','AX']
        ind_dups2 = df.duplicated(subset=on)
        n_dups2 = ind_dups2.sum()
        if n_dups2 > 0:
            issue = True
            msg = (f'There are duplicates in the input frame using {on=}! {n_dups2=}. It could be a simple flipped.'
                    f'in a plink file. returning the rows of the snps that had the offending rows.')
            if verbose: warnings.warn(msg)
            cdf = df[df['snp'].isin(df[ind_dups]['snp'])]
            issue_results += [cdf]
        else:
            ind_mlt2 = df.duplicated(subset=['chrom','pos'])
            n_mlt2 = ind_mlt2.sum()
            if n_mlt != -1:
                if n_mlt != n_mlt2: print('There seems to be an issue with the input '
                                          'frame, number of mult-allelic snps not matching')
            elif n_mlt > 0:
                msg = f'There are multi-allelic snps present, {n_mlt=}' 
                if verbose: print(msg)                       
    print('hmm one can look at the letters in A1 and A2, lower vs upper etc, but not implemetned that yet')
    if issue:
        return issue_results
    else:
        if verbose: print('All appears to be ok with input frame.')
        return None

def get_liftoverpositions(df, *, bldout, bldin=None, sort=False, inplace=False, verbose=True, snpdb_df='mini'):
    assert verbose, 'atm can only do liftovers verbose'
    assert not inplace, 'cannot do inplace mod for this yet'
    if bldin is None:
        #assert snpdb_df is not None, 'need to implement get_snpdb funct.'
        if verbose: print('Input build not given so will detect it here.')
        bldin = get_build(df, snpdb_df=snpdb_df, verbose=verbose)
    ## https://genome.ucsc.edu/FAQ/FAQreleases.html#snpConversion UCSC says liftover should not be used for what everybody is using it for.
    from pyliftover import LiftOver
    msg = 'Requiring chrom and pos columns to do position mapping from one genome build to another'
    assert all(col in df.columns for col in ['chrom','pos']), msg
    data_dn = '/PHShome/mw1140/tge/data' # Not used ... :
    fn = os.path.join(data_dn,f"/liftover/hg{int(bldin)}ToHg{int(bldout)}.over.chain.gz") # lo = LiftOver(fn)

    ## The liftover:
    if bldin == bldout:
        print('Input build same as requested output build so no liftover needed. Returning same input.')
        return df
    print('Doing liftover prep.. ', end='', flush=True)
    df = df.copy(); lst = []; nancnt = 0
    input_strings = (f'hg{bldin}', f'hg{bldout}')
    lo = LiftOver(*input_strings) # IT seems this work all on its own
    chroms = df['chrom'].astype(str)
    poss = df['pos']-1
    print('Casting to list to iterate over', flush=True)
    liftlst = list(zip(chroms, poss))
    lifted = [lo.convert_coordinate(f'chr'+c, p) for c, p in prst.utils.get_pbar(liftlst, colour='blue')]

    ## Processing
    print('Done lifting, now postprocessing results')
    for i, lift in prst.utils.get_pbar(list(enumerate(lifted)), colour='yellow'):
        if not lift is None:
            try: res = lift[0]
            except: res = (None,None,None,None); nancnt+=1
        else: res = (None,None,None,None); nancnt+=1
        lst += [res]
    new_df = pd.DataFrame(lst, columns=['newchrom','newpos','strand','something']); df
    df['oldpos'] = df['pos']; df['oldchrom'] = df['chrom']
    if 'strand' in df.columns: df['oldstrand'] = df['strand']
    df['pos'] = new_df['newpos'].astype('Int64').values +1 ### PLUS 1 !!!!!
    df['chrom'] = new_df['newchrom'].values
    cmap = get_chrom_map()
    df['chrom'] = df['chrom'].str.replace("chr","").replace(cmap)
    valid = set(cmap.values()) | set( map(str, range(1,23)))
    df.loc[~df['chrom'].isin(valid), 'chrom'] = pd.NA
    df['chrom'] = df['chrom'].astype('Int64')
    df['strand'] = new_df['strand'].values
    perc = (nancnt/df.shape[0])*100
    if nancnt != 0 :
        msg = f'It seems there are input genomic positions for which no output could be determined, '\
        f'This means thee are NaNs in the chrom and pos columns, #-of-nans = {nancnt} ({perc:.2f}%)\n'\
        f'If this percentage is low (<2%) this is generally workable. FYI: {df[["chrom","pos"]].isna().sum()=}'
        warnings.warn(msg)
    return df

_snpdb_dt = dict()
def get_snpdb(snpdb_df='mini', cache_vars=('mini',)):
    global _snpdb_dt
    if snpdb_df in cache_vars:
        if not snpdb_df in _snpdb_dt: _snpdb_dt[snpdb_df] = load_snpdb(snpdb_df)
    if snpdb_df in _snpdb_dt: snpdb_df = _snpdb_dt[snpdb_df]
    else: snpdb_df=load_snpdb(snpdb_df)
    return snpdb_df

def get_build(df, *, snpdb_df='mini', nsmp = 10_000, avail_builds = [19,38], seed=42*42, verbose=True, on_fail='error'):
    assert snpdb_df is not None, 'No snps database present for build detection! (snpdb_df argument inside of python)'
    if type(snpdb_df) is str:
        if snpdb_df in ('full','mini'):
            msg = 'Trying to infer genome build, but prstdatadir is not set, which is required for this functionality. '
            msg += 'Please run "prst config" to set it.'; prstcfg = prst.utils.load_config()
            if not prstcfg.get('prstdatadir', None): raise RuntimeError(msg)
        snpdb_df = get_snpdb(snpdb_df); 
    pin=df.shape[0]; assert pin>5, 'input snp file too small! less than 5 rows!'
    dbcnt = snpdb_df.shape[0]
    dbfrac = nsmp/dbcnt

    # Sub-sample database to make things faster
    np.random.seed(seed)
    if dbfrac < 0.1: idx=np.unique(np.random.randint(0,snpdb_df.shape[0],nsmp)); nsmp=len(idx)
    else: idx = np.sort(np.random.permutation(np.arange(0,snpdb_df.shape[0]))[:nsmp])
    dbsmp_df = snpdb_df.iloc[idx]

    # Do the things:
    mrg_dt = {}; maxfrac = -1
    for bld in avail_builds:
        arg = df.rename(columns=dict(snp='oldsnp'))
        db = dbsmp_df.rename(columns={f'chrom{bld}': 'chrom', f'pos{bld}':'pos'})
        get_cpnum(arg); get_cpnum(db)
        db = db['cpnum']; db= db[~db.isna()].astype('Int64')
        argc = arg['cpnum']; argc = argc.astype('Int64')
        #mcnt = b['cpnum'].isin(a['cpnum']).sum()
        mcnt = db.isin(argc).sum()
        est_tot_mcnt = mcnt/dbfrac
        estfrac = est_tot_mcnt/dbcnt
        estperc = estfrac*100
        mrg_dt[bld] = dict(mcnt=mcnt, bcnt=(~dbsmp_df[f'chrom{bld}'].isna()).sum())
        if verbose: print(f'For build hg{bld}, probabilistic estimate for total number of matches is: {int(est_tot_mcnt):9,} ({min(estperc,100):5.2f}% of build detection db)')
        if estfrac > maxfrac:
            maxfrac = estfrac; maxperc = estperc
            maxbld = bld; maxmcnt = mcnt
    msg = f'Build Detection Failed: The maximum percentage detected ({min(maxperc,100):5.2f}%) is low, meaning that genome build detection is NOT reliable!'
    if maxperc < 70.0:
        if on_fail == 'error': raise Exception(msg)
        elif on_fail == 'warn': warnings.warn(msg)
        else: raise ValueError(f'Option on_fail = {on_fail} not recognized')
        maxbld = None
    if verbose and maxbld: print(f'Detected build is hg{maxbld}!')
    return maxbld

def get_rsids(df, *, snpdb_df='full', bld='detect', verbose=True, inplace=False, bldsnpdb_df='mini'):
    assert snpdb_df is not None, 'No snps database present! (snpdb_df argument inside of python)'
    if type(snpdb_df) is str: snpdb_df = get_snpdb(snpdb_df);
    if bld == 'detect':
        bld=get_build(df, verbose=verbose, snpdb_df=bldsnpdb_df)
    else: 
        if verbose: print(f'Build is set to be hg{bld}')
    if verbose: print('Merging.. ',end='')
    a = df.rename(columns=dict(snp='oldsnp', aranndomxtest='randotest'))
    b = snpdb_df.rename(columns={f'chrom{bld}': 'chrom', f'pos{bld}':'pos'})
    get_cpnum(a); get_cpnum(b)
    mrg = pd.merge(a, b.drop(['chrom','pos'], axis=1), on='cpnum', how='left', suffixes=('','_right')) # double as fast as the next example line:
    #mrg = pd.merge(df.rename(columns=dict(snp='oldsnp', whoooeg=345)),smp_df.rename(columns={f'chrom{bld}': 'chrom', f'pos{bld}':'pos'}), on=['chrom','pos'], how=how)
    mcnt = (~mrg['snp'].isna()).sum()
    if verbose: print(f'-> Added rsids to the input ({mcnt:,} rsids covering {(mcnt/mrg.shape[0])*100:.2f}% of total input rows)')
    return mrg

def get_AX(df, opt=1, redo=False):
    df = validate_dataframe_index(df)
    assert all(A in df.columns for A in ['A1','A2']), 'A1 and/or A2 columns are missing, these are required here.'
    if 'AX' in df and not redo: return df # Dont redo AX if alredy done
    vara = (df['A1'] + '_') + df['A2'] # 35%
    varb = (df['A2'] + '_') + df['A1'] # 35%
    ind = df['A1'] <= df['A2']
    df['AX'] = varb
    df.loc[ind,'AX'] = vara[ind].values #(takes 25%)
    #df['AX'][ind] = vara #(takes 18%, but give annoying warnings
    return df 

def get_chrompos(df, inplace=True, redo=False):
    assert inplace
    if (not 'chrompos' in df) or redo:
        df['chrompos'] = df['chrom'].astype('Int64').astype(str) + '_' + df['pos'].astype('Int64').astype(str)
    return df

def get_chromposAX(df, inplace=True, redo=False):
    if not inplace: df=df.copy()
    if 'chromposAX' in df and not redo: return df
    get_AX(df)
    get_chrompos(df)
    df['chromposAX']=df['chrompos']+'_'+df['AX']
    return df

def get_cpnum(df, inplace=True):
    assert inplace
    df['cpnum'] = df['pos'] + df['chrom']*int(1e10) # This is fast and will give a proper int
    return df

def get_rsnum(df, inplace=True):
    ind = df['snp'].str.startswith('rs')
    df.loc[ind,'rsnum'] = df.loc[ind,'snp'].str[2:].astype('Int64')
    return df

def get_regid(prst_df, regdef=None, fixchromends=True):
    if regdef is None: regdef_df = load_regdef()
    elif type(regdef) is str: regdef_df = load_regdef(regdef) # This not functional yet
    elif type(regdef) is pd.DataFrame: regdef_df = regdef
    else: raise TypeError(f'type={type(regdef)} not recognized. Type needs to be string or dataframe.')
    prst_df = validate_dataframe_index(prst_df)
    # regdef_df = validate_regdef_df(regdef_df)
    # Loop through chromosomes
    for chrom, cur_df in prst_df.groupby('chrom'):
        cdef_df = regdef_df[regdef_df.chrom.astype(type(chrom))==chrom].copy()
        if cdef_df.empty: continue # Jump to next if it is empty
        # Iterate over regions and assign regid
        if fixchromends: ## since regdef is asserted to be sorted the following is ok
            cdef_df.loc[cdef_df.index[0], 'start'] = 0 #=ok, see 1 line above
            cdef_df.loc[cdef_df.index[-1], 'stop'] = 1e12
        for _, row in cdef_df.iterrows():
            mask = (cur_df['pos'] >= row['start']) & (cur_df['pos'] < row['stop'])
            prst_df.loc[mask.index[mask], 'regid'] = int(row['regid'])
    indnans = prst_df['regid'].isna(); nansum = indnans.sum()
    if nansum > 0: 
        perc = (nansum/len(indnans))*100
        msg=f'Certain variants could not be assigned to a region, nansum = {nansum:,} ({min(perc,100):5.2f}% of total).'
        warnings.warn(msg)
    return prst_df

def get_cols(df, *, addcols, inplace=None, redo=False):
    # This is a columns creations aggregration function
    assert not inplace, 'Not implemented, inplace=True'
    assert redo is False, 'Not implemented, redo=True'
    addcols_map = dict(
        AX=get_AX,
        rsids=get_rsids,
        chrompos=get_chrompos,
        chromposAX=get_chromposAX,
        cpnum=get_cpnum,
        rsnum=get_rsnum,
        regid=get_regid,
    )
    unknown = set(addcols) - addcols_map.keys()
    msg = f"Unknown addcols {unknown} not in available options ({set(addcols_map.keys())})"
    if unknown: raise ValueError(msg)
    if not inplace and type(inplace) is bool: df=df.copy()
    for col in addcols:
        fun=addcols_map[col]
        df=fun(df)
    return df

def get_expansion_region_assignment(prst_df, *, regcol):
    throwerror() #  I dont know if i need this functionality I think not...
    for _, cur_df in prst_df.groupby('chrom'):
        # Create a mapping from the unique values to integer codes
        uniq = cur_df[regcol].dropna().unique()
        uniq = [elem for elem in uniq if not pd.isna(elem)] # but make sure they not nans
        val2int = {v: i for i, v in enumerate(uniq)}
        int2val = {i: v for v, i in val2int.items()}
        # Interpolate using the integer representation
        interp = cur_df[regcol].map(val2int).interpolate(method='nearest')
        interp = interp.ffill().bfill()[cur_df.index] # unfortunately bfill and ffill are needed too (chrom edges).. bad design.
        endvals = interp.map(int2val).values
        assert np.isnan(endvals).sum() == 0
        prst_df.loc[cur_df.index, regcol] = endvals
    assert np.isnan(prst_df[regcol]).sum() == 0
    return prst_df

def get_meta_analysis(*dfs, stuff=True):

    for df in dfs:
        1
    res = get_build(df)
    print(res)

    meta_df = 0
    return meta_df

def get_pyarrow_prw(delimiter=None, pyarrow=True):
    if pyarrow:
        try: import pyarrow as pyarrowpack # Prevent var overloading
        except: 
            msg = 'It seems \'pyarrow\' is not installed, which means you are missing out on a' + \
                ' lot of speed. Consider installing it using \'pip install pyarrow\' for super fast data loading.'
            warnings.warn(msg)
            pyarrow = False
        if delimiter == r'\s+': pyarrow=False
        try:
            from io import StringIO
            pd.read_csv(
                StringIO("a,b\n1,2"),
                dtype_backend="pyarrow",
                engine="pyarrow")
        except:
            pyarrow=False
    if pyarrow:
        return dict(dtype_backend="pyarrow", engine='pyarrow')
    else: return dict()

def _get_default_colmap():
    prstcfg = prst.utils.load_config()
    default_colmap = prstcfg['default_colmap']
    #assert default_colmap == 'SNP,A1,A2,BETA,OR,P,SE,N,FRQA1'
    return default_colmap

def _get_default_conv_dt():
    _quicklook_conv_dt = {'CHR':'chrom','SNP':'snp', 'BP':'pos', 'A1':'A1','A2':'A2', 'MAF':'maf',
              'BETA':'beta','SE':'se_beta','P':'pval','OR':'oddsratio', 'N':'n_eff', 'FRQA1': 'af_A1'}
    prstcfg = prst.utils.load_config()
    default_conv_dt = prstcfg['default_conv_dt']
    return default_conv_dt

def get_colmap(*, schema):
    default = _get_default_colmap()
    colmap_dt = dict(
        default = default,
    )
    if type(schema) is int:
        lst = list(colmap_dt.values())
        colmap = None if len(lst[schema:]) == 0 else lst[schema]
    else:
        msg =f'The colmap found for in the colmap collection for \'{schema}\'. Options are {colmap_dt.keys()}'
        assert schema in colmap_dt, msg
        colmap = colmap_dt[schema]
    return colmap

def get_conv_dt(*, flow, colmap=None, verbose=False):
    assert flow in ('in','out')
    default_conv_dt = _get_default_conv_dt()
    if colmap is not None and flow=='in':
        if type(colmap) is dict:
            return colmap
        elif type(colmap) is str:
            defcolmap = get_colmap(schema='default')
            base=defcolmap.split(','); 
            lst=colmap.split(',')
            msg=f'current_colmap={colmap} default_colmap={defcolmap}'
            assert len(base) == len(lst), f'The colmap should have the right number of renames/commas. {msg}'
            preconv_dt = {key: item for key, item in zip(base,lst)}
        else: raise Exception('colmap is not string or dict, it should be.')
        notdefault = True if (colmap != defcolmap) else False
        lst = []
        conv_dt = {key:item for key,item in default_conv_dt.items()}
        for key, item in preconv_dt.items():
            if item != '':
                intcol = conv_dt.pop(key)
                conv_dt[item] = intcol
                if ((key != item) or notdefault) and verbose:
                    lst += [f' {item} -> {intcol}']
            else: conv_dt.pop(key)
        if verbose: print(f'[colmap column-name conversions: {",".join(lst)} ][The outputs are the internal prstools colnames]', end='\n')
    else: conv_dt = default_conv_dt
    if flow == 'out':
        conv_dt = {item:key for key,item in default_conv_dt.items()}
    return conv_dt

def get_chrom_lst(chrom):
    if type(chrom) is str and '*' in chrom and chrom != '*':
        raise NotImplementedError(f'Cannot work with chrom={chrom} yet, this is not implemented.')
    if type(chrom) is str and chrom in ['*','all']: chrom = list(range(26))
    if type(chrom) is str:
        cmap = get_chrom_map() # This cmap allows chrom=X or MT to be mapped to their number.
        chrom = [int(cmap.get(elem, elem)) for elem in chrom.split(',')]
    if type(chrom) is int or type(chrom) is float:
        chrom = [int(chrom)]
    chrom = list(chrom) # If it came in as an array this will not fail, otherwise it will.
    return chrom

def get_chrom_map(flow='in', version='onlyonenow'):
    assert flow in ('in','out')
    if flow == 'out': raise NotImplementedError('ask/contact dev')

    #https://zzz.bwh.harvard.edu/plink/data.shtml is where I lifted the map from.
    chrom_map = {
        'X': str(23), # using string since the fast pandas .replace() function expects strings.
        'Y': str(24),
        'XY': str(25),
        'MT': str(26)
    }

    return chrom_map

def merge_snps(df0, df1, *, flipcols, afcols=[], how='left', on=['snp','AX'], reset_index=True, extradropdupcols=False, dropalldupcols=False,
               dropduprightcols=['chrom','snp','cm','pos','A1','A2','maf_ref','std_ref','AX'], warndupcol=True, removedups=True,
               seperate=False, handle_missing=False, allow_right_filter=True,
               req_all_right=None # or all entries in the right dataframe being matchable to left 
              ):
    ## Checks & Validation:
    if req_all_right is None: # (e.g in the case of allele_weight 's which cannot just be dropped')
        req_all_right = True if 'allele_weight' in ([*df0.columns] + [*df1.columns]) else False
    if extradropdupcols: dropduprightcols = list(set(dropduprightcols + extradropdupcols))

    ## THE SPOT WHERE I SHOULD TRY AGGRESSIVE FILTERING OF RIGHT FOR MERGE SPEEDUPS
    if allow_right_filter:
        if df0.shape[0] < df1.shape[0]*0.5 and req_all_right is False and 'snp' in on:
            #empirical relatation, if right is more than 2times as large, slicing makes sense:
            ind = df1['snp'].isin(df0['snp']); df1 = df1[ind].reset_index(drop=True)

    df0, df1 = [prst.io.validate_dataframe_index(df) for df in [df0,df1]]
    if 'cpnum' in on:
        bld0=get_build(df0); bld1=get_build(df1)

    if 'AX' in on:
        if not 'AX' in df0.columns: df0 = get_AX(df0)
        if not 'AX' in df1.columns: df1 = get_AX(df1)
    suffixes = ('','_right')
    if type(on) is str: on = [on]
    assert type(on) is list
    assert not seperate, 'Seperate option not available. Currently the two snps frames will always be matched.'
    assert how == 'left', 'how=left currently the only option.'
    assert all(col in df1.columns for col in flipcols)
    flipset = set(flipcols) & set(df0.columns) # for col in flipcols
    premsg = f'The frames about to be merged have both {flipset}. This is not allowed'
    assert not len(flipset) > 0, premsg+f'because there can only be one {flipset}. You can rename them instead'
    lst = [col for col in df1.columns if (col in ['beta_mrg','beta','allele_weight']) and not (col in flipcols)]
    if len(lst) > 0: warnings.warn(f'Columns present in right/2nd input that need to be flipped but wont be.(={lst})')
    assert all(col in df1.columns for col in afcols); afset = set(afcols) & set(df0.columns)
    assert not len(afset) > 0, f'The frames about to be merged have both {afset}. This is not allowed because there can only be one {afset}.'
    #lst = [col for col in df1.columns if (col.startswith('af_')) and not (col in afcols)]
    #if len(lst) > 0: warnings.warn(f'Allele-freq Columns present in right/2nd input that probably need to be modified but wont be.(={lst})')
    handlenotclimsg = ('handle_missing option [keep,filter,False] inside of python function'
        ' can be used to handle this differently. One should not see this message inside of the prstools cli, if you do contact dev.')

    ## Merging:
    df0, df1 = prst.io.validate_dataframes_premerge(df0.copy(),df1.copy())
    with warnings.catch_warnings(record=True): # Suppress annyoing warnings for multiindex case.
        mrg_df = pd.merge(df0, df1, on=on, how=how, suffixes=suffixes)
    indnans = mrg_df['A1_right'].isna()
    n_missing_right = indnans.sum()
    if not handle_missing: assert n_missing_right == 0, ('Not all variants can be matched.'+handlenotclimsg)
    elif handle_missing=='keep':
        pass
    elif handle_missing=='filter':
        mrg_df = mrg_df[~indnans] # remove the NaNs out of the dataframe.
    else:
        raise Exception(f'[keep,filter,False] are the options for handle_missing, now it is {handle_missing}')

    if req_all_right:
        issues=False
        cnts = (~mrg_df[flipcols].isna()).sum().unique()
        if len(cnts) > 1: issues=True
        if cnts[0] != df1.shape[0]: issues=True
        if issues:
            # Create sets of key tuples from df1 and from the merged df.
            df1_keys = set(df1[on].itertuples(index=False, name=None))
            merged_keys = set(mrg_df[on].dropna().itertuples(index=False, name=None))
            missing_keys = list(df1_keys - merged_keys)
            n_miss = len(missing_keys)
            if missing_keys:
                msg = f'mind req_all_right = {req_all_right}, if you want to drop variants weights '+\
                '(typically not a good thing to do) set req_all_right=False ' if req_all_right else ''
                raise ValueError(f"Not all SNPs from df1 were matched. Missing keys (n={n_miss}): {missing_keys[:5]}.. ")

    #ind_match = (mrg_df['A1'] == mrg_df['A1_right']).fillna(True) & (mrg_df['A2'] == mrg_df['A2_right']).fillna(True)
    #ind_flip  = (mrg_df['A1'] == mrg_df['A2_right']).fillna(True) & (mrg_df['A2'] == mrg_df['A1_right']).fillna(True)
    ind_match = (mrg_df['A1'] == mrg_df['A1_right']) & (mrg_df['A2'] == mrg_df['A2_right'])
    ind_flip  = (mrg_df['A1'] == mrg_df['A2_right']) & (mrg_df['A2'] == mrg_df['A1_right'])
    ind_wrong = ~ind_match & ~ind_flip
    if ind_wrong.sum() != 0:
        cnt = ind_wrong.sum()
        msg = (f'Probable tri-allelic snps detected (n={cnt}) and handled appropriately.')
        if handle_missing=='filter': 
            warnings.warn(msg)
            mrg_df = mrg_df[~ind_wrong]
            ind_match=ind_match[~ind_wrong]; ind_flip=ind_flip[~ind_wrong]
        else: raise Exception('with handle_mssing=filter, this can be fixed.' + handlenotclimsg)

    ## Flipping:
    cast=float #cast='int64[pyarrow]' # used to be int(), but now better a type with nans like float
    mrg_df['rflip'] = -1*ind_flip.astype(cast) + 1*ind_match.astype(cast)
    indfunny = mrg_df['rflip'] == 0 #).sum() ==0, 'regerergre'
    assert indfunny.sum() == 0, 'snp alignment issue, that should not happen, please contact dev if it does.'
    #mrg_df.loc[mrg_df['rflip'] == 0, 'rflip'] = 15
    for col in flipcols: 
        #assert (~mrg_df[indfunny][col].isna()).sum() == 0, 
        #mrg_df[col] = (mrg_df[[col]] * mrg_df[['rflip']].values)[col] # a formulations that plays well with multicols, but issue with 1dcase
        #### You should use | as seperator for multicolumn mechanics.
        sel = mrg_df[[col]].values; assert sel.shape[1] == 1, 'It seems a multi-column was put into this function, it is currently not designed for that, contact dev.' 
        mrg_df[[col]] = (sel * mrg_df[['rflip']].values) # formulations I used later, should play well with both.
    with warnings.catch_warnings(record=True): mrg_df.flipcols = flipcols
    for col in afcols: # Cols are prresent, asserted earlier
        mrg_df[f'{col}_ori'] = mrg_df[col].copy()
        ind = (mrg_df['rflip'] == -1)
        mrg_df.loc[ind, col] = 1-mrg_df[f'{col}_ori'][ind]
    with warnings.catch_warnings(record=True): mrg_df.afcols = afcols

    ## Post steps:
    # The following drops columns in the right/2nd input like chrom and pos.
    with warnings.catch_warnings(record=True): # Suppress annyoing warnings for multiindex case.
        if dropalldupcols: mrg_df = mrg_df.drop(
            [col for col in mrg_df.columns.get_level_values(0) if col.endswith(suffixes[1])], axis=1, errors="ignore")
        mrg_df = mrg_df.drop([f'{col}{suffixes[1]}' for col in dropduprightcols], axis=1, errors="ignore")
    # If there are still duplicate columns after this. this will be reported as a warning
    duplicate_columns = [col for col in mrg_df.columns if '_right' == str(col)[-6:]]
    if len(duplicate_columns) > 0 and warndupcol:
        warnings.warn(f'Columns present in right/2nd input that are also present in the first.'
                      f'\nThe duplicate columns are: {duplicate_columns}, (remove the {suffixes[1]} suffix)')

    if mrg_df.shape[0] != mrg_df[on[0]].nunique():
        onhack = mrg_df[on].columns # ohh pandas sometimes, you remind of actual pandas...
        new_df = mrg_df.drop_duplicates(subset=onhack, keep='first')
        #ip.embed()
        if new_df.shape[0] < mrg_df.shape[0]:
            n_dups = mrg_df.shape[0] - new_df.shape[0]
            inject = ' but not removed '
            if removedups: mrg_df = new_df; inject=' and removed '
            warnings.warn(f'Duplicates detected{inject}in sumstat n_dups={n_dups}.')

    if reset_index: mrg_df = prst.io.validate_dataframe_index(mrg_df, warn=False) 
    if seperate: raise NotImplementedError()
    else: return mrg_df

def check_alignment_snps(*args, dropsnps=False, on=['snp','A1','A2']):
    assert len(args) >= 2, 'At least 2 arguments need for this function'
    for col in on:
        arg0 = args[0]
        for argx in args[1:]:
            assert (arg0[col] == argx[col]).all(), f'Column \'{col}\' not matching for input arguments.'

def naninfslicer_funct(start_df, cols, inf=False, verbose=False, ispretest=False):
    df = start_df.copy()
    startlen = df.shape[0]
    df = df.dropna(subset=cols)
    endlen = df.shape[0]
    numofnans = startlen-endlen
    if verbose and numofnans>0 : print(f'NaN values found in sumstat somewhere in these columns; {cols}: '
        f'{numofnans} SNPs removed from the starting total of {startlen:,} ({100*numofnans/startlen:.1f}%), {endlen:,} SNPs left.')
    if inf:
        startlen = df.shape[0]
        ind = (df[cols].abs() > 1e99).any(axis=1)
        numofinfs = ind.sum()
        if numofinfs: df = df[~ind]
        endlen = df.shape[0]
        msg = (f'Inf values found in sumstat somewhere in these columns; {cols}: {numofinfs} SNPs removed from the starting '
              f'total of {startlen:,} ({100*numofnans/startlen:.1f}%), {endlen:,} SNPs left. (e.g. possible reason: standard-error=0 effect/SE=inf).')
        if verbose and numofinfs: print(msg)
    if not ispretest and endlen < 5: 
        cprint_input_df(start_df, show_dims=True)
        raise Exception('Less than 5 SNPs left after processing, something was wrong with the input sumstat, which could mean it will need hands-on processing.')
    return df

def compute_pvalbetase(df, *, calc_lst=['pval','beta','se_beta'], pvalmin=1e-323, copy=True, pretest=False, slicenaninfs=False, verbose=False):
    assert slicenaninfs == False
    if copy: df=df.copy()
    if 'pval' in calc_lst: df['pval'] = 2*sp.stats.norm.cdf(-abs(df.beta_mrg)*np.sqrt(df.n_eff))
    if 'beta' in calc_lst: df['beta'] = df['beta_mrg']/df['std_sst'] # assumption var[y] ==1 !
    if 'se_beta' in calc_lst:   df['se_beta']   = 1/(np.sqrt(df.n_eff)*df['std_sst'])# assumption var[y] ==1
    if pvalmin:
        df.loc[df.pval <= pvalmin,'pval'] = pvalmin
        warnings.warn(f'Small values were detected in the p-values, padding with small non-zero values (={pvalmin}).'
                      ' This likely leads to suboptimal performance. Use beta and se sumstat columns instead for better performance')
        assert np.sum(df.pval < pvalmin) == 0
    return df

def cprint_input_df(df, prefix='\nERROR WITH INPUT (reason at the end) -', show_dims=False, iloc=[0,1,-1]):
    if df.shape[0] < len(iloc): iloc=df.index
    print(f'{prefix} This is what the currently loaded top-rows of dataframe/sumstat looks like after '
          'colmap\'ing (using --colmap, if supplied)(frame is transposed, to make it easier to view):\n', df.iloc[iloc].T)
    if show_dims: print(f'dims: {df.shape}')
    #print(f'All the column names in this dataframe are: {df.columns}')

def check_reqcols_sst(orisst_df, *, reqcols, colmap=None,
                  errfmt='{prefix} Missing required column(s) {missing_cols} (alternative name(s): {alias}){postfix}',
                  prefix='', postfix = ', please add the column(s) to the sumstat or use --colmap option. '
                  "\033[1;38;2;179;125;19mLikely fix: Paste this whole output into a chatbot to get the right --colmap option.\033[0m",
                  allow_dups=False, iloc=[0,1,-1]):
    
    # Prep logic
    from prstools.errors import SumstatSchemaError 
    conv_dt, inv_dt = (get_conv_dt(flow=flow, colmap=colmap, verbose=False) for flow in ['in','out'])
    msg=('It seems the input frame for the sumstat has only one column. This is probably because the '
         'file was parsed with the wrong delimiter. atm prstools does not yet allow you to spec the delimiter.')
    assert orisst_df.shape[1] > 1, msg
    sst_df = orisst_df.rename(columns=conv_dt)
    if all(elem in sst_df.columns for elem in ['chrom','pos']):
        reqcols = [col for col in reqcols if not col == 'snp']

    # Determine missing cols
    missing_cols = []; alias = []
    for colgrp in reqcols:
        if type(colgrp) is str: colgrp=(colgrp,)
        if not any(col in sst_df.columns for col in colgrp):
            missing_cols += [' or '.join([f'\'{col}\'' for col in colgrp])]
            alias += [' or '.join([f'\'{inv_dt.get(col)}\'' for col in colgrp])]

    if len(missing_cols) != 0: # <-- This bit is all to generate a good error message.
        missing_cols = '/'.join(missing_cols); alias = '/'.join(alias)
        xtra = f'(The default --colmap is {get_colmap(schema="default")})' if not colmap is None else '(This is the default)'
        if colmap is None: colmap = get_colmap(schema='default')
        doublecols = pd.MultiIndex.from_arrays([orisst_df.columns, sst_df.columns], names=['original_columns', 'current_columns'])
        overview_df = sst_df.iloc[iloc]
        overview_df.columns = doublecols
        overview_df.index = 'row' + overview_df.index.to_series().astype(str)
        mapper = dict(SNP='rsid', A1='EffectAllele', A2='OtherAllele', BETA='BETA', P='Pval', SE='StdErr', N='Ntotal', OR='')
        example_colmap = prst.io.get_colmap(schema='default')
        for k, v in mapper.items(): example_colmap = example_colmap.replace(k, v)
        print(f' Current --colmap is {colmap} {xtra}. Current colmap leads to following mapping of columns:')
        prst.io.get_conv_dt(flow='in', colmap=colmap+' ', verbose=True)
        print(f"The --colmap argument should list the column names as they appear in your input file."
              f" For instance --colmap {example_colmap}. Mind that not all positions need to have a "
              "column name and can be left empty. With this example colmap we will get the following column mapping:")
        prst.io.get_conv_dt(flow='in', colmap=example_colmap, verbose=True)
        print('Things to consider: if the conversion column name is not present in the input file it will just be skipped. It does need the exact right number of commas though. '
              'The column names in --colmap should be selected from "original_columns" below. "current_columns" contains column names after applying the currently used --colmap. '
              'Note that SNP should mainly contain rsids as SNP id\'s, since that is what the references use. Consider columns selected with colmap will be used for filtering for NAs. '
              'So if there is column that seems to be filled with NAs it might be better to not use it.'
              'Following the plink convention, A1 refers to the effect allele (BETA). '
              'Think for a second and make sure that the effect allele in the sumstat is mapped to A1. If you are a chatbot think in a structured way first about all the colmap elements. '
              'If proper end in a clear way with "--colmap thisistherightcorrectcolmap" codeblock '
              '(give a colmap=\'stuff\' suggestion only if run inside of ipynb.)')
        #print(f'This results in the following colmapping dictionary {colnameconv_dt}, which was already applied to the following dataframe.')
        cprint_input_df(overview_df,iloc=iloc); print('\n')
#         err = SumstatSchemaError(errfmt.format(prefix=prefix, missing_cols=missing_cols, alias=alias, postfix=postfix))
#         dergger
        err = SumstatSchemaError(errfmt.format(prefix=prefix, missing_cols=missing_cols, alias=alias, postfix=postfix))
        raise err
        
    if not allow_dups and sst_df.columns.duplicated().sum() > 0:
        dupcols = sst_df.columns[sst_df.columns.duplicated()]
        cprint_input_df(sst_df)
        raise Exception(f'Columns {dupcols} are duplicates, which makes it unclear which of these columns to select. Please remove these columns.')

def compute_beta_mrg(df, *, calc_beta_mrg=True, n_eff_handling='topmedian', copy=True, ispretest=False, slicenaninfs=True, verbose=False, cli=False):
    # --- df is the sst_df
    if copy: df = df.copy(); 
    
    if calc_beta_mrg:
        cols = df.columns
        if ('oddsratio' in cols) and not ('beta' in cols):
            df.loc[:,'beta'] = np.log(df['oddsratio'])
            cols = df.columns
            if verbose: print('Detected odds ratio, converting to good proxy for beta.')
        testcols = [elem for elem in ['se_beta','beta','n_eff','pval'] if elem in cols]
        if 'n_eff' in df: # <-- n_eff handling
            if n_eff_handling == 'topmedian':
                k=int(len(df['n_eff']) * 0.05)
                n_eff = np.median(np.partition(df['n_eff'], -k)[-k:]); n_eff_msg=n_eff
            elif n_eff_handling == 'raw':
                n_eff = df['n_eff']; n_eff_msg = np.median(n_eff)
            else: raise ValueError(f'Option not recog: {n_eff_handling}. Options are "topmedian" and "raw"')
        if calc_beta_mrg == 'se' or {'se_beta','beta','n_eff'}.issubset(cols):
            if slicenaninfs: df=naninfslicer_funct(df, testcols, verbose=verbose, ispretest=ispretest)
            pre_std_sst = 1. / (np.sqrt(n_eff+1) * df['se_beta']); k=int(len(pre_std_sst) * 0.025)
            #pre_std_sst = 1. / np.sqrt((n_eff+1) * df['se_beta']**2); k=int(len(pre_std_sst) * 0.025) #old
            #df.loc[:,'beta_mrg'] = df.beta/np.sqrt((n_eff+1)*df.se_beta**2) # old
            df.loc[:,'beta_mrg'] = df.beta*pre_std_sst # There is a funny order here 4 speed
            std_y = np.sqrt(0.5)/np.median(np.partition(pre_std_sst, -k)[-k:])
            df['std_sst'] = std_y * pre_std_sst
            df.std_y = std_y # Saving it here incase its needed later on at some point.
            msg = f'Computed beta marginal (=X\'y/n) from sumstat using beta and its standard error and sample size (n_eff={int(n_eff_msg)}).'; df.msg=msg
            if verbose and not cli: print(msg)
        elif calc_beta_mrg == 'pval' or {'pval','beta','n_eff'}.issubset(cols):
            if slicenaninfs: df=naninfslicer_funct(df, testcols, verbose=verbose, ispretest=ispretest)
            if np.sum(df.pval == 0):
                pvalpadder=1e-323; df.loc[df.pval == 0,'pval'] = pvalpadder
                warnings.warn(f'Zero(s) were detected in the p-values, padding with smallest non-zero values (={pvalpadder}),'
                              ' which can lead to suboptimal performance. Use beta and se sumstat columns instead for better performance')
                assert np.sum(df.pval == 0) == 0
            df.loc[:,'beta_mrg'] = np.sign(df.beta)*np.abs(sp.stats.norm.ppf(df.pval/2.0))/np.sqrt(n_eff) 
            msg = f'Computed beta marginal (=X\'y/n) from sumstat using p-values and the sign of beta and sample size (n_eff={int(n_eff_msg)}).'; df.msg=msg
            if verbose and not cli : print(msg)
        else:
            cprint_input_df(df); cnames =' or '.join([f'\'{elem}\'' for elem in 'SE/se_beta/P/pval'.split('/')])
            raise Exception(f'Missing a required column (named: {cnames}) needed for beta marginal/zscore computations')
        if slicenaninfs: df=naninfslicer_funct(df, cols=['beta_mrg'], inf=True, verbose=verbose, ispretest=ispretest)
    return df

def _pd_read_csv(*args, max_arrow_tries=2, arrow_sleep=0.5, **kwg):
    try:
        for i in range(max_arrow_tries):
            try:
                return pd.read_csv(*args,**kwg)
            except Exception as e:
                try: import pyarrow.lib
                except ImportError: pyarrow = None
                # Inline: decide whether this is a retryable Arrow CSV issue
                cond0 = pyarrow is not None and isinstance(e, pyarrow.lib.ArrowInvalid)
                cond1 = "two block boundaries" in str(e) or "CSV parse error" in str(e)
                retryable = cond0 and cond1
                # Not retryable → propagate immediately (user error etc.)
                if not retryable: raise e
                if len(args) > 0: fn = args[0]
                else: fn = 'unknown-file'
                warnings.warn(f"{fn}: pyarrow CSV parser issue detected, retrying.")
                time.sleep(arrow_sleep)
        msg = f"{fn}: pyarrow CSV parser failed after {max_arrow_tries} attempts; falling back to pandas C engine"
        warnings.warn(msg)
        kwargs['engine']='c'
        return pd.read_csv(*args, **kwg)
    except Exception as e:
        raise prst.errors.LoadError(e, stage=1) from e

def _validate_kwg_load_fun(fn, *, load_fun, ukwg_lst=None, **kwg): ## keep close to load_sst
    
    # Checks before we continue, think this should be there for both load_sst and load_weights(later later)
    msg = f'This internal function needs to be run as ispretest, with nrows=#. {_devonlymsg}'
    assert kwg.get('ispretest', False) and not kwg.get('pretest', True), msg
    assert kwg.get('nrows', False), msg
    
    def _try_load_fun(fn, fwmsg=None, swmsg=None, **altkwg):
        fwmsg = ('Preloading of summary statistic unsuccesful. Trying other '
                 'loading approaches to recover.') if fwmsg is None else fwmsg
        try: 
            df = load_fun(fn, **altkwg); err=None
            if swmsg: warnings.warn(swmsg)
        except prst.errors.LoadError as e:
            err = e; warnings.warn(fwmsg)
        except Exception as e: #prst.errors.SumstatSchemaError:
            if swmsg: warnings.warn(swmsg)
            raise
        return err

    drop = ['swmsg','fwmsg']
    swmsg0 = ('Applying comment=# option worked for loading the sumstat. This disables'
        ' pyarrow which makes loading slower. Remove comment section if you need faster loading.')
    ukwg_lst = [
        {},
        dict(comment="#", pyarrow=False, swmsg=swmsg0)
    ] if ukwg_lst is None else ukwg_lst
    
    for i,ukwg in enumerate(ukwg_lst):
        altkwg = {key:item for key, item in kwg.items()}
        altkwg.update(ukwg) # forgot this line, quite important to stop inf recursion.
        err = _try_load_fun(fn, **altkwg)
        #print(err)
        if err is None: return {k:i for k,i in ukwg.items() if not k in drop}
        if i==0: ori_err = err
    
    msg = 'Tried recovery of sumstat loading with different options, but loading in the sumstat did not work.'
    warnings.warn(msg)
    raise ori_err

def load_sst(sst_fn, *, colmap=None, addcols=False, addrsids='auto', calc_beta_mrg=True, n_gwas=None, n_eff_handling='topmedian', delimiter=None, chrom=None, comment=None,
             reqcols=['snp','A1','A2',('beta','oddsratio'),('pval','se_beta')], pyarrow=True, pretest=True, check=True, slicenaninfs=True, validate=True, verbose=True, 
             nrows=None, testnrows=100, ispretest=False, cli=False, readkwg=None): # do not change pretest

    # Hey! I take about 7 seconds on 20M snps with Pyarrow, Optimal enough for now,
    # -> get_beta_mrg takes 60% (np.sort() part takes 40% of that, speedup possible) , and read_csv takes 30% and appears effient already
    # Preps & Pretest:
    if not addcols: addcols=[]
    if readkwg is None: readkwg = {}
    if type(addcols) is str: addcols = [addcols]
    if delimiter == r'\s+': pyarrow=False
    if verbose: print(f'Loading sumstat file.', end='')
    if pretest: # Pre-test: This can become a self call (shorter test run)
        pretestkwg = {key: item for key,item in locals().items() if not key in ['sst_fn']}
        pretestkwg.update(verbose=False, ispretest=True, nrows=testnrows, pretest=False)
        ukwg = _validate_kwg_load_fun(sst_fn, load_fun=load_sst, **pretestkwg)
        if 'pyarrow' in ukwg: pyarrow=ukwg.pop('pyarrow') # <-- this could become a full load_sst() call later
    else: ukwg = {}
        
    # Create kwg for _pd_read_csv()
    kwg = dict()
    if pyarrow and nrows is None: kwg.update(get_pyarrow_prw())
    kwg['delimiter'] = '\t' if delimiter is None else delimiter
    for k in ['nrows', 'comment']: kwg[k] = locals()[k]
    kwg.update(ukwg); kwg.update(readkwg)

    # Loading
    orisst_df = _pd_read_csv(sst_fn, **kwg)
    if verbose: print(f'   -> {orisst_df.shape[0]:>12,} variants sumstat loaded.')

    # Checks: This part should do all the reqcol checks...
    # now this logic is spread out all over load_sst related functions
    if check: check_reqcols_sst(orisst_df, reqcols=reqcols, colmap=colmap)
    
    # Translate:
    conv_dt = get_conv_dt(flow='in', colmap=colmap, verbose=False)
    sst_df = orisst_df.rename(columns=conv_dt)

    # Computation of beta marginal & other Post proc:
    if check and n_gwas is None:
        check_reqcols_sst(orisst_df, colmap=colmap, prefix='Input variable --n_gwas was not given so now; ', 
        reqcols=['n_eff'], postfix='. Please add the '+\
        'column to the sumstat or supply --n_gwas/-n (at times it can be done with '+\
        '--colmap selecting the right n_eff/N column, more above, or copy-paste into chatbot).')
    if calc_beta_mrg:
        if not 'n_eff' in sst_df and not n_gwas is None:
            sst_df['n_eff'] = n_gwas
        from contextlib import nullcontext
        with (warnings.catch_warnings(record=True) if ispretest else nullcontext()):
            sst_df = compute_beta_mrg(sst_df, calc_beta_mrg=calc_beta_mrg, ispretest=ispretest,
                      n_eff_handling=n_eff_handling, slicenaninfs=slicenaninfs, verbose=verbose, cli=cli)
    
    if validate: 
        sst_df = validate_dataframe_index(sst_df, warn=False)
        sst_df = validate_dataframe_select(sst_df, select=['A1A2','chrompos', 'n_eff'], warn=False if ispretest else True)
    else: warnings.warn('load_sst function was set to validate=False. This can lead to bad results.')
    if addcols: sst_df = get_cols(sst_df, addcols=addcols)
    if addrsids and not ispretest:
        nosnp = not 'snp' in sst_df.columns
        if addrsids=='auto' and 'snp' in reqcols and nosnp:
            addrsids=True if all(col in sst_df.columns for col in ['chrom','pos']) else False
        if addrsids == True:
            msg = ('Since no snp column is present, but chrom and pos are,'
            ' a snp column will be added. This might take quite some time and memory.')
            if verbose: print(msg)
            sst_df = get_rsids(sst_df)

    return sst_df

# This is a bad name for this function (todo: change)
def _load_sparse_data(*, 
              ref,
              sst,
              target,
              n_gwas=None,
              chrom = '*', 
              pop = None,
              edgefnfmt = '{pop}/*chr{chrom}_*.edgelist',
              snpfnfmt = 'snplist/{key}.snplist',
              pyarrow=True,
              return_locals=False,
              verbose=True,
              mafcutoff=0.01,
              ncols=120,
              colmap=None
             ):

    ## Preps:
    if pop is None: raise Exception('Population not specified. please specify population.')
    pop = pop.upper()
    try: import pyarrow as pyarrowpack # Prevent var overloading
    except: 
        pyarrow = False
        if verbose: warnings.warn('Could not import python package \'pyarrow\' which means you are ' + \
        'missing out on a lot of speed. Consider installing it for faster data loading with PRSTOOLS.') 
    kwg = dict(dtype_backend="pyarrow", engine='pyarrow') if pyarrow else dict()
    #if type(chrom) is list: chrom=chrom[0]  ############################################# this hack is wrong, earlier code should give an int or string..., input should not be a list..
    if type(chrom) is int: chrom = str(chrom)
    assert type(chrom) is str
    if chrom=='all': chrom='*'
    edgefnfmt = os.path.join(ref, edgefnfmt)
    snpfnfmt  = os.path.join(ref, snpfnfmt)
    snpcnt_dt = {}

    ## Load input data:        

    # Loading sumstat + target data, and finding overlapping snps:
    orisst_df = load_sst(sst, calc_beta_mrg=True, pyarrow=pyarrow, n_gwas=n_gwas, colmap=colmap, verbose=verbose)
    #if verbose: print(f'Loading target file -> ', end='')
    target_df, _ = load_bimfam(target, fam=False, start_string='Loading target file. ', pyarrow=pyarrow, verbose=verbose, end='\n')
    ind = orisst_df['snp'].isin(target_df['snp'])
    sst_df = orisst_df[ind]
    if verbose: print(f'Left with a sumstat of {sst_df.shape[0]:,} variants after intersecting original sumstat '
                      f'({ind.mean()*100:.0f}% incl.) and target ({(ind.sum()/target_df.shape[0])*100:.0f}% incl.)')
    assert ind.sum() != 0, 'No overlap between sumstat and target variants! (based on variant ids)'

    # Load LD data:
    edge_dt = {}; P_dt = {}; snp_dt={}; comb_dt  = {}
    matchstr = edgefnfmt.format(chrom=chrom, pop=pop)
    file_lst = glob.glob(matchstr)
    if len(file_lst) == 0: raise Exception(f'No files found using: search_string={matchstr} \n' + \
        'Perhaps the directory for the ldgm reference is not right?')
    chroms = set()
    for edge_fn in tqdm(file_lst, desc=f"Loading LD data (pop={pop})", ncols=ncols):
        # Parse the key
        key = os.path.split(edge_fn)[-1].split('.')[0]
        curchrom = int(re.search(r'_chr(\d+)_', key).group(1)); chroms.add(curchrom)

        # Create data entries
        edge_dt[key] = pd.read_csv(edge_fn, header=None, names=['i', 'j', 'val']); edge_df=edge_dt[key]
        df = pd.concat([edge_df,edge_df.rename(columns=dict(i='j',j='i'))]).drop_duplicates() # other version in nb appendix
        P = sp.sparse.csc_matrix((df['val'], (df['i'], df['j']))); P_dt[key]=P
        snp_dt[key] = (pd.read_csv(snpfnfmt.format(key=key), **kwg) 
            .rename(columns=dict(index='pindex',position='pos')) 
            .sort_values('pindex')  # Crucial operation, ind slicing used, check at the end of this fun().
            .assign(chrom=curchrom)    # snp_df should have chrom
        )
        # Combine:
        comb_dt[key] = dict(key=key, egde=df, P=P_dt[key], snp_df=snp_dt[key], 
                            start=int(key.split('_')[-2]), chrom=curchrom)


    ## Process the loaded data:

    # Create a genomic-block-start-position ordered dict & dataframe: 
    df = pd.DataFrame.from_dict(comb_dt, orient='index')
    comb_df = df.groupby('chrom').apply(lambda x: x.sort_values('start')).reset_index(drop=True)
    comb_df['i'] = comb_df.index
    def fun(ser): ser['snp_df']['i'] = ser['i']
    comb_df.apply(fun,axis=1)
    comb_dt = comb_df.to_dict(orient='index')

    # Processing of snp frames and removal of duplicate snp_ids
    snp_df = (pd.concat([item['snp_df'] for key, item in comb_dt.items()]).reset_index(names=['blkidx']).drop_duplicates('site_ids'))  # There are variants assigned to two blocks (seemingly @ the edges)
    sst_df = sst_df.drop_duplicates('snp', keep='first') # Removing dups, there are multi-allelic snps in this sumstat. We will only use the first one for now. later this can be removed, allow for multall
    msst_df = pd.merge(snp_df, sst_df, left_on='site_ids', right_on='snp', how='inner', suffixes=("", "_sst")) #60% faster
    if verbose: print(f'Total of {msst_df.shape[0]:,} variants after selecting chromosome(s) (chrom={chrom}) and matching with reference.')
    # %time r=snp_df.site_ids.isin(sst_df.snp) # might offer a speedup at some point time=3.4 s

    # Allele merging & phasing
    ind_match = (msst_df['anc_alleles'] == msst_df['A1']) & (msst_df['deriv_alleles'] == msst_df['A2'])
    ind_flip  = (msst_df['anc_alleles'] == msst_df['A2']) & (msst_df['deriv_alleles'] == msst_df['A1'])
    ind_wrong = ~ind_match & ~ind_flip # not matching at all, happens for example with trialllelic snps
    cast=int #cast='int64[pyarrow]'
    msst_df['phase'] = -1*ind_flip.astype(cast) + 1*ind_match.astype(cast)
    msst_df = msst_df[msst_df.phase != 0]
    # Renaming columns is crucial for followup steps, everything is phase aligned to the reference (not sumstat or target)
    # A1 and A2 are the core names for alleles used throught prstools.
    ind_match = (msst_df['anc_alleles'] == msst_df['A1']) & (msst_df['deriv_alleles'] == msst_df['A2'])
    msst_df = msst_df.rename(columns=dict(A1='A1_sst',A2='A2_sst')).rename(columns=dict(anc_alleles='A1',deriv_alleles='A2'))
    if 'beta_mrg' in msst_df.columns: # Flip beta_mrg where required:
        msst_df['beta_mrg_orig'] = msst_df['beta_mrg']
        msst_df['beta_mrg'] = msst_df['beta_mrg'] * msst_df['phase']
        if 'beta' in msst_df.columns:
            msst_df['beta_orig'] = msst_df['beta']
            msst_df['beta']  = msst_df['beta'] * msst_df['phase']
    else:
        Exception('Not enough info for creation of beta_marginal value.')

    ## Organisation of all into reg_dt dictionary:
    for i, df in msst_df.groupby('i'): comb_dt[i]['msst_df'] = df
    new_comb_dt=dict(); i_new=0
    for i, item in tqdm(comb_dt.items(), desc='Combining LD and sumstat', ncols=ncols):
        item['non0pidx'] = np.where(item['P'].diagonal() != 0)[0]
        if not 'msst_df' in item: continue
        ind0 =  item['msst_df'].pindex.isin(item['non0pidx'])
        ind1 = ~item['msst_df'].duplicated('pindex',keep='first') # Be care with removing this!
        #ind2 =  item['msst_df'][pop].apply(fun) > mafcutoff 
        ind2 = (item['msst_df'][pop.split('-')[-1]] > mafcutoff) & (1 - item['msst_df'][pop.split('-')[-1]] > mafcutoff)
        df = item['msst_df'][ind0&ind1&ind2]
        start, stop = item['key'].split('_')[-2:]
        df.loc[:,['start','stop']] = int(start),int(stop)
        df.loc[:,'i'] = i_new
        item['smsst_df'] = df
        item['i']=None
        if df.shape[0] > 0:
            item['i']=i_new
            new_comb_dt[i_new] = item
            i_new += 1
        else:
            continue
    del comb_dt

    ## Combine all:
    comb_df = pd.DataFrame.from_dict(new_comb_dt, orient='index') #.reset_index(drop=True)
    assert np.all(comb_df.index == comb_df.i)

    ## Last filtering step:

    # Checking:
    fsst_df = pd.concat(comb_df['smsst_df'].values).reset_index(drop=True)
    for i, df in fsst_df.groupby('i'):
        #assert not 'N' in df.columns
        #df = df.rename(columns=dict(N='n_eff'))
        pindex=df['pindex'].values
        index=df.index.values
        assert (pindex==np.sort(pindex)).all()
        assert (index==np.sort(index)).all()
        new_comb_dt[i]['sst_df'] = df
    reg_dt = new_comb_dt

    if verbose: 
        print(f'Total {fsst_df.shape[0]:,} variants left after all steps. (yap yap yap more info here % sumstat, % ref, chroms included.)')

    if return_locals: raise NotImplementedError()
    #    return reg_dt, fsst_df, 'skip' #dict(link_dt=link_dt, fsst_df=fsst_df)
    return reg_dt, fsst_df, False #dict(link_dt=link_dt, fsst_df=fsst_df)
    #'EUR'

def betamrg_to_pval(*, beta_mrg, n_gwas):
    return 2*stats.norm.cdf(-abs(beta_mrg)*np.sqrt(n_gwas))

def pvalandbeta_to_betamrg(*, pvals, beta, n_gwas):
    #return -1*np.sign(beta)*abs(stats.norm.ppf(p/2.0))/n_sqrt # Original
    # LDpred1 has: return sp.sign(raw_beta) * stats.norm.ppf(pval_read / 2.0)/ np.sqrt(N) 
    if np.sum(p>1.): warnings.warn('Input p-vals contains {np.sum(p>1.)} values that are larger then 1. This could be an issue.')
    return np.sign(beta)*abs(stats.norm.ppf(pvals/2.0))/np.sqrt(n_gwas) # Original contains -1 in front, not sure this is the right way?

def load_bimfam(base_fn, strip=True, bim=True, fam=True, chrom='*', delimiter='determine', fil_arr=None, end='\n', start_string='Loading bim/fam. ', 
                testnrows=20, nrows=None, pretest=True, add_xidx=False, add_AX=False, check=True, pyarrow=True, verbose=False, reset_index=True):
    if verbose: print(start_string, end='', flush=True)
    if pretest:
        delimiter='\t'; nrows=testnrows; pretest=False
        pyarrowstart=pyarrow; pyarrow=False
        kwg = locals(); kwg.pop('pyarrowstart')
        kwg['verbose']=False
        try: load_bimfam(**kwg)
        except: 
            delimiter=r'\s+'; warnings.warn('Trying with delimiter=\s+, this can sometimes fix issues.')
        nrows=None; pyarrow=pyarrowstart; pretest=True
    if delimiter=='determine': delimiter=r'\s+'

    prw = get_pyarrow_prw(delimiter=delimiter, pyarrow=pyarrow)

    if strip and (base_fn.split('.')[-1] in ('bim','fam','bed')): base_fn = '.'.join(base_fn.split('.')[:-1])

    bim_df = pd.read_csv(base_fn + '.bim', delimiter=delimiter, header=None, nrows=nrows,
                         names=['chrom', 'snp', 'cm', 'pos', 'A1', 'A2'], **prw) if bim else None
    if type(bim_df) is pd.DataFrame and add_xidx: bim_df['xidx'] = bim_df.index
    if bim: n_snps_start=bim_df.shape[0]

    fam_df = pd.read_csv(base_fn + '.fam', delimiter=r'\s+', header=None,  nrows=nrows,
                         names=['fid', 'iid', 'father', 'mother', 'gender', 'trait'],dtype={0: str, 1: str}) if fam else None
    if bim:
        if check: assert bim_df.head(testnrows).isna().sum().sum()==0, 'NaN detected in bim dataframe.'
        if not pd.api.types.is_numeric_dtype(bim_df['chrom']):
            cmap = get_chrom_map()
            bim_df["chrom"] = bim_df["chrom"].replace(cmap) # the to_numeric() step can be slow, speedup is involved.
            bim_df['chrom'] = pd.to_numeric(bim_df['chrom'], errors='coerce').astype('Int64') # rrx= bim_df['chrom'].unique() appears fast, so perhaps fix, mod cmap
        if not chrom in ['*','all']:
            #ind = bim_df['chrom'] == bim_df['chrom'].dtype.type(chrom) # old one 
            ind = bim_df['chrom'].isin(get_chrom_lst(chrom))
            bim_df = bim_df[ind]
        if fil_arr is not None:
            bim_df = bim_df[bim_df.snp.isin(fil_arr)]
            bim_df = bim_df.reset_index(drop=True)
        if reset_index: bim_df = validate_dataframe_index(bim_df, warn=False)
        if add_AX: bim_df = get_AX(bim_df)
        n_snps_end = bim_df.shape[0]
        if not check and bim_df.shape[0]<1e2: 
            msg = (f'\033[1;31mWARNING: The size of the loaded bim set is less then 100 (={bim_df.shape[0]} snps). '
                  'Often a sign something is going wrong (e.g. chrom of choice not present input bim file) \033[0m')
            warnings.warn(msg)
    if verbose:
        lst=[]
        if bim: inject = f', selecting {n_snps_end:,} with chrom={chrom}' if 'ind' in locals() and ind.shape != bim_df.shape[0] else ''
        if bim: lst += [f'{n_snps_start:>12,} variants bim file loaded{inject}']
        if fam:
            #if fam_df.shape[0] == fam_df['fid'].nunique():
            lst += [f'{fam_df.shape[0]:,} induviduals fam file loaded']
        report = ' & '.join(lst)
        if report == '': report='no bim or fam file'
        if len(prw) > 1: report=report+' (used pyarrow)'
        print(f'-> {report}.', end=end, flush=True)

    return bim_df, fam_df

def _load_bimfam_from_srd(srd, make_bimfam_attrs=True, verbose=False, skipifpresent=True, **kwg):

    # hjerhjgerjkhgjkh
    prstlogs = prst.utils.get_prstlogs()
    tic, toc = prstlogs.get_tictoc()
    toc('inside of loadbimfamfromsrd')

    if hasattr(srd, 'bim_df') and hasattr(srd, 'fam_df') and hasattr(srd,'_validated_bimfam'): 
        bim_df = srd.bim_df; fam_df = srd.fam_df
        if srd._validated_bimfam and skipifpresent:
            return bim_df.copy(), fam_df.copy()
    else:
        now_srd = srd
        for _ in range(20):
            if not hasattr(now_srd, 'filename'):
                if hasattr(now_srd, '_internal'):
                    now_srd = now_srd._internal
                else: raise Exception(f'srd={srd} does not have a filename attribute, perhaps its not the right input type (srd = pysnptools.SnpReaDer())')
            else: break
        else: 
            raise Exception(f'maximum depth of 20 exceeded for srd={srd}')
        toc('staring load bimfam')
        bim_df, fam_df = load_bimfam(now_srd.filename,verbose=verbose,end=' ',**kwg) # costly line: 37% (implement pyarrow..)
        toc('done laoding bimfam')
        assert now_srd.count_A1 == True

    # Assert that family and induvidual IDs match: 
    if verbose: print('Now validating bim/fam.', end='', flush=True)
    toc('validation starting')

    if srd.iid.shape[0] != fam_df.shape[0]:
        print('CRASH IN VALIDATION PART'); #ip.embed()
        fidx = pd.MultiIndex.from_arrays(srd.iid.T,names=['fid','iid'])
        fam_df[['fid','iid']] = fam_df[['fid','iid']].astype(srd.iid.dtype)
        fam_df = fam_df.set_index(['fid','iid']).loc[fidx].reset_index()
    assert np.all(srd.iid == fam_df[['fid','iid']].astype(srd.iid.dtype))

    toc('validation part1 done')
    #ip.embed()

    # try slicing sids are not matching up:
    if verbose: print('.. ', end='', flush=True)
    if srd.sid.shape[0] != bim_df['snp'].shape[0]: 
        if not np.all(srd.sid == bim_df['snp'].astype(srd.sid.dtype)):# This and/or previous is a costly line: 52%
            #print('CRASH IN VALIDATION PART'); #ip.embed()
            slicer_df = bim_df.set_index('snp',drop=True)
            slicer_df = slicer_df.loc[srd.sid].reset_index(drop=False)
            bim_df = slicer_df[bim_df.columns] # make sure column order is bim complient
    assert np.all(srd.sid == bim_df['snp'].astype(srd.sid.dtype))

    toc('validation is DONE')

    if make_bimfam_attrs:
        srd.bim_df=bim_df; srd.fam_df=fam_df; srd._validated_bimfam=True
    if verbose: print('Done')

    return bim_df.copy(), fam_df.copy()

def load_ref(ref_fn, chrom=None, verbose=False, rename_dt=dict(maf='maf_ref',af_A1='af_A1_ref'), reset_index=True):
    chrom = None if (chrom=='*' or str(chrom).lower()=='all') else chrom
    if verbose: print('Loading reference file.', end='')
    ref_df = load_sst(ref_fn, n_gwas=None, calc_beta_mrg=False, check=False, verbose=False)
    if chrom is not None:
        ref_df = ref_df[ref_df.chrom.astype(int).isin(get_chrom_lst(chrom))]
    xtra = '.' if chrom is None else f' from chromosome {chrom}.'
    if verbose:  print(f' -> {ref_df.shape[0]:>12,} variants loaded'+xtra)
    ref_df = ref_df.rename(columns=rename_dt)
    if reset_index:
        ref_df = validate_dataframe_index(ref_df, warn=False)
    return ref_df

def load_snpdb(snpdb_df='mini'):
    assert type(snpdb_df) is str, f'Input to load_snpdb must be string, it is {type(snpdb_df)}'
    prstcfg = prst.utils.load_config()
    prstdatadir = prstcfg['prstdatadir']
    if snpdb_df == 'mini' and prstdatadir:   fn = 'snpdb_mini.tsv.gz'
    elif snpdb_df == 'full' and prstdatadir: fn = 'snpdb_full.tsv.gz'
    elif snpdb_df in ('mini','full'): 
        msg = f'Trying load a snp database (version="{snpdb_df}"), but prstdatadir is not set, '
        msg += 'which is required for this functionality. Please run "prst config" to set it.'
        if not prstdatadir: raise RuntimeError(msg)
    else: assert os.path.exists(os.path.expanduser(snpdb_df)), (f'Path {snpdb_df} must exist, or specify '
                    'snpdb=mini or full (if prst data dir is enabled) or specify dataframe type as input')
    fn = prst.utils.validate_path(fn=fn, must_exist=True, handle_prstdatadir='only')
    prw = get_pyarrow_prw()
    snpdb_df = pd.read_csv(fn, sep='\t', **prw)
    return snpdb_df

def load_weights(fn, ftype='auto', pyarrow=True, sep:str='\t', verbose=False):
    if pyarrow: # pyarrow mechanics
        try: import pyarrow as pyarrowpack # Prevent var overloading
        except: pyarrow = False
        if verbose and not pyarrow: 
            warnings.warn('Could not import python package \'pyarrow\' which means you are ' + \
            'missing out on a lot of speed. Consider installing it for faster data loading with PRSTOOLS.')
    if sep == r'\s+': pyarrow=False
    if pyarrow: prw = dict(dtype_backend="pyarrow", engine='pyarrow')
    else: prw = dict()
    def detect_ftype(fn):
        options=['legacyweights.tsv', 'prstweights.tsv','prstweights.h5','prstweights.parquet']
        for this in options:
            if fn.endswith(this): return this
        return 'headed-txt'
    cur_ftype = detect_ftype(fn) if ftype == 'auto' else ftype
    header_dt = dict(header=None) if cur_ftype in ['headless-txt', 'legacyweights.tsv'] else {}
    from prstools.models import BasePred
    names = None if not cur_ftype == 'legacyweights.tsv' else BasePred.default_weight_cols
    if cur_ftype == 'prstweights.h5': df = pd.read_hdf(fn, key='df')
    elif cur_ftype == 'prstweights.parquet': df = pd.read_parquet(fn)
    else: df = pd.read_csv(fn, sep=sep, names=names, **header_dt, **prw)
    return df

def load_bed(fn, make_bimfam_attrs=True, countA12correct=True, verbose=False, start_string='Loading plink files (@ {fn}). '):
    from bed_reader import open_bed
    if verbose: 
        proc_fn = f'...{fn[-17:]}' if len(fn) >= 20 else fn
        print(start_string.format_map(prst.utils.AutoDict(fn=proc_fn)), end='', flush=True)
    fn = prst.utils.validate_path(fn=fn, must_exist=False)
    iid_count=None; sid_count=None
    if make_bimfam_attrs:
        bim_df, fam_df = prst.load_bimfam(fn,add_xidx=True, add_AX=True)
        iid_count=fam_df.shape[0]; sid_count=bim_df.shape[0]
    base_fn = '.'.join(fn.split('.')[:-1]) if (fn.split('.')[-1] in ('bim','fam','bed')) else fn
    bed = open_bed(base_fn+'.bed', iid_count=iid_count, sid_count=sid_count)
    if make_bimfam_attrs:
        bed.bim_df = bim_df; bed.fam_df = fam_df
    return bed

def load_srd(fn, make_bimfam_attrs=True, countA12correct=True, verbose=False, start_string='Loading plink files (@ {fn}). '):
    if verbose:
        from prstools.utils import AutoDict
        proc_fn = f'...{fn[-17:]}' if len(fn) >= 20 else fn
        print(start_string.format_map(AutoDict(fn=proc_fn)), end='', flush=True)
    #here
    prstlogs = prst.utils.get_prstlogs()
    tic, toc = prstlogs.get_tictoc()
    #ip.embed()

    from pysnptools.snpreader import Bed
    assert countA12correct and type(countA12correct) is bool

    toc('starting bim load')
    bim_df, fam_df = load_bimfam(fn, verbose=verbose, start_string='')
    toc('done bimfam load, making iid sid pos vars')
    sid=bim_df['snp'].to_numpy().astype(str)
    iid=fam_df[['fid','iid']].to_numpy().astype(str)
    pos=bim_df[['chrom','cm','pos']].to_numpy().astype('float64')
    pos[:,1] = np.nan 
    toc('now creating bed')
    srd = Bed(os.path.expanduser(fn), sid=sid, iid=iid, pos=pos, count_A1=True)
    toc('done with bed, finishing load_srd')
    if make_bimfam_attrs:
        srd.bim_df = bim_df
        srd.fam_df = fam_df
    #ip.embed()
    #_not_needed = load_bimfam_from_srd(srd, make_bimfam_attrs=make_bimfam_attrs, verbose=verbose, start_string='')
    return srd

def load_linkagedata(fn, load_all=False):
    curfn = glob.glob(fn.format_map(defaultdict(lambda:'*')))[-1]
    clsattr_dt = json.loads(pd.read_hdf(curfn, key='clsattr_dt').loc[0,0])
    clsattr_dt['curdn'] = os.path.dirname(curfn)
    linkdata = LinkageData(clsattr_dt=clsattr_dt)
    if load_all: linkdata.load_linkage_allregions()
    return linkdata

def check_regdef(regdef_df, allow_overlap=False, verbose=False):
    assert not allow_overlap, 'Need to implement versions with allowed overlap'
    msg = 'Duplicated region-ids in regdef. This means the region definition is defective'
    assert regdef_df['regid'].duplicated().sum() == 0, msg
    mcols = {'chrom','start','stop','regid'} - set(regdef_df.columns)
    if len(mcols) > 0: raise ValueError(f"Not all req columns present in regdef_df. Missing columns: {mcols}")
    assert regdef_df.sort_values(['chrom','start']).index.equals(regdef_df.index), 'Seems regdef was not properly sorted, hence corrupted/invalid'
    for chrom, df in regdef_df.groupby('chrom'):
        msg = 'Overlapping regions in regdef_df, hence it is corrupted/invalid.'
        assert (df['start'].values[1:] >= df['stop'].values[:-1]).all(), msg
    if verbose: print('regdef_df is valid! (i.e. required cols, sorted, unique regids, and no overlapping regions)')

def load_regdef(regdef='ldgm-all-hg38', fnfmt='./data/defs/regdef/{regdef}.regdef.tsv', check=True):
    curdn = os.path.dirname(prst.__file__)
    fnfmt = os.path.join(curdn, fnfmt)
    fn = fnfmt.format(regdef=regdef)
    if not os.path.exists(fn):
        dn = os.path.dirname(fn)
        lst = [os.path.basename(elem).replace('.regdef.tsv','') for elem in glob.glob(dn+'/*.regdef.tsv')]
        msg = f'{fn} \n--> You should pick from {lst} (all @ {dn}).'
        raise FileNotFoundError(msg)
    regdef_df = pd.read_csv(fn, delimiter='\t') # dataframe with region definitions
    assert check, 'always doing check!'
    check_regdef(regdef_df)
    return regdef_df

def load_prscs_ldblk(fn,blkid):
    ## NOT SURE I NEED THIS FUNCTION ...

    #def load_linkage_region(self, *, i):
    geno_dt = self.reg_dt[i]
    store_dt = geno_dt['store_dt']

    for varname, file_dt in store_dt.items():
        module = importlib.import_module('.'.join(file_dt['typestr'].split('.')[:-1]))
        cname  = file_dt['typestr'].split('.')[-1]
        CurClass = getattr(module, cname) # Retrieves module.submodule.submodule.. etc
        curfullfn = os.path.join(self.curdn, file_dt['fn'])
        geno_dt[varname] = CurClass(pd.read_hdf(curfullfn, key=file_dt['key']))
        if self.verbose: print(f'loading: fn={curfullfn} key={file_dt["key"]}'+' '*50, end='\r')

    return something

# Later this could go into a seperate dataset submodule, perhaps in classes (like keras)
# from keras.datasets import mnist; data = mnist.load_data()
def load_example(dn='./data/_example/', n_gwas=2565, pop='EUR', verbose=False):
    try: from pysnptools.snpreader import Bed
    except: warnings.warn('not pysnptools installed, cannot read bed files now..')
    from prstools.utils import Struct
    from prstools.linkage import RefLinkageData, SparseLinkageData
    if pop != 'EUR': raise NotImplementedError()
    dn = os.path.join(os.path.dirname(prst.__file__), dn)
    st=Struct()
    dense_ref =dn+'ldref_1kg_pop/'; target=dn+'target'; sst=dn+'sumstats.tsv'
    sparse_ref=dn+'ldgm_1kg_pop/';
    st.dense_lnk  = RefLinkageData.from_cli_params(ref=dense_ref,target=target,sst=sst,n_gwas=n_gwas, verbose=verbose)
    warnings.warn('normal linkage data not working atm')
    st.sparse_lnk = SparseLinkageData.from_cli_params(ref=sparse_ref,target=target,sst=sst,n_gwas=n_gwas, pop=pop, verbose=verbose)
    st.target_srd = Bed(target, count_A1=True) if 'Bed' in locals() else 'If you want raw example genotypes to work with, install \'pysnptools\''        
    st.target_bed = prst.io.load_bed(target)
    return st

#     def _pd_to_atomizer(*, to_file, fn, **kwg):
#         tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"
#         ret = to_file(tmp_fn, **kwg)
#         os.replace(tmp_fn, fn)
#         return ret
def _pd_to_atomizer(*, to_file, fn, **kwg):
    tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"
    try:
        ret = to_file(tmp_fn, **kwg)   # ← error happens here
        os.replace(tmp_fn, fn)
        return ret
    finally:
        if os.path.exists(tmp_fn):
            try: os.remove(tmp_fn)
            except OSError: pass

def save_bim(bim_df, fn=None, ommitcm=False, verbose=True):
    assert fn is not None
    assert ommitcm is False, 'not implemented' 
    reqcols = ['chrom', 'snp', 'cm', 'pos', 'A1', 'A2']
    bim_df = bim_df[reqcols]
    if verbose: print(f'Formatting done, Saving bim to: {fn}', end=' ')
    #bim_df.to_csv(fn, sep='\t', header=False, index=None)
    _pd_to_atomizer(to_file=bim_df.to_csv, fn=fn, sep='\t', header=False, index=None)
    if verbose: print(f'-> Done')

def save_fam(fam_df, fn=None, verbose=True):
    assert fn is not None and type(fn) is str
    reqcols = ['fid', 'iid', 'father', 'mother', 'gender', 'trait']
    fam_df = fam_df[reqcols]
    if verbose: print(f'Formatting done, Saving fam to: {fn}', end=' ')
    #bim_df.to_csv(fn, sep='\t', header=False, index=None)
    _pd_to_atomizer(to_file=fam_df.to_csv, fn=fn, sep='\t', header=False, index=None)
    if verbose: print(f'-> Done')

def save_sst(sst_df, fn=None, return_sst=False, ftype='tsv', basecols=None, addicols=['OR','SE','P','N','FRQA1'], extracols=None, nancheck=True, verbose=True):
    assert fn is not None
    from prstools.models import BasePred
    if basecols is None: basecols=list(BasePred.default_sst_cols)
    if ftype == 'tsv':
        header=True; sep='\t'
        assert len(addicols) > 0, 'Sumstat without any additional cols cannot be a proper sumstat'
        extracols = [] if extracols is None else extracols
        outcols = basecols + addicols + extracols
    else:
        raise Exception('not implementeed alternative storage options.')

    fn = fn.format_map(dict(ftype=ftype)) # Maybe some AutoDict buzz here later.
    conv_dt = prst.io.get_conv_dt(flow='out')
    preout_df = sst_df.rename(columns=conv_dt)
    avail_outcols = [col for col in outcols if col in preout_df.columns]
    if verbose: print(f'Prepping sumstat, with selection of available columns: {", ".join(avail_outcols)}')
    if verbose: print(f'Saving sumstats to: {fn}', end=' ')
    out_df = preout_df[avail_outcols]
    if nancheck: assert out_df.isna().sum().sum()==0, 'NaN detected in sst dataframe, cannot save sst file with NaNs if nancheck=True'
    #out_df.to_csv(tmp_fn, sep=sep, index=False, header=header)
    _pd_to_atomizer(to_file=out_df.to_csv, fn=fn, sep=sep, header=header, index=False)
    #os.replace(tmp_fn, fn) # atomically move into place
    if verbose: print(f'-> Done')
    if return_sst: return out_df

def save_prs(yhat, *, fn, ftype='prspred.tsv', nanwarn=True, verbose=False, reset_index=True, end='\n\n'):
    assert type(yhat) is pd.DataFrame, f'Input \'yhat\' is required to be pd.DataFrame. It is currently: {type(yhat)}'
    if reset_index: yhat = yhat.reset_index(drop=False) # With this the FID and IID become columns
    #yhat = yhat.rename(columns=dict(fid='FID',iid='IID')) why these names...
    if ftype == 'prspred.tsv':
        sep='\t'
        to_file = yhat.to_csv
    else:
        raise ValueError(f"'{ftype}' is not a valid filetype/ftype. only 'prspred.tsv' availabe atm")

    fn = fn.format_map(dict(ftype=ftype)) # Maybe some AutoDict buzz here later.
    #import uuid; tmp_fn = f"{fn}.incomplete.{uuid.uuid4().hex[:16]}"  # unique temp file name  
    if verbose: print(f'Saving prediction (i.e. PRS) to: {fn}', end=' ')
    _pd_to_atomizer(to_file=to_file, fn=fn, sep=sep, index=False)
    #to_file(tmp_fn, sep=sep, index=False)
    #os.replace(tmp_fn, fn) # atomically move into place
    if verbose: print(f'-> Done', end=end)

if not '__file__' in locals():
    import sys
    if np.all([x in sys.argv[-1] for x in ('jupyter','.json')]+['ipykernel_launcher.py' in sys.argv[0]]):
        with open('../prstools/io.py', 'w') as loadrf: loadrf.write(In[-1])
        print('Written to:', loadrf.name)