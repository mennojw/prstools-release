"""
PrivacyPreservingMetricsComputer
MultiPGSComputer

"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
# from functools import partial
try: from tqdm.auto import tqdm
except: tqdm = lambda x:x
from scipy import stats

def _ols(y, X):
    X = np.asarray(X, dtype=float)
    if X.ndim == 1: X = X[:,None]
    X = np.column_stack([np.ones(len(y)), X])

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta

    rank = np.linalg.matrix_rank(X)
    dof = len(y) - rank
    if dof <= 0: raise ValueError("Not enough individuals to fit regression.")

    sse = resid @ resid
    sst = ((y - y.mean())**2).sum()
    r2 = 1 - sse/sst if sst else np.nan

    s2 = sse / dof
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)) * s2)
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)

    return beta, se, p, r2
    
def evaluate(pheno_df, pred_df, covar_df=None, metrics=['R2'], verbose=True):
    """Evaluate one or more PRSs against one or more quantitative phenotypes."""
    assert all(m in ['R2'] for m in metrics), 'R2 only avail metric for now'
    phenocols = list(pheno_df.columns)
    prscols   = list(pred_df.columns)
    covcols   = [] if covar_df is None else list(covar_df.columns)

    if not phenocols: raise ValueError("No phenotype columns found.")
    if not prscols: raise ValueError("No PRS columns found.")

    n_aligned = len(pheno_df); results = []

    for pheno in phenocols:
        for prs in prscols:
            mask = pheno_df[pheno].notna() & pred_df[prs].notna()
            if covcols: mask &= covar_df[covcols].notna().all(axis=1)

            n = int(mask.sum())
            if n < 3: raise ValueError(f"Too few complete individuals for phenotype={pheno!r}, prs={prs!r} ({n=}).")

            y = pheno_df.loc[mask, pheno].astype(float).to_numpy()
            s = pred_df.loc[mask, prs].astype(float).to_numpy()
            C = covar_df.loc[mask, covcols].astype(float).to_numpy() if covcols else np.empty((n, 0))

            _, _, _, r2_base = _ols(y, C)
            beta, se, p, r2 = _ols(y, np.column_stack([C, s]))

            results.append(dict(phenotype=pheno, prs=prs, n_aligned=n_aligned, n=n, n_missing=n_aligned-n,
                                beta=beta[-1], se=se[-1], p=p[-1], r2=r2, r2_base=r2_base, delta_r2=r2-r2_base))

            if verbose:
                print(f"Evaluating {pheno} ~ {prs}: {n:,}/{n_aligned:,} complete, "
                      f"R2={r2:.4g}, delta-R2={r2-r2_base:.4g}")

    return pd.DataFrame(results)


class PrivacyPreservingMetricsComputer():
    
    def __init__(self, *, linkdata, brd, Bm, s=None, dtype='float32',
                 clear_linkage=False, pbar=tqdm, verbose=True):
        
        self.linkdata   = linkdata
        self.brd        = brd
        self.s          = s
        self.Bm         = Bm 
        self.dtype      = dtype
        self.clear_linkage = clear_linkage
        if not pbar: pbar = lambda x: x
        self.pbar       = pbar
        self.verbose    = verbose
            
    def evaluate(self, debug=False):
        
        # Load and init variables:
        linkdata = self.linkdata
        #linkdata = self.linkdata.init() # init, in case required.
        brd = self.brd; Bm = self.Bm
        if self.verbose & (self.s is None): print('Retrieving Standard Dev. var \'s\'')
        s = linkdata.get_s() if self.s is None else self.s
        assert (np.isnan(s).sum()+np.isinf(s).sum()) == 0 # add s =standard dev as argument with object creation if this line keeps failing
        bCb = 0.; BmBt = 0.
        info_dt = dict()

        # Cycle through the blocks:
        for i, geno_dt in self.pbar(linkdata.reg_dt.items()):
            #if self.verbose: print(f'PPB: Processing region {i}', end='\r')
            # Ready the LD:
            L = linkdata.get_left_linkage_region(i=i)
            D = linkdata.get_auto_linkage_region(i=i)
            R = linkdata.get_right_linkage_region(i=i)
            lr = linkdata.get_left_range_region(i=i)
            ar = linkdata.get_auto_range_region(i=i)
            rr = linkdata.get_right_range_region(i=i)

            # Ready The Weights:
            B_L = brd[:,lr[0]:lr[1]].read().val.astype(self.dtype).T
            B_D = brd[:,ar[0]:ar[1]].read(dtype=self.dtype).val.T
            B_R = brd[:,rr[0]:rr[1]].read(dtype=self.dtype).val.T
            B_L = s[lr[0]:lr[1]]*B_L
            B_D = s[ar[0]:ar[1]]*B_D
            B_R = s[rr[0]:rr[1]]*B_R

            # Do the computation:
            CB = L.dot(B_L) + D.dot(B_D) + R.dot(B_R)
            bCb += (B_D*CB).sum(axis=0)
            BmBt += (B_D.T.dot(Bm.iloc[ar[0]:ar[1],:])).T
            info_dt[i] = dict(shapeL=L.shape, shapeD=D.shape, shapeR=R.shape, 
                              lr=lr, ar=ar, rr=rr)

            # Pruning to minimize memory overhead:
            if (i > 0) and self.clear_linkage:
                linkdata.clear_linkage_region(i=i-1)
            #if (i > 38) & debug: break; return locals()

        # Complete resutls:
        linkdata.clear_all_xda()
        cols = brd.row.astype(str).flatten()
        bCb  = pd.DataFrame(bCb[np.newaxis,:], index=['bCb'], columns=cols)
        BmBt = pd.DataFrame(BmBt, index=Bm.columns, columns=cols)
        ppbr2_df = (BmBt**2)/bCb.loc['bCb']
        res_dt = dict(ppbr2_df=ppbr2_df, bCb=bCb, BmBt=BmBt, info_dt=info_dt, s=s)

        return res_dt
    
# locals_dt = dict()
class MultiPGSComputer():
    
    def __init__(self, *, brd, unscaled=True, verbose=False, dtype='float32', allow_nan=False, pbar=tqdm):
        self.brd   = brd
        if hasattr(brd, 'val'):
            assert np.sum(np.isnan(brd.val)) == 0 
        self.unscaled = unscaled
        self.verbose = verbose
        self.dtype  = dtype
        self.allow_nan = allow_nan
        self.pbar = pbar
        
    def predict(self, *, srd, prd=None, n_inchunk=1000, stansda=None):
        try:
            from pysnptools.standardizer import UnitTrained
        except e:
            raise Exception(e,'You should install pysnptools for this functionality to work.')
            
        # Load that PGS (& optionaly phenos)
        brd = self.brd
        Yhat = np.zeros((srd.shape[0], brd.shape[0]), dtype=self.dtype)        
        assert np.all(brd.col.astype(str) == srd.sid)
        if prd: 
            assert np.all(srd.iid == prd.iid)
            pda   = prd.read(dtype=self.dtype).standardize()
            Ytru  = pda.val
            Bm    = np.zeros((srd.shape[1], Ytru.shape[1])) + np.nan
        
        # Loop through Genome:
        stansda_lst = []; start=0
        for start in self.pbar(range(0, srd.shape[1], n_inchunk)):
            stop = min(start+n_inchunk, srd.shape[1])
            sda, stansda = srd[:,start:stop].read(dtype=self.dtype).standardize(return_trained=True)
            X = sda.val
            s = stansda.stats[:,1][:,np.newaxis]; s[np.isinf(s)] = 1
            L = brd[:,start:stop].read(dtype=self.dtype).val
            B = s*L.T # s seems of little effect on time here, projected loading takes time 200s for HM3 8K betas. for 10k induv.
            Yhat += X@B
            if prd: Bm[start:stop] = X.T@Ytru
            stansda_lst.append(stansda)
            
        if prd:    
            Bm = Bm/Ytru.shape[0]
            if not self.allow_nan: assert np.isnan(Bm).sum() == 0
            if not self.allow_nan: assert np.isnan(Bm).sum() == 0
            Bm = pd.DataFrame(Bm, index=srd.sid, columns=prd.col)
            Ytru = pd.DataFrame(Ytru, # Make Ytru a proper dataframe
                index=pd.MultiIndex.from_arrays(prd.iid.T, names=('fid','iid')),
                columns=prd.col)
        else:
            Ytru=None; Bm=None
            
        # Combine Standardizers:
        sid     = np.concatenate([stan.sid   for stan in stansda_lst])
        assert  np.unique(sid).shape[0] == sid.shape[0]
        stats   = np.concatenate([stan.stats for stan in stansda_lst])
        stansda = UnitTrained(sid, stats)   
        s = stansda.stats[:,1][:,np.newaxis]; s[np.isinf(s)] = 1
        
        # Create Yhat dataframe:
        Yhat  = pd.DataFrame(
            data    = Yhat, 
            index   = pd.MultiIndex.from_arrays(srd.iid.T, names=('fid','iid')),
            columns = self.brd.row.astype(str)
        ); assert Yhat.isna().sum().sum() == 0
        
        res_dt = dict(Yhat=Yhat, Bm=Bm, brd=brd, Ytru=Ytru, stansda=stansda, s=s)
        
        return locals()
    
    def run(self):
        pass
    
    def fit(self):
        pass

if '_isdevenv_prstools' in locals():
    if _isdevenv_prstools:
        with open('../prstools/scores.py', 'w') as f: f.write(In[-1])
        print('Written to:', f.name)