"""Numba-accelerated BKL kernel + multiprocessing over samples. Same physics/protocol as haah_tmem.py
(RG decoder readout every T_ec = e^{4 beta}/100). Shares the code construction, RG decoder and logical test
with haah_tmem.py; only the BKL inner loop is re-implemented.
Usage: python3 haah_tmem_fast.py <beta> <L,L,..> <nsamp> [tag] [s4,s2,s0] [seed_offset] [Ltqo|auto]"""
import numpy as np, sys, json, time, multiprocessing as mp
from numba import njit
sys.path.insert(0, __file__.rsplit('/',1)[0])
from haah_tmem import Haah

@njit(cache=True)
def _seed(s): np.random.seed(s)

@njit(cache=True)
def bkl_until(e,viol,nviol,q2s,s2q,rate_cls,t,t_stop):
    N=e.shape[0]
    while True:
        counts=np.zeros(5)
        for q in range(N): counts[nviol[q]]+=1
        R=0.0
        for c in range(5): R+=counts[c]*rate_cls[c]
        dt=np.random.exponential(1.0/R)
        if t+dt>=t_stop: return t_stop
        t+=dt
        u=np.random.random()*R; cls=0; acc=counts[0]*rate_cls[0]
        while u>acc and cls<4:
            cls+=1; acc+=counts[cls]*rate_cls[cls]
        j=int(np.random.random()*counts[cls]); q=-1; seen=0
        for qq in range(N):
            if nviol[qq]==cls:
                if seen==j: q=qq; break
                seen+=1
        if q<0: raise RuntimeError('BKL class selection failed (no qubit in chosen class)')
        e[q]^=1
        for i in range(4):
            s=q2s[q,i]; viol[s]^=1; d=1 if viol[s]==1 else -1
            for k in range(8): nviol[s2q[s,k]]+=d

def memory_time_fast(H,beta,seed,spec=(1.0,1.0,1.0),max_checks=10**7):
    _seed(seed); N=H.N
    q2s=np.array(H.q2s,dtype=np.int64); s2q=np.array(H.s2q,dtype=np.int64)
    e=np.zeros(N,np.int8); viol=np.zeros(H.L**3,np.int8); nviol=np.zeros(N,np.int8)
    s4,s2,s0=spec; rate_cls=np.array([s4*np.exp(-4*beta),s2*np.exp(-2*beta),s0,s2,s4])
    Tec=np.exp(4*beta)/100; t=0.0
    for n in range(1,max_checks):
        t=bkl_until(e,viol,nviol,q2s,s2q,rate_cls,t,n*Tec)
        if not viol.any():
            if e.any() and H.is_logical_error(e): return n*Tec,'logical',0
            continue                                   # empty syndrome: nothing to decode
        P=H.decode(viol)
        if P is None: return n*Tec,'decoder',int(viol.sum())
        if H.is_logical_error((e+P)%2): return n*Tec,'logical',int(viol.sum())
    return np.nan,'timeout',int(viol.sum())

_H=None
def _worker(args):
    beta,seed,spec=args; return memory_time_fast(_H,beta,seed,spec)

def run(L,beta,nsamp,spec=(1.0,1.0,1.0),nproc=16,seed_offset=0,Ltqo=None):
    """seed_offset: shifts the per-sample seeds 1000*L+i (published runs used offset 0). Ltqo: decoder box cap (default L-2)."""
    global _H; _H=Haah(L,Ltqo=Ltqo); t0=time.time()
    with mp.get_context('fork').Pool(nproc) as pool:
        res=pool.map(_worker,[(beta,seed_offset+1000*L+i,spec) for i in range(nsamp)])
    T=np.array([r[0] for r in res]); modes=[r[1] for r in res]; nd=[r[2] for r in res]; ok=~np.isnan(T); T=T[ok]
    Tec=np.exp(4*beta)/100
    out={"beta":beta,"spec":list(spec),"n":int(ok.sum()),"Tmem_mean":float(T.mean()),"Tmem_sem":float(T.std()/np.sqrt(len(T))),"Tec":Tec,
         "Tmem_over_Tec":float(T.mean()/Tec),"modes":{m:modes.count(m) for m in set(modes)},"defects_at_failure_mean":float(np.mean(nd)),"wall_s":time.time()-t0,"k":_H.k}
    print(f"L={L} beta={beta} spec={spec} k={_H.k} n={len(T)} Tmem={T.mean():.3g} ± {T.std()/np.sqrt(len(T)):.2g} Tmem/Tec={T.mean()/Tec:.1f} modes={out['modes']} <defects>={np.mean(nd):.1f} ({time.time()-t0:.0f}s)",flush=True)
    return out

if __name__=="__main__":
    beta=float(sys.argv[1]); Ls=[int(x) for x in sys.argv[2].split(',')]; nsamp=int(sys.argv[3])
    tag=sys.argv[4] if len(sys.argv)>4 else ""; spec=tuple(float(x) for x in sys.argv[5].split(",")) if len(sys.argv)>5 else (1.0,1.0,1.0)
    seed_offset=int(sys.argv[6]) if len(sys.argv)>6 else 0; Ltqo=None if len(sys.argv)<=7 or sys.argv[7]=='auto' else int(sys.argv[7])
    out={L:run(L,beta,nsamp,spec,seed_offset=seed_offset,Ltqo=Ltqo) for L in Ls}
    json.dump(out,open(__file__.rsplit('/',2)[0]+f"/results/tmem_beta{beta}{tag}.json","w"),indent=1)
