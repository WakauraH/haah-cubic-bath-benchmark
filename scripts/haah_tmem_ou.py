"""H2 (classical memory): diffusion-channel rate s0(t) = exp(xi(t) - sigma^2/2), xi = Ornstein-Uhlenbeck with
correlation time tau_c and stationary variance sigma^2, so E[s0]=1 (same mean rate and same Boltzmann
populations as the flat bath). Rates are piecewise-constant over dt_ou = tau_c/20 while dE=0 moves exist;
otherwise xi is advanced with the exact OU transition over the event gap.
Self-check: sigma -> 0 reproduces the flat-bath T_mem (L=7, beta=4.5: 7.1e7 +/- 0.8e7)."""
import numpy as np, sys, json, time, multiprocessing as mp
from numba import njit
sys.path.insert(0, __file__.rsplit('/',1)[0])
from haah_tmem import Haah
from haah_tmem_fast import _seed

@njit(cache=True)
def ou_jump(xi,dt,tau,sig2):
    a=np.exp(-dt/tau); return xi*a+np.sqrt(sig2*(1-a*a))*np.random.normal()

@njit(cache=True)
def bkl_until_ou(e,viol,nviol,q2s,s2q,base,t,t_stop,xi,tau,sig2):
    """base = rate classes with s0=1. Returns (t_stop, xi at t_stop)."""
    N=e.shape[0]; dt_ou=tau/20.0
    while True:
        counts=np.zeros(5)
        for q in range(N): counts[nviol[q]]+=1
        s0=np.exp(xi-0.5*sig2)
        R=counts[0]*base[0]+counts[1]*base[1]+counts[2]*base[2]*s0+counts[3]*base[3]+counts[4]*base[4]
        dt=np.random.exponential(1.0/R)
        if counts[2]>0 and dt>dt_ou:               # rate may change before the event: advance OU, no event
            if t+dt_ou>=t_stop: return t_stop, ou_jump(xi,t_stop-t,tau,sig2)
            t+=dt_ou; xi=ou_jump(xi,dt_ou,tau,sig2); continue
        if t+dt>=t_stop: return t_stop, ou_jump(xi,t_stop-t,tau,sig2)
        t+=dt; xi=ou_jump(xi,dt,tau,sig2)
        u=np.random.random()*R; acc=counts[0]*base[0]; cls=0
        w=np.array([base[0],base[1],base[2]*s0,base[3],base[4]])
        while u>acc and cls<4:
            cls+=1; acc+=counts[cls]*w[cls]
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

def memory_time_ou(H,beta,seed,tau,sigma,max_checks=10**7):
    _seed(seed); np.random.seed(seed); N=H.N   # seeds both the numba RNG and the Python RNG used for xi(0)
    q2s=np.array(H.q2s,dtype=np.int64); s2q=np.array(H.s2q,dtype=np.int64)
    e=np.zeros(N,np.int8); viol=np.zeros(H.L**3,np.int8); nviol=np.zeros(N,np.int8)
    base=np.array([np.exp(-4*beta),np.exp(-2*beta),1.0,1.0,1.0]); Tec=np.exp(4*beta)/100; t=0.0
    sig2=sigma**2; xi=np.random.normal()*sigma
    for n in range(1,max_checks):
        t,xi=bkl_until_ou(e,viol,nviol,q2s,s2q,base,t,n*Tec,xi,tau,sig2)
        if not viol.any():
            if e.any() and H.is_logical_error(e): return n*Tec,'logical',0
            continue
        P=H.decode(viol)
        if P is None: return n*Tec,'decoder',int(viol.sum())
        if H.is_logical_error((e+P)%2): return n*Tec,'logical',int(viol.sum())
    return np.nan,'timeout',int(viol.sum())

_H=None
def _worker(a):
    beta,seed,tau,sigma=a; return memory_time_ou(_H,beta,seed,tau,sigma)
def run(L,beta,nsamp,tau,sigma,nproc=16):
    global _H; _H=Haah(L); t0=time.time()
    with mp.get_context('fork').Pool(nproc) as pool:
        res=pool.map(_worker,[(beta,7000+1000*L+i,tau,sigma) for i in range(nsamp)])
    T=np.array([r[0] for r in res]); modes=[r[1] for r in res]; ok=~np.isnan(T); T=T[ok]; Tec=np.exp(4*beta)/100
    out={"beta":beta,"tau_c":tau,"sigma":sigma,"n":int(ok.sum()),"Tmem_mean":float(T.mean()),"Tmem_sem":float(T.std()/np.sqrt(len(T))),"Tmem_over_Tec":float(T.mean()/Tec),"modes":{m:modes.count(m) for m in set(modes)},"wall_s":time.time()-t0}
    print(f"L={L} beta={beta} tau_c={tau:g} sigma={sigma} n={len(T)} Tmem={T.mean():.3g} ± {T.std()/np.sqrt(len(T)):.2g} Tmem/Tec={T.mean()/Tec:.1f} modes={out['modes']} ({time.time()-t0:.0f}s)",flush=True)
    return out

if __name__=="__main__":
    if sys.argv[1]=="selftest":
        o=run(7,4.5,60,1.0,1e-6); assert abs(o['Tmem_mean']-7.1e7)<3*np.hypot(o['Tmem_sem'],8e6), o; print("selftest ok (sigma->0 reproduces flat bath)"); sys.exit()
    beta=float(sys.argv[1]); L=int(sys.argv[2]); nsamp=int(sys.argv[3]); sigma=float(sys.argv[4]); taus=[float(x) for x in sys.argv[5].split(',')]
    fn=__file__.rsplit('/',2)[0]+f"/results/h2ou_beta{beta}_L{L}_sigma{sigma}.json"
    try: out=json.load(open(fn))
    except FileNotFoundError: out={}
    for tau in taus: out[str(tau)]=run(L,beta,nsamp,tau,sigma); json.dump(out,open(fn,"w"),indent=1)   # merge + checkpoint per point
