"""Tier 1: Markovian Davies-rate syndrome dynamics of Haah's cubic code (X-error sector).

Model: H = -Delta * sum_c S^Z_c ; X flip on qubit q changes energy by dE = 2*Delta*(4 - 2*n_viol(q))
where n_viol(q) = # violated Z-stabilizers containing q (each qubit sits in 4 Z-stabilizers).
Rate (Glauber/detailed balance, flat spectral density): gamma(dE) = gamma0 / (1 + exp(beta*dE)).
Memory time: tau_mem = first t with P_success(t) < 1-1/e, P_success = prob. that
min-weight readout decoder returns the correct logical class.

Two solvers, same chain:
  exact : master equation on all 2^N Pauli frames (L=2, N=16) by uniformization (Appendix B.10 of the paper).
  mc    : Gillespie trajectories (scales to large L).
Self-check: exact vs mc agree at L=2.
"""
import numpy as np, itertools, json, sys, time
from scipy.sparse import csr_matrix
sys.path.insert(0, __file__.rsplit('/',1)[0])
from haah_check import build, rank_f2

def nullspace_f2(A):
    """Basis of {v : A v = 0 mod 2}."""
    A=A.copy()%2; m,n=A.shape; piv=[]; r=0
    for c in range(n):
        p=np.nonzero(A[r:,c])[0]
        if len(p)==0: continue
        p=p[0]+r; A[[r,p]]=A[[p,r]]
        for i in range(m):
            if i!=r and A[i,c]: A[i]^=A[r]
        piv.append(c); r+=1
        if r==m: break
    free=[c for c in range(n) if c not in piv]
    basis=[]
    for f in free:
        v=np.zeros(n,int); v[f]=1
        for i,c in enumerate(piv): v[c]=A[i,f]
        basis.append(v)
    return np.array(basis)

def logicals_Z(X,Z):
    """Logical Z operators: commute with X-stabilizers, independent modulo Z-stabilizer rowspace."""
    ker=nullspace_f2(X)                      # all Z-type ops commuting with X stabs
    base=rank_f2(Z); L=[]
    M=Z.copy()
    for v in ker:
        if rank_f2(np.vstack([M,v]))>rank_f2(M): M=np.vstack([M,v]); L.append(v)
    return np.array(L)

class HaahDavies:
    def __init__(self,L,beta,Delta=1.0,gamma0=1.0):
        self.L=L; self.beta=beta; self.Delta=Delta; self.g0=gamma0
        X,Z,N=build(L); self.N=N; self.Z=Z
        self.LZ=logicals_Z(X,Z); self.k=len(self.LZ)
        self.q2s=[np.nonzero(Z[:,q])[0] for q in range(N)]   # Z-stabs containing qubit q
        assert all(len(s)==4 for s in self.q2s)
        self.rate=np.array([gamma0/(1+np.exp(beta*2*Delta*(4-2*nv))) for nv in range(5)])  # index by n_viol
    # ---- decoder (lookup, small L only) ----
    def build_lookup(self):
        N=self.N; E=np.arange(2**N,dtype=np.int64)
        bits=((E[:,None]>>np.arange(N))&1).astype(np.int8)          # 2^N x N
        syn=(bits@self.Z.T)%2; cls=(bits@self.LZ.T)%2
        w=bits.sum(1)
        synid=syn@(1<<np.arange(self.Z.shape[0])); clsid=cls@(1<<np.arange(self.k))
        # min-weight representative per syndrome -> its class
        order=np.lexsort((w,synid)); first=np.unique(synid[order],return_index=True)[1]
        dec=dict(zip(synid[order][first],clsid[order][first]))
        self.success=np.array([dec[s]==c for s,c in zip(synid,clsid)],dtype=float)
        self.bits=bits; self.synid=synid
    # ---- exact master equation ----
    def exact(self,ts):
        """Uniformization: P(t)=sum_n Pois(n;Lam t) (I+W/Lam)^n p0, Lam = N*gamma0 >= max out-rate."""
        N=self.N; S=2**N; bits=self.bits
        nviol=((bits@self.Z.T)%2)
        rows=[];cols=[];vals=[]
        for q in range(N):
            nv=nviol[:,self.q2s[q]].sum(1)
            rows.append(np.arange(S)^(1<<q)); cols.append(np.arange(S)); vals.append(self.rate[nv])
        rows=np.concatenate(rows); cols=np.concatenate(cols); vals=np.concatenate(vals)
        W=csr_matrix((vals,(rows,cols)),shape=(S,S))
        out=np.asarray(W.sum(0)).ravel(); Lam=N*self.g0
        Q=W/Lam+csr_matrix((1-out/Lam,(np.arange(S),np.arange(S))),shape=(S,S))   # stochastic
        ts=np.asarray(ts); Lt=Lam*ts; nmax=int(Lt.max()+10*np.sqrt(Lt.max())+20)
        logw=-Lt[:,None]+np.arange(nmax+1)[None,:]*np.log(Lt[:,None])-np.cumsum(np.log(np.maximum(np.arange(nmax+1),1)))[None,:]
        w=np.exp(logw)                                   # Poisson weights, ts x (nmax+1)
        v=np.zeros(S); v[0]=1; acc=np.zeros(len(ts))
        for n in range(nmax+1):
            acc+=w[:,n]*(v@self.success)
            if w[:,n:].max()<1e-14 and n>Lt.max(): break
            v=Q@v
        return acc
    # ---- Gillespie ----
    def mc(self,ts,ntraj,seed=0):
        rng=np.random.default_rng(seed); N=self.N; succ=np.zeros(len(ts))
        Zc=csr_matrix(self.Z)
        for _ in range(ntraj):
            e=np.zeros(N,np.int8); viol=np.zeros(self.Z.shape[0],np.int8); t=0.0; i=0
            while i<len(ts):
                nv=np.array([viol[s].sum() for s in self.q2s]); r=self.rate[nv]; R=r.sum()
                dt=rng.exponential(1/R); 
                while i<len(ts) and ts[i]<t+dt:
                    succ[i]+=self._success(e); i+=1
                t+=dt; q=rng.choice(N,p=r/R); e[q]^=1; viol[self.q2s[q]]^=1
        return succ/ntraj
    def _success(self,e):
        idx=int(e.astype(np.int64)@(1<<np.arange(self.N)))
        return float(self.success[idx])

def tau_mem(ts,P,thr=1-1/np.e):
    i=np.argmax(P<thr); return float(ts[i]) if P[i]<thr else np.inf

if __name__=="__main__":
    out={}
    for beta in [0.5,0.75,1.0,1.25,1.5]:
        m=HaahDavies(2,beta); m.build_lookup()
        ts=np.geomspace(1e-2,3e2,50)*np.exp(8*(beta-1.0)) if beta<=1 else np.geomspace(1e-2,3e2,50)*min(np.exp(8*(beta-1.0)),3e2)
        t0=time.time(); Pe=m.exact(ts); te=time.time()-t0
        t0=time.time(); Pm=m.mc(ts,ntraj=800,seed=1); tm=time.time()-t0
        err=np.abs(Pe-Pm).max()
        out[beta]={"k":m.k,"tau_exact":tau_mem(ts,Pe),"tau_mc":tau_mem(ts,Pm),"max_abs_err":float(err),
                   "t_exact_s":te,"t_mc_s":tm,"ts":ts.tolist(),"P_exact":Pe.tolist(),"P_mc":Pm.tolist()}
        print(f"beta={beta} k={m.k} tau_exact={out[beta]['tau_exact']:.3g} tau_mc={out[beta]['tau_mc']:.3g} max|dP|={err:.3f} ({te:.1f}s / {tm:.1f}s)")
        assert err<0.06, "MC disagrees with exact master equation"
    json.dump(out,open(__file__.rsplit('/',2)[0]+"/results/tier1_validation.json","w"),indent=1)
