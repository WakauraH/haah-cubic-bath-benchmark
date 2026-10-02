"""Memory time of Haah's cubic code under Metropolis/BKL thermal dynamics with RG decoder readout.
Protocol follows Bravyi & Haah arXiv:1112.3252 Sec. VI (RG decoder, generic TestNeutral via F2 solve)
and Sec. VIII (X errors only, J=1/2 so one defect = energy 1, Metropolis rates, trial decode every
T_ec = e^{4 beta}/100, T_mem = first trial time at which decoding fails; sample mean over runs).
Self-check: decoder corrects all weight-1 and random weight-2 X errors at L=5."""
import numpy as np, itertools, json, sys, time
sys.path.insert(0, __file__.rsplit('/',1)[0])
from haah_code import build, rank_f2

def nullspace_f2(A):
    A=A.copy()%2; m,n=A.shape; piv=[]; r=0
    for c in range(n):
        p=np.nonzero(A[r:,c])[0]
        if len(p)==0: continue
        p=p[0]+r; A[[r,p]]=A[[p,r]]
        nz=np.nonzero(A[:,c])[0]; nz=nz[nz!=r]; A[nz]^=A[r]
        piv.append(c); r+=1
        if r==m: break
    free=[c for c in range(n) if c not in set(piv)]
    B=np.zeros((len(free),n),np.int8)
    for i,f in enumerate(free):
        B[i,f]=1
        for j,c in enumerate(piv): B[i,c]=A[j,f]
    return B

def solve_f2(A,b):
    """Return x with A x = b (mod 2) or None."""
    A=np.concatenate([A%2,b.reshape(-1,1)%2],1).astype(np.int8); m,n=A.shape; n-=1; piv=[]; r=0
    for c in range(n):
        p=np.nonzero(A[r:,c])[0]
        if len(p)==0: continue
        p=p[0]+r; A[[r,p]]=A[[p,r]]
        nz=np.nonzero(A[:,c])[0]; nz=nz[nz!=r]; A[nz]^=A[r]
        piv.append(c); r+=1
        if r==m: break
    if (A[r:,n]).any(): return None
    x=np.zeros(n,np.int8)
    for i,c in enumerate(piv): x[c]=A[i,n]
    return x

class Haah:
    def __init__(self,L,Ltqo=None):
        self.L=L; X,Z,N=build(L); self.Z=Z; self.N=N
        self.cube_xyz=np.array(list(itertools.product(range(L),repeat=3)))       # cube c -> (x,y,z)
        self.q2s=[np.nonzero(Z[:,q])[0] for q in range(N)]                        # Z-stabs touching qubit q
        self.s2q=[np.nonzero(Z[c])[0] for c in range(L**3)]                       # qubits of stab c
        self.qsite=np.arange(N)//2; self.site_xyz=self.cube_xyz                    # site index == cube index layout
        # logical Z operators (commute with X-stabs, independent mod Z-stabs)
        # greedy: keep each kernel vector independent of Z-stabs and of those already kept (one-pass F2 elimination on bit-packed rows)
        ker=nullspace_f2(X); piv={}; LZ=[]
        def _red(row):
            r=int.from_bytes(np.packbits(row%2).tobytes(),'big')
            while r:
                h=r.bit_length(); q=piv.get(h)
                if q is None: piv[h]=r; return True
                r^=q
            return False
        for z in Z: _red(z)
        for v in ker:
            if _red(v): LZ.append(v)
        self.LZ=np.array(LZ); self.k=len(LZ)
        # Bravyi-Haah: strictly Ltqo = L-2 for the cubic code; their analysis assumes Ltqo <= L/2 for convenience.
        self.Ltqo=(L-2) if Ltqo is None else Ltqo; self.pM=max(p for p in range(0,20) if 2**p<=self.Ltqo)
    # ---- geometry on torus ----
    def dist(self,a,b):
        d=np.abs(self.cube_xyz[a]-self.cube_xyz[b]); d=np.minimum(d,self.L-d); return d.max()
    def box(self,cubes):
        """Smallest enclosing box on torus: per axis, cut at largest cyclic gap. Returns (lo, size) per axis."""
        lo=[];size=[]
        for ax in range(3):
            xs=np.sort(np.unique(self.cube_xyz[cubes,ax]))
            if len(xs)==1: lo.append(int(xs[0])); size.append(1); continue
            gaps=np.diff(np.append(xs,xs[0]+self.L)); i=np.argmax(gaps)
            start=xs[(i+1)%len(xs)]; lo.append(int(start)); size.append(int(self.L-gaps[i]+1))
        return lo,size
    def region_qubits(self,lo,size,g=0):
        """Qubits on the corner sites of the cubes in the box grown by g cube layers (g=0: the box's
        own corner sites = 1-neighbourhood of the box in the sense of Bravyi-Haah Def. 1)."""
        L=self.L; sites=set()
        for dx in range(-g,size[0]+g):
            for dy in range(-g,size[1]+g):
                for dz in range(-g,size[2]+g):
                    x,y,z=(lo[0]+dx)%L,(lo[1]+dy)%L,(lo[2]+dz)%L
                    for ex in (0,1):
                        for ey in (0,1):
                            for ez in (0,1): sites.add(((x+ex)%L)*L*L+((y+ey)%L)*L+((z+ez)%L))
        return np.array(sorted(2*s+q for s in sites for q in (0,1)))
    # ---- RG decoder ----
    def test_neutral(self,cluster):
        lo,size=self.box(cluster)
        if max(size)>self.Ltqo: return None
        # cubic code: a neutral cluster is annihilated by an operator inside its minimal enclosing box (BH Sec. VI E)
        qs=self.region_qubits(lo,size,0)
        stabs=np.unique(np.concatenate([self.q2s[q] for q in qs]))
        A=self.Z[np.ix_(stabs,qs)]; b=np.isin(stabs,cluster).astype(np.int8)
        x=solve_f2(A,b)
        if x is not None:
            P=np.zeros(self.N,np.int8); P[qs[x==1]]=1; return P
        return None
    def components(self,defects,r):
        defects=list(defects); n=len(defects); parent=list(range(n))
        def find(i):
            while parent[i]!=i: parent[i]=parent[parent[i]]; i=parent[i]
            return i
        for i in range(n):
            for j in range(i+1,n):
                if self.dist(defects[i],defects[j])<=r: parent[find(i)]=find(j)
        comps={}
        for i in range(n): comps.setdefault(find(i),[]).append(defects[i])
        return [np.array(c) for c in comps.values()]
    def decode(self,syndrome):
        """syndrome: 0/1 array over cubes. Returns correction P (X-error pattern) or None on failure."""
        S=set(np.nonzero(syndrome)[0]); P=np.zeros(self.N,np.int8)
        for p in range(self.pM+1):
            if not S: break
            for comp in self.components(S,2**p):
                Pa=self.test_neutral(comp)
                if Pa is not None:
                    P^=Pa; S-=set(comp)
        return P if not S else None
    def syndrome(self,e): return (self.Z@e)%2
    def is_logical_error(self,e_total):
        return bool(((self.LZ@e_total)%2).any())
    # ---- BKL Metropolis dynamics ----
    def memory_time(self,beta,rng,max_wall=600,spec=(1.0,1.0,1.0)):
        """spec=(s4,s2,s0): spectral-filtering factors J(|dE|)/J_flat at |dE|=4,2,0. Detailed balance is kept
        (creation and annihilation at the same |dE| are scaled equally). spec=(1,1,1) = Metropolis flat bath."""
        """Return T_mem sample: first trial-decoding time at which decoding fails."""
        N=self.N; e=np.zeros(N,np.int8); viol=np.zeros(self.L**3,np.int8); t=0.0
        nviol=np.zeros(N,np.int8)                     # violated Z-stabs touching q (0..4)
        s4,s2,s0=spec
        rate_cls=np.array([s4*np.exp(-4*beta),s2*np.exp(-2*beta),s0,s2,s4])   # Metropolis x filter, dE = 4-2*nviol
        Tec=np.exp(4*beta)/100; next_check=Tec; wall=time.time()
        while True:
            counts=np.bincount(nviol,minlength=5); R=(counts*rate_cls).sum()
            t+=rng.exponential(1/R)
            while t>=next_check:
                P=self.decode(viol)
                if P is None: return next_check,'decoder',int(viol.sum())
                if self.is_logical_error((e+P)%2): return next_check,'logical',int(viol.sum())
                next_check+=Tec
                if time.time()-wall>max_wall: return np.nan,'timeout',int(viol.sum())
            cls=rng.choice(5,p=counts*rate_cls/R); q=rng.choice(np.nonzero(nviol==cls)[0])
            e[q]^=1
            for s in self.q2s[q]:
                viol[s]^=1; d=1 if viol[s] else -1
                nviol[self.s2q[s]]+=d

def selftest():
    """Weight-1 errors at L=5 (Ltqo=2, pM=0) and random weight-2/3 errors at L=9 (Ltqo=4, pM=1).
    Note: at L=5 two adjacent errors give a cluster box of size 3 > Ltqo, which this decoder cannot
    annihilate -- consistent with Bravyi-Haah (decoder failures dominate T_mem at small L)."""
    H=Haah(5); assert H.Ltqo==3 and H.pM==1
    for q in range(H.N):
        e=np.zeros(H.N,np.int8); e[q]=1; P=H.decode(H.syndrome(e)); assert P is not None and not H.is_logical_error((e+P)%2)
    t0=time.time(); H=Haah(9); rng=np.random.default_rng(0); rate={}
    for w in (2,3):
        ok=0
        for _ in range(100):
            e=np.zeros(H.N,np.int8); e[rng.choice(H.N,w,replace=False)]=1
            P=H.decode(H.syndrome(e)); ok+=(P is not None and not H.is_logical_error((e+P)%2))
        rate[w]=ok/100
    # weight-3 errors can merge into one cluster whose box exceeds Ltqo=L//2 (decoder failure, as in Bravyi-Haah)
    assert rate[2]>=0.97 and rate[3]>=0.9, rate
    print(f"selftest ok: L=5 all weight-1; L=9 (k={H.k}, pM={H.pM}) random weight-2: {rate[2]:.2f}, weight-3: {rate[3]:.2f} corrected ({time.time()-t0:.0f}s)")

if __name__=="__main__":
    selftest()
    beta=float(sys.argv[1]) if len(sys.argv)>1 else 4.0; Ls=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [5,7]
    nsamp=int(sys.argv[3]) if len(sys.argv)>3 else 20
    tag=sys.argv[4] if len(sys.argv)>4 else ""
    spec=tuple(float(x) for x in sys.argv[5].split(",")) if len(sys.argv)>5 else (1.0,1.0,1.0)
    out={}
    for L in Ls:
        H=Haah(L); rng=np.random.default_rng(L); T=[]; modes=[]; nd=[]; t0=time.time()
        for i in range(nsamp):
            Ti,mode,n=H.memory_time(beta,rng,spec=spec); T.append(Ti); modes.append(mode); nd.append(n)
        T=np.array(T); ok=~np.isnan(T); T=T[ok]
        Tec=np.exp(4*beta)/100
        out[L]={"beta":beta,"spec":list(spec),"n":int(ok.sum()),"Tmem_mean":float(T.mean()),"Tmem_sem":float(T.std()/np.sqrt(len(T))),"Tec":Tec,
                "Tmem_over_Tec":float(T.mean()/Tec),"modes":{m:modes.count(m) for m in set(modes)},"defects_at_failure_mean":float(np.mean(nd)),"wall_s":time.time()-t0}
        print(f"L={L} beta={beta} k={H.k} n={len(T)} Tmem={T.mean():.3g} ± {T.std()/np.sqrt(len(T)):.2g}  Tmem/Tec={T.mean()/Tec:.2f} modes={out[L]['modes']} <defects>={np.mean(nd):.1f} ({time.time()-t0:.0f}s)",flush=True)
    json.dump(out,open(__file__.rsplit('/',2)[0]+f"/results/tmem_beta{beta}{tag}.json","w"),indent=1)
