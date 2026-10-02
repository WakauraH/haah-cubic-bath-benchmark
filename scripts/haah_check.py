"""Haah cubic code #1 on L^3 periodic lattice, 2 qubits/site.
X-stab (f, g), Z-stab (gbar, fbar); f=1+x+y+z, g=1+xy+yz+zx (Haah 2011 polynomial form).
Self-check: all stabilizers commute; report rank and k = N - rank."""
import numpy as np, itertools, sys
def rank_f2(M):
    M=M.copy()%2; r=0
    for c in range(M.shape[1]):
        p=np.nonzero(M[r:,c])[0]
        if len(p)==0: continue
        p=p[0]+r; M[[r,p]]=M[[p,r]]
        for i in np.nonzero(M[:,c])[0]:
            if i!=r: M[i]^=M[r]
        r+=1
        if r==M.shape[0]: break
    return r
def build(L):
    N=2*L**3; idx=lambda x,y,z,q: 2*((x%L)*L*L+(y%L)*L+(z%L))+q
    f=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]; g=[(0,0,0),(1,1,0),(0,1,1),(1,0,1)]
    neg=lambda P:[(-a,-b,-c) for a,b,c in P]
    X=[];Z=[]
    for x,y,z in itertools.product(range(L),repeat=3):
        sx=np.zeros(N,int); sz=np.zeros(N,int)
        for d in f: sx[idx(x+d[0],y+d[1],z+d[2],0)]^=1
        for d in g: sx[idx(x+d[0],y+d[1],z+d[2],1)]^=1
        for d in neg(g): sz[idx(x+d[0],y+d[1],z+d[2],0)]^=1
        for d in neg(f): sz[idx(x+d[0],y+d[1],z+d[2],1)]^=1
        X.append(sx); Z.append(sz)
    return np.array(X),np.array(Z),N
for L in [2,3,4,5,6]:
    X,Z,N=build(L)
    comm=(X@Z.T)%2
    ok=not comm.any()
    r=rank_f2(X)+rank_f2(Z)
    print(f"L={L} N={N} stabilizers={2*L**3} commute={ok} rank={r} k={N-r}")
    assert ok
