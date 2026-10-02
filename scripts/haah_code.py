"""Haah's 3D cubic code (Code 1) built from the corner table of Bravyi & Haah, Fig. 1
(arXiv:1112.3252). Coordinates: x = left->right, y = front->back, z = bottom->top.
Corner (dx,dy,dz) -> two-qubit Pauli string; 'ZI' = Z on species 0, 'IZ' = Z on species 1.
Self-check: all stabilizers commute; k(L) matches Haah 2011 Table II
(k = 2[1 - 2q2 + 2^{r+1}(q2 + 12 q15 + 60 q63)]); e.g. k(2)=6, k(3)=2, k(4)=14, k(15)=50."""
import numpy as np, itertools
GZ={(0,0,0):'IZ',(1,0,0):'ZI',(0,0,1):'ZI',(1,0,1):'ZZ',(0,1,0):'II',(1,1,0):'IZ',(0,1,1):'IZ',(1,1,1):'ZI'}
GX={(0,0,0):'IX',(1,0,0):'XI',(0,0,1):'XI',(1,0,1):'II',(0,1,0):'XX',(1,1,0):'IX',(0,1,1):'IX',(1,1,1):'XI'}
def site(L,x,y,z): return (x%L)*L*L+(y%L)*L+(z%L)
def build(L):
    N=2*L**3; X=np.zeros((L**3,N),np.int8); Z=np.zeros((L**3,N),np.int8)
    for c,(x,y,z) in enumerate(itertools.product(range(L),repeat=3)):
        for (dx,dy,dz),p in GZ.items():
            s=site(L,x+dx,y+dy,z+dz)
            for q in range(2):
                if p[q]=='Z': Z[c,2*s+q]^=1
        for (dx,dy,dz),p in GX.items():
            s=site(L,x+dx,y+dy,z+dz)
            for q in range(2):
                if p[q]=='X': X[c,2*s+q]^=1
    return X,Z,N
def rank_f2(M):
    rows=[int(''.join(map(str,r)),2) for r in M]; r=0
    while rows:
        p=max(rows)
        if p==0: break
        r+=1; hb=p.bit_length()-1
        rows=[m^p if (m>>hb)&1 else m for m in rows if m!=p]
    return r
def k_table(L):
    q=lambda n: int(L%n==0); r=0
    while L%(2**(r+1))==0: r+=1
    return 2*(1-2*q(2)+2**(r+1)*(q(2)+12*q(15)+60*q(63)))
if __name__=="__main__":
    for L in [2,3,4,5,6,7,15]:
        X,Z,N=build(L); assert not ((X@Z.T)%2).any()
        k=N-rank_f2(X)-rank_f2(Z); print(f"L={L} N={N} k={k} table={k_table(L)}"); assert k==k_table(L)
    print("Code 1 verified against Haah 2011 Table II")
