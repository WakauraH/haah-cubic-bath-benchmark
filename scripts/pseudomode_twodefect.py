"""Gap item G-B: two defect qubits (gap Delta=4 each) coupled to (i) one shared pseudomode or (ii) two independent
pseudomodes, rotating-wave coupling g, mode frequency w0. Exact reduced Gibbs state of the full system.
Sector free energies F_ij = -(1/beta) ln Tr[Pi_ij e^{-beta H}] with Pi_ij projecting on defect occupations (i,j).
Reports eps_1 = F_10 - F_00 (single-defect free energy relative to Delta), U_12 = F_11 - F_10 - F_01 + F_00.
Additivity (uniform rescaling) requires U_12 = 0 and eps_1 = Delta + delta. Also reports the JC spectral shift
of the qubit-like one-excitation level for comparison."""
import numpy as np, itertools, json, warnings; warnings.filterwarnings('ignore')
Delta=4.0
def kron(*ops):
    out=np.array([[1.0]])
    for o in ops: out=np.kron(out,o)
    return out
def model(w0,g,nmax,shared):
    I2=np.eye(2); sz=np.diag([1.0,-1.0]); sp=np.array([[0,1],[0,0]],float); Pe=np.diag([1.0,0.0])
    a=np.diag(np.sqrt(np.arange(1,nmax+1)),1); Im=np.eye(nmax+1); N=a.T@a
    if shared:
        H=0.5*Delta*(kron(sz,I2,Im)+kron(I2,sz,Im))+w0*kron(I2,I2,N)+g*(kron(sp,I2,a)+kron(sp.T,I2,a.T)+kron(I2,sp,a)+kron(I2,sp.T,a.T))
        P={(i,j):kron(Pe if i else I2-Pe, Pe if j else I2-Pe, Im) for i in (0,1) for j in (0,1)}
    else:
        H=0.5*Delta*(kron(sz,I2,Im,Im)+kron(I2,sz,Im,Im))+w0*(kron(I2,I2,N,Im)+kron(I2,I2,Im,N))+g*(kron(sp,I2,a,Im)+kron(sp.T,I2,a.T,Im)+kron(I2,sp,Im,a)+kron(I2,sp.T,Im,a.T))
        P={(i,j):kron(Pe if i else I2-Pe, Pe if j else I2-Pe, Im, Im) for i in (0,1) for j in (0,1)}
    return H,P
def sector_free_energies(beta,w0,g,nmax=10,shared=True):
    H,P=model(w0,g,nmax,shared); E,V=np.linalg.eigh(H); w=np.exp(-beta*(E-E.min()))
    F={}
    for k,Pk in P.items():
        diag=np.real(np.einsum('ik,ij,jk->k',V.conj(),Pk,V)); Z=(w*diag).sum(); F[k]=-np.log(Z)/beta+E.min()
    return F
if __name__=="__main__":
    beta=4.5; out=[]
    print("beta=4.5, Delta=4.  eps_1 = F10-F00 (bare Delta=4), U12 = F11-F10-F01+F00, JC shift of qubit-like level = [sqrt(D^2+4g^2)-D]/2")
    for w0 in (2.0,6.0):
        for g in (0.3,0.6,1.0):
            D=Delta-w0; jc=(np.sqrt(D*D+4*g*g)-abs(D))/2*np.sign(D) if D!=0 else g
            for shared in (True,False):
                F=sector_free_energies(beta,w0,g,shared=shared)
                eps1=F[(1,0)]-F[(0,0)]; eps2=F[(0,1)]-F[(0,0)]; U=F[(1,1)]-F[(1,0)]-F[(0,1)]+F[(0,0)]
                # independent-mode reference for the shared case: shift of a single defect with its own mode
                out.append(dict(w0=w0,g=g,shared=shared,eps1=eps1,eps2=eps2,U12=U,jc_shift=jc))
                print(f"w0={w0} g={g} {'shared' if shared else 'indep '}: eps_1-Delta = {eps1-Delta:+.3f}  U_12 = {U:+.3f}   (JC level shift {jc:+.3f})")
    json.dump(out,open(__file__.rsplit('/',2)[0]+"/results/gap/B_twodefect_meanforce.json","w"),indent=1)
