"""Quantum H2 test: one defect degree of freedom (qubit, gap Delta = 4 = energy of the 4 stabilizers flipped
by a single X error in the cubic code) coupled to a structured bath represented by ONE damped pseudomode
(Tamascelli et al. 2018): H = (Delta/2) sz + w0 a'a + g sx (a + a'),  mode damped at rate kappa into a
thermal bath at inverse temperature beta (n_th = 1/(e^{beta w0}-1)). Effective spectral density seen by the
qubit is Lorentzian centred at w0 with width kappa; memory time tau_c = 1/kappa.
Questions: (1) steady-state excited population P_exc vs Boltzmann P_B = e^{-beta Delta}/(1+e^{-beta Delta});
(2) does the deviation depend on kappa (memory) or only on g, w0 (Hamiltonian of mean force)?
(3) Markov limit kappa >> g: relaxation rate -> golden rule 2 pi J(Delta) (n+1), J Lorentzian.
Exact: vectorised Liouvillian (dim 2(nmax+1))^2, null vector + spectrum. numpy only.
Self-check: at g=0.02, kappa=1, w0=Delta the relaxation rate matches 4 g^2/kappa within 5%, and
P_exc -> P_B within 1% as g -> 0."""
import numpy as np, json, sys, itertools, warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)   # spurious BLAS matmul warnings on this machine
Delta=4.0
def ops(nmax):
    a=np.diag(np.sqrt(np.arange(1,nmax+1)),1); I2=np.eye(2); Im=np.eye(nmax+1)
    sz=np.diag([1.0,-1.0]); sx=np.array([[0,1],[1,0]],float)
    return np.kron(sz,Im),np.kron(sx,Im),np.kron(I2,a),np.kron(I2,a.T@a)
def liouvillian(H,Ls):
    d=H.shape[0]; I=np.eye(d); L=-1j*(np.kron(I,H)-np.kron(H.T,I))
    for c in Ls:
        cd=c.conj().T; L+=np.kron(c.conj(),c)-0.5*np.kron(I,cd@c)-0.5*np.kron((cd@c).T,I)
    return L
def steady(L,d):
    A=L.copy(); b=np.zeros(d*d,complex); A[0,:]=np.eye(d).reshape(-1,order='F'); b[0]=1     # replace one eq by Tr rho = 1
    rho=np.linalg.solve(A,b).reshape(d,d,order='F'); return rho
def solve(beta,g,kappa,w0,nmax=12,rwa=True):
    """rwa=True: H_int = g (s+ a + s- a'), excitation-conserving (Jaynes-Cummings); artifact-free with a local
    thermal Lindblad on the mode. rwa=False keeps counter-rotating terms: the local damping then converts the
    ground state's virtual dressing (eps = g/(Delta+w0)) into real qubit excitations at rate kappa*eps^2."""
    SZ,SX,A,N=ops(nmax); sp=np.kron(np.array([[0,1],[0,0]],float),np.eye(nmax+1))   # s+ : |g>->|e> with sz=+1 excited
    Hint=g*(sp@A+sp.T@A.T) if rwa else g*SX@(A+A.T)
    H=0.5*Delta*SZ+w0*N+Hint
    nth=1/np.expm1(beta*w0)
    Ls=[np.sqrt(kappa*(nth+1))*A, np.sqrt(kappa*nth)*A.T]
    d=H.shape[0]; L=liouvillian(H,Ls); rho=steady(L,d)
    Pexc=np.real(np.trace(rho@((np.eye(d)+SZ)/2)))            # sz=+1 = excited (energy +Delta/2)
    # population relaxation rate: slowest Liouvillian mode that actually carries weight in P_exc(t)
    # starting from the qubit excited x thermal mode
    ev,V=np.linalg.eig(L); O=((np.eye(d)+SZ)/2).reshape(-1,order='F')
    nb=nth**np.arange(nmax+1)/(nth+1)**(np.arange(nmax+1)+1); rho0=np.kron(np.diag([1.0,0.0]),np.diag(nb))
    c=np.linalg.solve(V,rho0.reshape(-1,order='F')); w=np.abs(c*(O@V))
    sel=[k for k in range(len(ev)) if w[k]>1e-9*w.max() and abs(ev[k])>1e-10]
    rate=min(-ev[k].real for k in sel)
    # global Gibbs (mean force) reference and bare Boltzmann
    E,V=np.linalg.eigh(H); w=np.exp(-beta*(E-E.min())); w/=w.sum(); rhoG=(V*w)@V.conj().T
    PG=np.real(np.trace(rhoG@((np.eye(d)+SZ)/2))); PB=np.exp(-beta*Delta)/(1+np.exp(-beta*Delta))
    return dict(Pexc=Pexc,P_meanforce=PG,P_boltz=PB,rate=rate,nth=nth,rwa=rwa)
def golden(g,kappa,w0,nth): J=g*g*kappa/((Delta-w0)**2+kappa**2/4); return J*(nth+1)+J*nth
if __name__=="__main__":
    r=solve(4.5,0.02,1.0,4.0); gr=golden(0.02,1.0,4.0,r['nth'])
    assert abs(r['rate']/gr-1)<0.05, (r['rate'],gr); assert abs(r['Pexc']/r['P_boltz']-1)<0.01, r
    print(f"selftest ok (RWA, resonant): rate {r['rate']:.3e} vs golden rule {gr:.3e}; P_exc/P_B = {r['Pexc']/r['P_boltz']:.4f}")
    rc=solve(4.5,0.02,1.0,4.0,rwa=False); eps2=(0.02/8)**2
    print(f"counter-rotating + local thermal damping ARTIFACT: P_exc/P_B = {rc['Pexc']/rc['P_boltz']:.3g}; predicted kappa*eps^2/Gamma / P_B = {1.0*eps2/rc['rate']/rc['P_boltz']:.3g}")
    out={"dynamic_resonant":[], "meanforce_static":[]}
    # (1) dynamics: memory time 1/kappa and coupling g at resonance (detailed balance exact at the qubit gap)
    for beta in (2.0,4.5):
        for g in (0.05,0.2,0.5,1.0):
            for kappa in (0.03,0.1,0.3,1.0,3.0,10.0):
                r=solve(beta,g,kappa,Delta); r.update(beta=beta,w0=Delta,g=g,kappa=kappa); out["dynamic_resonant"].append(r)
                print(f"[dyn] beta={beta} g={g} kappa={kappa:5}: P_exc/P_B={r['Pexc']/r['P_boltz']:.4f}  meanforce/P_B={r['P_meanforce']/r['P_boltz']:.4f}  rate={r['rate']:.2e} (golden {golden(g,kappa,Delta,r['nth']):.2e})")
    # (2) static mean-force (exact reduced Gibbs state of qubit+mode, RWA): detuning-induced gap renormalisation
    for beta in (2.0,4.5):
        for w0 in (1.0,2.0,3.0,4.0,5.0,6.0,8.0):
            for g in (0.1,0.3,0.6,1.0):
                SZ,SX,A,N=ops(12); sp=np.kron(np.array([[0,1],[0,0]],float),np.eye(13)); H=0.5*Delta*SZ+w0*N+g*(sp@A+sp.T@A.T)
                E,V=np.linalg.eigh(H); w=np.exp(-beta*(E-E.min())); w/=w.sum(); d=H.shape[0]
                PG=np.real(np.trace(((V*w)@V.conj().T)@((np.eye(d)+SZ)/2))); PB=np.exp(-beta*Delta)/(1+np.exp(-beta*Delta))
                # dressed defect population: thermal weight of eigenstates whose qubit-excited character exceeds 1/2
                qexc=np.real(np.einsum('ik,ij,jk->k',V.conj(),(np.eye(d)+SZ)/2,V)); Pdr=w[qexc>0.5].sum()
                lamb=g*g/(Delta-w0) if w0!=Delta else 0.0     # 2nd-order shift of the qubit-excited level (RWA)
                out["meanforce_static"].append(dict(beta=beta,w0=w0,g=g,P_meanforce=PG,P_dressed=Pdr,P_boltz=PB,ratio=PG/PB,ratio_dressed=Pdr/PB,lamb_pred=np.exp(-beta*lamb)))
                print(f"[mf] beta={beta} w0={w0} g={g}: bare P_mf/P_B={PG/PB:.3g}  dressed P_dr/P_B={Pdr/PB:.3f}  Lamb prediction e^(-beta g^2/(Delta-w0))={np.exp(-beta*lamb):.3f}")
    json.dump(out,open(__file__.rsplit('/',2)[0]+"/results/pseudomode_gap.json","w"),indent=1)
