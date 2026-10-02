"""Emits results/gap/pseudomode_uniqueness.json: over the resonant RWA grid of pseudomode_gap.py (beta in {2,4.5},
g in {0.05,0.2,0.5,1.0}, kappa in {0.03,...,10}), the number of Liouvillian eigenvalues with |lambda|<1e-10 (must be 1)
and the smallest nonzero |Re lambda| (spectral gap), documenting the numerical uniqueness check of Appendix B.6 step 5."""
import sys,os,json,warnings,numpy as np; warnings.filterwarnings('ignore'); sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pseudomode_gap as pm
out=[]
for beta in (2.0,4.5):
    for g in (0.05,0.2,0.5,1.0):
        for kappa in (0.03,0.1,0.3,1.0,3.0,10.0):
            SZ,SX,A,N=pm.ops(12); sp=np.kron(np.array([[0,1],[0,0]],float),np.eye(13)); H=0.5*pm.Delta*SZ+4.0*N+g*(sp@A+sp.T@A.T)
            nth=1/np.expm1(beta*4.0); L=pm.liouvillian(H,[np.sqrt(kappa*(nth+1))*A,np.sqrt(kappa*nth)*A.T]); ev=np.linalg.eigvals(L)
            nzero=int((np.abs(ev)<1e-10).sum()); gap=float(np.min(-ev.real[np.abs(ev)>=1e-10]))
            out.append(dict(beta=beta,g=g,kappa=kappa,n_zero_modes=nzero,spectral_gap=gap))
json.dump(out,open(os.path.dirname(os.path.abspath(__file__))+'/../results/gap/pseudomode_uniqueness.json','w'),indent=1)
print('zero modes per point:',sorted({o['n_zero_modes'] for o in out}),' min spectral gap:',min(o['spectral_gap'] for o in out))
