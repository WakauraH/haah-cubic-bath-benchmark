"""Generates results/gap/pseudomode_cutoff_sweep.json: oscillator-cutoff convergence of the one-defect (P_exc/P_B, rate at
g=1, kappa=0.03, beta=4.5, resonant) and two-defect (eps1-Delta, U12 at w0=2, g=1, beta=4.5, shared mode) pseudomode models."""
import sys,os,json,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import pseudomode_gap as pm, pseudomode_twodefect as td
out={'one_defect_cutoff':{},'two_defect_cutoff':{}}
for nmax in (8,12,20,28):
    r=pm.solve(4.5,1.0,0.03,4.0,nmax=nmax); out['one_defect_cutoff'][nmax]={'Pexc_over_PB':r['Pexc']/r['P_boltz'],'rate':r['rate']}
for nmax in (6,10,14,20):
    F=td.sector_free_energies(4.5,2.0,1.0,nmax=nmax,shared=True); out['two_defect_cutoff'][nmax]={'eps1_minus_Delta':F[(1,0)]-F[(0,0)]-4,'U12':F[(1,1)]-F[(1,0)]-F[(0,1)]+F[(0,0)]}
json.dump(out,open(os.path.dirname(os.path.abspath(__file__))+'/../results/gap/pseudomode_cutoff_sweep.json','w'),indent=1,default=float); print('sweep written')
