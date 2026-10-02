"""Canonical analysis script (Protocol F). Regenerates the derived numbers quoted in the main text and Table I-III from results/*.json
(Appendix A validation numbers come from results/tier1_validation.json and the *.log files listed in the manifest; the pseudomode
and two-defect tables from scripts/pseudomode_gap.py and scripts/pseudomode_twodefect.py):
size exponents (weighted log-log LS, window 5-11 primary, 7-13 alternate, weights (T/SEM)^2, error scaled by
sqrt(max(1,chi2/dof)), window systematic |p1-p2|/2), Lambda_passive with numerator+denominator SEM, three-point
channel sensitivities over s in {0.1,1,10}, OU gains with denominator SEM, log-normal hazard averages (annealed and
quenched, pure power law and measured curve), quenched CV, creation-limited crossover. Output: results/derived_numbers.json"""
import json,glob,os,numpy as np
R=os.path.dirname(os.path.abspath(__file__))+'/../results/'
def load(f): return json.load(open(R+f))
def wls(d,Ls):
    Ls=np.array(Ls); T=np.array([d[str(L)]['Tmem_mean'] for L in Ls]); E=np.array([d[str(L)]['Tmem_sem'] for L in Ls]); s=E/T
    X=np.vstack([np.ones(len(Ls)),np.log(Ls)]).T; W=np.diag(1/s**2); cov=np.linalg.inv(X.T@W@X); b=cov@X.T@W@np.log(T)
    r=(np.log(T)-X@b)/s; chi2dof=(r**2).sum()/(len(Ls)-2); e=np.sqrt(cov[1,1])
    return dict(p=b[1],err=e,chi2dof=chi2dof,err_pdg=e*np.sqrt(max(1,chi2dof)),residuals=r.tolist())
out={'fits':{},'lambda':{},'sensitivity':{},'ou':{},'hazard':{},'crossover':{}}
series={'flat_b4.5':'tmem_beta4.5_fast_flat.json','flat_b5.0':'tmem_beta5.0_fast_flat.json','s4_0.01':'tmem_beta4.5_fast_s4_0.01.json','s4_10':'tmem_beta4.5_fast_s4_10.json','s2_0.1':'tmem_beta4.5_fast_s2_0.1.json','s0_0.1':'tmem_beta4.5_fast_s0_0.1.json','s0_0.01':'tmem_beta4.5_fast_s0_0.01.json','bp4.70':'tmem_beta4.7_fast_flat_lamb.json','bp4.875':'tmem_beta4.875_fast_flat_lamb.json','bp5.06':'tmem_beta5.06_fast_flat_lamb.json'}
for k,f in series.items():
    d=load(f); Ls=sorted(int(x) for x in d); w1=[L for L in (5,7,9,11) if L in Ls]; w2=[L for L in (7,9,11,13) if L in Ls]
    f1=wls(d,w1); f2=wls(d,w2) if len(w2)>=3 else None
    out['fits'][k]={'w5_11':f1,'w7_13':f2,'window_syst':(abs(f1['p']-f2['p'])/2 if f2 else None)}
    T=[d[str(L)]['Tmem_mean'] for L in Ls]; E=[d[str(L)]['Tmem_sem'] for L in Ls]
    out['lambda'][k]={f"{Ls[i]}->{Ls[i+1]}":(T[i+1]/T[i], T[i+1]/T[i]*np.hypot(E[i+1]/T[i+1],E[i]/T[i])) for i in range(len(Ls)-1)}
flat=load('tmem_beta4.5_fast_flat.json')
runs={os.path.basename(f)[len('tmem_beta4.5_fast_'):-5]:load(os.path.basename(f)) for f in glob.glob(R+'tmem_beta4.5_fast_*.json')}
for L in (9,13):
    for ch,idx in (('s4',0),('s2',1),('s0',2)):
        pts=sorted((r[str(L)]['spec'][idx],r[str(L)]['Tmem_mean']) for r in runs.values() if str(L) in r and all(r[str(L)]['spec'][j]==1.0 for j in range(3) if j!=idx) and r[str(L)]['spec'][idx] in (0.1,1.0,10.0))
        x=np.log([p[0] for p in pts]); y=np.log([p[1] for p in pts]); out['sensitivity'][f'{ch}_L{L}']=float(np.polyfit(x,y,1)[0])
F9=flat['9']
for sig in ('1.0','2.0'):
    try: d=load(f'h2ou_beta4.5_L9_sigma{sig}.json')
    except FileNotFoundError: continue
    for k in sorted(d,key=float):
        g=d[k]['Tmem_mean']/F9['Tmem_mean']; eg=g*np.hypot(d[k]['Tmem_sem']/d[k]['Tmem_mean'],F9['Tmem_sem']/F9['Tmem_mean']); cv=d[k]['Tmem_sem']*np.sqrt(d[k]['n'])/d[k]['Tmem_mean']
        out['ou'][f'sigma{sig}_tau{float(k):g}']={'gain':g,'err':eg,'cv':cv,'n':d[k]['n']}
a=0.6
for sig in (1.0,2.0):
    out['hazard'][f'sigma{sig}']={'annealed_powerlaw':float(np.exp(a*(1-a)*sig**2/2)),'quenched_powerlaw':float(np.exp(a*(1+a)*sig**2/2)),'quenched_cv':float(np.sqrt(2*np.exp(a*a*sig*sig)-1))}
    # measured hazard curve: log-log interpolation through s0 in {0.01,0.1,1,10} at L=9, power-law tails with the end slopes
    xs=np.log([0.001,0.01,0.1,1,10]); hs=np.log([F9['Tmem_mean']/runs[t]['9']['Tmem_mean'] for t in ('s0_0.001','s0_0.01','s0_0.1')]+[1.0,F9['Tmem_mean']/runs['s0_10']['9']['Tmem_mean']])
    rng=np.random.default_rng(0); xi=rng.normal(0,sig,2_000_000)-sig**2/2
    y=np.interp(xi,xs,hs); y=np.where(xi<xs[0],hs[0]+(xi-xs[0])*(hs[1]-hs[0])/(xs[1]-xs[0]),y); y=np.where(xi>xs[-1],hs[-1]+(xi-xs[-1])*(hs[-1]-hs[-2])/(xs[-1]-xs[-2]),y)
    out['hazard'][f'sigma{sig}']['annealed_measured']=float(1/np.exp(y).mean()); out['hazard'][f'sigma{sig}']['quenched_measured']=float(np.exp(-y).mean())
    xs4,hs4=xs[1:],hs[1:]; y4=np.interp(xi,xs4,hs4); y4=np.where(xi<xs4[0],hs4[0]+(xi-xs4[0])*(hs4[1]-hs4[0])/(xs4[1]-xs4[0]),y4); y4=np.where(xi>xs4[-1],hs4[-1]+(xi-xs4[-1])*(hs4[-1]-hs4[-2])/(xs4[-1]-xs4[-2]),y4)
    out['hazard'][f'sigma{sig}']['quenched_measured_without_1e-3']=float(np.exp(-y4).mean()); out['hazard'][f'sigma{sig}']['annealed_measured_without_1e-3']=float(1/np.exp(y4).mean())
N=2*9**3; Tf=F9['Tmem_mean']; out['crossover']={'s4_star':1/(N*np.exp(-18)*Tf)}
for tag,s4 in (('3e-4',3e-4),('1e-4',1e-4)):
    c=load(f'tmem_beta4.5_gapC_s4_{tag}.json')['9']; out['crossover'][f's4_{tag}']={'Tmem':c['Tmem_mean'],'sem':c['Tmem_sem'],'T_first':1/(N*s4*np.exp(-18)),'additive_pred':Tf+1/(N*s4*np.exp(-18))}
# aggregate CV over fixed-rate series; local hazard slopes; per-L OU gains with PDG-scaled slope; CV replica study; constants
cvs=[]
for f in glob.glob(R+'tmem_beta*_fast_*.json'):
    if 'fastcheck' in f or 'reviewA' in f: continue
    for L,r in load(os.path.basename(f)).items(): cvs.append(r['Tmem_sem']*np.sqrt(r['n'])/r['Tmem_mean'])
out['cv_fixed_rate']={'mean':float(np.mean(cvs)),'std':float(np.std(cvs)),'n_points':len(cvs)}
out['hazard']['local_slopes_L9']=dict(zip(['1e-3..1e-2','1e-2..1e-1','1e-1..1','1..10'],(np.diff(hs)/np.diff(xs)).tolist()))
Ls=[5,7,9,11,13]; gs=[];es=[]
for L in Ls:
    d=load(f'h2ou_beta4.5_L{L}_sigma2.0.json'); k=[k for k in d if float(k)==100.0][0]; fl=flat[str(L)]
    g=d[k]['Tmem_mean']/fl['Tmem_mean']; gs.append(g); es.append(g*np.hypot(d[k]['Tmem_sem']/d[k]['Tmem_mean'],fl['Tmem_sem']/fl['Tmem_mean']))
gs=np.array(gs);es=np.array(es); w=(gs/es)**2; X=np.vstack([np.ones(5),np.log(Ls)]).T; W=np.diag(w); cov=np.linalg.inv(X.T@W@X); b=cov@X.T@W@np.log(gs)
chi=((np.log(gs)-X@b)**2*w).sum()/3
out['ou']['per_L_sigma2_tau100']={'gains':gs.tolist(),'errs':es.tolist(),'slope':float(b[1]),'slope_err_pdg':float(np.sqrt(cov[1,1])*np.sqrt(max(1,chi))),'chi2dof':float(chi)}
# paired comparison tau_c=1 vs 100 (same flat denominator cancels): numerator errors only
d2=load('h2ou_beta4.5_L9_sigma2.0.json'); a1,a2=d2['1.0'],d2['100.0']
out['ou']['tau1_vs_tau100_sigma']=float((a2['Tmem_mean']-a1['Tmem_mean'])/np.hypot(a1['Tmem_sem'],a2['Tmem_sem']))
# CV replica study: sample CV at n=200 under the fully quenched power-law model and under an exponential (annealed) law
rng=np.random.default_rng(1); a=0.6; sig=2.0; n=200; reps=400000
xi=rng.normal(0,sig,(reps,n))-sig*sig/2; m=np.exp(-a*xi); T=rng.exponential(m); cv_q=T.std(1,ddof=0)/T.mean(1)
Te=rng.exponential(1.0,(reps,n)); cv_a=Te.std(1,ddof=0)/Te.mean(1)
cv_meas=out['ou']['sigma2.0_tau1e+10']['cv']
Te100=rng.exponential(1.0,(reps//4,100)); cv_a100=Te100.std(1,ddof=0)/Te100.mean(1)
out['ou']['cv_sampling_sd_exponential_n100']=float(cv_a100.std())
out['ou']['cv_replica_n200']={'quenched_median':float(np.median(cv_q)),'quenched_5_95':[float(np.percentile(cv_q,5)),float(np.percentile(cv_q,95))],'P_cv_le_meas_quenched':float((cv_q<=cv_meas).mean()),'annealed_median':float(np.median(cv_a)),'P_cv_ge_meas_annealed':float((cv_a>=cv_meas).mean()),'cv_measured':cv_meas}
# extension runs (2026-09-17): flat bath beyond L*, s0=1e-3 size series, quenched OU at other L and sigma
ext={}
fl={**flat}
for L in (17,19):
    try: fl.update(load(f'tmem_beta4.5_fast_flat_L{L}.json'))
    except FileNotFoundError: pass
Lx=sorted(int(x) for x in fl); ext['flat_Tmem']={L:(fl[str(L)]['Tmem_mean'],fl[str(L)]['Tmem_sem']) for L in Lx}
ext['flat_lambda']={f"{Lx[i]}->{Lx[i+1]}":(fl[str(Lx[i+1])]['Tmem_mean']/fl[str(Lx[i])]['Tmem_mean'],fl[str(Lx[i+1])]['Tmem_mean']/fl[str(Lx[i])]['Tmem_mean']*np.hypot(fl[str(Lx[i+1])]['Tmem_sem']/fl[str(Lx[i+1])]['Tmem_mean'],fl[str(Lx[i])]['Tmem_sem']/fl[str(Lx[i])]['Tmem_mean'])) for i in range(len(Lx)-1)}
big=[L for L in Lx if L>=13]
if len(big)>=2:   # plateau test: weighted mean of T(L>=13) and chi2 about it
    T=np.array([fl[str(L)]['Tmem_mean'] for L in big]); E=np.array([fl[str(L)]['Tmem_sem'] for L in big]); w=1/E**2; m=(w*T).sum()/w.sum()
    ext['flat_plateau_L>=13']={'Ls':big,'mean':float(m),'err':float(w.sum()**-0.5),'chi2dof':float((((T-m)/E)**2).sum()/(len(big)-1))}
s3={**load('tmem_beta4.5_fast_s0_0.001_allL.json'),**load('tmem_beta4.5_fast_s0_0.001.json')}
ext['s0_1e-3_gain']={int(L):(s3[L]['Tmem_mean']/flat[L]['Tmem_mean'],s3[L]['Tmem_mean']/flat[L]['Tmem_mean']*np.hypot(s3[L]['Tmem_sem']/s3[L]['Tmem_mean'],flat[L]['Tmem_sem']/flat[L]['Tmem_mean'])) for L in sorted(s3,key=int)}
f1=wls(s3,[5,7,9,11]); f2=wls(s3,[7,9,11,13]); ext['s0_1e-3_fit']={'w5_11':f1,'w7_13':f2,'window_syst':abs(f1['p']-f2['p'])/2}
ext['quenched']={}
for L,sig in ((7,2.0),(9,2.0),(11,2.0),(9,1.0)):
    q=load(f'h2ou_beta4.5_L{L}_sigma{sig}.json')['10000000000.0']; F=flat[str(L)]; g=q['Tmem_mean']/F['Tmem_mean']
    pred=float(np.exp(a*(1+a)*sig**2/2)); eg=g*np.hypot(q['Tmem_sem']/q['Tmem_mean'],F['Tmem_sem']/F['Tmem_mean'])
    ext['quenched'][f'L{L}_sigma{sig}']={'gain':g,'err':eg,'pred_powerlaw':pred,'pull':(g-pred)/eg,'cv':q['Tmem_sem']*np.sqrt(q['n'])/q['Tmem_mean'],'n':q['n']}
pl=np.array([v['pull'] for v in ext['quenched'].values()]); ext['quenched_chi2']={'chi2':float((pl**2).sum()),'npts':len(pl)}
s1={**load('tmem_beta4.5_fast_s0_0.01.json'),**load('tmem_beta4.5_fast_s0_0.01_L17.json'),**load('tmem_beta4.5_fast_s0_0.01_L19.json')}
ext['s0_1e-2_beyond_Lstar']={'gain':{L:(s1[str(L)]['Tmem_mean']/fl[str(L)]['Tmem_mean'],s1[str(L)]['Tmem_mean']/fl[str(L)]['Tmem_mean']*np.hypot(s1[str(L)]['Tmem_sem']/s1[str(L)]['Tmem_mean'],fl[str(L)]['Tmem_sem']/fl[str(L)]['Tmem_mean'])) for L in (13,17,19)},
 'ratio':{f'{a}->{b}':(s1[str(b)]['Tmem_mean']/s1[str(a)]['Tmem_mean'],s1[str(b)]['Tmem_mean']/s1[str(a)]['Tmem_mean']*np.hypot(s1[str(b)]['Tmem_sem']/s1[str(b)]['Tmem_mean'],s1[str(a)]['Tmem_sem']/s1[str(a)]['Tmem_mean'])) for a,b in ((13,17),(17,19))},
 'Tmem':{L:(s1[str(L)]['Tmem_mean'],s1[str(L)]['Tmem_sem'],s1[str(L)]['defects_at_failure_mean']) for L in (13,17,19)}}
out['extension']=ext; print('extension',json.dumps(ext,indent=1,default=float))
out['constants']={'N_L13':2*13**3,'qubits_per_step_L11_13':2*13**3-2*11**3,'qubits_per_step_d7_9':(2*9*9-1)-(2*7*7-1),'Lstar':{b:float(np.exp(0.78*b-0.87)) for b in (4.5,4.7,4.875,5.0,5.06)},'BH_exponent':{b:2.93*b-10.52 for b in (4.5,4.7,4.875,5.0,5.06)},'T13_over_T5':{k:float(load(f)['13']['Tmem_mean']/load(f)['5']['Tmem_mean']) for k,f in series.items() if '13' in load(f)}}
out['constants']['Tmem_over_Tec']={k:{L:r['Tmem_over_Tec'] for L,r in load(f).items()} for k,f in series.items()}
b40=load('tmem_beta4.0.json'); out['constants']['beta4.0_Tmem_over_Tec']={L:r['Tmem_mean']/(np.exp(16)/100) for L,r in b40.items()}
json.dump(out,open(R+'derived_numbers.json','w'),indent=1,default=float)
print('cv_fixed_rate',out['cv_fixed_rate']); print('local slopes',out['hazard']['local_slopes_L9']); print('per-L OU',out['ou']['per_L_sigma2_tau100']['slope'],out['ou']['per_L_sigma2_tau100']['slope_err_pdg']); print('tau1 vs tau100 sigma',out['ou']['tau1_vs_tau100_sigma']); print('cv replica',out['ou']['cv_replica_n200'])
for k,v in out['fits'].items(): print(f"{k:10s} p(5-11)={v['w5_11']['p']:.2f}±{v['w5_11']['err_pdg']:.2f} chi2/dof={v['w5_11']['chi2dof']:.2f}"+(f"  p(7-13)={v['w7_13']['p']:.2f}±{v['w7_13']['err_pdg']:.2f}  syst={v['window_syst']:.2f}" if v['w7_13'] else ""))
print('sensitivity',{k:round(v,2) for k,v in out['sensitivity'].items()}); print('hazard',out['hazard']); print('crossover',out['crossover'])
