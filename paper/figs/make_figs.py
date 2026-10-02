"""Regenerate all manuscript figures as vector PDFs with clean legends. Column width 3.4 in, full width 7.0 in, 8 pt."""
import json, glob, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':8,'axes.labelsize':8,'legend.fontsize':6.5,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':0.6,'lines.markersize':3.5,'pdf.fonttype':42})
R='results/'; out='paper/figs/'
def load(f): return json.load(open(R+f))
def series(d): Ls=np.array(sorted(int(k) for k in d)); return Ls,np.array([d[str(L)]['Tmem_mean'] for L in Ls]),np.array([d[str(L)]['Tmem_sem'] for L in Ls])
BH=lambda b:2.93*b-10.52; Lstar=lambda b:np.exp(0.78*b-0.87)
# ---- Fig 1 baseline
fig,ax=plt.subplots(figsize=(3.4,2.6))
for f,b,c in [("tmem_beta4.5_fast_flat.json",4.5,'C0'),("tmem_beta5.0_fast_flat.json",5.0,'C3')]:
    d=load(f)
    if b==4.5: d={**d,**load('tmem_beta4.5_fast_flat_L17.json'),**load('tmem_beta4.5_fast_flat_L19.json')}   # sizes beyond L*
    Ls,T,E=series(d); ax.errorbar(Ls,T,E,fmt='o',color=c,capsize=2,label=fr'$\beta={b}$')
    Lf=np.linspace(4.5,14); ax.plot(Lf,T[0]*(Lf/Ls[0])**BH(b),'--',color=c,lw=1,label=fr'$L^{{2.93\beta-10.52}}=L^{{{BH(b):.2f}}}$'); ax.axvline(Lstar(b),ls=':',color=c,lw=0.8)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xticks([5,7,9,11,13,17,19]); ax.set_xticklabels(['5','7','9','11','13','17','19']); ax.minorticks_off()
ax.set_xlabel(r'linear size $L$'); ax.set_ylabel(r'$T_{\rm mem}$ (units of $1/\gamma_0$)'); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(out+'fig1_baseline.pdf'); plt.close()
# ---- Fig 2 channels
runs={os.path.basename(f)[len('tmem_beta4.5_fast_'):-5]:load(os.path.basename(f)) for f in glob.glob(R+'tmem_beta4.5_fast_*.json')}
for tag in ('s4_3e-4','s4_1e-4'): runs[tag]=load(f'tmem_beta4.5_gapC_{tag}.json')
lab={'flat':'flat','s4_0.01':r'$s_4=0.01$','s4_0.1':r'$s_4=0.1$','s4_10':r'$s_4=10$','s2_0.1':r'$s_2=0.1$','s2_10':r'$s_2=10$','s0_0.1':r'$s_0=0.1$','s0_0.01':r'$s_0=0.01$','s0_10':r'$s_0=10$'}
sty={'flat':('k','o'),'s4_0.01':('C0','v'),'s4_0.1':('C0','^'),'s4_10':('C0','s'),'s2_0.1':('C1','^'),'s2_10':('C1','s'),'s0_0.1':('C2','^'),'s0_0.01':('C2','v'),'s0_10':('C2','s')}
runs['s0_0.001']={**runs.pop('s0_0.001_allL'),**runs['s0_0.001']}; runs.pop('flat_L17',None); runs.pop('flat_L19',None)   # merge the s0=1e-3 size series; flat L>13 is shown in Fig. 1 only
lab['s0_0.001']=r'$s_0=0.001$'; sty['s0_0.001']=('C2','D')
runs['s0_0.01']={**runs['s0_0.01'],**runs.pop('s0_0.01_L17'),**runs.pop('s0_0.01_L19')}   # filtered bath beyond L*
fig,ax=plt.subplots(1,2,figsize=(7.0,2.7))
for tag in ['flat','s4_0.01','s4_0.1','s4_10','s2_0.1','s2_10','s0_0.1','s0_0.01','s0_0.001','s0_10']:
    Ls,T,E=series(runs[tag]); c,m=sty[tag]; ax[0].errorbar(Ls,T,E,fmt=m+'-',color=c,lw=0.8,capsize=1.5,label=lab[tag])
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xticks([5,7,9,11,13,17,19]); ax[0].set_xticklabels(['5','7','9','11','13','17','19']); ax[0].minorticks_off(); ax[0].set_xlabel(r'$L$'); ax[0].set_ylabel(r'$T_{\rm mem}$'); ax[0].set_ylim(2e5,4e10); ax[0].legend(ncol=4,frameon=False,loc='lower center',fontsize=6.5,columnspacing=0.8,handletextpad=0.3)
flat=runs['flat']['9']['Tmem_mean']
for ch,idx,c,name in ((r'$s_4$: pair creation, $|\Delta E|=4$',0,'C0','s4'),(r'$s_2$: barrier, $|\Delta E|=2$',1,'C1','s2'),(r'$s_0$: diffusion, $\Delta E=0$',2,'C2','s0')):
    pts=[]
    for tag,d in runs.items():
        sp=d['9']['spec']
        if all(sp[j]==1.0 for j in range(3) if j!=idx): pts.append((sp[idx],d['9']['Tmem_mean'],d['9']['Tmem_sem']))
    pts=np.array(sorted(pts)); ax[1].errorbar(pts[:,0],pts[:,1],pts[:,2],fmt='o-',color=c,lw=0.8,capsize=1.5,label=ch)
x=np.array([1e-3,10]); ax[1].plot(x,flat/x,'k:',lw=0.8,label=r'$\propto 1/s$'); ax[1].axhline(flat,color='gray',ls='--',lw=0.8,label='flat bath')
xs4=np.logspace(-4.3,1,50); ax[1].plot(xs4,flat+1/(1458*xs4*np.exp(-18)),'-',color='C0',lw=0.7,alpha=0.6,label=r'$T_{\rm flat}+T_{\rm first}(s_4)$'); ax[1].axvline(3.8e-4,color='C0',ls=':',lw=0.7)
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlim(5e-5,20); ax[1].set_xlabel(r'filter factor $s=J(|\Delta E|)/J_{\rm flat}$'); ax[1].set_ylabel(r'$T_{\rm mem}(L=9)$'); ax[1].legend(frameon=False)
for a,l in zip(ax,'ab'): a.text(-0.18,1.02,f'({l})',transform=a.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(out+'fig2_channels.pdf'); plt.close()
# ---- Fig 3 OU
fig,ax=plt.subplots(figsize=(3.4,2.6)); a=0.6
for f,c in [('h2ou_beta4.5_L9_sigma1.0.json','C0'),('h2ou_beta4.5_L9_sigma2.0.json','C1')]:
    d=load(f); taus=np.array(sorted(float(k) for k in d)); key=lambda t:[k for k in d if float(k)==t][0]; sig=d[key(taus[0])]['sigma']
    T=np.array([d[key(t)]['Tmem_mean'] for t in taus]); E=np.array([d[key(t)]['Tmem_sem'] for t in taus])
    Ef=load('tmem_beta4.5_fast_flat.json')['9']['Tmem_sem']; Eg=(T/flat)*np.hypot(E/T,Ef/flat)
    ax.errorbar(taus,T/flat,Eg,fmt='o-',color=c,lw=0.8,capsize=2,label=fr'$\sigma={sig:g}$'); ax.axhline(np.exp(a*(1-a)*sig**2/2),ls='--',color=c,lw=0.8,label=fr'annealed $1/\langle s_0^{{0.6}}\rangle$, $\sigma={sig:g}$')
    if sig==2.0: ax.axhline(np.exp(a*(1+a)*sig**2/2),ls=':',color=c,lw=0.8,label=r'quenched $\langle s_0^{-0.6}\rangle$, $\sigma=2$'); ax.axvline(1.2e8,color='gray',ls=':',lw=0.7)
ax.axhline(1,color='k',lw=0.6); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\tau_c$ (hop times); dotted vertical line: $T_{\rm mem}^{\rm flat}$'); ax.set_ylabel(r'$T_{\rm mem}/T_{\rm mem}^{\rm flat}$  ($L=9$, $\beta=4.5$)'); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(out+'fig3_ou.pdf'); plt.close()
# ---- Fig 4 pseudomode
P=load('pseudomode_gap.json'); fig,ax=plt.subplots(1,2,figsize=(7.0,2.7))
for beta,ls in ((2.0,'--'),(4.5,'-')):
    for g,c in zip(sorted({r['g'] for r in P['dynamic_resonant']}),['C0','C1','C2','C3']):
        rs=sorted([r for r in P['dynamic_resonant'] if r['beta']==beta and r['g']==g],key=lambda r:r['kappa'])
        ax[0].plot([r['kappa'] for r in rs],[r['Pexc']/r['P_boltz'] for r in rs],'o'+ls,color=c,lw=0.8,label=fr'$\beta={beta:g}$, $g={g:g}$')
    for g,c in zip(sorted({r['g'] for r in P['meanforce_static']}),['C0','C1','C2','C3']):
        rs=sorted([r for r in P['meanforce_static'] if r['beta']==beta and r['g']==g],key=lambda r:r['w0'])
        ax[1].plot([r['w0'] for r in rs],[r['ratio'] for r in rs],'o'+ls,color=c,lw=0.8,label=fr'$\beta={beta:g}$, $g={g:g}$')
ax[0].axhline(1,color='k',lw=0.6); ax[0].set_xscale('log'); ax[0].set_ylim(0.98,1.02); ax[0].set_xlabel(r'$\kappa=1/\tau_c$  ($\omega_0=\Delta$)'); ax[0].set_ylabel(r'$P_{\rm exc}/P_{\rm Boltzmann}$'); ax[0].legend(ncol=2,frameon=False)
ax[1].axhline(1,color='k',lw=0.6); ax[1].axvline(4,color='gray',ls=':',lw=0.8); ax[1].set_yscale('log'); ax[1].set_xlabel(r'bath peak $\omega_0$  ($\Delta=4$)'); ax[1].set_ylabel(r'$\mathrm{Tr}[\Pi_e\,\rho_{\rm Gibbs}]\,/\,P_{\rm Boltzmann}$'); ax[1].legend(ncol=2,frameon=False)
for a,l in zip(ax,'ab'): a.text(-0.18,1.02,f'({l})',transform=a.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(out+'fig4_pseudomode.pdf'); plt.close()
# ---- Fig 5 Lamb lever
files={4.5:'tmem_beta4.5_fast_flat.json',4.7:'tmem_beta4.7_fast_flat_lamb.json',4.875:'tmem_beta4.875_fast_flat_lamb.json',5.06:'tmem_beta5.06_fast_flat_lamb.json'}
fig,ax=plt.subplots(1,2,figsize=(7.0,2.7))
for (b,f),c in zip(sorted(files.items()),['C0','C1','C2','C3']):
    Ls,T,E=series(load(f)); delta=(b/4.5-1)*4; ax[0].errorbar(Ls,T,E,fmt='o-',color=c,lw=0.8,capsize=1.5,label=fr"$\beta'={b}$ ($\delta={delta:.2f}$)"); ax[0].axvline(Lstar(b),ls=':',color=c,lw=0.7)
    lam=T[1:]/T[:-1]; err=lam*np.hypot(E[1:]/T[1:],E[:-1]/T[:-1]); ax[1].errorbar(Ls[:-1]+1,lam,err,fmt='o-',color=c,lw=0.8,capsize=1.5,label=fr"$\beta'={b}$")
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xticks([5,7,9,11,13]); ax[0].set_xticklabels(['5','7','9','11','13']); ax[0].minorticks_off(); ax[0].set_xlabel(r'$L$'); ax[0].set_ylabel(r'$T_{\rm mem}$'); ax[0].legend(frameon=False)
ax[1].axhline(2.14,color='k',ls='--',lw=0.8,label=r'$\Lambda_{\rm active}=2.14$'); ax[1].set_xlabel(r'$L$ (midpoint of $L\to L+2$)'); ax[1].set_ylabel(r'$\Lambda_{\rm passive}=T_{\rm mem}(L+2)/T_{\rm mem}(L)$'); ax[1].legend(frameon=False)
for a,l in zip(ax,'ab'): a.text(-0.18,1.02,f'({l})',transform=a.transAxes,fontweight='bold')
fig.tight_layout(); fig.savefig(out+'fig5_lamb.pdf'); plt.close(); print('figs ok')
