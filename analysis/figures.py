from pathlib import Path
import os
os.environ['MPLCONFIGDIR']=str(Path(__file__).parent/'mplcache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
from scipy.stats import pearsonr
B=Path(__file__).parent;O=B/'results';F=B/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':350})
def save(fig,name):
 fig.savefig(F/(name+'.png'),bbox_inches='tight',pad_inches=.035);fig.savefig(F/(name+'.pdf'),bbox_inches='tight',pad_inches=.035);plt.close(fig)
m=pd.read_csv(O/'mouse_genomewide.csv');s=pd.read_csv(O/'signature.csv');cv=pd.read_csv(O/'mouse_crossvalidation.csv')
fig,axs=plt.subplots(2,1,figsize=(3.48,4.1),gridspec_kw={'height_ratios':[1.1,1]});ax=axs[0];ax.scatter(m.atlas_effect,m.control_effect,s=2,c='#b9c4cc',alpha=.35,rasterized=True)
for direction,col in [(1,'#176a86'),(-1,'#a24d27')]:
 z=s[s.direction==direction];ax.scatter(z.atlas_effect,z.control_effect,s=9,c=col,label='Resilient-enriched' if direction==1 else 'Other-type-enriched')
ax.plot([-3,6],[-3,6],color='#777777',ls=':',lw=.7);ax.set(xlabel='Atlas association (log2(CPM + 1) difference)',ylabel='Control association',xlim=(-3,6.2),ylim=(-3,6.2));ax.set_title('A  Baseline associations across collections',loc='left',fontweight='bold',fontsize=8);ax.text(.04,.92,f'r = {pearsonr(m.atlas_effect,m.control_effect).statistic:.3f}',transform=ax.transAxes);ax.legend(frameon=False,fontsize=6,loc='lower right')
ax=axs[1]
for field,col,label in [('heldout_score','#176a86','50-gene signature'),('known_marker_score','#a24d27','Six established markers')]:
 yy=cv.published_resilient.to_numpy();vv=cv[field].to_numpy();rank=np.argsort(-vv);tpr=np.r_[0,np.cumsum(yy[rank])/yy.sum()];fpr=np.r_[0,np.cumsum(~yy[rank])/(~yy).sum()];auc=np.trapezoid(tpr,fpr);ax.plot(fpr,tpr,color=col,lw=1.2,label=f'{label}: AUC {auc:.3f}')
ax.plot([0,1],[0,1],color='#777777',ls=':',lw=.7);ax.set(xlabel='False-positive rate',ylabel='True-positive rate',xlim=(0,1),ylim=(0,1.04));ax.set_title('B  Leave one mouse subtype out',loc='left',fontweight='bold',fontsize=8);ax.legend(frameon=False,fontsize=6,loc='lower right');fig.tight_layout(pad=.55,h_pad=1.3);save(fig,'figure1_mouse')
h=pd.read_csv(O/'human_signature_results.csv');order=['MG_ON','MG_OFF','PG_ON','PG_OFF','ipRGC','HRGC7','HRGC8','HRGC9','HRGC10','HRGC11','HRGC12'];labels={'MG_ON':'ON midget','MG_OFF':'OFF midget','PG_ON':'ON parasol','PG_OFF':'OFF parasol'}
fig,ax=plt.subplots(figsize=(3.48,3.2))
for ref,col,offset in [('all','#176a86',-.13),('nonmidget','#a24d27',.13)]:
 for i,t in enumerate(order):
  z=h[(h.reference==ref)&(h.subtype==t)]
  if len(z):
   v=z.iloc[0];ax.errorbar(v.score,i+offset,xerr=[[v.score-v.ci_low],[v.ci_high-v.score]],fmt='o',color=col,ms=3,capsize=2,lw=.8,markerfacecolor=col if v.q<.05 else 'white')
 ax.plot([],[],color=col,marker='o',ms=3,label='All other RGCs' if ref=='all' else 'Other nonmidget RGCs')
ax.axvline(0,color='#777777',ls=':',lw=.8);ax.set_yticks(range(len(order)),[labels.get(t,t) for t in order]);ax.invert_yaxis();ax.set_xlabel('Mean donor signature contrast (95% bootstrap interval)');ax.legend(frameon=False,fontsize=6,loc='lower left');ax.set_ylim(11.1,-.6);ax.set_xlim(-1.55,.75);fig.tight_layout(pad=.5);save(fig,'figure2_human')
genes=['PTPRT','MAN1A1','IGFBP5','RASGRP1','GLDN','PARM1','LPAR1'];mouse=['Ptprt','Man1a','Igfbp5','Rasgrp1','Gldn','Parm1','Lpar1'];types=['PG_ON','PG_OFF','ipRGC','HRGC7','HRGC8','HRGC9'];hh=pd.read_csv(O/'human_candidate_results.csv');a=hh[hh.reference=='nonmidget'].pivot(index='gene',columns='subtype',values='effect').reindex(index=genes,columns=types);q=hh[hh.reference=='nonmidget'].pivot(index='gene',columns='subtype',values='q').reindex(index=genes,columns=types)
fig,axs=plt.subplots(2,1,figsize=(3.48,4.0),gridspec_kw={'height_ratios':[1.2,1]});ax=axs[0];im=ax.imshow(a.values,cmap='RdBu_r',vmin=-3.5,vmax=3.5,aspect='auto');ax.set_yticks(range(len(genes)),genes);ax.set_xticks(range(len(types)),['PG ON','PG OFF','ipRGC','HRGC7','HRGC8','HRGC9'],rotation=30,ha='right');ax.set_title('A  Human candidate expression contrasts',loc='left',fontweight='bold',fontsize=8)
for i in range(len(genes)):
 for j in range(len(types)):
  if q.iloc[i,j]<.05:ax.text(j,i,'•',ha='center',va='center',fontsize=8,color='black' if abs(a.iloc[i,j])<2 else 'white')
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.03);cb.ax.tick_params(labelsize=6);cb.set_label('log2(CPM + 1) difference',fontsize=6)
inj=pd.read_csv(O/'mouse_within_type_injury.csv');ax=axs[1]
for i,(g,mg) in enumerate(zip(genes,mouse)):
 vals=inj[(inj.mouse==mg)&(inj.strata==2)].injury_difference.to_numpy();jitter=np.linspace(-.17,.17,len(vals));ax.scatter(i+jitter,vals,s=4,alpha=.6,color='#687a86');ax.plot([i-.2,i+.2],[np.median(vals)]*2,color='#a24d27',lw=1.7)
ax.axhline(0,color='#777777',ls=':',lw=.8);ax.set_xticks(range(len(genes)),genes,rotation=45,ha='right');ax.set_ylabel('14-day minus control');ax.set_title('B  Mouse injury changes within 23 subtypes',loc='left',fontweight='bold',fontsize=8);fig.tight_layout(pad=.55,h_pad=1.5);save(fig,'figure3_candidates')
print('figures complete')
