from pathlib import Path
import json,warnings
import numpy as np,pandas as pd
from scipy import stats
from discover import bh
B=Path(__file__).parent;O=B/'results';rng=np.random.default_rng(20260914)
def run():
 z=np.load(O/'human_pseudobulk.npz');Y=z['counts'];G=z['genes'];m=pd.read_csv(O/'human_groups.csv');valid=~m.donor_id.str.startswith('pooled');L=Y.sum(axis=1);lookup=pd.Series(np.arange(len(G)),index=G)
 sig=pd.read_csv(O/'signature.csv');gi=lookup[sig.human].to_numpy();up=gi[sig.direction==1];dn=gi[sig.direction==-1]
 discovery=pd.read_csv(O/'mouse_genomewide.csv');cand=discovery[(discovery.q<.05)&(discovery.atlas_effect>1)&(discovery.control_effect>0)&(discovery.min_batch_loo_effect>0)].copy();ci=lookup[cand.human].to_numpy();cand.to_csv(O/'mouse_candidate_set.csv',index=False)
 universe=lookup[discovery.human].to_numpy();meanexp=np.log2(1+1e6*Y.sum(0)/L.sum());bins=pd.qcut(meanexp[universe],10,labels=False,duplicates='drop');bm=dict(zip(universe,bins));pools={b:universe[bins==b] for b in np.unique(bins)}
 rp=np.array([[rng.choice(pools[bm[g]]) for g in up] for _ in range(2000)]);rn=np.array([[rng.choice(pools[bm[g]]) for g in dn] for _ in range(2000)])
 rows=[];gene_rows=[];donor_rows=[];meanrows=[]
 # Keep every subtype's cells in the reference. Only target eligibility is thresholded.
 samples=list(m.loc[valid].groupby(['donor_id','sampleid']).groups.items())
 for ref in ['all','nonmidget']:
  for t in sorted(m.author_cell_type.unique()):
   if ref=='nonmidget' and t.startswith('MG_'):continue
   vals={};studies={}
   for (donor,sample),ii in samples:
    ii=np.array(list(ii));aa=ii[m.author_cell_type.iloc[ii].to_numpy()==t];bb=ii[m.author_cell_type.iloc[ii].to_numpy()!=t]
    if ref=='nonmidget':bb=bb[~m.author_cell_type.iloc[bb].str.startswith('MG_').to_numpy()]
    if not len(aa) or not len(bb) or m.n_cells.iloc[aa].sum()<20 or m.n_cells.iloc[bb].sum()<20:continue
    delta=np.log2(1+1e6*Y[aa].sum(0)/L[aa].sum())-np.log2(1+1e6*Y[bb].sum(0)/L[bb].sum());vals.setdefault(donor,[]).append(delta);studies[donor]=str(m.study.iloc[aa[0]])
   if len(vals)<8:continue
   donors=list(vals);D=np.stack([np.mean(vals[d],axis=0) for d in donors]);score=D[:,up].mean(1)-D[:,dn].mean(1);mu=D.mean(0);meanrows.append({'reference':ref,'subtype':t,'donors':len(D)});np.save(O/f'human_delta_{ref}_{t}.npy',D.astype(np.float32))
   for d,s in zip(donors,score):donor_rows.append({'reference':ref,'subtype':t,'donor':d,'study':studies[d],'score':s})
   boot=rng.choice(score,(5000,len(score)),replace=True).mean(1);lo,hi=np.quantile(boot,[.025,.975]);p=stats.wilcoxon(score,method='auto').pvalue if np.any(score) else 1.;null=mu[rp].mean(1)-mu[rn].mean(1);pnull=(1+(np.abs(null)>=abs(score.mean())).sum())/2001
   rows.append({'reference':ref,'subtype':t,'donors':len(D),'score':score.mean(),'ci_low':lo,'ci_high':hi,'positive_fraction':(score>0).mean(),'p':p,'matched_gene_p':pnull,'matched_gene_null_mean':null.mean(),'loo_min':np.min((score.sum()-score)/(len(score)-1)),'loo_max':np.max((score.sum()-score)/(len(score)-1))})
   with warnings.catch_warnings():
    warnings.simplefilter('ignore');ps=stats.wilcoxon(D[:,ci],axis=0,method='approx').pvalue
   for k,(ix,g) in enumerate(zip(ci,cand.human)):
    v=D[:,ix];gene_rows.append({'reference':ref,'subtype':t,'gene':g,'donors':len(D),'effect':v.mean(),'positive_fraction':(v>0).mean(),'p':ps[k],'loo_min':np.min((v.sum()-v)/(len(v)-1))})
 out=pd.DataFrame(rows);out['q']=bh(out.p);out['matched_gene_q']=bh(out.matched_gene_p);out.to_csv(O/'human_signature_results.csv',index=False)
 gg=pd.DataFrame(gene_rows);gg['q']=bh(gg.p);gg.to_csv(O/'human_candidate_results.csv',index=False);pd.DataFrame(donor_rows).to_csv(O/'human_signature_donors.csv',index=False)
 print('mouse candidates',len(cand),'human gene tests',len(gg),flush=True);print(out.to_string(index=False),flush=True)
if __name__=='__main__':run()
