from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import stats
from discover import bh
B=Path(__file__).parent;O=B/'results';rng=np.random.default_rng(20260918)
def main():
 z=np.load(O/'human_pseudobulk.npz');Y=z['counts'];G=z['genes'];L=Y.sum(1);m=pd.read_csv(O/'human_groups.csv');sig=pd.read_csv(O/'signature.csv');lookup=pd.Series(np.arange(len(G)),index=G);gi=lookup[sig.human].values;direction=sig.direction.values
 known={'IGF1','OPN4','SPP1','TIMP2','NDNF','PRPH'};retain=~sig.human.isin(known).to_numpy();rows=[]
 # Additional analyses use the exact primary donor contrasts already computed.
 donors=pd.read_csv(O/'human_signature_donors.csv')
 for (ref,t),v in donors.groupby(['reference','subtype'],sort=False):
  D=np.load(O/f'human_delta_{ref}_{t}.npy');s=D[:,gi[retain & (direction==1)]].mean(1)-D[:,gi[retain & (direction==-1)]].mean(1)
  row={'reference':ref,'subtype':t,'donors':len(s),'score_without_known_markers':s.mean(),'p_without_known_markers':stats.wilcoxon(s).pvalue}
  los=[]
  for st in v.study.unique():
   a=v.loc[v.study!=st,'score'].to_numpy()
   if len(a)>=3:los.append(a.mean())
  row.update({'studies':v.study.nunique(),'leave_study_out_min':min(los) if los else np.nan,'leave_study_out_max':max(los) if los else np.nan});rows.append(row)
 ss=pd.DataFrame(rows);ss['q_without_known_markers']=bh(ss.p_without_known_markers);ss.to_csv(O/'human_sensitivity.csv',index=False)
 # Vary the cell threshold without changing the learned gene panel.
 threshold_rows=[];samples=list(m.loc[~m.donor_id.str.startswith('pooled')].groupby(['donor_id','sampleid']).groups.items());YY=Y[:,gi]
 for threshold in [10,30]:
  for ref in ['all','nonmidget']:
   for t in sorted(m.author_cell_type.unique()):
    if ref=='nonmidget' and t.startswith('MG_'):continue
    vals={}
    for (donor,sample),ii in samples:
     ii=np.array(list(ii));aa=ii[m.author_cell_type.iloc[ii].to_numpy()==t];bb=ii[m.author_cell_type.iloc[ii].to_numpy()!=t]
     if ref=='nonmidget':bb=bb[~m.author_cell_type.iloc[bb].str.startswith('MG_').to_numpy()]
     if not len(aa) or not len(bb) or m.n_cells.iloc[aa].sum()<threshold or m.n_cells.iloc[bb].sum()<threshold:continue
     d=np.log2(1+1e6*YY[aa].sum(0)/L[aa].sum())-np.log2(1+1e6*YY[bb].sum(0)/L[bb].sum());vals.setdefault(donor,[]).append(d[direction==1].mean()-d[direction==-1].mean())
    if len(vals)>=8:
     s=np.array([np.mean(a) for a in vals.values()]);threshold_rows.append({'threshold':threshold,'reference':ref,'subtype':t,'donors':len(s),'score':s.mean(),'p':stats.wilcoxon(s).pvalue})
 tt=pd.DataFrame(threshold_rows);tt['q']=bh(tt.p);tt.to_csv(O/'human_cell_thresholds.csv',index=False)
 pairs=[]
 for a,b in [('MG_ON','MG_OFF'),('PG_ON','PG_OFF'),('ipRGC','PG_ON')]:
  vals={}
  for (donor,sample),ii in samples:
   ii=np.array(list(ii));ia=ii[m.author_cell_type.iloc[ii].to_numpy()==a];ib=ii[m.author_cell_type.iloc[ii].to_numpy()==b]
   if not len(ia) or not len(ib) or min(m.n_cells.iloc[ia].sum(),m.n_cells.iloc[ib].sum())<20:continue
   d=np.log2(1+1e6*YY[ia].sum(0)/L[ia].sum())-np.log2(1+1e6*YY[ib].sum(0)/L[ib].sum());vals.setdefault(donor,[]).append(d[direction==1].mean()-d[direction==-1].mean())
  s=np.array([np.mean(x) for x in vals.values()]);boot=rng.choice(s,(5000,len(s)),replace=True).mean(1);lo,hi=np.quantile(boot,[.025,.975]);pairs.append({'a':a,'b':b,'donors':len(s),'effect':s.mean(),'ci_low':lo,'ci_high':hi,'p':stats.wilcoxon(s).pvalue})
 pp=pd.DataFrame(pairs);pp['q']=bh(pp.p);pp.to_csv(O/'human_direct_pairs.csv',index=False)
 # Gene panel size is fixed before human analysis. A family-held-out mouse check addresses identity leakage.
 z=np.load(O/'mouse_type_expression.npz');E=z['atlas'];y=z['resilient'];types=z['types'];cv=[]
 families={'ipRGC':[22,31,33,40,43],'alpha':[41,42,43,45],'ON_DS':[10]}
 for family,held in families.items():
  test=np.isin(types,held);train=~test;X=E[train];lab=y[train];r=(X[lab].mean(0)-X[~lab].mean(0))/(X.std(0,ddof=1)+1);pos=np.argsort(r)[-25:];neg=np.argsort(r)[:25];V=(E-X.mean(0))/(X.std(0,ddof=1)+1);s=V[:,pos].mean(1)-V[:,neg].mean(1)
  for i in np.flatnonzero(test):cv.append({'family':family,'type':types[i],'resilient':bool(y[i]),'score':s[i],'percentile_vs_training':100*(s[train]<s[i]).mean()})
 pd.DataFrame(cv).to_csv(O/'mouse_family_holdout.csv',index=False)
 print(pp.to_string(index=False));print(ss.to_string(index=False));print(pd.DataFrame(cv).to_string(index=False))
if __name__=='__main__':main()
