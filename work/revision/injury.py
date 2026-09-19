from pathlib import Path
import numpy as np,pandas as pd
B=Path(__file__).parent;O=B/'results'
def main():
 sig=pd.read_csv(O/'signature.csv');cand=pd.read_csv(O/'mouse_candidate_set.csv');wanted=list(dict.fromkeys(list(sig.mouse)+list(cand.mouse)));profiles={};rows=[]
 for cond in ['control','14d']:
  m=pd.read_csv(O/f'{cond}_pseudobulk_groups.csv');z=np.load(O/f'{cond}_pseudobulk.npz');Y=z['counts'];G=z['genes'];lookup=pd.Series(np.arange(len(G)),index=G);ind=lookup[wanted].values
  # Only Vglut2 background. Match CD45-depleted and CD90-only purification strata.
  for t in range(1,46):
   for stratum in ['CD45','CD90only']:
    sel=m.library.str.contains('Vg') & (m.type==t);sel &= m.library.str.contains('CD45') if stratum=='CD45' else ~m.library.str.contains('CD45')
    if m.loc[sel,'cells'].sum()<5:continue
    c=Y[sel].sum(0);profiles[(cond,t,stratum)]=(np.log2(1+1e6*c[ind]/c.sum()),int(m.loc[sel,'cells'].sum()))
 for t in range(1,46):
  diffs=[];ns=[]
  for s in ['CD45','CD90only']:
   if ('control',t,s) in profiles and ('14d',t,s) in profiles:
    a,na=profiles[('control',t,s)];b,nb=profiles[('14d',t,s)];diffs.append(b-a);ns.append((na,nb))
  if not diffs:continue
  delta=np.mean(diffs,axis=0)
  for k,g in enumerate(wanted):rows.append({'type':t,'mouse':g,'strata':len(diffs),'control_cells':sum(a for a,b in ns),'injured_cells':sum(b for a,b in ns),'injury_difference':delta[k]})
 out=pd.DataFrame(rows);out.to_csv(O/'mouse_within_type_injury.csv',index=False)
 print('types with two matched strata',out[out.strata==2].type.nunique())
 for g in ['Ptprt','Man1a','Igfbp5','Rasgrp1','Gldn','Parm1','Lpar1']:
  a=out[(out.mouse==g)&(out.strata==2)];print(g,'types',len(a),'median',a.injury_difference.median(),'positive',int((a.injury_difference>0).sum()))
if __name__=='__main__':main()
