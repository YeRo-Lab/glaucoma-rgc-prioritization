from pathlib import Path
import sys,json,gzip,csv,hashlib
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
from scipy.stats import wilcoxon
BASE=Path(__file__).resolve().parent
OUT=BASE/'results';OUT.mkdir(exist_ok=True)
CAND=['SPP1','UCN','CRHBP','TIMP2','MMP12','MMP9','NDNF','PRPH','CRHR1']
MARK=['OPN4','NPPB','NEUROD2','SATB2','TAC1','RBPMS','SLC17A6','TBR1','TPBG','FABP4','CHRNA2','BNC2']
GENES=CAND+MARK

def vec(g):
 if isinstance(g,h5py.Group):
  cats=vec(g['categories']);codes=g['codes'][:];return np.array([cats[x] if x>=0 else 'NA' for x in codes])
 return g.asstr()[:] if g.dtype.kind in 'OSU' else g[:]

def human():
 f=h5py.File(BASE/'human_all_rgc.h5ad','r');obs=pd.DataFrame({k:vec(f['obs'][k]) for k in ['author_cell_type','donor_id','sampleid','tissue','disease','study_name']})
 raw=f['raw/X'];names=vec(f['raw/var/feature_name']);idx=[np.flatnonzero(names==g)[0] for g in GENES];n=len(obs);gnum=len(names)
 y=np.zeros((n,len(GENES)));lib=np.zeros(n);detect=np.zeros_like(y,dtype=bool)
 p=raw['indptr'][:]
 for a in range(0,n,5000):
  b=min(a+5000,n);l,r=p[a],p[b];x=csr_matrix((raw['data'][l:r],raw['indices'][l:r],p[a:b+1]-l),shape=(b-a,gnum));assert np.all(x.data>=0) and np.all(x.data==np.floor(x.data));lib[a:b]=np.asarray(x.sum(axis=1)).ravel();y[a:b]=x[:,idx].toarray()
  if a%50000==0:print('human rows',a,flush=True)
 assert (lib>0).all();detect=y>0
 np.savez_compressed(OUT/'human_candidate_counts.npz',counts=y,library_size=lib,genes=np.array(GENES))
 obs.to_csv(OUT/'human_metadata.csv',index=False)
 summary={'cells':n,'genes':gnum,'types':obs.author_cell_type.value_counts().to_dict(),'donor_records':obs.donor_id.nunique(),'pooled_cells':int(obs.donor_id.str.startswith('pooled').sum()),'unpooled_donors':obs.loc[~obs.donor_id.str.startswith('pooled'),'donor_id'].nunique(),'disease':obs.disease.value_counts().to_dict(),'samples':obs.sampleid.nunique()}
 (OUT/'human_summary.json').write_text(json.dumps(summary,indent=2))
 rows=[]
 for t,ii in obs.groupby('author_cell_type',observed=True).indices.items():
  for k,g in enumerate(GENES):rows.append({'subtype':t,'gene':g,'n_cells':len(ii),'detection_pct':100*detect[ii,k].mean(),'pooled_cpm':1e6*y[ii,k].sum()/lib[ii].sum()})
 pd.DataFrame(rows).to_csv(OUT/'human_detection.csv',index=False)
 # Preserve sample strata (region, library preparation) before equal weighting donors.
 paired=[]
 valid=~obs.donor_id.str.startswith('pooled')
 for (donor,sample),ii in obs.loc[valid].groupby(['donor_id','sampleid'],observed=True).groups.items():
  ii=np.array(list(ii));tt=obs.author_cell_type.iloc[ii].values
  for t in sorted(set(tt)):
   aa=ii[tt==t];bb=ii[tt!=t]
   if len(aa)<10 or len(bb)<10:continue
   delta=np.log2(1+1e6*y[aa].sum(axis=0)/lib[aa].sum())-np.log2(1+1e6*y[bb].sum(axis=0)/lib[bb].sum())
   for k,g in enumerate(CAND):paired.append({'donor':donor,'sample':sample,'tissue':obs.tissue.iloc[ii[0]],'subtype':t,'gene':g,'n_target':len(aa),'n_other':len(bb),'delta':delta[k]})
 pp=pd.DataFrame(paired);pp.to_csv(OUT/'human_sample_contrasts.csv',index=False)
 for threshold in [10,20,30]:
  z=pp[(pp.n_target>=threshold)&(pp.n_other>=threshold)];dd=z.groupby(['subtype','gene','donor']).delta.mean().reset_index();rows=[];rng=np.random.default_rng(20260913)
  for (t,g),v in dd.groupby(['subtype','gene']):
   d=v.delta.to_numpy();n_d=len(d)
   if n_d<8:continue
   ci=np.quantile(rng.choice(d,(2000,n_d),replace=True).mean(axis=1),[.025,.975]);pval=wilcoxon(d,alternative='two-sided',method='auto').pvalue if np.any(d!=0) else 1.
   rows.append({'subtype':t,'gene':g,'donors':n_d,'mean_log2_cpm1_difference':d.mean(),'ci_low':ci[0],'ci_high':ci[1],'positive_fraction':(d>0).mean(),'p':pval,'min_leave_one_donor_out_mean':min((d.sum()-d)/(n_d-1))})
  rr=pd.DataFrame(rows);order=np.argsort(rr.p.values);q=np.minimum.accumulate((rr.p.values[order]*len(rr)/np.arange(1,len(rr)+1))[::-1])[::-1];rr['q']=np.nan;rr.loc[order,'q']=np.minimum(q,1);rr.to_csv(OUT/f'human_contrasts_min{threshold}.csv',index=False)
 print(summary,flush=True)

def mouse():
 allrows=[];details=[]
 for cond,file in [('control','mouse_control.csv.gz'),('14d','mouse_14d.csv.gz')]:
  with gzip.open(BASE/file,'rt') as f:
   cells=next(csv.reader([f.readline()]))[1:];samples=np.array([x.split('_')[0] for x in cells]);uniq=sorted(set(samples));si=[np.flatnonzero(samples==s) for s in uniq];lib=np.zeros(len(cells));nf=np.zeros(len(cells));mt=np.zeros(len(cells));selected={};ng=0
   for line in f:
    name,values=line.split(',',1);name=name.strip('"');a=np.fromstring(values,sep=',');assert len(a)==len(cells);lib+=a;nf+=(a>0);ng+=1
    if name.lower().startswith('mt-'):mt+=a
    if name.upper() in GENES:selected[name.upper()]=a
  qc=(lib>=1000)&(nf>=500)&(mt/np.maximum(lib,1)<=.2)&((selected['RBPMS']>0)|(selected['SLC17A6']>0))
  np.savez_compressed(OUT/f'mouse_{cond}_candidate_counts.npz',cells=np.array(cells),counts=np.stack([selected.get(g,np.zeros(len(cells))) for g in GENES],axis=1),genes=np.array(GENES),library_size=lib,n_genes=nf,mt_counts=mt,qc=qc)
  for s,ii in zip(uniq,si):
   n_raw=len(ii);ii=ii[qc[ii]]
   details.append({'condition':cond,'library':s,'cells':len(ii),'raw_cells':n_raw,'total_counts':int(lib[ii].sum()),'genes':ng})
   for g,a in selected.items():allrows.append({'condition':cond,'library':s,'gene':g,'cells':len(ii),'cpm':1e6*a[ii].sum()/lib[ii].sum(),'detection_pct':100*(a[ii]>0).mean()})
  print(cond,len(cells),len(uniq),ng,flush=True)
 pd.DataFrame(allrows).to_csv(OUT/'mouse_library_expression.csv',index=False);pd.DataFrame(details).to_csv(OUT/'mouse_library_summary.csv',index=False)
if __name__=='__main__':
 if sys.argv[1]=='human':human()
 else:mouse()
