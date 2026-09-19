from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import sparse,stats
B=Path(__file__).parent;O=B/'results';RES={10,22,31,33,40,42,43}
def bh(p):
 p=np.nan_to_num(p,nan=1.);ix=np.argsort(p);q=np.empty(len(p));q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def aggregate(cond):
 m=pd.read_csv(O/f'mouse_{cond}_metadata.csv');genes=np.load(O/f'mouse_{cond}_genes.npy');x=sparse.load_npz(O/f'mouse_{cond}.npz');keep=m.subtype.notna() & ~m.subtype.isin(['46','Unassigned']);m=m.loc[keep].copy();x=x[:,keep.to_numpy()]
 batchfile='mouse_atlas_batch.tsv' if cond=='atlas' else 'mouse_crush_batch.tsv';ba=pd.read_csv(B/batchfile,sep='\t').set_index('NAME').iloc[:,0];m['batch']=m.cell.map(ba);m['type']=m.subtype.str.split('_').str[0].astype(int)
 unit='batch' if cond=='atlas' else 'library';group=m.groupby([unit,'type'],sort=True);idx=group.ngroup().to_numpy();meta=group.size().reset_index(name='cells');a=sparse.csr_matrix((np.ones(len(m)),(np.arange(len(m)),idx)),shape=(len(m),len(meta)));y=np.asarray((x@a).todense()).T;lib=y.sum(axis=1);cpm=1e6*y/lib[:,None];e=np.log2(1+cpm)
 meta.to_csv(O/f'{cond}_pseudobulk_groups.csv',index=False);np.savez_compressed(O/f'{cond}_pseudobulk.npz',counts=y,genes=genes,logcpm=e,cpm=cpm)
 return meta,genes,e,cpm
def run():
 om=pd.read_csv(B/'orthology.tsv',sep='\t',low_memory=False);pairs=[]
 for _,v in om.groupby('DB Class Key'):
  m=v.loc[v['NCBI Taxon ID']==10090,'Symbol'].unique();h=v.loc[v['NCBI Taxon ID']==9606,'Symbol'].unique()
  if len(m)==len(h)==1:pairs.append((m[0],h[0]))
 ort=pd.DataFrame(pairs,columns=['mouse','human']).drop_duplicates();ort=ort[~ort.mouse.duplicated(False)&~ort.human.duplicated(False)];ort.to_csv(O/'one_to_one_orthologs.csv',index=False)
 ma,ga,ea,ca=aggregate('atlas');mc,gc,ec,cc=aggregate('control');mi,gi,ei,ci=aggregate('14d')
 hg=np.load(O/'human_pseudobulk.npz')['genes'];hg=hg[~pd.Series(hg).duplicated(False).to_numpy()];ort=ort[ort.mouse.isin(ga)&ort.mouse.isin(gc)&ort.human.isin(hg)];gidx=pd.Series(np.arange(len(ga)),index=ga);cidx=pd.Series(np.arange(len(gc)),index=gc);aidx=gidx[ort.mouse].values;bidx=cidx[ort.mouse].values
 types=np.arange(1,46);E=[];C=[];CP=[]
 for t in types:
  ii=(ma.type==t)&(ma.cells>=5);jj=(mc.type==t)&(mc.cells>=5)
  E.append(ea[ii][:,aidx].mean(axis=0));CP.append(ca[ii][:,aidx].mean(axis=0));C.append(ec[jj][:,bidx].mean(axis=0) if jj.any() else np.full(len(bidx),np.nan))
 E=np.array(E);C=np.array(C);CP=np.array(CP);y=np.isin(types,list(RES));ok=(CP>=1).sum(axis=0)>=5;ok &=np.isfinite(E).all(axis=0)
 E=E[:,ok];C=C[:,ok];ort=ort.iloc[np.flatnonzero(ok)].reset_index(drop=True)
 ef=E[y].mean(axis=0)-E[~y].mean(axis=0);cf=np.nanmean(C[y],axis=0)-np.nanmean(C[~y],axis=0);u,p=stats.mannwhitneyu(E[y],E[~y],axis=0,alternative='two-sided',method='asymptotic');q=bh(p)
 result=ort.copy();result['atlas_effect']=ef;result['control_effect']=cf;result['auc']=u/(y.sum()*(~y).sum());result['p']=p;result['q']=q
 # Leave one atlas biological batch out; do not treat its sequencing channels as independent replicates.
 loo=[]
 for batch in sorted(ma.batch.unique()):
  Z=np.array([ea[(ma.type==t)&(ma.cells>=5)&(ma.batch!=batch)][:,aidx][:,ok].mean(axis=0) for t in types]);loo.append(Z[y].mean(axis=0)-Z[~y].mean(axis=0))
 result['min_batch_loo_effect']=np.min(loo,axis=0);result['max_batch_loo_effect']=np.max(loo,axis=0)
 result.to_csv(O/'mouse_genomewide.csv',index=False);np.savez_compressed(O/'mouse_type_expression.npz',atlas=E,control=C,types=types,resilient=y,mouse=result.mouse.to_numpy(dtype=str),human=result.human.to_numpy(dtype=str))
 # Fixed-size signature uses atlas alone. Control data provide separate baseline reproducibility.
 rank=ef/(E.std(axis=0,ddof=1)+1);up=np.argsort(rank)[-25:][::-1];down=np.argsort(rank)[:25];sel=np.r_[up,down];sig=result.iloc[sel].copy();sig['direction']=np.r_[np.ones(25),-np.ones(25)];sig.to_csv(O/'signature.csv',index=False)
 scores=[];known=['Igf1','Opn4','Spp1','Timp2','Ndnf','Prph'];knownidx=np.flatnonzero(result.mouse.isin(known));ks=[]
 for i in range(45):
  train=np.arange(45)!=i;X=E[train];lab=y[train];r=(X[lab].mean(0)-X[~lab].mean(0))/(X.std(0,ddof=1)+1);pos=np.argsort(r)[-25:];neg=np.argsort(r)[:25];z=(E[i]-X.mean(0))/(X.std(0,ddof=1)+1);scores.append(z[pos].mean()-z[neg].mean());ks.append(z[knownidx].mean())
 def auc(s):return float(stats.mannwhitneyu(np.array(s)[y],np.array(s)[~y]).statistic/(7*38))
 pd.DataFrame({'type':types,'published_resilient':y,'heldout_score':scores,'known_marker_score':ks}).to_csv(O/'mouse_crossvalidation.csv',index=False)
 print('genes',len(result),'q05',int((q<.05).sum()),'resup',int(((q<.05)&(ef>0)).sum()),'cv auc',auc(scores),'known auc',auc(ks),flush=True)
 print(sig[['mouse','atlas_effect','control_effect','q','direction']].to_string(index=False),flush=True)
 (O/'discovery_summary.json').write_text(json.dumps({'genes_tested':len(result),'q05':int((q<.05).sum()),'cv_auc':auc(scores),'known_auc':auc(ks),'atlas_cells':int(ma.cells.sum()),'control_cells':int(mc.cells.sum()),'injured_cells':int(mi.cells.sum()),'resilient_types':sorted(RES)},indent=2))
if __name__=='__main__':run()
