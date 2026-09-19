from pathlib import Path
import sys,json,gzip,csv
import numpy as np,pandas as pd,h5py
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from analyze import vec
B=Path(__file__).resolve().parents[1];O=B/'revision/results'
def human():
 f=h5py.File(B/'human_all_rgc.h5ad','r');obs=pd.DataFrame({k:vec(f['obs'][k]) for k in ['author_cell_type','donor_id','sampleid','tissue','study_name']})
 keys=['donor_id','sampleid','author_cell_type'];grp=obs.groupby(keys,sort=True).ngroup().to_numpy();meta=obs.groupby(keys,sort=True).agg(n_cells=('tissue','size'),tissue=('tissue','first'),study=('study_name','first')).reset_index()
 raw=f['raw/X'];names=vec(f['raw/var/feature_name']);p=raw['indptr'][:];N=len(obs);G=len(names);K=len(meta);y=np.zeros((K,G),dtype=np.float64);det=np.zeros((K,G),dtype=np.int32)
 for a in range(0,N,4000):
  b=min(a+4000,N);l,r=p[a],p[b];x=sparse.csr_matrix((raw['data'][l:r],raw['indices'][l:r],p[a:b+1]-l),shape=(b-a,G));m=sparse.csr_matrix((np.ones(b-a),(grp[a:b],np.arange(b-a))),shape=(K,b-a));y+=(m@x).toarray();x.data=np.ones_like(x.data);det+=(m@x).toarray().astype(np.int32)
  if a%40000==0:print('human',a,N,flush=True)
 np.savez_compressed(O/'human_pseudobulk.npz',counts=y,detected=det,genes=names);meta.to_csv(O/'human_groups.csv',index=False);print('human groups',K,flush=True)
def mouse():
 for cond in ['control','14d']:
  z=np.load(B/f'results/mouse_{cond}_candidate_counts.npz');qc=z['qc'];rows=[];names=[]
  with gzip.open(B/f'mouse_{cond}.csv.gz','rt') as f:
   cells=np.array(next(csv.reader([f.readline()]))[1:]);assert np.array_equal(cells,z['cells'])
   for line in f:
    g,s=line.split(',',1);x=np.fromstring(s,sep=',');rows.append(sparse.csr_matrix(x[qc][None,:],dtype=np.float32));names.append(g.strip('"'))
  mat=sparse.vstack(rows,format='csr');sparse.save_npz(O/f'mouse_{cond}.npz',mat);pd.DataFrame({'cell':cells[qc],'library':[x.split('_')[0] for x in cells[qc]]}).to_csv(O/f'mouse_{cond}_cells.csv',index=False);np.save(O/f'mouse_{cond}_genes.npy',np.array(names));print(cond,mat.shape,flush=True)
if __name__=='__main__':globals()[sys.argv[1]]()
