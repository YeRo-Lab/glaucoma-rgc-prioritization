from pathlib import Path
import gzip,csv
import numpy as np,pandas as pd
from scipy import sparse
B=Path(__file__).parent;O=B/'results'
ann=pd.read_csv(B/'mouse_atlas_cluster.tsv',sep='\t').set_index('NAME').Cluster
with gzip.open(B/'mouse_atlas.csv.gz','rt') as f:
 cells=np.array(next(csv.reader([f.readline()]))[1:]);labels=ann.reindex(cells);keep=labels.notna().to_numpy() & labels.ne('46').to_numpy() & labels.ne('Unassigned').to_numpy();rows=[];genes=[]
 for line in f:
  g,s=line.split(',',1);v=np.fromstring(s,sep=',');rows.append(sparse.csr_matrix(v[keep][None,:],dtype=np.float32));genes.append(g.strip('"'))
  if len(genes)%10000==0:print('atlas genes',len(genes),flush=True)
x=sparse.vstack(rows,format='csr');lib=np.asarray(x.sum(axis=0)).ravel();nf=np.asarray((x>0).sum(axis=0)).ravel();genes=np.array(genes);mt=np.asarray(x[np.char.startswith(np.char.lower(genes),'mt-')].sum(axis=0)).ravel();ident=np.asarray(x[np.isin(genes,['Rbpms','Slc17a6'])].sum(axis=0)).ravel();qc=(lib>=1000)&(nf>=500)&(mt/np.maximum(lib,1)<=.2)&(ident>0)
sparse.save_npz(O/'mouse_atlas.npz',x[:,qc]);np.save(O/'mouse_atlas_genes.npy',genes);cc=cells[keep][qc];pd.DataFrame({'cell':cc,'library':[s.split('_')[0] for s in cc],'subtype':labels.to_numpy()[keep][qc]}).to_csv(O/'mouse_atlas_metadata.csv',index=False);print('atlas raw',len(cells),'annotated',keep.sum(),'qc',qc.sum(),flush=True)
