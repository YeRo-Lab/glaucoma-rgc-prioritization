from pathlib import Path
import pandas as pd
B=Path(__file__).parent;O=B/'results'
a=pd.read_csv(B/'mouse_cluster.tsv',sep='\t').set_index('NAME').Cluster
assert a.index.is_unique
for c in ['control','14d']:
 m=pd.read_csv(O/f'mouse_{c}_cells.csv');m['subtype']=m.cell.map(a)
 old=O/f'mouse_{c}_metadata.csv'
 if old.exists():pd.testing.assert_frame_equal(m,pd.read_csv(old))
 m.to_csv(old,index=False)
print('Cell-to-subtype annotation joins verified')
