from pathlib import Path
import numpy as np,pandas as pd,json
from scipy import sparse,stats
from discover import bh
B=Path(__file__).parent;O=B/'results';checks={}
z=np.load(O/'human_pseudobulk.npz');m=pd.read_csv(O/'human_groups.csv');checks['human_nuclei']=int(m.n_cells.sum());checks['human_features']=len(z['genes']);checks['unpooled_donors']=m.loc[~m.donor_id.str.startswith('pooled'),'donor_id'].nunique();assert checks['human_nuclei']==372816 and checks['human_features']==35475 and checks['unpooled_donors']==94
for c in ['control','14d','atlas']:
 x=sparse.load_npz(O/f'mouse_{c}.npz');mm=pd.read_csv(O/f'mouse_{c}_metadata.csv');assert x.shape[1]==len(mm);assert mm.cell.is_unique
 checks[c+'_qc_cells']=len(mm)
 if c!='atlas':
  old=np.load(B.parent/f'results/mouse_{c}_candidate_counts.npz');checks[c+'_original_qc']=int(old['qc'].sum());assert old['qc'].sum()==len(mm)
g=pd.read_csv(O/'mouse_genomewide.csv');assert len(g)==11169 and g.mouse.is_unique and g.human.is_unique;assert np.allclose(g.q,bh(g.p));checks['mouse_q05']=int((g.q<.05).sum())
s=pd.read_csv(O/'signature.csv');assert len(s)==50 and (s.direction==1).sum()==25;checks['signature_concordant']=int(((s.atlas_effect*s.control_effect)>0).sum());checks['known_markers_in_signature']=s.loc[s.human.isin(['IGF1','OPN4','SPP1','TIMP2','NDNF','PRPH']),'human'].tolist()
h=pd.read_csv(O/'human_signature_results.csv');d=pd.read_csv(O/'human_signature_donors.csv');assert len(h)==20 and np.allclose(h.q,bh(h.p));means=d.groupby(['reference','subtype']).score.mean();assert np.allclose(h.score,[means.loc[r.reference,r.subtype] for r in h.itertuples()]);assert not d.donor.str.startswith('pooled').any()
a=pd.read_csv(O/'human_candidate_results.csv');assert len(a)==880 and np.allclose(a.q,bh(a.p));checks['human_tests']=len(a)
w=a.pivot(index=['gene','subtype'],columns='reference',values=['effect','q','positive_fraction']).dropna();keep=(w.effect.min(axis=1)>.5)&(w.q.max(axis=1)<.05)&(w.positive_fraction.min(axis=1)>=.75);checks['both_reference_candidate_pairs']=int(keep.sum());w.loc[keep].to_csv(O/'robust_candidates_verified.csv');assert keep.sum()==33
checks['baseline_effect_correlation']=float(stats.pearsonr(g.atlas_effect,g.control_effect).statistic)
(O/'verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
