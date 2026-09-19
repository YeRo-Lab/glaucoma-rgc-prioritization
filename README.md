# Human RGC Subtype and Candidate Target Prioritization in Glaucoma

Author: Yeganeh Madadi, Computer Science Department, Appalachian State University, NC, USA.

This repository contains the computational analysis underlying the manuscript. It contains numerical results, the 50-gene panel, the 44-gene mouse candidate set, 33 qualifying human gene–subtype pairs, figures, source annotations, analysis scripts, source checksums, and a verification report. Raw matrices are not bundled because the human file alone is 7.34 GB.

## Scope and interpretation

The mouse source is Tran et al. 2019: adult atlas GSE133382 and control/day-14 ONC GSE137398, within GSE137400. Adding the adult atlas expands analysis within the same agreed Tran study; Rheaume and human organoid datasets are not used. Labels are inherited from Tran. The human source is Li/Rui Chen 2026, CELLxGENE collection https://cellxgene.cziscience.com/collections/4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e and dataset 95ce06a0-7de9-4ab2-bcd1-a6625f55e7e3. Exact downloadable files and hashes appear in source_manifest.json.

Human results concern healthy retinal expression. They do not identify experimentally confirmed human glaucoma-resilient or susceptible types. The 33 pairs are expression associations, not 33 novel biomarkers. Statistical associations do not establish causal neuroprotection. Both primary and exploratory findings are included; no preregistration is claimed. The mouse 50-gene score does not outperform established markers on internal assessment. Literature novelty is not guaranteed by this reanalysis.

## Main tables in analysis/results (inside analysis_package.zip)

- mouse_genomewide.csv: all 11,169 tested orthologs; atlas/control effects, subtype-level p and BH q, leave-batch-out effects. Seven resilient versus 38 other types; related types are not independent animals.
- signature.csv: positive and negative 25-gene panels. Not all 50 genes are individually significant.
- mouse_candidate_set.csv: 44 genes meeting mouse criteria.
- human_signature_results.csv and human_signature_donors.csv: 20 comparisons, donor scores, bootstrap intervals, signed-rank tests and expression-matched competitive tests.
- human_candidate_results.csv: 880 tests with one BH correction family.
- robust_candidates.csv: 33 qualifying gene–subtype pairs, requiring both references; illustrative manuscript table is a subset.
- human_direct_pairs.csv: three exploratory shared-sample comparisons, BH correction across three.
- human_sensitivity.csv: remove four known panel markers; omit contributing studies. A study label is taken from each donor's contributing sample metadata.
- human_cell_thresholds.csv: 10/30-cell sensitivity; BH across all eligible threshold comparisons together.
- mouse_family_holdout.csv: internal family holdout, not external validation.
- mouse_within_type_injury.csv: descriptive within-type changes. Manuscript uses rows with strata == 2 (23 types). Single-stratum rows are included for transparency but not in Fig. 3B.
- verification.json: independent consistency checks on totals, annotation joins, correction families, and candidate counts.

Effects use differences in log2(CPM + 1), not log fold changes. Bootstrap intervals describe mean scores; signed-rank tests need not agree with those intervals. All-other and nonmidget references use the target's own sample; eligible samples are averaged within donor. Human pooled donor records are excluded. Missing rare-type results indicate insufficient donor coverage, not biological absence.

## Reproduction

Use Python and versions in analysis/requirements.txt. At the extracted package root, retain the analysis/ directory structure. Download raw matrices to destinations in source_manifest.json and verify SHA-256. Annotation and orthology snapshots are included. A full rerun needs several GB of memory, substantial disk space, and processing time. Human data were accessed September 13; additional atlas/annotation/orthology resources were accessed September 14, 2026. Analyses were finalized September 18.

First extract `analysis_package.zip` into a separate working directory. It includes the result tables and figures as well as the same scripts displayed in this repository. Run the following commands from that extracted directory:

```sh
python analyze.py mouse
python analysis/prepare.py mouse
python analysis/annotate.py
python analysis/atlas.py
python analysis/prepare.py human
python analysis/discover.py
python analysis/human.py
python analysis/sensitivity.py
python analysis/injury.py
python analysis/verify.py
python analysis/figures.py
```

The first script is retained solely to reproduce full-matrix mouse QC and cache masks. Its optional legacy human analysis is not part of this manuscript and should not be run. prepare.py imports its HDF5 helper. The manuscript's baseline mouse discovery uses adult biological batches, not sequencing channels as independent biological replicates. verify.py regenerates robust_candidates_verified.csv with the same 33 qualifying pairs; robust_candidates.csv is a more readable merged version sorted by minimum effect.

All numerical stages above were run on the source data; annotation and result verification passed. The complete sequence was not repeated in a newly isolated software environment. Raw caches are intentionally omitted. Supplied result tables permit checking reported numbers without downloading raw data; plotting uses the supplied CSV files.


## Availability and validation

Raw sequencing matrices are not redistributed. Source URLs and SHA-256 checksums are in `analysis/source_manifest.json`. Download instructions and the analysis sequence above preserve the original directory structure. The archive `analysis_package.zip` contains the same analysis files and results for convenient download. A clean-environment end-to-end rerun has not been performed; this is stated separately from the completed numerical analyses.
