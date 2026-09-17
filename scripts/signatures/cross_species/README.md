# Cross-Species Signature Analyses

This directory contains scripts that test a signature defined in one species within the single-cell object from the other species.

| Direction | Status |
| --- | --- |
| `mouse_to_human/` | Contains the cluster 2 IFN/ISG responder-vs-non-responder interrogation in the human object. |
| `human_to_mouse/` | No final uploaded script was present for this direction at the time of repository cleanup. |

`translate_mouse_signatures_to_human.R` is an optional transparency helper that converts the final mouse signature gene vectors to human orthologs. It is not required to run the manuscript analyses; its `orthogene` dependency is included in `envs/renv.lock`.

Human scripts read the local object path from `HUMAN_SEURAT_RDS`, defaulting to `data/human_seurat.rds`. Human objects and generated outputs are local derived data and are not tracked in git.
