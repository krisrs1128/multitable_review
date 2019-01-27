# A Review of Multitable Methods

This is code to accompany *Multitable Methods for Microbiome Data Integration* by Sankaran and Holmes.

To rerun any particular method, edit the associated `read_data()` command into `read_data(simulate=TRUE)`, and put simulated data [`sim.rda`](https://drive.google.com/file/d/1wZpOIofWjzAApIrIG54-hqqjqkDVzZJK/view?usp=sharing) in a `data/` subdirectory of the main repo. To run all methods, just use `bash run_all.sh` in the main repository.

Alternatively, you can pull a docker image with the associated data and all packages pre-installed by visiting [this page](https://hub.docker.com/r/krisrs1128/multitable_review/).

Implementations of specific methods are given in the `dimension_red` (CCA, CCpnA, CoIA, LDA, PCA, PCA-IV, and PMD) and `supervised` (Graph-Fused and Multitask Lasso) folders.
