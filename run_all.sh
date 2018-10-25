mkdir -p chapter/figure/CCA
mkdir -p chapter/figure/ccpna
mkdir -p chapter/figure/coia
mkdir -p chapter/figure/graph_lasso
mkdir -p chapter/figure/lda_cca
mkdir -p chapter/figure/pca
mkdir -p chapter/figure/pca_iv
mkdir -p chapter/figure/pmd
mkdir -p chapter/figure/spls

cd dimension_red
Rscript pca.R
Rscript pmd.R
Rscript cca.R
Rscript ccpna.R
Rscript coia.R
Rscript lda_cca.R
Rscript pca_iv.R
Rscript illustration.R
Rscript illustration_pmd.R

cd ../supervised/
Rscript graph_lasso.R
Rscript multitask_lasso.R
Rscript spls.R
