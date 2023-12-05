# install CMS classifier
"$CONDA_PREFIX/bin/R" -e 'devtools::install_github("Sage-Bionetworks/CMSclassifier", build_vignettes = FALSE, quiet=TRUE, upgrade="never")'
