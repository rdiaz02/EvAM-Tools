Five example data sets for upload. ovx and ov and ov3 are from the Oncotree package[1] , using either the original names of genes or a transformation or a much shorter one-letter name for each gene (all work, though, even if for safety we say no gene names should start with a number, and only the ones with short gene names can, of course, be displayed nicely in histograms in output plots).

The second is a tiny data to easily do things such as make a whole gene of 0 presence by removing a genotype, etc.

The remaining two are single column and a 0 column file, for testing correct behavior.





[1] https://cran.r-project.org/web/packages/Oncotree/index.html

@Manual{,
    title = {Oncotree: Estimating Oncogenetic Trees},
    author = {Aniko Szabo and Lisa Pappas},
    year = {2022},
    note = {R package version 0.3.4},
    url = {https://CRAN.R-project.org/package=Oncotree},
  }
