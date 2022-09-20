Example data sets for upload to the web app.

ovx, ov2, and ov3 are from the Oncotree package[1] , using either the original names of genes or a transformation or a much shorter one-letter name for each gene (all work, though, even if for safety we say no gene names should start with a number, and only the ones with short gene names can, of course, be displayed nicely in histograms in output plots).

tinydata.csv is a tiny data set to easily do things such as make a whole gene of 0 presence by removing a genotype, etc.

d1col.csv and d0col.csv are single column and a 0 column file, for testing correct behavior.

BRCA_ba_s.csv is a data set for  basal-like subtypes from [2, 3], originally from [4] (see Supplementary File S5\_Text, \url{https://doi.org/10.1371/journal.pcbi.1007246.s007} of [5] for full details about data origins and preprocessing)


When downloading the data, *make sure you download the raw data*. You can do this as follows:

- clone the repository (or get a zip file, and uncompress);

- click on the file, and then click on "Raw"; now, use your web browser to download/save the file.




[1] https://cran.r-project.org/web/packages/Oncotree/index.html

@Manual{,
    title = {Oncotree: Estimating Oncogenetic Trees},
    author = {Aniko Szabo and Lisa Pappas},
    year = {2022},
    note = {R package version 0.3.4},
    url = {https://CRAN.R-project.org/package=Oncotree},
  }


[2] @article{,
  title = {The {{cBio Cancer Genomics Portal}}: {{An Open Platform}} for {{Exploring Multidimensional Cancer Genomics Data}}: {{Figure}} 1.},
  author = {Cerami, Ethan and Gao, JianJiong and Dogrusoz, Ugur and Gross, Benjamin E and Sumer, Selcuk Onur and Aksoy, B"ulent Arman and Jacobsen, Anders and Byrne, Caitlin J and Heuer, Michael L and Larsson, Erik and Antipin, Yevgeniy and Reva, Boris and Goldberg, Arthur P and Sander, Chris and Schultz, Nikolaus},
  year = {2012},
  month = may,
  journal = {Cancer Discovery},
  volume = {2},
  number = {5},
  pages = {401--404},

}

[3] @article{, 
  title = {Integrative Analysis of Complex Cancer Genomics and Clinical Profiles Using the {{cBioPortal}}.},
  author = {Gao, JianJiong and Aksoy, B"ulent Arman and Dogrusoz, Ugur and Dresdner, Gideon and Gross, Benjamin and Sumer, S Onur and Sun, Yichao and Jacobsen, Anders and Sinha, Rileen and Larsson, Erik and Cerami, Ethan and Sander, Chris and Schultz, Nikolaus},
  year = {2013},
  month = apr,
  journal = {Science Signaling},
  volume = {6},
  number = {269},
  pages = {pl1},
  doi = {10.1126/scisignal.2004088},
}

[4] @article{,
  title = {Comprehensive Molecular Characterization of Human Colon and Rectal Cancer},
  author = {{Cancer Genome Atlas Research Network}},
  year = {2012},
  month = jul,
  journal = {Nature},
  volume = {487},
  number = {7407},
  pages = {330--337},
  issn = {1476-4687},
  doi = {10.1038/nature11252},
  url = {https://www.nature.com/articles/nature11252},

}

[5] @article{, 
  title = {Every Which Way? {{On}} Predicting Tumor Evolution Using Cancer Progression Models ({{bioRxiv}})},
  shorttitle = {Every Which Way?},
  author = {{Diaz-Uriarte}, Ramon and Vasallo, Claudia},
  year = {2019},
  journal = {bioRxiv},
  pages = {371039},
  doi = {10.1101/371039},
  url = {https://www.biorxiv.org/content/10.1101/371039v5},
}
