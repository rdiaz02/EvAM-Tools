# EvAM Tools

## EvAM Tools: purpose
Tools for evolutionary accumulation, or event accumulation, models. For now, this refers to "Cancer Progression Models" (CPM) but these are not limited to cancer.

We provide an R package, evamtools, and a shiny (https://shiny.rstudio.com/) app that allows to:
  * Run state-of-the-art CPM methods, including Conjuntive Bayesian Networks (CBN), Oncogenetic trees (OT), Mutual Hazard Networks (MHN), and Hidden Extended Suppes-Bayes Causal Networks (H-ESBCNs ---PMCE---), with a single function call.
  * From the fitted models, represent, graphically, the fitted models (DAGs of restrictions or matrix of hazards, as appropriate), the transition matrices and transition rate matrices (where appropriate) between genotypes.
  * Using the shiny app, easily visualize the effects of changes in genotype composition on the fitted models by entering user-defined cross-sectional data using a GUI.

<!-- Remember the long name of output, as seen in server.R, around line 1344 -->
<!--  choiceNames =  c("Transition rates", -->
<!--                                          "Genotype transitions counts", -->
<!--                                          "Genotype frequencies", -->
<!--                                          "Transition probabilities", -->
<!--                                          "Lambdas/probabilities", -->
<!--                                          "Time-discretized transition matrix"), -->
<!--                         choiceValues =  c("trans_rate_mat", -->
<!--                                           "genotype_transitions", -->
<!--                                           "freqs", -->
<!--                                           "trans_mat", -->
<!--                                           "lambdas", -->
<!--                                           "td_trans_mat"), -->


### References ###

#### OT ####


- Szabo, A., & Boucher, K. M. (2008). Oncogenetic Trees. In W. Tan, & L. Hanin
  (Eds.), Handbook of Cancer Models with Applications (pp. 1–24). : World
  Scientific.


- Desper, R., Jiang, F., Kallioniemi, O. P., Moch, H., Papadimitriou, C. H., &
  Sch\"affer, A A (1999). Inferring tree models for oncogenesis from comparative
  genome hybridization data. J Comput Biol, 6(1), 37–51.


#### CBN ####


- Gerstung, M., Baudis, M., Moch, H., & Beerenwinkel, N. (2009). Quantifying
  cancer progression with conjunctive Bayesian networks. Bioinformatics, 25(21),
  2809–2815. http://dx.doi.org/10.1093/bioinformatics/btp505


- Gerstung, M., Eriksson, N., Lin, J., Vogelstein, B., & Beerenwinkel,
  N. (2011). The Temporal Order of Genetic and Pathway Alterations in
  Tumorigenesis. PLoS ONE, 6(11),
  27136. http://dx.doi.org/10.1371/journal.pone.0027136


- Montazeri, H., Kuipers, J., Kouyos, R., B\"oni, J\"urg, Yerly, S., Klimkait,
  T., Aubert, V., … (2016). Large-scale inference of conjunctive Bayesian
  networks. Bioinformatics, 32(17),
  727–735. http://dx.doi.org/10.1093/bioinformatics/btw459


#### MHN ####


- Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer
  progression using Mutual Hazard Networks. Bioinformatics, 36(1),
  241–249. http://dx.doi.org/10.1093/bioinformatics/btz513

#### HESBCN (PMCE) ####


- Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., &
  Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer
  evolution with high prognostic power. Bioinformatics, 38(3),
  754–762. http://dx.doi.org/10.1093/bioinformatics/btab717

 (About terminology: we will often refer to HESBCN, as that is the program we use, as shown here: https://github.com/danro9685/HESBCN. H-ESBCN is part of the PMCE procedure).


#### Conditional prediction of genotypes and probabilities of paths from CPMs ####

- Diaz-Colunga}, J., & Diaz-Uriarte, Ramon (2021). Conditional prediction of
  consecutive tumor evolution using cancer progression models: What genotype
  comes next? PLOS Computational Biology, 17(12),
  1009055. http://dx.doi.org/10.1371/journal.pcbi.1009055


- Diaz-Uriarte}, R., & Vasallo, C. (2019). Every which way? On predicting tumor
  evolution using cancer progression models. PLOS Computational Biology, 15(8),
  1007246. http://dx.doi.org/10.1371/journal.pcbi.1007246






## Copyright and origin of files ##

### Copyright and origin of files under evamtools/R ###

- All files under evamtools/R are copyright Pablo Herrera Nieto and Ramon Diaz-Uriarte (and released under the GPL v3 license) except for the following:

- File HESBCN__import.hesbcn.R: 
   This file contains function import.hesbcn (with a minor modification to return "Best Lambdas").
   
   Code from https://github.com/BIMIB-DISCo/PMCE/blob/main/Utilities/R/utils.R
   Commit 5578c79 from 2021-09-29

   License: Apache License 2.0, which is compatible with the GPL 3 used by the rest of this project.
   
   Author of code: from commit history, most likely Daniele Ramazzotti (danro9685)

   Authors of project: F. Angaroni, K. Chen, C. Damiani, G. Caravagna, A. Graudenzi, D. Ramazotti.
   
  Paper:  Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., &
  Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer
  evolution with high prognostic power. Bioinformatics, 38(3): 754-762. http://dx.doi.org/10.1093/bioinformatics/btab717


- Files MHN__*.R: MHN__UtilityFunctions.R, MHN__RegularizedOptimization.R, MHN__ModelConstruction.R, MHN__Likelihood.R, MHN__InlineFunctions.R,,  MHN__ExampleApplications.R

  Files obtained from https://github.com/RudiSchill/MHN
  Commit 49a8cc0 from 2018-08-16
  We have added the "MHN__" and made minor modifications to conform to usage within an R package. We have moved the inline C code to MHN.c and done the rest of the scaffolding for it to be used from the R package.
  
  License: no license information available in the repository nor the files.
  
  Author of code: from commit history, most likely Rudolf Schill.
  
  Authors of paper/project: Schill, R., Solbrig, S., Wettig, T., & Spang, R.
  
  Paper: Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1),  241–249. http://dx.doi.org/10.1093/bioinformatics/btz513

- This information is also provided under evamtools/inst/COPYRIGHTS and in the header of the files themselves, as comments.
  

- The authors of the above code have been added to the DESCRIPTION file, under "Auhor".


- Some of the code for transition rate matrices, input from different CPMs, tests, etc (authored by Ramon Diaz-Uriarte, released under the GPL-3) has been previously used in Diaz-Uriarte and Vasallo, 2019 (all code available from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007246#sec019, files S1 Dataset and S2 Dataset) and Diaz-Colunga and Diaz-Uriarte, 2021 (repository: https://github.com/rdiaz02/what_genotype_next).

    
### ct-cbn: ###

This repository includes code by:

- Niko Beerenwinkel, Moritz Gerstung, and Seth Sullivant. This is file ct-cbn-0.1.04b-RDU.tar.gz, a
   modification by RDU of the ct-cbn code from
   https://bsse.ethz.ch/cbg/software/ct-cbn.html.
   - The modification involves a minor bug fix (which, however, could be related to non-identifiability) as well as producing output with lambdas and likelihood from the initial run and each of the iterations.
  


## How to install to run just the R package? ##

  * Install CBN
    - Use the file ct-cbn-0.1.04b-RDU.tar.gz.
    - Uncompress. Then the usual ./configure, make.
    - Put the h-cbn and ct-cbn executables in your $PATH.
  * Install HESBCN
    - From  https://github.com/danro9685/HESBCN . <!-- For proper compilation you should modify line 23 in *Makefile* to include *LDLIBS = -lgsl -lm -lgslcblas*. You can also change compilation from gcc-5 to something a bit more up to date, like gcc-10.  -->
	When compiled add the h-esbcn executable to your $PATH.
  * Install the evamtools package
     - Make sure you have the require dependencies and imports, as listed in DESCRIPTION: igraph, OncoSimulR, stringr, Matrix, parallel, Oncotree , gtools , plot.matrix , DT, shinyjs, shiny
     - Build (R CMD build evamtools) and install (R CMD INSTALL evamtools_x.y.z.tar.gz, with x.y.z replaced by the current version number). File `build-test.sh` builds, tests, and installs the package (and takes care of the version number).
      
	  
<!-- In order to run all CPMs you will need to clone their repositories: -->
<!-- include github repo, and instructions to get them running -->
<!-- * DBN: https://github.com/phillipnicol/OncoBN. -->
<!-- FIXME: commit hash and date when we got the package. -->
<!-- * HyperTraps: https://github.com/sgreenbury/HyperTraPS. Follow the instruction there to install the package through conda  -->
<!-- and get a working enviroment with HyperTraPS installed. Add both the `bin` and the `src/python` to your `$PATH`. Also, add the `bin` folder to your $PATH. It is also convenient to run `echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc`. -->

	  
## How to install to also run the Shiny app? [FIXME: Pablo writes this](#) ##

  * Install the evamtools R package as explained [above](#how-to-install-to-run-just-the-r-package).
  * Then ... [FIXME: Pablo writes this](#)



## How to run the shiny app locally  [FIXME: Pablo writes this](#) ##

## How to create a Docker image  [FIXME: Pablo writes this](#) ##

### How to update the Docker image if you change the code  
[FIXME: Pablo writes this](#) ###

## How to run the shiny app in a local intranet  [FIXME: Pablo writes this](#) ##

### From the Docker image ###

### Without the Docker image ###



## Main files and directories

### Dockerfile    [FIXME: Pablo writes this](#)  
### docker   [FIXME: Pablo writes this](#) 
### evamtools
The R package itself with standard organization. Directories and files under inst:
  * shiny-examples: [FIXME: Pablo writes this](#) Explicar los subdirectorios y ficheros principales.
  * miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex: explanation of using OncoSimulR to check transition matrices, the equivalence of lambdas to terms in fitness expressions, and interpretation of the lambdas for HESBCN with OR and XOR.
  * miscell/examples: examples referred to from other files (for example, from the former tex file).
  * miscell/tests-sample_genotypes_from_trm: output of tests that were run to verify the code for sampling genotypes from the transition rate matrices. We compared the output of our code with that from the code of the original authors (MHN, MCCBN) for a large set of cases.
      
	Note that the R package uses testthat to test our R code. Those tests will run automatically with the usual procedures from testthat or while doing `R CMD check`. For example, we check that transition rate matrices and transition probability matrices give identical results when compared to finding them via OncoSimulR (file test.OT-CBN-trans-mat-against-oncosimul.R and test.HESBCN-trans-mat-against-oncosimul.R) and against hand-computed examples (file test.trans-rates-f-graphs.R and test.HESBCN-transition-rate-matrices.R). The tests in evamtools/tests/testthat are separate from the tests under  inst/miscell/tests-sample_genotypes_from_trm  
	




[FIXME: Pablo writes this: Selenium tests of the server](#)


