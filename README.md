# EvAM Tools #
<!-- Create talk with markdown-toc, from node. -->
<!-- /markdown-toc ~/Proyectos/EvAM-Tools/README.md -i -->

<!-- toc -->

- [EvAM Tools: purpose](#evam-tools-purpose)
- [Copyright and origin of files](#copyright-and-origin-of-files)
  * [Copyright and origin of files under evamtools/R](#copyright-and-origin-of-files-under-evamtoolsr)
  * [ct-cbn](#ct-cbn)
- [How to install the R package](#how-to-install-the-r-package)
- [How to run the shiny app locally](#how-to-run-the-shiny-app-locally)
- [How to run from Docker images](#how-to-run-from-docker-images)
  * [Using the downloaded Docker image to run the Shiny app](#using-the-downloaded-docker-image-to-run-the-shiny-app)
  * [Using the downloaded Docker image to run the evamtools package from R](#using-the-downloaded-docker-image-to-run-the-evamtools-package-from-r)
- [How to create a Docker image](#how-to-create-a-docker-image)
        * [Creating the image fails because of no internet connection from the container](#creating-the-image-fails-because-of-no-internet-connection-from-the-container)
  * [How to update the Docker image if you change the code](#how-to-update-the-docker-image-if-you-change-the-code)
- [How to run the shiny app in a local intranet](#how-to-run-the-shiny-app-in-a-local-intranet)
  * [From the Docker image](#from-the-docker-image)
  * [How to run the Shiny app without the Docker image](#how-to-run-the-shiny-app-without-the-docker-image)
- [Main files and directories](#main-files-and-directories)
  * [Dockerfile](#dockerfile)
  * [docker](#docker)
  * [evamtools](#evamtools)
  * [Selenium tests of the server](#selenium-tests-of-the-server)
      - [Running the Selenium tests](#running-the-selenium-tests)
- [References](#references)
  * [OT](#ot)
  * [CBN](#cbn)
  * [MHN](#mhn)
  * [HESBCN (PMCE)](#hesbcn-pmce)
  * [OncoBN, DBN](#oncobn-dbn)
  * [Conditional prediction of genotypes and probabilities of paths from CPMs](#conditional-prediction-of-genotypes-and-probabilities-of-paths-from-cpms)

<!-- tocstop -->

## EvAM Tools: purpose
Tools for evolutionary accumulation, or event accumulation, models. For now, this refers to "Cancer Progression Models" (CPM) but these are not limited to cancer.

We provide an R package, evamtools, and a Shiny (https://shiny.rstudio.com/) app that allows to:
  * Run state-of-the-art CPM methods, including Conjuntive Bayesian Networks (CBN ---and their MC-CBN version), Oncogenetic trees (OT), Mutual Hazard Networks (MHN), Hidden Extended Suppes-Bayes Causal Networks (H-ESBCNs ---PMCE---), and Disjunctive Bayesian Networks (DBN) with a single function call.
  * From the fitted models, represent, graphically, the fitted models (DAGs of restrictions or matrix of hazards, as appropriate), the transition matrices and transition rate matrices (where appropriate) between genotypes.
  * Using the shiny app, easily visualize the effects of changes in genotype composition on the fitted models by entering user-defined cross-sectional data using a GUI.


For easier use, we provide links to Docker images that you can download and run, as well as instructions on how to build Docker images.


---
---


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

    
### ct-cbn ###

   This repository includes ct-​cbn-0.1.04b, from
   https://bsse.ethz.ch/cbg/software/ct-cbn.html, whose authors are Niko Beerenwinkel, Moritz Gerstung, and Seth Sullivant. It is released under the GNU GPL.
   
   The code included in this repo is file ct-cbn-0.1.04b-RDU.tar.gz, a
   modification by RDU of the above code that includes: a minor bug fix (which, however, could be related to non-identifiability); output with lambdas and likelihood from the initial run and each of the iterations.
  
   (For references about CBN see [References](#references)). 


---
---
## Installing and running ##

### Overview

You can install:
  * The [package and its dependencies] (#how-to-install-the-r-package)
  * A Docker image
	
You can run:
  * The package from R and, if you want, the shiny app 
  * The package in an RStudio session from the Docker image
  * The Shiny app from a Docker image
  
You can also build your own Docker image


### How to install the R package 

  * Install CBN
    - Use the file ct-cbn-0.1.04b-RDU.tar.gz.
    - Uncompress the directory. Then the usual configure, make dance. (Go inside the uncompressed directory and type `./configure`; when finished, type `make`).
    - Put the `h-cbn` and `ct-cbn` executables in your $PATH.
	
  * Install HESBCN
    - Clone the repository from  https://github.com/danro9685/HESBCN . 
	- Go inside that directory, and do the usual configure, make dance.
	- When finished, add the `h-esbcn` executable to your $PATH.
	
  * Install OncoBN
    - OncoBN is available from https://github.com/phillipnicol/OncoBN
    - But this should work: start R, install the devtools package if you don't have it, and then issue `devtools::install_github("phillipnicol/OncoBN")`.
	
  * Install MC-CBN
    - Go to https://github.com/cbg-ethz/MC-CBN and follow the installation instructions: https://github.com/cbg-ethz/MC-CBN#installation
    - If you have all dependencies installed, and the correct version of libboost, then from R this should work: `install.packages("https://github.com/cbg-ethz/MC-CBN/releases/download/v2.1.0/mccbn_2.1.0.tar.gz", repos=NULL)`
    - If it fails, make sure to review the installation instructions and then install manually: https://github.com/cbg-ethz/MC-CBN#installation-from-source 
  
  * Install the evamtools package
     - Make sure you have the require dependencies and imports, as listed in DESCRIPTION: igraph, OncoSimulR, stringr, Matrix, parallel, Oncotree , gtools , plot.matrix , DT, shinyjs, shiny, RhpcBLASctl.
         - Note that we list, as imports, OncoBN, mccbn. You need those (from above). 
     - Build (R CMD build evamtools) and install (R CMD INSTALL evamtools_x.y.z.tar.gz, with x.y.z replaced by the current version number). File `build-test.sh` builds, tests, and installs the package (and takes care of the version number).
      
	  

	  

## How to run the shiny app locally ##

Once the package is installed, to launch the shiny app open an R terminal and type

```
library(evamtools)
runShiny()
```


## How to run from Docker images ##
We provide two docker images, one for running the Shiny app, and another with  RStudio to run the evamtools package directly.  They are available from **FIXME**.

To use them:

### Using the downloaded Docker image to run the Shiny app ###
**FIXME**.

### Using the downloaded Docker image to run the evamtools package from R ###
**FIXME**.



## How to create a Docker image ##

From the `EvAM-Tools` directory run: 

```
docker build --tag evamtools .
``` 

to build the docker image. Of course, you need to have Docker installed for that to work. (Details about Docker are available here: https://docs.docker.com/get-docker/ .
Details about R with Docker and Rocker project here: https://www.rocker-project.org/ **


**FIXME: Pablo completes this**

- Does that need to be run as root?
- Are you sure this is using R-4.1.2? I think we might need something else:
    - https://github.com/rocker-org/rocker-versioned2
    - Please, give, explicitly, the *current version of R being used***	


** Creating the image fails because of no internet connection from the container**
Creating the above image requires installing R packages and that might fail because the Docker container cannot connect with the internet. The following might help: https://superuser.com/a/1582710 , https://superuser.com/a/1619378 . 



### How to update the Docker image if you change the code ###
Just as above. But this time it should run much faster, because many steps will be skipped.


---


## How to run the shiny app in a local intranet  ##

### From the Docker image ###

Once you have the Docker image built run the following command to run the image connecting port 3000 of the computer with port 3000 of the container

```
docker run -p 3000:3000 evamtools ##
```


**FIXME: Pablo writes this**

- Is this run as root or a user? Explain the pros and cons. If possible, set up to run as user, not root.



### How to run the Shiny app without the Docker image ###

To run the shiny app you may want to change the port (right now it run in 3000). This can be done by modifying line 27 in `evamtools/R/runShiny.R`, the one with
`shiny::runApp(appDir, port = 3000, host = "0.0.0.0", display.mode = "normal")`. Then, the app has to be launched as explained above and the corresponding port in the server has to be open to make it visible. 


---
---



## Main files and directories

### Dockerfile ###  

The Dockerfile includes all the information to create a container with all dependencies. It first uses a default image that includes the latest R version. Then install all R dependencies. Finally it also deals with the installation of third party code. 

### docker ###
**FIXME: Pablo writes this**


### evamtools
The R package itself with standard organization. Directories and files under inst:
  * shiny-examples: code for the shiny application. The application consists on two main files: `server.R` (that controls the logic) and `ui.R` (includes all the interface). There are two additional directories: `assets` (with files for the landing page) and `test_shiny_app` (with Selenium tests for the app and testing related files ---see [Selenium tests of the server](#selenium-tests-of-the-server)).
  * miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex: explanation of using OncoSimulR to check transition matrices for OT, CBN, OncoBN, and HESBCN, the equivalence of lambdas to terms in fitness expressions, interpretation of the lambdas for HESBCN with OR and XOR.
  * miscell/examples: examples referred to from other files (for example, from the former tex file).
  * miscell/tests-sample_genotypes_from_trm: output of tests that were run to verify the code for sampling genotypes from the transition rate matrices. We compared the output of our code with that from the code of the original authors (MHN, MCCBN) for a large set of cases.
      
	Note that the R package uses testthat to test our R code. Those tests will run automatically with the usual procedures from testthat or while doing `R CMD check`. For example, we check that transition rate matrices and transition probability matrices give identical results when compared to finding them via OncoSimulR (file test.OT-CBN-trans-mat-against-oncosimul.R and test.HESBCN-trans-mat-against-oncosimul.R) and against hand-computed examples (file test.trans-rates-f-graphs.R and test.HESBCN-transition-rate-matrices.R). The tests in evamtools/tests/testthat are separate from the tests under  inst/miscell/tests-sample_genotypes_from_trm  

### Selenium tests of the server ###
  
 The shiny web app also include its own set of tests that are run with Selenium. They are found under `inst/shiny_examples/evamtools/test_shiny_app/test.shiny.py.`
  There are test of the basic functionlity of the web page (navigation, loading files...), working with csd, DAG and matrix inputs (loading examples, modifying values, changes gene names and number), and checking interface behaviour in the results' tab.

##### Running the Selenium tests ####

Before launching the test you have to install a web driver. For testing, I use Chrome. To download the Chrome web driver go to [https://chromedriver.chromium.org/downloads]/(https://chromedriver.chromium.org/downloads) and select the version supporting your browser. Once the web driver is downloaded, extract it and make it available from $PATH.


To run selenium test first the shiny serve has to be running. 
Then, in a separate console, go to
`evamtools/inst/shiny-examples/evamtools/test_shiny_app` and type:

```bash
## To run all tests
python test.shiny.py 
## or
## To run a particular test
python test.shiny.py className.testName 
```

If you change the server where the shiny app is running, then the driver has to be changed. This is done in line 14 of `evamtools/inst/shiny-examples/evamtools/test_shiny_app/test.shiny.py` (the one with  `self.driver.get("http://127.0.0.1:3000/")`). Testing the public server directly can also be considered once it is up.


Running all test takes around 8 minutes, writes considerable amount of temporary files (more than 1Gb) and makes the terminal useless (windows will be continuously opening). Also, the tests are designed to run in a screen with 1377x768 resolution; they have adapted to bigger screens, but not fully tested in those settings. 


---
---



## References ##

### OT ###


- Szabo, A., & Boucher, K. M. (2008). Oncogenetic Trees. In W. Tan, & L. Hanin
  (Eds.), Handbook of Cancer Models with Applications (pp. 1–24). : World
  Scientific.


- Desper, R., Jiang, F., Kallioniemi, O. P., Moch, H., Papadimitriou, C. H., &
  Sch\"affer, A A (1999). Inferring tree models for oncogenesis from comparative
  genome hybridization data. J Comput Biol, 6(1), 37–51.


### CBN ###

- Beerenwinkel, N., & Sullivant, S. (2009). Markov models for accumulating
  mutations. Biometrika, 96(3), 645.

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


### MHN ###


- Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer
  progression using Mutual Hazard Networks. Bioinformatics, 36(1),
  241–249. http://dx.doi.org/10.1093/bioinformatics/btz513

### HESBCN (PMCE) ###


- Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., &
  Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer
  evolution with high prognostic power. Bioinformatics, 38(3),
  754–762. http://dx.doi.org/10.1093/bioinformatics/btab717

   (About terminology: we will often refer to HESBCN, as that is the program we use, as shown here: https://github.com/danro9685/HESBCN. H-ESBCN is part of the PMCE procedure).

### OncoBN, DBN ###

- Nicol, P. B., Coombes, K. R., Deaver, C., Chkrebtii, O., Paul, S., Toland,
  A. E., & Asiaee, A. (2021). Oncogenetic network estimation with disjunctive
  Bayesian networks. Computational and Systems Oncology, 1(2),
  1027. http://dx.doi.org/10.1002/cso2.1027




### Conditional prediction of genotypes and probabilities of paths from CPMs ###

- Diaz-Colunga}, J., & Diaz-Uriarte, Ramon (2021). Conditional prediction of
  consecutive tumor evolution using cancer progression models: What genotype
  comes next? PLOS Computational Biology, 17(12),
  1009055. http://dx.doi.org/10.1371/journal.pcbi.1009055


- Diaz-Uriarte}, R., & Vasallo, C. (2019). Every which way? On predicting tumor
  evolution using cancer progression models. PLOS Computational Biology, 15(8),
  1007246. http://dx.doi.org/10.1371/journal.pcbi.1007246

