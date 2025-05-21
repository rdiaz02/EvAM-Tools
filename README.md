# EvAM Tools #

<!-- Simpler: use markdown-toc emacs package start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
<!-- Create toc with markdown-toc, from node. -->
<!-- markdown-toc README.md -i -->

<!-- toc -->

- [EvAM Tools: purpose](#evam-tools-purpose)
  * [Some examples of use](#some-examples-of-use)
- [Copyright and origin of files](#copyright-and-origin-of-files)
  * [Copyright and origin of files under evamtools/R](#copyright-and-origin-of-files-under-evamtoolsr)
  * [ct-cbn](#ct-cbn)
- [Installing and running](#installing-and-running)
  * [Overview of options](#overview-of-options)
  * [How to install the R package](#how-to-install-the-r-package)
  * [Docker images](#docker-images)
  * [How to run the R package and the shiny app locally without Docker](#how-to-run-the-r-package-and-the-shiny-app-locally-without-docker)
  * [How to run the R package from the Docker image](#how-to-run-the-r-package-from-the-docker-image)
  * [How to run the Shiny app from the Docker image](#how-to-run-the-shiny-app-from-the-docker-image)
- [Main files and directories](#main-files-and-directories)
  * [Dockerfiles](#dockerfiles)
  * [evamtools](#evamtools)
- [References](#references)
  * [Overview paper on CPMs and EvAMs](#overview-paper-on-cpms-and-evams)
  * [OT](#ot)
  * [CBN and MCCBN](#cbn-and-mccbn)
  * [MHN](#mhn)
  * [H-ESBCN (PMCE)](#h-esbcn-pmce)
  * [OncoBN (DBN)](#oncobn-dbn)
  * [HyperTraPS-CT](#hypertraps-ct)
  * [BML](#bml)
  * [Conditional prediction of genotypes and probabilities of paths from CPMs](#conditional-prediction-of-genotypes-and-probabilities-of-paths-from-cpms)
- [Additional documentation](#additional-documentation)
- [Citing EvAM-Tools](#citing-evam-tools)
- [Funding](#funding)
<!-- tocstop -->

## EvAM Tools: purpose
Tools for evolutionary accumulation models (EvAMs). We use code from what are usually referred to as "Cancer Progression Models" (CPM) but these are not limited to cancer (the key idea is that events are gained one by one, but not lost).

We provide an R package, evamtools, that can also launch a GUI as a Shiny (https://shiny.rstudio.com/) web app (running on https://www.iib.uam.es/evamtools/) that allows you to:

  * Run state-of-the-art EvAM methods, including Conjuntive Bayesian Networks (CBN ---and their MC-CBN version---), Oncogenetic trees (OT), Mutual Hazard Networks (MHN), Hidden Extended Suppes-Bayes Causal Networks (H-ESBCNs ---PMCE---), Disjunctive Bayesian Networks (DBN, from the OncoBN package), Hypercubic transition path sampling (HyperTraPS-CT), and Bayesian Mutational Landscape (BML), with a single function call.
  * From the fitted models, represent, graphically, the fitted models (DAGs of restrictions or matrix of hazards, as appropriate), the transition matrices and transition rate matrices (where appropriate) between genotypes and show frequencies of genotypes sampled from the fitted models.
  * Using the shiny app, easily visualize the effects of changes in genotype composition on the fitted models by entering user-defined cross-sectional data using a GUI.


For easier use, we provide links to Docker images that you can download and run, as well as instructions on how to build Docker images. You can also run the Shiny app on our servers: https://www.iib.uam.es/evamtools/ .


If you use the package or the web app, please **cite** the Bioinformatics paper:
Diaz-Uriarte, R & Herrera-Nieto, P. 2022. EvAM-Tools: tools for evolutionary accumulation and cancer progression models. _Bioinformatics_, 38 (24): 5457-5459.  https://doi.org/10.1093/bioinformatics/btac710 .


<b>Funding:</b> Supported by grant PID2019-111256RB-I00 funded by MCIN/AEI/10.13039/501100011033 and Comunidad de Madrid's PEJ-2019-AI/BMD-13961 to R. Diaz-Uriarte.
<p align="center">
<img src="micin-aei.png" alt="micin-aei logo" width="200"/>
</p>

---

### Some examples of use

What can EvAM-Tools be used for?

   * To understand EvAMs and CPMs and how different inputs (which can be easily modified interactively in the web app) affect the fitted models. For example, change the genotype frequencies according to sensible models of dependencies and run the EvAMs.

   * To understand what different models imply about how the cross-sectional data looks like. Create the dependency structures (DAGs or MHN log-Theta matrix) and generate data from them, possibly playing with the amount of noise. This does not even require to run the methods themselves.

   * As a research tool to analyze cross-sectional data with state-of-the-art EvAMs.

   * As a research tool in methodological work. For example:
       * We can examine how well a method can recover the true structure when the data fulfills the assumptions of a method. We would generate data under a particular model and see if the method that implements that model can recover the true structure under different sample sizes. (The web app only allows for playing with this; for serious work one would use the package itself and build code using the package functions).
       * We can examine how a give method works, and what type of inferences it performs, when data are generated under the model of another method. For example, what is the output from MHN if the data are really coming from an H-ESBCN model?
       * We can compare new methods against existing state-of-the-art ones (using different data sources ---simulated from EvAM-Tools or imported from other sources).



The web app in particular should allow for easy exploration of EvAM models without the need for installing any software. Details about the use of the web app are provided in https://www.iib.uam.es/evamtools/#input. The following figure provides an overview of the workflow with the web app:

![Overview EvAM-Tools web app](./evamtools/inst/shiny-examples/evamtools/assets/figure-overview.png?raw=true)

As can be seen, the web app workflows encompass different major functionalities and use cases, mainly:


1. Inference of EvAMs from user uploaded from a file.

2. Exploration of the inferences that different EvAM methods yield from manually constructed synthetic data.

3. Construction of EvAM models (DAGs and rates/probabilities or MHN models) and simulation of synthetic data from them.

    3.1. Examination of the consequences of different models and their parameters on the simulated data.

    3.2. Analysis of the data simulated under one model with methods that have different models (e.g., data simulated from CBN analyzed with OT and OncoBN).

	3.3. Analysis of the data simulated under model after manual modification of specific genotype frequencies (e.g., data simulated under CBN but where, prior to analysis, we remove all observations with the WT genotype and the genotype with all loci mutated).


Furthermore, note that in all cases, when data are analyzed, in addition to returning the fitted models, the web app returns the analysis of the EvAMs in terms of their predictions such as predicted genotype frequencies and transition probabilities between genotypes. Some of the above use cases (e.g., exhaustive examination of the performance of one method under data generated by a very different model) are probably better served using the R package directly; still, the web app makes it particularly easy to gain intuition about what these methods do with different data sets.


---
---


## Copyright and origin of files ##

### Copyright and origin of files under evamtools/R ###

- All files under evamtools/R are copyright Pablo Herrera Nieto, Ramon Diaz-Uriarte, and Javier Pérez de Lema Díez (and released under the GNU Affero General Public License (AGPL) v3 license) except for the following:

- File HESBCN__import.hesbcn.R:
   This file contains function import.hesbcn (with a minor modification to return "Best Lambdas").

   Code from https://github.com/BIMIB-DISCo/PMCE/blob/main/Utilities/R/utils.R .
   Commit 5578c79 from 2021-09-29.

   License: Apache License 2.0, which can be combined with software under the  AGPL 3 as used by the rest of this project.

   Author of code: from commit history, most likely Daniele Ramazzotti (danro9685)

   Authors of project: F. Angaroni, K. Chen, C. Damiani, G. Caravagna, A. Graudenzi, D. Ramazotti.

  Paper:  Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., &
  Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer
  evolution with high prognostic power. Bioinformatics, 38(3): 754-762. http://dx.doi.org/10.1093/bioinformatics/btab717


- Files MHN__*.R: MHN__UtilityFunctions.R, MHN__RegularizedOptimization.R, MHN__ModelConstruction.R, MHN__Likelihood.R, MHN__InlineFunctions.R,  MHN__ExampleApplications.R

  Files obtained from https://github.com/RudiSchill/MHN .
  Commit 49a8cc0 from 2018-08-16 (updated to reflect explicit MIT license on 2022-04-04). We have added the "MHN__" and made minor modifications to conform to usage within an R package. We have moved the inline C code in InlineFunctions.R (now MHN__InlineFuctions.R) to MHN.c and done the rest of the scaffolding for it to be used from the R package.

  License: MIT, which can be combined with software under the  AGPL 3 as used by the rest of this project.

  Author of code: Rudolf Schill (inferred from commit history).

  Authors of paper/project: Schill, R., Solbrig, S., Wettig, T., & Spang, R.

  Paper: Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1),  241–249. http://dx.doi.org/10.1093/bioinformatics/btz513


- File MCCBN__random_poset.R: This file contains function random_poset and additional required functions for random_poset to work.

  File obtained from file common.r from https://github.com/cbg-ethz/MC-CBN , commit a6eceb8 from 2022-02-05.

  License:  GNU General Public License v2.0, which can be combined with software under the AGPL 3 as used by the rest of this project.

  Author of code:  from comments and commit history (e.g., d2d8dfd, f08add1) most likely Hesam Montazeri.

  Authors of project: Hesam Montazeri, Susana Posada-Cespedes.


- This information is also provided under evamtools/inst/COPYRIGHTS and in the header of the files themselves, as comments.


- The authors of the above code have been added to the DESCRIPTION file, under "Auhor".


- Some of the code for transition rate matrices, input from different EvAMs, tests, etc (authored by Ramon Diaz-Uriarte, released under the GPL-3) has been previously used in Diaz-Uriarte and Vasallo, 2019 (all code available from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007246#sec019, files S1 Dataset and S2 Dataset) and Diaz-Colunga and Diaz-Uriarte, 2021 (repository: https://github.com/rdiaz02/what_genotype_next).


### ct-cbn ###

   This repository includes ct-​cbn-0.1.04b, from
   https://bsse.ethz.ch/cbg/software/ct-cbn.html, whose authors are Niko Beerenwinkel, Moritz Gerstung, and Seth Sullivant. It is released under the GNU GPL ("either version 2 of the License, or (at your option) any later version"). The GPL v3 can be combined with software under the AGPL v3, as used by the rest of this project.

   The code included in this repository is file ct-cbn-0.1.04b-RDU.tar.gz, a
   modification by RDU of the code in ct-​cbn-0.1.04b that includes: a minor bug fix (which, however, could be related to non-identifiability) and output with lambdas and likelihood from the initial run and each of the iterations.

   (For references about CBN see [References](#references)).


---
---
## Installing and running ##

EvAM-Tools is available as an R package, evamtools, that can launch a GUI as a Shiny web app.


### Overview of options

If you just want to run the Shiny app:
  * Go to http://iib.uam.es/evamtools .


You can install on your machines:
  * The [package and its dependencies](#how-to-install-the-r-package)
  * A [Docker image](#docker-images)

You can run on your machines:
  * The [R package, which includes the Shiny app](#how-to-run-the-r-package-and-the-shiny-app-locally-without-docker)
  * The [R package in an RStudio session from the Docker image](#how-to-run-the-r-package-from-the-docker-image)
  * The [Shiny app from a Docker image](#how-to-run-the-shiny-app-from-the-docker-image)

You can also build your own Docker image and you might want to run the Shiny app in a local intranet, possibly after modifying some settings. See the FAQ in the [EvAM-Tools: methods' details and FAQ](https://rdiaz02.github.io/EvAM-Tools/pdfs/evamtools_methods_details_faq.pdf) for details.


### How to install the R package

  * Install CBN
    - Use the file ct-cbn-0.1.04b-RDU.tar.gz.
    - Uncompress the directory. Then the usual configure, make dance. (Go inside the uncompressed directory and type `./configure`; when finished, type `make`).
    - Put the `h-cbn` and `ct-cbn` executables in your $PATH.

  * Install H-ESBCN
    - Clone the repository from  https://github.com/danro9685/HESBCN .
	- Go inside that directory, and type make.
	- When finished, add the `h-esbcn` executable to your $PATH.

  * Install OncoBN
    - OncoBN is available from https://github.com/phillipnicol/OncoBN
    - But this should work: start R, install the remotes package if you don't have it, and then issue `remotes::install_github("phillipnicol/OncoBN")`.

  *	Install HyperTraPS-CT
    - HyperTraPS-CT is available from https://github.com/StochasticBiology/hypertrapsct . <!-- We are using the bioconductor branch. Recall the hypertraps-ct repo used the bioconductor branch -->
	- This is what we do: start R, install the remotes package if you don't have it, and then issue `remotes::install_github("StochasticBiology/hypertraps-ct", ref = "bioconductor")`

  *	Install BML
    - An R wrapper for BML is available from https://github.com/Deschain/BML
	- This is what we do: start R, install the remotes package if you don't have it, and then issue `remotes::install_github("Deschain/BML")`

  * Install MC-CBN: this is optional.
    - Installing MC-CBN used to be complicated, because it required old versions of Boost (see https://github.com/cbg-ethz/MC-CBN/issues/5). This is no longer the case (see github issue) as of 2022-12-12, but MC-CBN is optional for EvAM-Tools. Of course, if you do not have MC-CBN, you will not be able to run MC-CBN. (MC-CBN is included in the methods available both from the web app and the Docker images).
    - Go to https://github.com/cbg-ethz/MC-CBN and follow the installation instructions: https://github.com/cbg-ethz/MC-CBN#installation
    <!-- - If you have all dependencies installed, and the correct version of libboost, then from R this should work: `install.packages("https://github.com/cbg-ethz/MC-CBN/releases/download/v2.1.0/mccbn_2.1.0.tar.gz", repos=NULL)` -->
    - Review the installation instructions and then install manually: https://github.com/cbg-ethz/MC-CBN#installation-from-source .
      (We suggest that at least for now you install manually, which is what we do in the Dockefiles, instead of doing `install.packages("https://github.com/cbg-ethz/MC-CBN/releases/download/v2.1.0/mccbn_2.1.0.tar.gz", repos=NULL)` because file `mccbn_2.1.0.tar.gz` is from December 2020, and thus it does not incorporate several bug fixes, including changes to the NAMESPACE).

  * Install the evamtools package
     - Make sure you have the required dependencies and imports, as listed in the DESCRIPTION file: igraph, OncoSimulR, stringr, Matrix, parallel, Oncotree , gtools , plot.matrix , DT, shinyjs, shiny, RhpcBLASctl, Rlinsolve, fastmatrix, graph, Rgraphviz, R.utils.
         - Note that we list, as imports, OncoBN, hypertracpst and BML. You need those (from above). Under Enhances we list mccbn: install it if you want to use it (see above).
     - Build (R CMD build evamtools) and install (R CMD INSTALL evamtools_x.y.z.tar.gz, with x.y.z replaced by the current version number). File `build-test.sh` builds, tests, and installs the package (and takes care of the version number).

	   Testing is, by default, parallelized and will use all CPUs except 1 (up to 20, the number of test files): the package includes over 1400 tests, with a test coverage of more than 90%. If you want to use fewer CPUs modify variable `TESTTHAT_CPUS` in script `build-test.sh` (see also https://testthat.r-lib.org/articles/parallel.html).


### Docker images
We provide two Docker images, one for running the Shiny app, and another with  RStudio to run the evamtools package directly.  They are available from
https://hub.docker.com/r/rdiaz02/evamshiny
and
https://hub.docker.com/r/rdiaz02/evamrstudio ; the first for running the Shiny app, the second for using the package from RStudio. Pull the one you need (`docker pull rdiaz02/evamshiny` or `docker pull rdiaz02/evamrstudio`). (When creating those, I first create a common one, with R and dependencies, evam-common-0, then one with the package  evam-common ---that depends on evam-common-0--; evamshiny and evamrstudio use evam-common, but you should be able to pull just the evamshiny or the evamrstudio images, independently.)

<!-- (Note: as of 2024-10-22, the docker image for using the package from RStudio is using R-4.4.1 and OncoSimulR 4.7.8, whereas the one for running the Shiny app is still using R-4.2.2 and OncoSimulR 4.0. I have updated the RStudio image because of internal lab needs ---access to the latest functionality in OncoSimulR--- but have not updated the Shiny app image since it is working and thus there is no need for updates). -->



Details about Docker are available here: https://docs.docker.com/get-docker/ .
Details about R with Docker and Rocker project here: https://www.rocker-project.org/ . Our images are based on the <!-- r-ver (https://hub.docker.com/r/rocker/r-ver) and --> rstudio (https://hub.docker.com/r/rocker/rstudio) Docker images from the Rocker project (https://www.rocker-project.org/).


More details about building and modifying the Docker images are provided in the FAQ available at [EvAM-Tools: method's details and FAQ](https://rdiaz02.github.io/EvAM-Tools/pdfs/evamtools_methods_details_faq.pdf) .


### How to run the R package and the shiny app locally without Docker

Install the [dependencies first](#how-to-install-the-r-package) and then the package. Once the package is installed, if you want to run the Shiny app open an R terminal and type

```
library(evamtools)
runShiny()
```
(If you do not want to run the Shiny app, do not issue `runShiny`).


### How to run the R package from the Docker image


Same as https://hub.docker.com/r/rocker/rstudio  (and see further options there):

```
docker run --rm -p 8787:8787 -e PASSWORD=yourpasswordhere rdiaz02/evamrstudio
```

Go to `localhost:8787` and log in. The password will the password you set above. The username:  If you are running rootless, the username has to be "root" (if you are running rootless, and as explained in https://rocker-project.org/images/versioned/rstudio.html , "the privileges of the root user in the container are already bounded by your regular user permissions in the host system."). Otherwise, the username should be  "rstudio". See https://rocker-project.org/images/versioned/rstudio.html for further options.

If you want to easily share data between your local file system and the Docker image you can do

```
docker run --rm -p 8787:8787 -v $HOME/tmp/rst:/home/rstudio -e PASSWORD=yourpasswordhere rdiaz02/evamrstudio
```
which will allow you to use your local `~/tmp/rst` to read from/write to the RStudio container (if using Windows you could write `-v /c/Users/someuser/somedirectory:/home/rstudio`); see additional documentation in https://www.rocker-project.org/use/shared_volumes/. See also https://docs.docker.com/storage/bind-mounts/.



Some additional notes:
   - If you get errors such as "docker: Error response from daemon: driver failed programming external connectivity on" you might want to restart the docker service.
   - We do not show Docker commands with `sudo`: it is possible, and generally preferable, to run docker without sudo; look a the Docker documentation in https://docs.docker.com/engine/security/rootless/ .

### How to run the Shiny app from the Docker image

```
docker run -d -p 4080:3000 --memory="2g" --name EVAM1 rdiaz02/evamshiny
```

This runs the `evamshiny` Docker image, mapping port 3000 of the container to port 4080 of the host  (so if you want to use the usual port 80, write 80 instead of 4080). You can use whatever you want instead of "EVAM1"; it is just a name to make other operations  (like stopping the container) simpler. In this example, we also limit the maximum memory to 2 GB.

This is a *non-interactive run*, and we use the "-d" or "--detach" option, so it runs in detached mode. You can point your browser to 0.0.0.0:4080 and the Shiny app should be there.  (If you want to keep the container running, you might want to add `tail -f /dev/null` to the above command).

(If you launch it this way, you can launch an arbitrary number of containers. For example, launch 15 different ones by specifying 15 different ports and 15 different names, and use HAProxy, https://www.haproxy.org/, to load-balance them ---you will want to use "sticky sessions").

---
---



## Main files and directories

### Dockerfiles ###

The Dockerfiles (Dockerfile-evam-shiny, Dockerfile-evam-rstudio) include all the information to create the containers with all dependencies.


### evamtools
The R package itself with standard organization. Directories and files under inst:
  * shiny-examples: code for the shiny application. The application consists on two main files: `server.R` (that controls the logic) and `ui.R` (includes all the interface). There are two additional directories: `assets` (with files for the landing page) and `test_shiny_app`.

  * miscell/Using_OncoSimulR_to_get_accessible_genotypes_trans_mats.tex: explanation of using OncoSimulR to check transition matrices for OT, CBN, OncoBN, and HESBCN, the equivalence of lambdas to terms in fitness expressions.
  * miscell/examples: examples referred to from other files (for example, from the former tex file).
  * miscell/tests-sample_genotypes_from_trm: output of tests that were run to verify the code for sampling genotypes from the transition rate matrices. <!-- We compared the output of our code with that from the code of the original authors (MHN, MCCBN) for a large set of cases. -->

	Note that the R package uses testthat to test our R code. Those tests will run automatically with the usual procedures from testthat or while doing `R CMD check`. <!-- For example, we check that transition rate matrices and transition probability matrices give identical results when compared to finding them via OncoSimulR (file test.OT-CBN-trans-mat-against-oncosimul.R and test.HESBCN-trans-mat-against-oncosimul.R) and against hand-computed examples (file test.trans-rates-f-graphs.R and test.HESBCN-transition-rate-matrices.R). --> The tests in evamtools/tests/testthat are separate from the tests under  inst/miscell/tests-sample_genotypes_from_trm



---
---



## References ##

### Overview paper on CPMs and EvAMs ###

- Diaz-Uriarte, R., & Johnston, I. G. (2025). A picture guide to cancer progression and evolutionary accumulation models: Systematic critique, plausible interpretations, and alternative uses. Ieee Access, 13, 62306–62340. https://doi.org/10.1109/ACCESS.2025.3558392



### OT ###


- Szabo, A., & Boucher, K. M. (2008). Oncogenetic Trees. In W. Tan, & L. Hanin
  (Eds.), Handbook of Cancer Models with Applications (pp. 1–24). : World
  Scientific. https://doi.org/10.1142/9789812779489_0001


- Desper, R., Jiang, F., Kallioniemi, O. P., Moch, H., Papadimitriou, C. H., &
  Sch\"affer, A A (1999). Inferring tree models for oncogenesis from comparative
  genome hybridization data. J Comput Biol, 6(1), 37–51.


### CBN and MCCBN ###

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

### H-ESBCN (PMCE) ###


- Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., &
  Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer
  evolution with high prognostic power. Bioinformatics, 38(3),
  754–762. http://dx.doi.org/10.1093/bioinformatics/btab717

   (About terminology: we will often refer to H-ESBCN, as that is the program we use, as shown here: https://github.com/danro9685/HESBCN. H-ESBCN is part of the PMCE procedure).

### OncoBN (DBN) ###

- Nicol, P. B., Coombes, K. R., Deaver, C., Chkrebtii, O., Paul, S., Toland,
  A. E., & Asiaee, A. (2021). Oncogenetic network estimation with disjunctive
  Bayesian networks. Computational and Systems Oncology, 1(2),
  1027. http://dx.doi.org/10.1002/cso2.1027


### HyperTraPS-CT ###

- Aga, O. N. L., Brun, M., Dauda, K. A., Diaz-Uriarte, R., Giannakis,
  K., & Johnston, I. G. (2024). Hypertraps-ct: Inference and prediction
  for accumulation pathways with flexible data and model structures.
  Plos Computational Biology, 20(9), e1012393.
  https://doi.org/10.1371/journal.pcbi.1012393



### BML ###

- Misra, N., Szczurek, E., & Vingron, M. (2014). Inferring the paths of somatic evolution in cancer. Bioinformatics (Oxford, England), 30(17), 2456–2463. https://doi.org/10.1093/bioinformatics/btu319


### Conditional prediction of genotypes and probabilities of paths from CPMs ###

- Hosseini, S., Diaz-Uriarte, Ramon, Markowetz, F., & Beerenwinkel,
  N. (2019). Estimating the predictability of cancer evolution. Bioinformatics,
  35(14), 389–397. http://dx.doi.org/10.1093/bioinformatics/btz332

- Diaz-Uriarte, R., & Vasallo, C. (2019). Every which way? On predicting tumor
  evolution using cancer progression models. PLOS Computational Biology, 15(8),
  1007246. http://dx.doi.org/10.1371/journal.pcbi.1007246

- Diaz-Colunga, J., & Diaz-Uriarte, R. (2021). Conditional prediction of
  consecutive tumor evolution using cancer progression models: What genotype
  comes next? PLOS Computational Biology, 17(12),
  1009055. http://dx.doi.org/10.1371/journal.pcbi.1009055



## Additional documentation
   Additional documentation is available from: https://rdiaz02.github.io/EvAM-Tools/

   The Bioinformatics paper "EvAM-Tools: tools for evolutionary accumulation and cancer progression models" is available at
  https://doi.org/10.1093/bioinformatics/btac710 .

   A video presentation of EvAM-Tools at the BCB Hub on 7 May 2025: https://youtu.be/WkqCfMFiJ3E


## Citing EvAM-Tools
If you use the package or the web app, please **cite** the Bioinformatics paper:

Diaz-Uriarte, R & Herrera-Nieto, P. 2022. EvAM-Tools: tools for evolutionary accumulation and cancer progression models. _Bioinformatics_.  https://doi.org/10.1093/bioinformatics/btac710 .



Ideally, also provide a link to the web app itself, https://iib.uam.es/evamtools , or the code repository, https://github.com/rdiaz02/EvAM-Tools.



   <!-- The preprint is available from bioRxiv at  https://doi.org/10.1101/2022.07.05.498481 . -->

   <!-- (or the direct bioRxiv link: https://www.biorxiv.org/content/10.1101/2022.07.05.498481v1 ) <!-- -\->, which includes additional examples ( https://www.biorxiv.org/content/10.1101/2022.07.05.498481v1.supplementary-material ) -->


## Funding
Supported by grant PID2019-111256RB-I00 funded by MCIN/AEI/10.13039/501100011033 and Comunidad de Madrid's PEJ-2019-AI/BMD-13961 to R. Diaz-Uriarte.
<p align="center">
<img src="micin-aei.png" alt="micin-aei logo" width="200"/>
</p>
