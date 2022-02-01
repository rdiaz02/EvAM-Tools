# License and copyright



FIXME: add this as file COPYRIGHTS under /inst
# Copyright and origin of files under evamtools/R

- All files under evamtools/R are copyright Pablo Herrera Nieto and Ramon Diaz-Uriarte, except for the following:

- File HESBCN__import.hesbcn.R: 
   This file contains function import.hesbcn.
   
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
  

# FIXME: need to update this file

- Add a file of copyrights (e.g., like in
  https://cran.r-project.org/web/packages/RcppEigen/COPYRIGHTS), for all that we
  have not authored, that lists
    - original location of the file, date we copied it, committ (if it exists),
      authors, copyright and license, and location in our repo
	  
- Add that information here too, at the bottom	

- In text above, wherever it is needed, clarify that original license (e.g., Apache 2, or
  whatever) is compatible with our GPL v3 license.

    
## FIXME : changes in other files
- In files we include from other authors, if it is not already there, add a
  header in the file itself that lists: original location of the file, date we
  copied it, committ (if it exists), authors, copyright and original license.
    - This would apply to files from MHN, for example, if we include them.
	
- DESCRIPTION of R package
    - remember to update author list if needed
	- remember to add COPYRIGHTS file, if needed (see above) or symlink it.Yes,
      it is the same file, but we want to make sure that if only the R package
      subdirectory is grabbed, that file is present.




All original code here is released under the GPL v3 license.


This repository includes code by:

- Niko Beerenwinkel, Moritz Gerstung, and Seth Sullivant: the file
   ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood.tar.gz, a
   modification of the ct-cbn code from
   https://bsse.ethz.ch/cbg/software/ct-cbn.html

- Schill, R., Solbrig, S., Wettig, T., & Spang, R: the files under MHN,
  downloaded from https://github.com/RudiSchill/MHN

# How to install?
In order to run all CPMs you will need to clone their repositories:

include github repo, and instructions to get them running

* DBN: https://github.com/phillipnicol/OncoBN.

FIXME: commit hash and date when we got the package.

* HyperTraps: https://github.com/sgreenbury/HyperTraPS. Follow the instruction there to install the package through conda 
and get a working enviroment with HyperTraPS installed. Add both the `bin` and the `src/python` to your `$PATH`. Also, add the `bin` folder to your $PATH. It is also convenient to run `echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc`.

* HESBCN: https://github.com/danro9685/HESBCN . For proper compilation you should modify line 23 in *Makefile* to include *LDLIBS = -lgsl -lm -lgslcblas*. You can also change compilation from gcc-5 to something a bit more up to date, like gcc-10. Finally, add the h-esbcn executable to your $PATH.


# Directory code_from_what_genotype_next 

## Origin of the files  
  - I copy the code from "what_genotype_next": https://github.com/rdiaz02/what_genotype_next
    - Directories:
      - what_genotype_next/run-biol-examples is examples
      - MHN under run-biol-examples is, well, MHN under examples
  - So, yes, cp. No symlinks.
    - Recall in that code there were symlinks:
      - what_genotype_next/run-biol-examples/schill-trans-mat.R
      - what_genotype_next/run-biol-examples/MHN
    - The directory MHN contains the code from
    https://github.com/RudiSchill/MHN.
    
      - I've added it here after downloading.
      - There is no author nor license information for that code. The code has
      been used in the paper

      Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer
      progression using Mutual Hazard Networks. Bioinformatics, 36(1),
      241–249. https://doi.org/10.1093/bioinformatics/btz513

      so I assume their authors to be the authors of the paper.
      
  - Sequence of cp done:
    - cd ./what_genotype_next
    - cp -a ~/tmp/MHN .
    - rm -r -f ./MHN/.git
    - cp ../../what_genotype_next/run-biol-examples/*.R .

  - I add the file
    ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood.tar.gz
    - This has the ct-cbn code from
    https://bsse.ethz.ch/cbg/software/ct-cbn.html
    plus a couple of bug fixes by myself (RDU) plus outputing the final lambda.
    - The authors and copyright holders of this code are Niko Beerenwinkel,
    Moritz Gerstung, and Seth Sullivant.
    - Their code is released under the GNU GPL.

    - You need to install that code and place the h-cbn and ct-cbn executables
      (or links to them) in you path.  Uncompress and the usual dance
      (./configure, make). Beware: if you have stale files around (from runs of make with other versions of gcc, for example), make might not work; the simplest fix is to checkout a new copy from the repo, extract, and repeat the ./configure, make, dance again.


# Other optional code
  - MC-CBN: https://github.com/cbg-ethz/MC-CBN
      - For now, the default is not to use MC-CBN. If you want to use it, you'll
      need to install the MC-CBN package. Note that you need older versions of
      libboost: https://github.com/cbg-ethz/MC-CBN/issues/5
	  
	  - Note also that to run MC-CBN you need to unset the environment variable `_R_CHECK_LENGTH_1_LOGIC2_` , because of a problem explained here:
	   https://github.com/cbg-ethz/MC-CBN/issues/9

	  

# Tests
FIXME write this better

We use testthat to test our R code. Those tests will run automatically with the usual procedures from testthat or while doing `R CMD check`. For example, we check that transition rate matrices and transition probability matrices from OT and CBN give identical results when compared to finding them via OncoSimulR (file test.OT-CBN-trans-mat-against-oncosimul.R) and against hand-computed examples (file test.trans-rates-f-graphs.R).


In addition, under evamtools/inst/miscell/tests-sample_genotypes_from_trm we have the output of tests that were run to verify the code for sampling genotypes from the transition rate matrices. We compared the output of our code with that from the code of the original authors (MHN, MCCBN) for a large set of cases.


