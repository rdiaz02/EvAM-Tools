# License and copyright

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
* HyperTraps: https://github.com/sgreenbury/HyperTraPS. Follow the instruction there to install the package through conda 
and get a working enviroment with HyperTraPS installed. Add both the `bin` and the `src/python` to your `$PATH`. Also, add the `bin` folder to your $PATH. 
* HESBCN: https://github.com/danro9685/HESBCN. For proper compilation you should modify line 23 in *Makefile* to include *LDLIBS = -lgsl -lm -lgslcblas*. You can also change compilation from gcc-5 to something a bit more up to date, like gcc-10. Finally, add this folder to your $PATH.

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
      (or links to them) in you path. Uncompress and the usual dance
      (./configure, make)


# Other optional code
  - MC-CBN: https://github.com/cbg-ethz/MC-CBN
      - For now, the default is not to use MC-CBN. If you want to use it, you'll
      need to install the MC-CBN package. Note that you need older versions of
      libboost: https://github.com/cbg-ethz/MC-CBN/issues/5
