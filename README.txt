* example directory

** origin of the files  
  - I copy the code from "what_genotype_next": https://github.com/rdiaz02/what_genotype_next
    - Directories:
      - what_genotype_next/run-biol-examples is examples
      - MHN under run-biol-examples is, well, MHN under examples
  - So, yes, cp. No symlinks.
    - Recall in that code there were symlinks:
      - what_genotype_next/run-biol-examples/schill-trans-mat.R
      - what_genotype_next/run-biol-examples/MHN
    - Recall you need to download the Schill code from  https://github.com/RudiSchill/MHN.
      - I've added it here after downloading.
  - Sequence of cp done:
    - cd ./examples
    - cp -a ~/tmp/MHN .
    - rm -r -f ./MHN/.git
    - cp ../../what_genotype_next/run-biol-examples/*.R .

    
