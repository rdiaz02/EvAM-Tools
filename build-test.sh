#!/bin/bash
rm -r -f evamtools.Rcheck
V_PKG=$(cat ./evamtools/DESCRIPTION | grep Version | cut -d' ' -f2)
rm evamtools_$V_PKG.tar.gz
rm -r -f ./evamtools/inst/miscell/auto
R_ENVIRON_USER=./Renviron.based_on_bioc R --no-site-file --no-init-file CMD build evamtools
R_ENVIRON_USER=./Renviron.based_on_bioc R --no-site-file --no-init-file CMD INSTALL --install-tests evamtools_$V_PKG.tar.gz
R_ENVIRON_USER=./Renviron.based_on_bioc R --no-site-file --no-init-file CMD check evamtools_$V_PKG.tar.gz

