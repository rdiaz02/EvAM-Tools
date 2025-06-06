#!/bin/bash
rm -r -f ./evamtools/tests/testthat/_snaps
rm -r -f evamtools.Rcheck
V_PKG=$(cat ./evamtools/DESCRIPTION | grep Version | cut -d' ' -f2)
rm evamtools_$V_PKG.tar.gz
rm -r -f ./evamtools/inst/miscell/auto
R_ENVIRON_USER=./Renviron.based_on_bioc R --no-site-file --no-init-file CMD build evamtools
R_ENVIRON_USER=./Renviron.based_on_bioc R --no-site-file --no-init-file CMD INSTALL --install-tests evamtools_$V_PKG.tar.gz

## For parallel testing. All except 1
export TESTTHAT_CPUS=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | tail -n 1)
R_ENVIRON_USER=./Renviron.based_on_bioc R --no-site-file --no-init-file CMD check evamtools_$V_PKG.tar.gz
## Use also R-devel
R_ENVIRON_USER=./Renviron.based_on_bioc R-devel --no-site-file --no-init-file CMD check evamtools_$V_PKG.tar.gz
