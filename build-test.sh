rm -r -f evamtools.Rcheck
V_PKG=$(cat ./evamtools/DESCRIPTION | grep Version | cut -d' ' -f2)
rm evamtools_$V_PKG.tar.gz
R_ENVIRON_USER=./Renviron.based_on_bioc R CMD build evamtools
R_ENVIRON_USER=./Renviron.based_on_bioc R CMD check evamtools_$V_PKG.tar.gz
R_ENVIRON_USER=./Renviron.based_on_bioc R CMD INSTALL evamtools_$V_PKG.tar.gz

