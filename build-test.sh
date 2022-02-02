# export R_LIBS_SITE="/usr/local/lib/R-devel/lib/R/library:~/R/x86_64-pc-linux-gnu-library/3.7"

rm -r evamtools_0.0.0.93.tar.gz
rm -revamtools.Rcheck

R_ENVIRON_USER=~/.Renviron.bioc R CMD build evamtools

R_ENVIRON_USER=~/.Renviron.bioc R CMD check evamtools_0.0.0.93.tar.gz

R_ENVIRON_USER=~/.Renviron.bioc R CMD INSTALL evamtools_0.0.0.93.tar.gz

## FIXME:
##    FIXME:   R version should be at least 4.1.2
##    FIXME:   but do we really need to export R_LIBS_SITE?
