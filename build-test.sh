export R_LIBS_SITE="/usr/local/lib/R-devel/lib/R/library:~/R/x86_64-pc-linux-gnu-library/3.7"

rm -r guloMAM_0.0.0.9015.tar.gz
rm -r guloMAM.Rcheck

R_ENVIRON_USER=~/.Renviron.bioc R CMD build guloMAM

R_ENVIRON_USER=~/.Renviron.bioc R CMD check guloMAM_0.0.0.9015.tar.gz
