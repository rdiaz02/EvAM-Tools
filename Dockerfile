FROM rocker/verse 

RUN apt-get update && \
    apt-get install -f -y software-properties-common && \
    rm -rf /var/lib/apt/lists/* && \
    add-apt-repository universe && \
    apt-get update && \
    apt-get install -y libboost-all-dev && \
    apt-get install -y libgsl-dev 

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('TRONCO')"
RUN R -e "BiocManager::install('OncoSimulR')"
RUN R -e "library(devtools);devtools::install_github('phillipnicol/OncoBN')"
RUN R -e "install.packages('relations')"
RUN R -e "install.packages('https://github.com/cbg-ethz/MC-CBN/releases/download/v2.1.0/mccbn_2.1.0.tar.gz', repos=NULL)"
RUN R -e "install.packages('RhpcBLASctl')"
RUN R -e "install.packages('igraph')"
RUN R -e "install.packages('inline')"
RUN R -e "install.packages('Oncotree')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('pryr')"
RUN R -e "install.packages('readr')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('foreach')"
RUN R -e "install.packages('stringr')"
RUN R -e "install.packages('testthat')"
RUN R -e "install.packages('plot.matrix')"
RUN R -e "install.packages('httr')"
RUN R -e "install.packages('openssl')"
RUN R -e "install.packages('xml2')"
RUN R -e "install.packages('usethids')"
RUN R -e "install.packages('credentials')"
RUN R -e "install.packages('roxygen2')"
RUN R -e "install.packages('optparse')"
RUN R -e "install.packages('rversions')"
RUN R -e "install.packages('imager')"
RUN R -e "install.packages('DT')"
RUN R -e "install.packages('shinyjs')"
RUN R -e "install.packages('markdown')"

COPY . /app/

RUN mkdir -p /app/Sources/

#Install cbn
RUN cd /app/guloMAM/R/ && \  
    cp ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood.tar.gz /app/Sources && \
    cd /app/Sources && \
    tar -xvzf ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood.tar.gz && \
    cd ct-cbn-0.1.04b-with-rdu-bug-fix-write-lambda-likelihood && \
    ./configure && \
    make install && \
    cp -t /usr/local/bin/ ct-cbn h-cbn 

#Install HESBCN
RUN cd /app/Sources && \
    git clone https://github.com/danro9685/HESBCN && \
    cd HESBCN && \ 
    sed -i "s/gcc-5/gcc-9/g" Makefile && \
    sed -i  "s/^LDLIBS = $/LDLIBS = -lgsl -lm -lgslcblas/g" Makefile && \ 
    make && \
    cp h-esbcn /usr/local/bin/ 

#Install HyperTraPS

#Create folder to share with the outside
RUN mkdir /app/outside

WORKDIR /app/guloMAM

# CMD ["Rscript", "insts/miscell/examples/simple_test.R"]
# CMD ["Rscript", "docker/runFromImage.R"]
ENTRYPOINT  ["../docker/runFromImage.R", "-f"]
CMD ["NULL"]

EXPOSE 3000