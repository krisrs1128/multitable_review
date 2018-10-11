FROM r-base:3.5.1

RUN apt-get update
RUN apt-get install -y libssl-dev
RUN apt-get install -y libcurl4-gnutls-dev
RUN apt-get install -y libxml2-dev

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('PMA')"
RUN Rscript -e "install.packages('ade4')"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('expm')"
RUN Rscript -e "install.packages('forcats')"
RUN Rscript -e "install.packages('ggrepel')"
RUN Rscript -e "install.packages('glasso')"
RUN Rscript -e "install.packages('glmnet')"
RUN Rscript -e "install.packages('plot3D')"
RUN Rscript -e "install.packages('reshape2')"
RUN Rscript -e "install.packages('rstan')"
RUN Rscript -e "install.packages('spls')"
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('vegan')"
RUN Rscript -e "install.packages('viridis')"
RUN Rscript -e "devtools::install_github('krisrs1128/gflasso')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('phyloseq')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('DESeq2')"

COPY data/sim.rda /home/data/sim.rda
COPY dimension_red /home/dimension_red/
COPY supervised /home/supervised/
COPY run_all.sh /home/run_all.sh