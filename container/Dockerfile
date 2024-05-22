FROM ghcr.io/nextflow-io/training:latest

USER root

## installing R from Scratch

RUN gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -
RUN apt update
RUN apt install r-base r-base-dev -y

## installing RStudio

RUN apt-get update && \
    apt-get install -y gdebi-core && \
    wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.09.1-494-amd64.deb && \
    gdebi -n rstudio-server-2023.09.1-494-amd64.deb && \
    rm rstudio-server-2023.09.1-494-amd64.deb

RUN chmod -R 777 /usr/local/lib/R/site-library

RUN apt-get install libharfbuzz-dev libfribidi-dev -y

USER gitpod

RUN echo "gitpod:pass" | sudo chpasswd

WORKDIR /home/gitpod

RUN Rscript -e "install.packages(c(\
'BiocManager', \
'kableExtra', \
'remotes', \
'ggnewscale', \
'tidyverse'\
), \
repos = 'https://cloud.r-project.org', \
dependencies = TRUE\
)"

RUN Rscript -e "remotes::install_github('YuLab-SMU/ggtree')"

RUN Rscript -e "BiocManager::install(c(\
'Gviz', \
'VariantAnnotation', \
'GenomicFeatures', \
'rtracklayer', \
'Biostrings', \
'ggtree', \
'DESeq2', \
'tximport', \
'tximeta', \
'msa', \
'seqinr', \
'fgsea', \
'GOSemSim', \
'qvalue', \
'reshape2', \
'downloader', \
'enrichplot', \
'plotly'\
))"


COPY install_packages.R install_packages.R
RUN wget https://www.bioconductor.org/packages/3.12/data/annotation/src/contrib/org.Hs.eg.db_3.12.0.tar.gz
RUN wget https://bioconductor.statistik.tu-dortmund.de/packages/3.12/bioc/src/contrib/AnnotationDbi_1.52.0.tar.gz
RUN wget https://bioconductor.statistik.tu-dortmund.de/packages/3.12/bioc/src/contrib/clusterProfiler_3.18.1.tar.gz
RUN wget https://bioconductor.org/packages/3.12/data/annotation/src/contrib/DO.db_2.9.tar.gz
RUN wget https://bioconductor.statistik.tu-dortmund.de/packages/3.12/bioc/src/contrib/DOSE_3.16.0.tar.gz
RUN wget https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz
RUN Rscript install_packages.R
RUN rm *.tar.gz
RUN rm install_packages.R
