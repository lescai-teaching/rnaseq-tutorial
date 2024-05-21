FROM ghcr.io/lescai-teaching/rstudio-docker:latest as build

FROM ghcr.io/nextflow-io/training:latest as main

RUN /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
RUN brew install R

ENV NB_USER=gitpod

COPY --from=build /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=build /usr/local/lib/R/library /usr/local/lib/R/library

USER root

RUN apt-get update && \
    apt-get install -y gdebi-core && \
    wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.09.1-494-amd64.deb && \
    gdebi -n rstudio-server-2023.09.1-494-amd64.deb && \
    rm rstudio-server-2023.09.1-494-amd64.deb

COPY .Rprofile /opt/conda/lib/R/.Rprofile

USER ${NB_USER}

COPY .Rprofile /home/${NB_USER}/.Rprofile

RUN echo "${USER}:pass" | sudo chpasswd

WORKDIR /home/${NB_USER}
