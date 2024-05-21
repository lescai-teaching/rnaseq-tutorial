FROM ghcr.io/lescai-teaching/rstudio-docker:latest as build

FROM ghcr.io/nextflow-io/training:latest as main

RUN echo "\n" | /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
RUN eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
RUN /home/linuxbrew/.linuxbrew/bin/brew install R

COPY --from=build /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=build /usr/local/lib/R/library /usr/local/lib/R/library

USER root

RUN apt-get update && \
    apt-get install -y gdebi-core && \
    wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.09.1-494-amd64.deb && \
    gdebi -n rstudio-server-2023.09.1-494-amd64.deb && \
    rm rstudio-server-2023.09.1-494-amd64.deb

COPY .Rprofile /opt/conda/lib/R/.Rprofile

USER gitpod

COPY .Rprofile /home/gitpod/.Rprofile

RUN echo "gitpod:pass" | sudo chpasswd

WORKDIR /home/gitpod
