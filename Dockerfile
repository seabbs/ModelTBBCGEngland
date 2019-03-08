## Start with the tidyverse docker image
FROM seabbs/tidyverse-gpu

MAINTAINER "Sam Abbott" sam.abbott@bristol.ac.uk

ADD . /home/rstudio/ModelTBBCGEngland

WORKDIR /home/rstudio

## Get thrust
RUN wget https://github.com/thrust/thrust/archive/1.9.4.zip \
    && unzip 1.9.4.zip \
    && mv thrust-1.9.4 thrust \
    && mv thrust /usr/local/include

## Set env for Libbi
ENV PERL_MM_USE_DEFAULT=1

## Get LibBi (and deps)
RUN git clone https://github.com/lawmurray/LibBi.git \
  && cd LibBi \
  && sudo apt-get install -y \
    libblas-dev \
    liblapack-dev \
    libqrupdate-dev \
    libboost-all-dev \
    libgsl0-dev \
    libnetcdf-dev \
    autoconf \
    automake \
  && sudo cpan .

### Get R package libraries
RUN apt-get install -y \
     libnetcdf-dev \
     && apt-get clean

WORKDIR /home/rstudio/ModelTBBCGEngland

RUN Rscript -e 'devtools::install_dev_deps()'
