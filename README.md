
Model Tuberculosis, including BCG vaccination, in England
=========================================================

`ModelTBBCGEngland` is an R package that contains models, summarised data, and fitting scripts for a model of TB, including BCG vaccination, fit to routine surveillance data in England.

Installation
------------

You can install `ModelTBBCGEngland` directly from github with:

``` r
# install.packages("devtools")
## Set github username
devtools::install_github("seabbs/ModelTBBCGEngland", username = 'seabbs')
```

### Using docker

This analysis was developed in a docker container based on the [tidyverse](https://hub.docker.com/r/rocker/tidyverse/) docker image. To run the docker image run:

``` bash
docker run -d -p 8787:8787 -e USER=ModelTBBCGEngland -e PASSWORD=ModelTBBCGEngland \
    --mount type=bind,source=$(pwd)/results/modeltbbcgengland,target=/home/rstudio/ModelTBBCGEngland/vignettes/results \
    --mount type=bind,source=$(pwd)/data/tb_data,target=/home/rstudio/ModelTBBCGEngland/data-raw/tb_data \
    --name modeltransvsdirect seabbs/modeltbbcgengland
```

This command expects both a data and results folder to be available to be mounted. Create them using the following (Mac/Linux).

``` bash
mkdir -p results/modeltbbcgengland
mkdir -p data/tb_data
```

The rstudio client can be found on port `:8787` at your local machines ip. The default username:password is seabbs:seabbs, set the user with `-e USER=username`, and the password with `- e PASSWORD=newpasswordhere`.

To run a plain R terminal use:

``` bash
docker run --rm -it --user seabbs modeltbbcgengland /usr/bin/R
```

To run a plain bash session:

``` bash
docker run --rm -it --user seabbs modeltbbcgengland /bin/bash
```
