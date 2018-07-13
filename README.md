
Model Tuberculosis, including BCG vaccination, in England
=========================================================

`ModelTBBCGEngland` is an R package that contains models, summarised data, and fitting scripts for a model of TB, including BCG vaccination, fit to routine surveillance data in England.

Installation
------------

You can install `ModelTBBCGEngland` directly from github with:

``` r
# install.packages("devtools")
## Set github username
## Apply for PAT code and store as GITHUB_PAT in .Reviron
devtools::install_github("seabbs/ModelTBBCGEngland", username = 'seabbs')
```

### Using docker

This analysis was developed in a docker container based on the [tidyverse](https://hub.docker.com/r/rocker/tidyverse/) docker image. To run the docker image run:

``` bash
docker run -d -p 8787:8787 -e USER=ModelTBBCGEngland -e PASSWORD=ModelTBBCGEngland --name modeltransvsdirect seabbs/modeltbbcgengland
```

The rstudio client can be found on port `:8787` at your local machines ip. The default username:password is seabbs:seabbs, set the user with `-e USER=username`, and the password with `- e PASSWORD=newpasswordhere`. If the raw data required for this analysis is available on your computer, and regenerate of the model data is required, then use the following to mount it into the docker container `--mount type=bind,source=$(pwd)/data/tb_data,target=/home/rstudio/ModelTBBCGEngland/data-raw/tb_data`.

To run a plain R terminal use:

``` bash
docker run --rm -it --user seabbs modeltbbcgengland /usr/bin/R
```

To run a plain bash session:

``` bash
docker run --rm -it --user seabbs modeltbbcgengland /bin/bash
```
