
# Model Tuberculosis, including BCG vaccination, in England

`ModelTBBCGEngland` is an R package that contains models, summarised
data, and fitting scripts for a model of TB, including BCG vaccination,
fit to routine surveillance data in England.

## Installation

You can install `ModelTBBCGEngland` directly from github with:

``` r
# install.packages("devtools")
devtools::install_github("seabbs/ModelTBBCGEngland")
```

### Model development

For a discussion of the model development process please see
[here](https://www.samabbott.co.uk/thesis/8-model-development.html). For
a discussion of the model fitting pipeline please see
[here](https://www.samabbott.co.uk/thesis/9-model-fitting.html).

### Model fitting

Run the following code in the terminal after installing and building the
package and opening the Rstudio project. The results (and logs) will be
found at
`./vignettes/results`.

``` bash
nohup Rscript inst/scripts/run-scenarios.R --sample_priors --fit --reports
```

To fit only the main scenario add `--scenario baseline`. If interested
in rapid exploration or if using limited compute resources experiment
with the `--calib_run` flag which uses point estimates for the initial
conditions and process noise - reducing the number of particles
required. Alternatively explore the other options available using
`--help`. An example of a fitting run used for exploration is
below.

``` bash
nohup Rscript inst/scripts/run-scenarios.R --sample_priors --fit --reports --cores 16 --calib_run --scenario baseline
```

For more interactive testing see `./tests/manual-tests/test-rbi-fit.R`.
This script was used during pipeline development.

### Using docker

This analysis was developed in a docker container based on the
[tidyverse](https://hub.docker.com/r/rocker/tidyverse/) docker image. To
run the docker image
run:

``` bash
docker run -d -p 8787:8787 -e USER=ModelTBBCGEngland -e PASSWORD=ModelTBBCGEngland \
    --mount type=bind,source=$(pwd)/results/modeltbbcgengland,target=/home/rstudio/ModelTBBCGEngland/vignettes/results \
    --mount type=bind,source=$(pwd)/data/tb_data,target=/home/rstudio/ModelTBBCGEngland/data-raw/tb_data \
    --name modeltransvsdirect seabbs/modeltbbcgengland
```

This command expects both a data and results folder to be available to
be mounted. Create them using the following (Mac/Linux).

``` bash
mkdir -p results/modeltbbcgengland
mkdir -p data/tb_data
```

The rstudio client can be found on port `:8787` at your local machines
ip. The default username:password is seabbs:seabbs, set the user with
`-e USER=username`, and the password with `- e
PASSWORD=newpasswordhere`.

To run a plain R terminal use:

``` bash
docker run --rm -it --user seabbs modeltbbcgengland /usr/bin/R
```

To run a plain bash session:

``` bash
docker run --rm -it --user seabbs modeltbbcgengland /bin/bash
```
