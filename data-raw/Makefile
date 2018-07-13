default: all
RDIR = .

all: data-proc/births.csv data-proc/england_demographics.csv data-proc/incidence.csv data-proc/mortality_rates.csv data-proc/sputum_by_age.csv data-proc/latent_parameters.csv data-proc/contacts.csv data-proc/eff_bcg.csv data-proc/disease_params_table.csv data-proc/demographic_params_table.csv

## Extract and format births in England
data-proc/births.csv: births_in_england.Rmd
			Rscript -e 'rmarkdown::render("births_in_england.Rmd", c("html_document"), output_dir = "docs", output_options = list())'



## Extract and format england population demographics by modelled age groups
data-proc/england_demographics.csv: birth_status_demographics.Rmd
			Rscript -e 'rmarkdown::render("birth_status_demographics.Rmd", c("html_document"), output_dir = "docs", output_options = list())'



## Extract and format ETS data into incidence data
data-proc/incidence.csv: clean_ets.Rmd
			Rscript -e 'rmarkdown::render("clean_ets.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


## Make parameter definitions into a table for the disease model
data-proc/disease_params_table.csv: disease_params_table.Rmd
			Rscript -e 'rmarkdown::render("disease_params_table.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


## Make parameter definitions into a table for the demographic model
data-proc/demographic_params_table.csv: demographic_params_table.Rmd
			Rscript -e 'rmarkdown::render("demographic_params_table.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


## Make sources definitions into a table
data-proc/sources_table.csv: sources_table.Rmd
			Rscript -e 'rmarkdown::render("sources_table.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


## Make mortality rates in England
data-proc/mortality_rates.csv: mortality_rate_england.Rmd
			Rscript -e 'rmarkdown::render("mortality_rate_england.Rmd", c("html_document"), output_dir = "docs", output_options = list())'

## Make ETS priors
data-proc/sputum_by_age.csv: ets_priors.Rmd
			Rscript -e 'rmarkdown::render("ets_priors.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


## Make Latent parameters distributions
data-proc/latent_parameters.csv: norm_dist_for_latent_params.Rmd
			Rscript -e 'rmarkdown::render("norm_dist_for_latent_params.Rmd", c("html_document"), output_dir = "docs", output_options = list())'

## Make contact matrix
data-proc/contacts.csv: generate_contact_matrix.Rmd
			Rscript -e 'rmarkdown::render("generate_contact_matrix.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


## Make BCG effectiveness
data-proc/eff_bcg.csv: BCG_eff_active_tb.Rmd
			Rscript -e 'rmarkdown::render("BCG_eff_active_tb.Rmd", c("html_document"), output_dir = "docs", output_options = list())'


