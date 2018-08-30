library('rbi')
library('rbi.helpers')
library('ModelTBBCGEngland')

## Should particles be adapted
adapt_part <- FALSE
adapt_prop <- TRUE
sample_post <- TRUE
verbose <- TRUE
## Need to preload input
input <- setup_model_input(run_time = 73, time_scale_numeric = 1)

model_file <- system.file(package="ModelTBBCGEngland", "bi/BaseLineModel.bi") # get full file name from package
SIRmodel <- bi_model(model_file) # load model

if (adapt_part) {
  SIRmodel <- everything_from_model(SIRmodel)
}

SIRdata <- bi_generate_dataset(SIRmodel, end_time=73, noutputs=6, seed=12345678, input = input)

bi_prior <- sample(proposal="prior", SIRmodel, nsamples = 1000, end_time = 73,
                   nparticles = 8, obs = SIRdata, 
                   input = input, seed=1234,
                   options = list(with="transform-initial-to-param"), verbose = verbose)

if (adapt_part) {
adapted <- adapt_particles(bi_prior, min = 4, max = 16)

adapted$options$nparticles
}else{
  adapted <- bi_prior
}

if (adapt_prop) {
  adapted <- adapt_proposal(adapted, min=0.05, max=0.4, adapt = "both")
  
  get_block(adapted$model, "proposal_parameter")
}

if (sample_post) {
  posterior <- sample(adapted, nsamples=1000, sample_obs=TRUE)
  
  plot(posterior)
}
