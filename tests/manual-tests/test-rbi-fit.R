library('rbi')
library('rbi.helpers')


model_file <- system.file(package="ModelTBBCGEngland", "bi/BaseLineModel.bi") # get full file name from package
SIRmodel <- bi_model(model_file) # load model

SIRmodel <- rbi::replace_all(SIRmodel, "\\(has_output = 0, has_input = 0\\)", "")
SIRmodel <- rbi::replace_all(SIRmodel, "\\(has_input = 0, has_output = 0\\)", "")
SIRmodel <- rbi::replace_all(SIRmodel, "\\(has_output = 0\\)", "")

SIRdata <- bi_generate_dataset(SIRmodel, end_time=73, noutputs=6, seed=12345678, input = input)

bi_prior <- sample(proposal="prior", SIRmodel, nsamples=1000, end_time=73, nparticles=4, obs=SIRdata, 
                   input = input, seed=1234,
                   options = list(with="transform-initial-to-param"))

adapted <- adapt_particles(bi_prior, min = 4, max = 128)

adapted$options$nparticles

adapted <- adapt_proposal(adapted, min=0.05, max=0.4, adapt = "both")

get_block(adapted$model, "proposal_parameter")


posterior <- sample(adapted, nsamples=500, sample_obs=TRUE)

plot(posterior)