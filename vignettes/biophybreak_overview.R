## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  out.width = "100%",
  comment = "#>"
)

## ----packages-----------------------------------------------------------------
library(biophybreak)
library(ggplot2)

seed <- 4321
set.seed(seed)

## ----simulate_data------------------------------------------------------------
#number of individuals
nInd <- 5

#time betweeen infection and sampling is gamma distributed with mean 2 and sd 1.5
infection.ages <- rgamma(nInd, shape = (16/9), rate = (8/9))

#type of distribution for generation function
#(relative chance of transmission based on how long an individual has been infected)
gen.dist.type <- "function"
gen.dist <- stepfun(x = c(0, 0.5), y = c(0, 3, 1)) #three times more likely to transmit in first 6 months

#make a transmission tree
tt <- rtrans_tree(nInd = nInd, 
                  infection.ages = infection.ages, 
                  gen.dist.type = "function", 
                  gen.dist = stepfun(x = c(0, 0.5), y = c(0, 3, 1))) 

#within host model to use
a <- 5 #effective population size
b <- 5 #effective population growth rate per generation (1.5 days)

#mutation model to use
mut.model = "pol" #can be "pol" or "env"

#generate phylogeny and create phybreakdata object
#sim_data <- sim.coal.phybreak(tt = tt, 
#                              a = rep(a, nInd), 
#                              b = rep(b, nInd), 
#                              mut_model = mut.model) 

#load pregenerated sequences for this example
data("example_seqs4321")

#generate phylogeny and create phybreakdata object
sim_data <- sim.coal.phybreak(tt = tt, 
                              a = rep(a, nInd), 
                              b = rep(b, nInd), 
                              user_seqs = example_seqs4321$seqs,
                              user_seq_hosts = example_seqs4321$seq_hosts,
                              user_seq_names = example_seqs4321$seq_names)

#plot simulated transmission history and phylogeny
plotPhyloTrans(sim_data)

#load parameters for multiple biomarker model
data("MBM_pars")

#simluate biomarker values from infection ages
biomarkers <- sim.biomarkers(t.inf = infection.ages, 
                             mub = MBM_pars$mub, #means of MBM parameters
                             Sigmab = MBM_pars$Sigmab, #covariance matrix of MBM parameters
                             sigmae = MBM_pars$sigmae) #standard deviations of gaussian noise for each biomarker

## -----------------------------------------------------------------------------
infection.dists <- mbm.predict(nInds = 5,
                               BED = biomarkers[,1], #values from BED assay (OD-n))
                               LAg = biomarkers[,2], #value from LAg-Avidity assay
                               CD4 = biomarkers[,3], #CD$+ T cell count
                               pol = biomarkers[,4], #fraction of polymorphic sites in pol gene
                               pol2 = biomarkers[,5], #NGS diversity in pol gene
                               mub = MBM_pars$mub, 
                               Sigmab = MBM_pars$Sigmab, 
                               sigmae = MBM_pars$sigmae,
                               n.adapt = 10000, n.burn = 50000, n.iter = 100000, output.raw = TRUE)

## -----------------------------------------------------------------------------
data <- phybreakdata(sequences = sim_data$sequences, #sequences
                     sample.times = sim_data$sample.times, #times of samples
                     host.names = sim_data$sample.hosts) #which host each sample was taken from

## -----------------------------------------------------------------------------
#set up phybreak object to run MCMC
MCMCstate <- phybreak(dataset = data, 
                      wh.model = "linear", #linear effective population size growth
                      wh.bottleneck = "wide", #wide transmission bottleneck
                      prior.wh.level.mean = 0.1, #pretty conservative; corresponds to a mean of about 24 for a
                      est.sample.mean = FALSE, #do not estimate the sample mean (because we are using the nonparametric distributions)
                      sample.nonpar = TRUE, #we are using the nonparametric distributions from the MBM
                      sample.density = infection.dists$pdf.num[order(tt$t_sam)], #distrubutions from MBM, ordered by sampling time
                      gen.type = gen.dist.type, #kind of generation function
                      gen.density = gen.dist, #generation function
                      post.sam.trans.rate = rep(1, nInd), #no penalty for transmission after sampling
                      est.gen.mean = FALSE, #do not estimate gen mean (because function is provided)
                      use.pml = TRUE, #use the pml function in phangorn for the sequence likelihood with GTR model
                      mut.model = mut.model)

## -----------------------------------------------------------------------------
MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 10000)

## -----------------------------------------------------------------------------
MCMCstate <- sample_phybreak(MCMCstate, nsample = 25000)

## ----effective_size-----------------------------------------------------------
ESS <- ESS(MCMCstate)
barplot(log10(ESS))
abline(h = log10(200))

## ----mpc_tree-----------------------------------------------------------------
plotPhyloTrans(MCMCstate, plot.which = "mpc")

## -----------------------------------------------------------------------------
#since we have simulated data, we can find the accuracy
accuracy <- phybreak.accuracy(phybreak.true = sim_data, MCMCstate = MCMCstate)

#without the simulated data, we can find the posterior pro
post.prob <- phybreak.infector.posts(MCMCstate = MCMCstate)

## -----------------------------------------------------------------------------
#if true history is known
plot(phybreak.plot.posteriors(accuracy$post_prob_df))

#if true history is unknown
plot(phybreak.plot.posteriors(post.prob$post_prob_df))

