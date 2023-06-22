#' @title Simulate biomarker values 
#' @description Function to simulate biomarkers
#' This version uses the actual CD4 values, not the square root transformed values of previous versions
#' @param t.inf A vector or list of the true infection ages for each individual. 
#' If given as a list, multiple infection ages can used per individual, and the output will be returned as a list.
#' @param mub A vector of length 12 for the means of the random effects parameters for the Multiple Biomarker Model.
#' @param Sigmab A 12x12 matrix for the covariance matrix of the random effects parameters for the Multiple Biomarker Model
#' @param sigmae A vector for the variances (not standard deviations) of the error paramemeters in the Multiple Biomarker Model
#' @param which.biomarkers A vector of length 5 or matrix of size (number of individuals)x5 for the biomarkers to simulate for each individual.
#'   If 1, the corresponding biomarker will be simulated. If 0, the corresponding biomarker will not be simulated.
#'   The order of the biomarkers is: BED, LAg, CD4, pol, and pol2.
#'   If a vector is given, the same set of biomarkers will be used for each individual (although the actual biomarker values are gemerated independently)
#' @param seed The RNG seed to use.
#' 
#' @return A Matrix of simulated biomarkers from the supplied model parameters. The order of the biomarkers is: BED, LAg, CD4, pol, and pol2.
#' @export

sim.biomarkers <- function(t.inf, 
                           mub = NULL, Sigmab = NULL, sigmae = NULL, 
                           which.biomarkers = matrix(1, nrow = length(t.inf), ncol = 5), 
                           seed = sample(1e9, 1)){
  #set RNG seed 
  set.seed(seed)
  
  #load trained parameters if not specified
  if(is.null(mub) && is.null(Sigmab) && is.null(sigmae)){
    data("MBM_pars")
    mub <- MBM_pars$mub
    Sigmab <- MBM_pars$Sigmab
    sigmae <- MBM_pars$sigmae
  } else if(is.null(mub) || is.null(Sigmab) || is.null(sigmae)){
    stop("mub, Sigmab, and sigmae must either all be specified or all null")
  }
  
  #number of individuals to simulate
  nInds <- length(t.inf)
  
  #make which.biomarkers into matrix if it is given as a vector
  if(length(which.biomarkers == 5) & is.null(dim(which.biomarkers))){
    which.biomarkers <- matrix(rep(which.biomarkers, nInds), nrow = nInds, byrow = TRUE)
  }
  #declare matrix for individual curve coefficients
  b_sim <- matrix(nrow = nInds, ncol = 12)
  #declare vectors for biomarker values
  #bed_sim <- rep(-1, nInd)
  #lag_sim <- rep(-1, nInd)
  #cd4_sim <- rep(-1, nInd)
  #pol_sim <- rep(-1, nInd)
  #pol2_sim <- rep(-1, nInd)
  
  if(!is.list(t.inf)){
    t.inf <- as.list(t.inf)
    return.list <- FALSE
  } else{
    return.list <- TRUE
  }
  
  #if(is.list(t_inf)){
  biomarkers <- list()
  for(i in 1:nInds){
    biomarkers[[i]] <- matrix(-1, nrow = length(t.inf[[i]]), ncol = 5)
  }
  #} else{
  #  biomarkers <- matrix(-1, nrow = nInds, ncol = 5)
  #}
  
  #simulate biomarker values
  for(i in 1:nInds){
    #reroll if any values are negative
    while(any(biomarkers[[i]][which(!is.na(biomarkers[[i]]))] < 0)){ #are any of the used biomarkers negative?
      b_sim[i,] <- mvtnorm::rmvnorm(1, mub, Sigmab)
      
      #BED
      biomarkers[[i]][,1] <- b_sim[i,1] + (b_sim[i,2] - b_sim[i])*exp(-exp(b_sim[i,3])*t.inf[[i]]) + rnorm(length(t.inf[[i]]), 0, sqrt(sigmae[1]))
      #LAg
      biomarkers[[i]][,2] <- b_sim[i,8]+(b_sim[i,9]-b_sim[i,8])*exp(-exp(b_sim[i,10])*t.inf[[i]]) + rnorm(length(t.inf[[i]]), 0, sqrt(sigmae[2]))
      #CD4
      biomarkers[[i]][,3] <- ((b_sim[i,4]+b_sim[i,5]*t.inf[[i]] + rnorm(length(t.inf[[i]]), 0, sqrt(sigmae[3])))*24)^2
      #pol
      biomarkers[[i]][,4] <- (b_sim[i,6]+b_sim[i,7]*t.inf[[i]] + rnorm(length(t.inf[[i]]), 0, sqrt(sigmae[4])))/200
      #pol2
      biomarkers[[i]][,5] <- (b_sim[i,11]+b_sim[i,12]*t.inf[[i]] + rnorm(length(t.inf[[i]]), 0, sqrt(sigmae[5])))/200
      
      #make it so unused biomarkers are NA
      biomarkers[[i]][,which(which.biomarkers[i,] == 0)] <- NA
    }
  }
  
  #add names to biomarker matrices
  for(i in 1:length(biomarkers)){
    colnames(biomarkers[[i]]) <- c("BED", "LAg", "CD4", "pol", "pol2")
  }
  
  if(return.list){
    return(biomarkers)
  } else{
    biomarkers_matrix <- biomarkers[[1]]
    for(i in 2:length(biomarkers)){
      biomarkers_matrix <- rbind2(biomarkers_matrix, biomarkers[[i]])
    }
    return(biomarkers_matrix)
  }
}

#' @title Predict infection times with MBM
#' @description Function to use the multiple biomarker model for infection time prediction
#'This version uses the actual CD4 values, not the square root transformed values of previous versions
#' @param BED A vector of BED values, one for each individual. Can be omitted.
#' @param LAg A vector of LAg values, one for each individual. Can be omitted.
#' @param CD4 A vector of CD4 values, one for each individual. These are the actual CD4+ T-cell concentrations, not the square-root transformed ones used in previous versions. 
#'   Can be omitted.
#' @param pol A vector for the proportion of polymorphic sites in the HIV polymerase gene, one for each individual.
#'   Can be omitted.
#' @param pol2 A vector of values for the diversity in the HIV polymerase gene as obtained from NGS, one for each individual. Can be omitted.
#' @param prev.neg.time The amount of time in between an individual's positive HIV test and a previous negative test. 
#'   Use Inf or NA if no previous negative test is available.
#' @param t.BED.delay A vector or matrix for the difference in time (in years (days/365.25)) between the dates of the positive test and when BED was sampled.
#' @param t.LAg.delay A vector or matrix for the difference in time (in years (days/365.25)) between the dates of the positive test and when LAg was sampled.
#' @param t.CD4.delay A vector or matrix for the difference in time (in years (days/365.25)) between the dates of the positive test and when CD4 was sampled.
#' @param t.pol.delay A vector or matrix for the difference in time (in years (days/365.25)) between the dates of the positive test and when pol was sampled.
#' @param t.pol2.delay A vector or matrix for the difference in time (in years (days/365.25)) between the dates of the positive test and when pol2 was sampled.
#' @param mub A vector of length 12 for the means of the random effects parameters for the Multiple Biomarker Model.
#' @param Sigmab A 12x12 matrix for the covariance matrix of the random effects parameters for the Multiple Biomarker Model
#' @param sigmae A vector for the variances (not standard deviations) of the error paramemeters in the Multiple Biomarker Model
#' @param n.adapt The number of iterations to adapt the rjags MCMC model
#' @param n.burn The number of iterations of burn-in for the rjags MCMC model
#' @param n.iter The number of iterations of samples for the rjags MCMC model.
#' @param prior.type The type of prior used for the time between infection and diagnosis in the multiple biomarker model.
#' 1 indicates a gamma distribution with mean `inf.mean` years and standard deviation `inf.sd` years.
#' 2 indicates a continuous uniform distribution with minimum 0 and maximum 12 years.
#' 3 indicates a distribution for an individual that is HIV positive but has not developed AIDS symptoms, 
#' assuming that AIDS symptoms develop after a length of time according to a gamma distribution with shape 3.349 and rate 0.327
#' 4 indicates a user-supplied distribution from the parameter user.prior.pdf
#' @param inf.mean The mean of the gamma distribution used for the prior distribution when there is no previous negative test. Default value is 2 years.
#' @param inf.sd The mean of the gamma distribution used for the prior distribution when there is no previous negative test. Default value is 1.5 years.
#' @param max.seroconvert.delay The maximum reasonable about of time that someone can be infected before an HIV test would be positive.
# '@param u1.pdf A pdf to use as the prior distribution for the time between infection and diagnosis in the multiple biomarker model.
#' It must be a list with x and y components corresponding to the time before diagnosis and probability density.
#' @param output.raw If TRUE, the output will include all MCMC samples for the infection ages.
#' @param seed The RNG seed to use.
#' @return The probability density of the infection times for each patient. 
#'   "pdf" is a continuous function and "pdf.num" is a density class object as produced by the density function.
#'   If output.raw is TRUE, "raw" is the matrix of MCMC samples of the infection ages for all individuals.
#' @export

mbm.predict <- function(BED = rep(NA, 1), 
                        LAg = rep(NA, 1), 
                        CD4 = rep(NA, 1), 
                        pol = rep(NA, 1), 
                        pol2 = rep(NA, 1),
                        prev.neg.time = NULL, 
                        t.BED.delay = rep(0, length(BED)),
                        t.LAg.delay = rep(0, length(LAg)),
                        t.CD4.delay = rep(0, length(CD4)),
                        t.pol.delay = rep(0, length(pol)),
                        t.pol2.delay = rep(0, length(pol2)),
                        mub = NULL, 
                        Sigmab = NULL, 
                        sigmae = NULL,
                        n.adapt = 1000, n.burn = 1000, n.iter = 1000, 
                        prior.type = NULL,
                        inf.mean = 2, inf.sd = 1.5, 
                        max.seroconvert.delay = 2/12,
                        u1.pdf = NULL,
                        output.raw = FALSE,
                        seed = sample(1e9, 1),
                        ...){
  #set RNG seed 
  set.seed(seed)
  
  #load trained parameters if not specified
  if(is.null(mub) && is.null(Sigmab) && is.null(sigmae)){
    data("MBM_pars")
    mub <- MBM_pars$mub
    Sigmab <- MBM_pars$Sigmab
    sigmae <- MBM_pars$sigmae
  } else if(is.null(mub) || is.null(Sigmab) || is.null(sigmae)){
    stop("mub, Sigmab, and sigmae must either all be specified or all null")
  }
  
  #convert legacy input
  addn.input <- list(...)
  if("t.sam.delay" %in% names(addn.input)) t.pol.delay <- addn.input$t.sam.delay
  
  #put into matrix format for JAGS model if it is not already
  if(!is.matrix(CD4)){
    cd4.test <- sqrt(matrix(CD4, ncol = 1))
  } else{
    cd4.test <- sqrt(CD4)
  }
  if(!is.matrix(pol)){
    pol.test <- matrix(pol, ncol = 1)
  } else{
    pol.test <- pol
  }
  if(!is.matrix(BED)){
    bed.test <- matrix(BED, ncol = 1)
  } else{
    bed.test <- BED
  }
  if(!is.matrix(LAg)){
    lag.test <- matrix(LAg, ncol = 1)
  } else{
    lag.test <- LAg
  }
  if(!is.matrix(pol2)){
    pol2.test <- matrix(pol2, ncol = 1)
  } else{
    pol2.test <- pol2
  }
  if(!is.matrix(t.pol.delay)){
    t.pol.delay.test <- matrix(t.pol.delay, ncol = dim(pol.test)[2])
  } else{
    t.pol.delay.test <- t.pol.delay
  }
  if(!is.matrix(t.CD4.delay)){
    t.cd4.delay.test <- matrix(t.CD4.delay, ncol = dim(cd4.test)[2])
  } else{
    t.cd4.delay.test <- t.CD4.delay
  }
  if(!is.matrix(t.BED.delay)){
    t.bed.delay.test <- matrix(t.BED.delay, ncol = dim(bed.test)[2])
  } else{
    t.bed.delay.test <- t.BED.delay
  }
  if(!is.matrix(t.LAg.delay)){
    t.lag.delay.test <- matrix(t.LAg.delay, ncol = dim(lag.test)[2])
  } else{
    t.lag.delay.test <- t.LAg.delay
  }
  if(!is.matrix(t.pol2.delay)){
    t.pol2.delay.test <- matrix(t.pol2.delay, ncol = dim(pol2.test)[2])
  } else{
    t.pol2.delay.test <- t.pol2.delay
  }
  
  nInds <- dim(bed.test)[1]
  nSamples <- dim(bed.test)[2]
  
  #set previous negative times to infinity if they are not provided
  if(is.null(prev.neg.time)){
    prev.neg.time <- rep(Inf, nInds)
  }
  
  #check that all biomarker vectors have the same length
  if(dim(lag.test)[1] != nInds || dim(cd4.test)[1] != nInds || dim(pol.test)[1] != nInds || 
     dim(pol2.test)[1] != nInds || dim(t.pol.delay.test)[1] != nInds || dim(t.cd4.delay.test)[1] != nInds || 
     dim(t.bed.delay.test)[1] != nInds || dim(t.lag.delay.test)[1] != nInds || dim(t.pol2.delay.test)[1] != nInds || 
     length(prev.neg.time) != nInds){
    stop("BED, LAg, CD4, pol, pol2, t.pol.delay, t.CD4.delay, t.BED.delay, t.LAg.delay, t.pol2.delay,
         and prev.neg.test must all have the same number of individuals")
  }
  
  #find number of non NA samples per individual (return 1 if only sample is NA)
  mbed <- max(rowSums(!is.na(bed.test)), 1)
  mlag <- max(rowSums(!is.na(lag.test)), 1)
  mcd4 <- max(rowSums(!is.na(cd4.test)), 1)
  mpol <- max(rowSums(!is.na(pol.test)), 1)
  mpol2 <- max(rowSums(!is.na(pol2.test)), 1)
  mt.bed.delay <- rowSums(!is.na(t.bed.delay.test))
  mt.lag.delay <- rowSums(!is.na(t.lag.delay.test))
  mt.cd4.delay <- rowSums(!is.na(t.cd4.delay.test))
  mt.pol.delay <- rowSums(!is.na(t.pol.delay.test))
  mt.pol2.delay <- rowSums(!is.na(t.pol2.delay.test))
  
  mp <- apply(cbind(mbed, mlag, mcd4, mpol, mpol2, 
                    mt.bed.delay, mt.lag.delay, mt.cd4.delay, mt.pol.delay, mt.pol2.delay), 1, max)
  
  #check to make sure output request is valid
  #if(!(output == "continuous" | output == "numeric" | output == "both")){
  #  stop("Invalid output type. Must be continuous, numeric, or both.")
  #}
  
  #parse prior distribution type
  if(is.null(prior.type)){
    prior.type <- rep(1L, nInds)
  } else if(any(!(prior.type %in% c(1,2,3,4)))){
    stop("prior.type values must be 1, 2, 3, or 4")
  } else if(length(prior.type) == 1){
    prior.type <- rep(prior.type, nInds)
  } else if(length(prior.type) != nInds){
    stop("Length of prior.type must be equal to the number of individuals")
  }
  #use one-hot encoding to set up a switch for the prior
  priors <- matrix(0, nrow = nInds, ncol = 4)
  for(i in seq_len(nInds)){
    priors[i,prior.type[i]] <- 1
  }
  
  #modify priors based on previous negative test
  
  #change mean and sd to shape and rate
  inf.shape <- (inf.mean/inf.sd)^2
  inf.rate <- inf.shape/inf.mean
  
  #cutoff value for when to start shrinking mean and sd in case of a previous negative test
  cutoff <- qgamma(.95, shape = inf.shape, rate = inf.rate)
  
  #maximum reasonable time until seroconversion
  #max.seroconvert.delay <- 2/12 #2 months
  
  #maximum reasonable time until and individual is uninfected
  max.uninf.time <- prev.neg.time + max.seroconvert.delay
  
  #mean and sd for individual specific truncated distributions
  trunc.mean <- rep(inf.mean, nInds)
  trunc.sd <- rep(inf.sd, nInds)
  
  for(i in 1:nInds){
    if(!is.na(prev.neg.time[i])){
      if(max.uninf.time[i] < cutoff){
        trunc.mean[i] <- inf.mean*(max.uninf.time[i]/cutoff)
        trunc.sd[i] <- inf.sd*(max.uninf.time[i]/cutoff)
      }
    }
  }
  #change mean and sd to shape and rate
  trunc.shape <- (trunc.mean/trunc.sd)^2
  trunc.rate <- trunc.shape/trunc.mean
  
  #change NAs in max.uninf.time to Inf
  max.uninf.time[is.na(max.uninf.time)] <- Inf
  
  #create uninformative prior distribution for time between infection and diagnosis
  #for patients HIV positive but without AIDS
  #"survival gamma"
  sg.shape <- 3.349
  sg.rate <- .3270 #gamma parameters for the amount of time to develop AIDS
  sg.step.size <- .05
  ts <- seq(0, 50, by = sg.step.size)
  sg_pdf <- (1 - pgamma(ts, shape = sg.shape, rate = sg.rate))/(sg.shape/sg.rate)
  
  sg_cdf <- rep(0, length = length(sg_pdf))
  for(i in 2:length(sg_cdf)){
    sg_cdf[i] <- sg_cdf[i-1] + mean(sg_pdf[i-1],sg_pdf[i])*sg.step.size #trapezoid rule
  }
  sg_cdf <- sg_cdf/max(sg_cdf)
  
  #find truncation value for icdf sampling
  sg_cdf_fn <- approxfun(ts, sg_cdf, yright = 1) #function
  sg_upper <- sg_cdf_fn(max.uninf.time)
  
  #process user-input prior distribution
  u1_pdf <- u1.pdf
  if(!is.null(u1_pdf)){
    #check to make sure it has x and y components
    if(!is.numeric(u1_pdf$x) || !is.numeric(u1_pdf$y)){
      stop("If provided, u1.pdf must have x and y components")
    } else if(length(u1_pdf$x) != length(u1_pdf$y)){
      stop("x and y components of u1.pdf must be the same length")
    }
  } else{
    #placeholder uniform distribution if not provided
    u1_pdf <- list(x = c(0, 10), y = c(1/10, 1/10))
  }
  
  #find cdf
  u1_cdf <- rep(0, length = length(u1_pdf$y))
  for(i in 2:length(u1_cdf)){
    u1_cdf[i] <- u1_cdf[i-1] + mean(u1_pdf$y[i-1],u1_pdf$y[i])*(u1_pdf$x[i]-u1_pdf$x[i-1]) #trapezoid rule
  }
  u1_cdf <- u1_cdf/max(u1_cdf)
  
  #collapse to unique cdf values
  u1_uniques <- sapply(unique(u1_cdf), 
                       FUN = function(unique_cdf_y, cdf_y) which(unique_cdf_y == cdf_y)[1], 
                       cdf_y = u1_cdf) #first of each value
  
  u1_cdf <- u1_cdf[u1_uniques]
  u1_t_ind <- u1_pdf$x[u1_uniques]
  
  #print(length(unique(u1_cdf)))
  #find truncation value for icdf sampling
  u1_cdf_fn <- approxfun(u1_t_ind, u1_cdf, yright = 1) #function
  u1_upper <- u1_cdf_fn(max.uninf.time)
  
  
  #add right values
  #u1_t_ind <- c(u1_t_ind, 100)
  #u1_cdf <- c(u1_cdf, 1)
  
  input_data <- list(np = nInds, 
                     mbed = mbed, mlag = mlag, mcd4 = mcd4, mpol = mpol, mpol2 = mpol2, 
                     mub = mub, Precb = solve(Sigmab), prece = 1/sigmae,
                     t_bed_delay = t.bed.delay.test, 
                     t_lag_delay = t.lag.delay.test, 
                     t_cd4_delay = t.cd4.delay.test,
                     t_pol_delay = t.pol.delay.test, 
                     t_pol2_delay = t.pol2.delay.test, 
                     bed_test = bed.test, 
                     lag_test = lag.test, 
                     cd4_test = cd4.test/24, 
                     pol_test = pol.test*200,
                     pol2_test = pol2.test*200, 
                     priors = priors,
                     prior_shape = trunc.shape, prior_rate = trunc.rate, neg_time = max.uninf.time,
                     sg_cdf = sg_cdf, sg_t_ind = ts, sg_upper = sg_upper, 
                     u1_cdf = u1_cdf, u1_t_ind = u1_t_ind, u1_upper = u1_upper)
  
  #prediction model
  mbm_predict <- "
  model{
    #####Prediction########
    #prior for age of infection
    for(i in 1:np){
      #t_pred[i] ~ dgamma(16/9, 0.75) #slightly higher mean
      #t_pred[i] ~ dgamma(16/9, 8/9) #corresponds to mean of 2 with stdev 1.5
      t_pred1[i] ~ dgamma(prior_shape[i], prior_rate[i]) T(,neg_time[i])
      t_pred2[i] ~ dunif(0,min(12, neg_time[i]))
      t_unif_sg[i] ~ dunif(0,sg_upper[i])
      t_pred3[i] <- interp.lin(t_unif_sg[i], sg_cdf, sg_t_ind) #use inverse transform sampling to generate values from desired prior
      t_unif_u1[i] ~ dunif(0,u1_upper[i]) #first user-input function
      t_pred4[i] <- interp.lin(t_unif_u1[i], u1_cdf, u1_t_ind) #use inverse transform sampling to generate values from desired prior
      t_pred[i] <- step(-0.5+priors[i,1])*t_pred1[i] + step(-0.5+priors[i,2])*t_pred2[i] +
                   step(-0.5+priors[i,3])*t_pred3[i] + step(-0.5+priors[i,4])*t_pred4[i]
    }
    
    #parameters for predictions
    for(i in 1:np) {
      b_pred[i,1:12] ~ dmnorm(mub[1:12], Precb[1:12,1:12])
    }
    
    #prediction from testing values
    for(i in 1:np){
      for(j in 1:mbed[i]){
        bed_test[i,j] ~ dnorm(mu_bed_pred[i,j],prece[1])
        mu_bed_pred[i,j] <- b_pred[i,1] + (b_pred[i,2] - b_pred[i,1])*exp(-exp(b_pred[i,3])*(t_pred[i] + t_bed_delay[i,j]))
      }
      for(j in 1:mlag[i]){
        lag_test[i,j] ~ dnorm(mu_lag_pred[i,j],prece[2])
        mu_lag_pred[i,j] <- b_pred[i,8] + (b_pred[i,9] - b_pred[i,8])*exp(-exp(b_pred[i,10])*(t_pred[i] + t_lag_delay[i,j]))
      }
      for(j in 1:mcd4[i]){
        cd4_test[i,j] ~ dnorm(mu_cd4_pred[i,j],prece[3])
        mu_cd4_pred[i,j]<- b_pred[i,4] + b_pred[i,5]*(t_pred[i] + t_cd4_delay[i,j])
      }
      for(j in 1:mpol[i]){
        pol_test[i,j] ~ dnorm(mu_pol_pred[i,j],prece[4])
        mu_pol_pred[i,j]<- b_pred[i,6] + b_pred[i,7]*(t_pred[i] + t_pol_delay[i,j])
      }
      for(j in 1:mpol2[i]){
        pol2_test[i,j] ~ dnorm(mu_pol2_pred[i,j],prece[5])
        mu_pol2_pred[i,j]<- b_pred[i,11] + b_pred[i,12]*(t_pred[i] + t_pol2_delay[i,j])
      }
    }
    
  }
  
  "
  
  model_jags <- rjags::jags.model(textConnection(mbm_predict), 
                                  data = input_data,
                                  #inits = list(mub = rep(0,12)),
                                  n.adapt = n.adapt)
  
  update(model_jags, n.burn)
  
  mbm_pred <- rjags::coda.samples(model_jags, 
                                  variable.names = c("t_pred"), 
                                  n.iter = n.iter)
  
  #uncorrected density
  t_dens_unc <- apply(mbm_pred[[1]],2,density, from = 0)
  
  #correct density for edge effects
  t_dens_sim <- lapply(t_dens_unc, FUN = function(t_dens) {
    t_dens$y <- t_dens$y/pnorm(t_dens$x, mean = 0, sd = t_dens$bw)
    return(t_dens)
  })
  
  #scale so it integrates to 1
  pdfs <- lapply(t_dens_sim, FUN = normalize.density)
  
  #extract continuous and numeric pdfs
  pdf <- lapply(pdfs, FUN = function(x){x$pdf})
  pdf.num <- lapply(pdfs, FUN = function(x){x$pdf.num})
  
  #find time between diagnosis and the first sequence
  diag.first.seq.delay <- apply(t.pol.delay.test, MARGIN = 1, FUN = min, na.rm = TRUE)
  
  #shift infection age distributions to be relative to first sequence time rather than diagnosis time
  pdf.num.seq <- mapply(FUN = diag.first.seq.shift, 
                        infection.age.dists.diag = pdf.num, diag.first.seq.delay = diag.first.seq.delay)
  
  #return(list(t_dens_unc, t_dens_sim, sample.pdf))
  if(output.raw == TRUE){
    return(list(pdf = pdf, pdf.num = pdf.num, pdf.num.seq = pdf.num.seq, raw = mbm_pred[[1]]))
  } else{
    return(list(pdf = pdf, pdf.num = pdf.num, pdf.num.seq = pdf.num.seq))
  }
}

normalize.density <- function(sample.density){
  density.length <- length(sample.density$x) #number of values in discrete density
  #make piecewise linear function for pdf
  sample.pdf <- approxfun(x = sample.density$x, y = sample.density$y, 
                          method = "linear", yleft = 0, yright = 0)
  tt <- seq(min(sample.density$x), max(sample.density$x), length.out = max(1024, density.length))
  sample.cdf_num <- rep(0, length(tt))
  for(i in 2:length(tt)){
    sample.cdf_num[i] <- sample.cdf_num[i-1] + integrate(f = sample.pdf, lower = tt[i-1], upper = tt[i])$value
  }
  #rescale functions so they integrate to 1
  integral <- max(sample.cdf_num)
  sample.pdf_num <- sample.density
  sample.pdf_num$y <- sample.pdf_num$y/integral
  sample.pdf <- approxfun(x = sample.pdf_num$x, y = sample.pdf_num$y, 
                          method = "linear", yleft = 0, yright = 0)
  #sample.cdf_num <- sample.cdf_num/integral
  #make into list with x and y values
  #sample.cdf_num <- list(x = tt, y = sample.cdf_num)
  #make piecewise linear function for cdf
  #sample.cdf <- approxfun(x = sample.cdf_num$x, y = sample.cdf_num$y, 
  #                        method = "linear", yleft = 0, yright = 1)
  #sample.dfs <- list(cdf_num = sample.cdf_num, cdf = sample.cdf, pdf = sample.pdf)
  #if(output == "continuous"){
  #  return(sample.pdf)
  #} else if(output == "numeric"){
  #  return(sample.pdf_num)
  #} else if(output == "both"){
  return(list(pdf = sample.pdf, pdf.num = sample.pdf_num))
  #}
}

#function to shift infection age distributions to be relative to first sequence time rather than diagnosis time
diag.first.seq.shift <- function(infection.age.dists.diag, diag.first.seq.delay){
  infection.age.dists.seq <- infection.age.dists.diag
  #shift distributions if needed
  if(diag.first.seq.delay > 0){
    infection.age.dists.seq$x <- c(0, diag.first.seq.delay-1e-10, diag.first.seq.delay + infection.age.dists.seq$x)
    infection.age.dists.seq$y <- c(1e-16, 1e-16, infection.age.dists.seq$y)
  }
  return(list(infection.age.dists.seq))
}