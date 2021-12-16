#' @title Simulate biomarker values 
#' @description Function to simulate biomarkers
#' This version uses the actual CD4 values, not the square root transformed values of previous versions
#' @param t.inf A vector of the true infection ages for each individual
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
                           mub, Sigmab, sigmae, 
                           which.biomarkers = matrix(1, nrow = length(t.inf), ncol = 5), 
                           seed = sample(1e9, 1)){
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
  
  biomarkers <- matrix(-1, nrow = nInds, ncol = 5)
  
  #simulate biomarker values
  for(i in 1:nInds){
    #reroll if any values are negative
    while(any(biomarkers[i,which(!is.na(biomarkers[i,]))] < 0)){ #are any of the used biomarkers negative?
      b_sim[i,] <- mvtnorm::rmvnorm(1, mub, Sigmab)
      
      #BED
      biomarkers[i,1] <- b_sim[i,1] + (b_sim[i,2] - b_sim[i])*exp(-exp(b_sim[i,3])*t.inf[i]) + rnorm(1, 0, sqrt(sigmae[1]))
      #LAg
      biomarkers[i,2] <- b_sim[i,8]+(b_sim[i,9]-b_sim[i,8])*exp(-exp(b_sim[i,10])*t.inf[i]) + rnorm(1, 0, sqrt(sigmae[2]))
      #CD4
      biomarkers[i,3] <- ((b_sim[i,4]+b_sim[i,5]*t.inf[i] + rnorm(1, 0, sqrt(sigmae[3])))*24)^2
      #pol
      biomarkers[i,4] <- (b_sim[i,6]+b_sim[i,7]*t.inf[i] + rnorm(1, 0, sqrt(sigmae[4])))/200
      #pol2
      biomarkers[i,5] <- (b_sim[i,11]+b_sim[i,12]*t.inf[i] + rnorm(1, 0, sqrt(sigmae[5])))/200
      
      #make it so unused biomarkers are NA
      biomarkers[i,which(which.biomarkers[i,] == 0)] <- NA
    }
  }
  return(biomarkers)
}

#' @title Predict infection times with MBM
#' @description Function to use the multiple biomarker model for infection time prediction
#'This version uses the actual CD4 values, not the square root transformed values of previous versions
#' @param BED A vector of BED values, one for each individual. Can be omitted.
#' @param LAg A vector of LAg values, one for each individual. Can be omitted.
#' @param CD4 A vector of CD4 values, one for each individual. These are the actual CD4+ T-cell concentrations, not the square-root transformed ones used in previous versions. 
#'   This value must be provided. 
#' @param pol A vector for the proportion of polymorphic sites in the HIV polymerase gene, one for each individual.
#'   This value must be provided.
#' @param pol2 A vector of values for the diversity in the HIV polymerase gene as obtained from NGS, one for each individual. Can be omitted.
#' @param prev.neg.time The amount of time in between an individual's positive HIV test and a previous negative test. 
#'   Use Inf or NA if no previous negative test is available.
#' @param t.sam.delay A vector for the difference in time (in years (days/365.25)) between the dates of the positive test and when the biomarkers aside from CD4 were sampled.
#' @param t.CD4.delay A vector for the difference in time (in years (days/365.25)) between the dates of the positive test and when CD4 was sampled.
#' @param mub A vector of length 12 for the means of the random effects parameters for the Multiple Biomarker Model.
#' @param Sigmab A 12x12 matrix for the covariance matrix of the random effects parameters for the Multiple Biomarker Model
#' @param sigmae A vector for the variances (not standard deviations) of the error paramemeters in the Multiple Biomarker Model
#' @param n.adapt The number of iterations to adapt the rjags MCMC model
#' @param n.burn The number of iterations of burn-in for the rjags MCMC model
#' @param n.iter The number of iterations of samples for the rjags MCMC model.
#' @param inf.mean The mean of the gamma distribution used for the prior distribution when there is no previous negative test. Default value is 2 years.
#' @param inf.sd The mean of the gamma distribution used for the prior distribution when there is no previous negative test. Default value is 1.5 years.
#' @param max.seroconvert.delay The maximum reasonable about of time that someone can be infected before an HIV test would be positive.
#' @param output.raw If TRUE, the output will include all MCMC samples for the infection ages.
#' @param seed The RNG seed to use.
#' 
#' @return The probability density of the infection times for each patient. 
#'   "pdf" is a continous function and "pdf.num" is a density class object as produced by the density function.
#'   If output.raw is TRUE, "raw" is the matrix of MCMC samples of the infection ages for all individuals.
#' @export

mbm.predict <- function(BED = rep(NA, length(CD4)), LAg = rep(NA, length(CD4)), CD4, pol, pol2 = rep(NA, length(CD4)),
                        prev.neg.time = NULL, 
                        t.sam.delay = rep(0, length(CD4)),
                        t.CD4.delay = rep(0, length(CD4)),
                        mub, 
                        Sigmab, 
                        sigmae,
                        n.adapt = 1000, n.burn = 1000, n.iter = 1000, 
                        inf.mean = 2, inf.sd = 1.5, 
                        max.seroconvert.delay = 2/12,
                        output.raw = FALSE,
                        seed = sample(1e9, 1)){
  #output can be "continuous", "numeric", or "both"
  
  #put into matrix format for JAGS model if it is not already
  if(!is.matrix(CD4)){
    cd4.test <- sqrt(matrix(CD4, ncol = 1))
  } else{
    cd4.test <- CD4
  }
  if(!is.matrix(pol)){
    pol.test <- matrix(pol, ncol = 1)
  } else{
    pol.test <- pol
  }
  if(!is.matrix(BED)){
    bed.test <- matrix(BED, ncol = dim(CD4)[2])
  } else{
    bed.test <- BED
  }
  if(!is.matrix(LAg)){
    lag.test <- matrix(LAg, ncol = dim(CD4)[2])
  } else{
    lag.test <- LAg
  }
  if(!is.matrix(pol2)){
    pol2.test <- matrix(pol2, ncol = dim(CD4)[2])
  } else{
    pol2.test <- pol2
  }
  if(!is.matrix(t.sam.delay)){
    t.sam.delay.test <- matrix(t.sam.delay, ncol = dim(CD4)[2])
  } else{
    t.sam.delay.test <- t.sam.delay
  }
  if(!is.matrix(t.CD4.delay)){
    t.CD4.delay.test <- matrix(t.CD4.delay, ncol = dim(CD4)[2])
  } else{
    t.CD4.delay.test <- t.CD4.delay
  }
  
  nInds <- dim(bed.test)[1]
  nSamples <- dim(bed.test)[2]
  
  #set previous negative times to infinity if they are not provided
  if(is.null(prev.neg.time)){
    prev.neg.time <- rep(Inf, nInds)
  }
  
  #check that all biomarker vectors have the same length
  if(dim(lag.test)[1] != nInds | dim(cd4.test)[1] != nInds | dim(pol.test)[1] != nInds | 
     dim(pol2.test)[1] != nInds | dim(t.sam.delay.test)[1] != nInds | dim(t.CD4.delay.test)[1] != nInds | 
     length(prev.neg.time) != nInds){
    stop("BED, LAg, CD4, pol, pol2, t.sam.delay, t.CD4.delay, and prev.neg.test must all have the same number of individuals")
  }
  #check that all biomarkers have the same number of samples
  if(dim(lag.test)[2] != nSamples | dim(cd4.test)[2] != nSamples | dim(pol.test)[2] != nSamples | 
     dim(pol2.test)[2] != nSamples | dim(t.sam.delay.test)[2] != nSamples | dim(t.CD4.delay.test)[2] != nSamples){
    stop("BED, LAg, CD4, pol, pol2, t.sam.delay, and t.CD4.delay 
         must all have the same number of samples (use NA for unavailable samples)")
  }
  
  #find number of non NA samples per individual
  mbed <- rowSums(!is.na(bed.test))
  mlag <- rowSums(!is.na(lag.test))
  mcd4 <- rowSums(!is.na(cd4.test))
  mpol <- rowSums(!is.na(pol.test))
  mpol2 <- rowSums(!is.na(pol2.test))
  mt.sam.delay <- rowSums(!is.na(t.sam.delay.test))
  mt.CD4.delay <- rowSums(!is.na(t.CD4.delay.test))
  
  mp <- apply(cbind(mbed, mlag, mcd4, mpol, mpol2, mt.sam.delay, mt.sam.delay, mt.CD4.delay), 1, max)
  
  #check to make sure output request is valid
  #if(!(output == "continuous" | output == "numeric" | output == "both")){
  #  stop("Invalid output type. Must be continuous, numeric, or both.")
  #}
  
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
  
  input_data <- list(np = nInds, mp = mp, 
                     mub = mub, Precb = solve(Sigmab), prece = 1/sigmae,
                     #t_test and tcd4 test are there in case the tests were at different times
                     t_test = t.sam.delay.test, tcd4_test = t.CD4.delay.test,
                     bed_test = bed.test, cd4_test = cd4.test/24, 
                     pol_test = pol.test*200,
                     lag_test = lag.test, pol2_test = pol2.test*200, 
                     prior_shape = trunc.shape, prior_rate = trunc.rate, neg_time = max.uninf.time)
  
  #prediction model
  mbm_predict <- "
  data{
    #scale stuff here?
    
    #make identity matrix
    for(i in 1:12){
      for(j in 1:12){
        eye[i,j] <- equals(i,j) #eye is the identity matrix
      }
    }
  }
  
  model{
    #####Prediction########
    #prior for age of infection
    for(i in 1:np){
      #t_pred[i] ~ dunif(0,12)
      #t_pred[i] ~ dgamma(16/9, 0.75) #slightly higher mean
      #t_pred[i] ~ dgamma(16/9, 8/9) #corresponds to mean of 2 with stdev 1.5
      t_pred[i] ~ dgamma(prior_shape[i], prior_rate[i]) T(,neg_time[i])
      #t_unif[i] ~ dunif(0,1)
      #t_pred[i] <- interp.lin(t_unif[i], sg_cdf, t_ind) #use inverse transform sampling to generate values from desired prior
    }
    
    #parameters for predictions
    for(i in 1:np) {
      b_pred[i,1:12] ~ dmnorm(mub[1:12], Precb[1:12,1:12])
    }
    
    #prediction from testing values
    for(i in 1:np){
      for(j in 1:mp[i]){
        bed_test[i,j] ~ dnorm(mu_bed_pred[i,j],prece[1])
        mu_bed_pred[i,j] <- b_pred[i,1] + (b_pred[i,2] - b_pred[i,1])*exp(-exp(b_pred[i,3])*(t_pred[i] + t_test[i,j]))
        
        lag_test[i,j] ~ dnorm(mu_lag_pred[i,j],prece[2])
        mu_lag_pred[i,j] <- b_pred[i,8] + (b_pred[i,9] - b_pred[i,8])*exp(-exp(b_pred[i,10])*(t_pred[i] + t_test[i,j]))
        
        cd4_test[i,j] ~ dnorm(mu_cd4_pred[i,j],prece[3])
        mu_cd4_pred[i,j]<- b_pred[i,4] + b_pred[i,5]*(t_pred[i] + tcd4_test[i,j])
        
        pol_test[i,j] ~ dnorm(mu_pol_pred[i,j],prece[4])
        mu_pol_pred[i,j]<- b_pred[i,6] + b_pred[i,7]*(t_pred[i] + t_test[i,j])
        
        pol2_test[i,j] ~ dnorm(mu_pol2_pred[i,j],prece[5])
        mu_pol2_pred[i,j]<- b_pred[i,11] + b_pred[i,12]*(t_pred[i] + t_test[i,j])
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
  
  #return(list(t_dens_unc, t_dens_sim, sample.pdf))
  if(output.raw == TRUE){
    return(list(pdf = pdf, pdf.num = pdf.num, raw = mbm_pred[[1]]))
  } else{
    return(list(pdf = pdf, pdf.num = pdf.num))
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