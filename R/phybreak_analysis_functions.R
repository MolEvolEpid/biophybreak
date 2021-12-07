#phybreak analysis functions

#' @title Find accuracy of phybreak inference
#' @description Function to calculate the accuracy when comparing to the true simulated values
#' @param phybreak.true A phybreakdata object with the true transmission history
#' @param MCMCstate A phybreak object containing the information from the MCMC inference
#' @param complete Whether or not the MCMCstate inference was done with every individual in the true history.
#' Set to FALSE if there is at least one missing.
#' @return A list containing the posterior support, accuracy, posterior increase factor, mean true posterior, 
#' a vector of whether each infector was inferred correctly, a vector of the posterior support the true infectors, 
#' the true infectors, and a dataframe with the posterior supports and true infectors.
#' @export
#' 
phybreak.accuracy <- function(phybreak.true, MCMCstate, complete = TRUE){
  #names of individuals
  ind_names_true <- names(phybreak.true$sim.infectors)
  potential_infectors_true <- c("index", ind_names_true)
  
  ind_names_infer <- unique(MCMCstate$d$hostnames)
  potential_infectors_infer <- c("index", ind_names_infer)
  
  nInd_true <- length(ind_names_true)
  nInd_infer <- length(ind_names_infer)
  
  if(complete == FALSE){
    true_infectors <- phybreak.true$sim.infectors
    unsampled <- setdiff(potential_infectors_true, potential_infectors_infer)
    #find individuals infected by unsampled individuals
    infectees_of_unsampled <- names(true_infectors)[which(true_infectors %in% unsampled)]
    #find individuals infected by the unsampled individuals
    while(length(which(true_infectors %in% unsampled)) > 0){
      unsam_inf_inds <- which(true_infectors %in% unsampled)
      true_infectors[unsam_inf_inds] <- true_infectors[true_infectors[unsam_inf_inds]]
    }
    #removed unsampled individuals
    true_infectors <- true_infectors[-which(names(true_infectors) %in% unsampled)]
    #switch from name to number
    infectees_of_unsampled <- which(names(true_infectors) %in% infectees_of_unsampled)
  } else{
    true_infectors <- phybreak.true$sim.infectors
  }
  
  #vector for correctly identified infectors
  #corr_infect <- logical(nInd_true)
  corr_infect <- logical(nInd_infer)
  #vector for posterior support of each true infector
  #post_true_each <- numeric(nInd_true)
  post_true_each <- numeric(nInd_infer)
  
  #find counts of each infector in all the samples
  post_counts <- apply(MCMCstate$s$infectors + 1, 1, tabulate, nbins = nInd_infer + 1) #"+1" must be used because zeros are not tabulated
  #this means that post_counts[1,i] is the number of times individual i was sampled as the index case 
  #and post_counts[j,i] is the number of samples where individual i was infected by individual j-1
  #name rows and columns
  colnames(post_counts) <- ind_names_infer
  rownames(post_counts) <- potential_infectors_infer
  nIters <- length(MCMCstate$s$infectors[1,])
  post_support <- post_counts/nIters #normalize by number of samples to get percentage support
  
  for(j in 1:length(true_infectors)){
    corr_infect[j] <- rownames(post_counts)[which.max(post_counts[,j])] == true_infectors[j] #is max post infector true infector?
    post_true_each[j] <- post_support[true_infectors[j], j] #posterior support for true infector
    #}
  }
  if(nInd_true == 2 & nInd_infer == 3){
    accuracy <- as.numeric(post_support[2] > .5) #only 0 or 1
    mean_enrich <- post_support[2,2]*3
  } else{
    accuracy <- sum(corr_infect)/nInd_infer
    mean_enrich <- sum(post_true_each) #because you divide by 1/nInd to get enrichment, then divide by nInd to get the mean
    mean_post <- mean(post_true_each)
    if(complete == FALSE){
      accuracy_unsampled_infector <- sum(corr_infect[infectees_of_unsampled])/length(infectees_of_unsampled)
      mean_enrich_unsampled_infector <- sum(post_true_each[infectees_of_unsampled])
      mean_post_unsampled_infector <- mean(post_true_each[infectees_of_unsampled])
    }
    
  }
  
  post_prob_df <- as.data.frame(t(post_support), stringsAsFactors = TRUE)
  post_prob_df$Individual <- rownames(post_prob_df)
  post_prob_df <- reshape2::melt(post_prob_df, id = "Individual")
  colnames(post_prob_df) <- c("Individual", "Infector", "Posterior.Support")
  
  true_inf_num <- rep(0, length(post_prob_df$Individual))
  for(i in 1:length(post_prob_df$Individual)){
    if(length(which(names(true_infectors) == post_prob_df$Individual[i])) > 0){
      if((post_prob_df$Infector[i] == true_infectors[which(names(true_infectors) == post_prob_df$Individual[i])])) true_inf_num[i] <- 1
    } else{true_inf_num[i] <- NA}
  }
  
  post_prob_df$True.Infector <- as.logical(true_inf_num)
  #set order of individuals and infectors
  post_prob_df$Individual <- factor(post_prob_df$Individual, levels = colnames(post_support))
  post_prob_df$Infector <- factor(post_prob_df$Infector, levels = rownames(post_support))
  
  post_prob_df <- post_prob_df[order(post_prob_df$Individual, post_prob_df$Infector),]
  if(isTRUE(complete)){
    return(list(post_support = post_support, 
                accuracy = accuracy, mean_enrich = mean_enrich, mean_post = mean_post,
                corr_infect = corr_infect, post_true = post_true_each,
                true_infectors = true_infectors,
                post_prob_df = post_prob_df))
  } else{
    return(list(post_support = post_support, 
                accuracy = accuracy, accuracy_unsampled_infector = accuracy_unsampled_infector, 
                mean_enrich = mean_enrich, mean_enrich_unsampled_infector = mean_enrich_unsampled_infector,
                mean_post = mean_post, mean_post_unsampled_infector = mean_post_unsampled_infector, 
                corr_infect = corr_infect, post_true = post_true_each,
                true_infectors = true_infectors,
                post_prob_df = post_prob_df))
  }
  
}

#' @title Find phybreak posterior supports
#' @description Function to compute the posterior supports when the true infectors are unknown, 
#' i.e., when using real data.
#' @param MCMCstate A phybreak object containing the information from the MCMC inference
#' @return A list containing a matrix and a dataframe of the posterior supports
#' @export
#' 
phybreak.infector.posts <- function(MCMCstate){
  #names of individuals
  ind_names_infer <- unique(MCMCstate$d$hostnames)
  potential_infectors_infer <- c("index", ind_names_infer)
  
  nInd_infer <- length(ind_names_infer)
  
  #find counts of each infector in all the samples
  post_counts <- apply(MCMCstate$s$infectors + 1, 1, tabulate, nbins = nInd_infer + 1) #"+1" must be used because zeros are not tabulated
  #this means that post_counts[1,i] is the number of times individual i was sampled as the index case 
  #and post_counts[j,i] is the number of samples where individual i was infected by individual j-1
  #name rows and columns
  colnames(post_counts) <- ind_names_infer
  rownames(post_counts) <- potential_infectors_infer
  nIters <- length(MCMCstate$s$infectors[1,])
  post_support <- post_counts/nIters #normalize by number of samples to get percentage support
  
  post_prob_df <- as.data.frame(t(post_support), stringsAsFactors = TRUE)
  post_prob_df$Individual <- rownames(post_prob_df)
  post_prob_df <- reshape2::melt(post_prob_df, id = "Individual")
  colnames(post_prob_df) <- c("Individual", "Infector", "Posterior.Support")
  
  #set order of individuals and infectors
  post_prob_df$Individual <- factor(post_prob_df$Individual, levels = colnames(post_support))
  post_prob_df$Infector <- factor(post_prob_df$Infector, levels = rownames(post_support))
  
  post_prob_df <- post_prob_df[order(post_prob_df$Individual, post_prob_df$Infector),]
  
  return(list(post_support = post_support, post_prob_df = post_prob_df))
}

#' @title Find phybreak posterior probabilities
#' @description Function to compute the posterior probabilities for each MCMC iteration
#' @param MCMCstate A phybreak object containing the information from the MCMC inference
#' @return A list of the posterior probability and likelihoods
#' @export
#' 
phybreak.posterior <- function(MCMCstate){
  iters <- length(MCMCstate$s$logLik)
  logp_mu <- log(1/MCMCstate$s$mu)
  logp_wh.0 <- dgamma(x = MCMCstate$s$wh.0, shape = MCMCstate$h$wh.0.sh, rate = MCMCstate$h$wh.0.sh/MCMCstate$h$wh.0.av, log = TRUE)
  logp_wh.s <- dgamma(x = MCMCstate$s$wh.s, shape = MCMCstate$h$wh.s.sh, rate = MCMCstate$h$wh.0.sh/MCMCstate$h$wh.s.av, log = TRUE)
  logp_theta <- logp_mu + logp_wh.0 + logp_wh.s
  logPost <- MCMCstate$s$logLik + logp_mu + logp_wh.0 + logp_wh.s
  return(list(logPost = logPost, logLik = MCMCstate$s$logLik, 
              logp_theta = logp_theta, logp_mu = logp_mu, logp_wh.0 = logp_wh.0, logp_wh.s = logp_wh.s))
}

log.mean.prob <- function(logLik, values = 1:length(logLik)){
  mean <- log(mean(exp(logLik[values]-min(logLik)))) + min(logLik)
  return(mean)
}

#' @title Plot phybreak traces
#' @description Function to plot the MCMC traces
#' @param MCMCstate A phybreak object containing the information from the MCMC inference
#' @export
#' 
phybreak.plot.traces <- function(MCMCstate){
  nInds <- length(MCMCstate$d$hostnames)
  plot(MCMCstate$s$mu, type = 'l', ylab = "mu")
  plot(MCMCstate$s$wh.0, type = 'l', ylab = "Within Host Level")
  plot(MCMCstate$s$wh.s, type = 'l', ylab = "Within Host Slope")
  for(i in 1:nInds){
    plot(MCMCstate$s$inftimes[i,], type = 'l', ylab = paste0("tinf.", MCMCstate$d$hostnames[i]))
  }
  for(i in 1:nInds){
    plot(MCMCstate$s$infectors[i,], type = 'l', 
         ylim = c(0, nInds), ylab = paste0("infector.", MCMCstate$d$hostnames[i]))
  }
  plot(MCMCstate$s$logLik, type = 'l', ylab = "logLik")
}

#' @title Barplot of posterior supports for infectors
#' @description Function to create a barplot of the posterior supports for each infector
#' @param post_prob_df A dataframe from either the phybreak.accuracy or phybreak.infector.posts function
#' @param treecolors Colors for each individual
#' @param unsampled Indices of individuals that were not considered sampled during inference if using simulated data
#' @param ylab Label of the y axis
#' ... Other parameters that can be passed to ggplot::labs
#' @export
#' 
phybreak.plot.posteriors <- function(post_prob_df, treecolors = NULL, unsampled = NULL, ylab = "Posterior Support", ...){
  if(is.null(treecolors)){
    #number of sampled and unsampled individuals (if any)
    nInd <- length(levels(post_prob_df$Individual)) + length(unsampled)
    treecolors <- hcl(unlist(sapply(1:floor(sqrt(nInd)) - 1, 
                                    function(xx) seq(xx, nInd - 1, 
                                                     floor(sqrt(nInd))))) * 360/nInd, 
                      c = 100, l = 65)
    if(length(unsampled) > 0){
      treecolors <- treecolors[-unsampled]
    }
  }
  #label for y axis
  if(!is.null(post_prob_df$True.Infector)){
    ggplot(data = post_prob_df, aes(x = Individual, y = Posterior.Support, fill = Infector, color = True.Infector)) + 
      geom_bar(position = "dodge", stat = "identity", width = 0.8) + 
      ylim(c(0,1)) + labs(...) + ylab(label = ylab) +
      scale_fill_manual(values = c("#888888", treecolors)) +
      scale_color_manual(values = c("#DDDDDD", "#000000")) +
      theme_bw()
  } else{
    ggplot(data = post_prob_df, aes(x = Individual, y = Posterior.Support, fill = Infector)) + 
      geom_bar(position = "dodge", stat = "identity", width = 0.8) + 
      ylim(c(0,1)) + labs(...) + ylab(label = ylab) +
      scale_fill_manual(values = c("#888888", treecolors)) +
      theme_bw()
  }
}