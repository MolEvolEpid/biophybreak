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
  
  #mean of the maximum posterior probability infector for each individual
  mean_max_post <- mean(apply(post_support, 2, max))
  
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
                accuracy = accuracy, mean_enrich = mean_enrich, mean_post = mean_post, mean_max_post = mean_max_post,
                corr_infect = corr_infect, post_true = post_true_each,
                true_infectors = true_infectors,
                post_prob_df = post_prob_df))
  } else{
    return(list(post_support = post_support, 
                accuracy = accuracy, accuracy_unsampled_infector = accuracy_unsampled_infector, 
                mean_enrich = mean_enrich, mean_enrich_unsampled_infector = mean_enrich_unsampled_infector,
                mean_post = mean_post, mean_post_unsampled_infector = mean_post_unsampled_infector, mean_max_post = mean_max_post,
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
  
  #mean of the maximum posterior probability infector for each individual
  mean_max_post <- mean(apply(post_support, 2, max))
  
  post_prob_df <- as.data.frame(t(post_support), stringsAsFactors = TRUE)
  post_prob_df$Individual <- rownames(post_prob_df)
  post_prob_df <- reshape2::melt(post_prob_df, id = "Individual")
  colnames(post_prob_df) <- c("Individual", "Infector", "Posterior.Support")
  
  #set order of individuals and infectors
  post_prob_df$Individual <- factor(post_prob_df$Individual, levels = colnames(post_support))
  post_prob_df$Infector <- factor(post_prob_df$Infector, levels = rownames(post_support))
  
  post_prob_df <- post_prob_df[order(post_prob_df$Individual, post_prob_df$Infector),]
  
  return(list(post_support = post_support, 
              mean_max_post = mean_max_post,
              post_prob_df = post_prob_df))
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
phybreak.plot.traces <- function(MCMCstate, main = ""){
  nInds <- length(unique(MCMCstate$d$hostnames))
  plot(MCMCstate$s$mu, type = 'l', ylab = "mu", main = main)
  plot(MCMCstate$s$wh.0, type = 'l', ylab = "Within Host Level", main = main)
  plot(MCMCstate$s$wh.s, type = 'l', ylab = "Within Host Slope", main = main)
  for(i in 1:nInds){
    plot(MCMCstate$s$inftimes[i,], type = 'l', ylab = paste0("tinf.", MCMCstate$d$hostnames[i]), main = main)
  }
  for(i in 1:nInds){
    plot(MCMCstate$s$infectors[i,], type = 'l', 
         ylim = c(0, nInds), ylab = paste0("infector.", MCMCstate$d$hostnames[i]), main = main)
  }
  plot(MCMCstate$s$logLik, type = 'l', ylab = "logLik", main = main)
}

#' @title Barplot of posterior supports for infectors
#' @description Function to create a barplot of the posterior supports for each infector
#' @param post_prob_df A dataframe from either the phybreak.accuracy or phybreak.infector.posts function
#' @param treecolors Colors for each individual
#' @param unsampled Indices of individuals that were not considered sampled during inference if using simulated data
#' @param xlab Label of the x axis
#' @param ylab Label of the y axis
#' @param angle Angle of x-axis labels
#' ... Other parameters that can be passed to ggplot::labs
#' @export
#' 
phybreak.plot.posteriors <- function(post_prob_df, treecolors = NULL, unsampled = NULL, 
                                     xlab = "Recipient", ylab = "Posterior Support", 
                                     angle = 90, ...){
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
    ggplot2::ggplot(data = post_prob_df, ggplot2::aes(x = Individual, y = Posterior.Support, fill = Infector, color = True.Infector)) + 
      ggplot2::geom_bar(position = "dodge", stat = "identity", width = 0.8) + 
      ggplot2::ylim(c(0,1)) + ggplot2::labs(...) + ggplot2::xlab(label = xlab) + ggplot2::ylab(label = ylab) +
      ggplot2::scale_fill_manual(values = c("#888888", treecolors)) +
      ggplot2::scale_color_manual(values = c("#DDDDDD", "#000000")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggtext::element_markdown(color = treecolors, angle = angle))
  } else{
    ggplot2::ggplot(data = post_prob_df, ggplot2::aes(x = Individual, y = Posterior.Support, fill = Infector)) + 
      ggplot2::geom_bar(position = "dodge", stat = "identity", width = 0.8) + 
      ggplot2::ylim(c(0,1)) + ggplot2::labs(...) + ggplot2::xlab(label = xlab) + ggplot2::ylab(label = ylab) +
      ggplot2::scale_fill_manual(values = c("#888888", treecolors)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggtext::element_markdown(color = treecolors, angle = angle))
  }
}

#' @title Plot MPC tree, infection age distributions, and infector posterior supports
#' @description Function to plot the maximum parent credibility tree, infection age distributions, 
#'   and posterior supports for infectors in one page
#' @param MCMCstate The output from sample_phybreak
#' @param phybreak.true The true transmission history (if using simulated data)
#' @param angle Angle of x-axis labels in the infector posterior density plot
#' @export
#' 
phybreak.plot.triple <- function(MCMCstate, phybreak.true = NULL,
                                 angle = 90){
  #number of individuals
  nInd <- length(unique(MCMCstate$d$hostnames))
  #find xmin appropriate from biomarker distributions
  time_min.1 <- mapply(MCMCstate$p$sample.pdf.nonpar, FUN = function(pdf.nonpar, sample.time){
    x <- rev(-pdf.nonpar$x) + sample.time
    y <- rev(pdf.nonpar$y)
    cdf <- find_cdf(x, y)
    #plot(x, y, type = 'l')
    #lines(x, cdf, type = 'l', lty = 2)
    min.eff <- min(x[cdf >= 0.1])
  }, sample.time = MCMCstate$d$sample.times[1:nInd])
  
  time_min.1.all <- min(time_min.1)
  
  #pdf(file = "3_test_margins.pdf", width = 5, height = 10)
  layout(matrix(1:3, nrow = 3))
  #plot MPC tree
  plotPhyloTrans(MCMCstate, plot.which = "mpc", 
                 xlim.adjust = c(time_min.1.all, 0),
                 xlab = "", mar = c(2,3.85,0,.2))
  #plot infection times
  treecolors <- hcl(unlist(sapply(1:floor(sqrt(nInd)) - 1, 
                                  function(xx) seq(xx, nInd - 1, 
                                                   floor(sqrt(nInd))))) * 360/nInd, 
                    c = 100, l = 65)
  xlims <- par("usr")[1:2]
  xmid <- mean(xlims)
  new_xlims <- (xlims-xmid)*.926+xmid
  new_xlims[2] <- diff(new_xlims) + new_xlims[1]
  ymax <- max(sapply(MCMCstate$p$sample.pdf.nonpar, FUN = function(dens){max(dens$y)}))
  par(mar = c(3, 2.6, 1, .2), mgp = c(1.5, 0.5, 0))
  plot(-MCMCstate$p$sample.pdf.nonpar[[1]]$x + MCMCstate$d$sample.times[1], 
       MCMCstate$p$sample.pdf.nonpar[[1]]$y, 
       type = 'l', col = treecolors[1], 
       xlim = new_xlims, 
       ylim = c(0, ymax),
       #main = paste0("Cluster "), 
       xlab = "Time (years after first diagnosis)", ylab = "Density",
       axes = FALSE)
  axis(1)
  axis(2)
  box(bty = 'l')
  for(j in 2:nInd){
    lines(-MCMCstate$p$sample.pdf.nonpar[[j]]$x + MCMCstate$d$sample.times[j], 
          MCMCstate$p$sample.pdf.nonpar[[j]]$y,
          type = 'l', col = treecolors[j])
  }
  if(class(phybreak.true) == "phybreakdata"){
    #plot true infection ages
    abline(v = phybreak.true$sim.infection.times, col = treecolors, lty = 2)
  }
  #plot infector posterior densities
  #find posteriors
  if(class(phybreak.true) == "phybreakdata"){
    post_prob_df <- phybreak.accuracy(phybreak.true = phybreak.true, MCMCstate = MCMCstate)$post_prob_df
  } else{
    post_prob_df <- phybreak.infector.posts(MCMCstate)$post_prob_df
  }
  plot.new()
  vps <- gridBase::baseViewports()
  grid::pushViewport(vps$figure)
  vp1 <-grid::plotViewport(c(0,0,0,0))
  posterior_dens <- phybreak.plot.posteriors(post_prob_df, 
                                             ylab = "Infector Posterior Support", angle = angle)
  print(posterior_dens, vp = vp1)
  grid::popViewport()
}

#get the cdf from the pdf for the infection age distribution
find_cdf <- function(x,y){
  #unnormalized pdf
  cdf_unnorm <- numeric(length = length(x))
  cdf_unnorm[1] <- 0
  x_diffs <- diff(x)
  for(i in 2:length(x)){ 
    cdf_unnorm[i] <- cdf_unnorm[i-1] + mean(y[i],y[i-1])*x_diffs[i-1] #linear approximation integral
  }
  cdf_norm <- cdf_unnorm/max(cdf_unnorm)
  return(cdf_norm)
}