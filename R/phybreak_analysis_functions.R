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
  max_posts <- apply(post_support, 2, max)
  mean_max_post <- mean(max_posts)
  
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
                max_posts = max_posts, mean_max_post = mean_max_post,
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
  max_posts <- apply(post_support, 2, max)
  mean_max_post <- mean(max_posts)
  
  post_prob_df <- as.data.frame(t(post_support), stringsAsFactors = TRUE)
  post_prob_df$Individual <- rownames(post_prob_df)
  post_prob_df <- reshape2::melt(post_prob_df, id = "Individual")
  colnames(post_prob_df) <- c("Individual", "Infector", "Posterior.Support")
  
  #set order of individuals and infectors
  post_prob_df$Individual <- factor(post_prob_df$Individual, levels = colnames(post_support))
  post_prob_df$Infector <- factor(post_prob_df$Infector, levels = rownames(post_support))
  
  post_prob_df <- post_prob_df[order(post_prob_df$Individual, post_prob_df$Infector),]
  
  return(list(post_support = post_support, 
              max_posts = max_posts,
              mean_max_post = mean_max_post,
              post_prob_df = post_prob_df))
}

#' @title Find MPC infectors
#' @description Function to find the maximum parent credibility infectors from a phybreak object
#' @param MCMCstate A phybreak object containing the information from the MCMC inference
#' @return The names and indices of the MPC infectors for each individual and the transmission heterogeneity
#' @export
mpcinfectors <- function(MCMCstate){
  samplenr <- .mpcinfector(MCMCstate, length(MCMCstate$s$logLik), TRUE, FALSE)
  mpcinfectors <- MCMCstate$s$infectors[,samplenr]
  nInd <- length(mpcinfectors)
  infectornames <- c("index", MCMCstate$d$hostnames[1:nInd])
  infectors <- infectornames[mpcinfectors+1] #+1 is to account for index
  names(infectors) <- MCMCstate$d$hostnames[1:nInd]
  #find transmission heterogeneity
  #(sorted and removed lowest number of transmissions since it must be 0)
  heterogeneity <- sd(sort(tabulate(mpcinfectors, nbins = length(mpcinfectors)))[-1]) #sorted and removed lowest number of transmissions
  return(list(names = infectors, indices = mpcinfectors, heterogeneity = heterogeneity))
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
#' @param ... Other parameters that can be passed to ggplot::labs
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
#' @param tree.col Colors to use for each individual
#' @param color.shuffle Changes to the default color ordering
#' @param ... Other parameters that can be passed to ggplot::labs
#' @export
#' 
phybreak.plot.triple <- function(MCMCstate, phybreak.true = NULL,
                                 angle = 90, tree.col = NULL, color.shuffle = NULL, ...){
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
  
  if(class(phybreak.true) != "phybreakdata"){
    layout(matrix(1:3, nrow = 3))
  } else{
    layout(matrix(c(2,1,3,4), nrow = 4))
  }
  
  #find colors
  if(is.null(tree.col)){
    treecolors <- hcl(unlist(sapply(1:floor(sqrt(nInd)) - 1, 
                                    function(xx) seq(xx, nInd - 1, 
                                                     floor(sqrt(nInd))))) * 360/nInd, 
                      c = 100, l = 65)
  } else{
    if(length(tree.col) != nInd){
      stop("tree.col must be the same length as the number of individuals")
    }
    treecolors <- tree.col
  }
  if(!is.null(color.shuffle)){
    if(length(color.shuffle) != nInd){
      stop("color.shuffle must be the same length as the number of individuals")
    }
    treecolors <- treecolors[color.shuffle]
  }
  
  #plot MPC tree
  plotPhyloTrans(MCMCstate, plot.which = "mpc", tree.col = treecolors,
                 xlim.adjust = c(time_min.1.all, 0),
                 xlab = "", mar = c(2,3.85,0,.2))
  #find x limits
  xlims <- par("usr")[1:2]
  xmid <- mean(xlims)
  new_xlims <- (xlims-xmid)*.926+xmid
  new_xlims[2] <- diff(new_xlims) + new_xlims[1]
  #plot true transmission history if it is known
  if(class(phybreak.true) == "phybreakdata"){
    #shift sampling and infection times if they are on a different scale as inferred
    shift <- MCMCstate$d$sample.times[1] - phybreak.true$sample.times[1]
    phybreak.true$sample.times <- phybreak.true$sample.times + shift
    phybreak.true$sim.infection.times <- phybreak.true$sim.infection.times + shift
    plotPhyloTrans(phybreak.true, tree.col = treecolors,
                   xlim.override = new_xlims, #may not work great if true infection times are older than predicted ones
                   xlab = "", mar = c(2,2.55,0,.2))
  }
  #plot infection times
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
                                             ylab = "Infector Posterior Support", angle = angle, treecolors = treecolors, ...)
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

#' @title Find transmission rates between individuals with certain labels
#' @description Function to find the expected number of transmission between and within individuals 
#' that belong to certain groups 
#' @param MCMCstate The output from sample_phybreak
#' @param labels The labels for each individual
#' @param infector.posterior.probabilities The matrix of posterior probabilities of the infectors, 
#' as supplied by phybreak.infector.posts(MCMCstate)$post_support
#' @param permute_test Whether or not to run a permutation test of how the individuals are labeled 
#' to find a null distribution for the label transmission rates
#' @param nPermute Number of permutations to use in the permutation test
#' @param percentiles Percentiles at which to find quantiles from the null distribution
#' @return A matrix containing the the probabilities for each individual to be infected by a sampled individual of a certain label
#' and a matrix with the number of transmissions between and within individuals with each label
#' @export
#' 
label.transmissions <- function(MCMCstate, labels, infector.posterior.probabilities = phybreak.infector.posts(MCMCstate)$post_support,
                                permute_test = FALSE, nPermute = 10000, percentiles = c(0.025, 0.975)){
  #make sure labels are a factor
  labels <- as.factor(labels)
  
  #levels of label factors
  label.levels <- levels(labels)
  
  #number of unique labels
  nLabels <- length(label.levels)
  
  #matrix for number of transmissions between pairs of labels
  label.transmissions <- matrix(0, nrow = nLabels, ncol = nLabels)
  
  #number of individuals
  nInd <- dim(infector.posterior.probabilities)[2]
  
  #labels of infectors for each individual
  infector.labels <- matrix(nrow = nInd, ncol = nLabels)
  for(i in seq_len(nInd)){
    for(j in seq_len(nLabels)){
      infector.labels[i,j] <- infector.posterior.probabilities[2:(nInd+1),i] %*% (labels == label.levels[j])
    }
  }
  
  for(i in seq_len(nInd)){
    label.transmissions[which(label.levels == labels[i]),] <- label.transmissions[which(label.levels == labels[i]),] + infector.labels[i,]
  }
  
  rownames(infector.labels) <- MCMCstate$d$hostnames[1:nInd]
  colnames(infector.labels) <- label.levels
  
  rownames(label.transmissions) <- label.levels
  colnames(label.transmissions) <- label.levels
  
  if(permute_test){
    #matrix for number of transmissions between pairs of labels
    label.transmissions.null <- array(dim = c(nLabels, nLabels, nPermute))
    #super slow but only takes a few seconds for 10k permutations
    for(i in seq_len(nPermute)){
      labels_permute <- sample(labels, replace = FALSE)
      label.transmissions.null[,,i] <- label.transmissions(MCMCstate = MCMCstate, 
                                                           labels = labels_permute,
                                                           infector.posterior.probabilities = infector.posterior.probabilities)[[2]]
    }
    
    #find quantiles for each entry of the matrix
    quantiles <- array(dim = c(nLabels, nLabels, length(percentiles)))
    for(i in seq_along(percentiles)){
      quantiles[,,i] <- apply(label.transmissions.null, MARGIN = c(1,2), FUN = quantile, probs = percentiles[i])
    }
    return(list(infector.labels = infector.labels, label.transmissions = label.transmissions,
                label.transmissions.null = label.transmissions.null, quantiles = quantiles))
  } else{
    return(list(infector.labels = infector.labels, label.transmissions = label.transmissions))
  }
}

#' @title Find transmission rates between individuals with certain labels over multiple clusters
#' @description Function to find the expected number of transmission between and within individuals 
#' that belong to certain groups over multiple clusters
#' @param output A list output as from run.biophybreak
#' @param df A dataframe as from prepare.HIV.data
#' @param labelname The name of the column in df to use as the label
#' @param permute_test Whether or not to run a permutation test of how the individuals are labeled 
#' to find a null distribution for the label transmission rates
#' @param nPermute Number of permutations to use in the permutation test
#' @param percentiles Percentiles at which to find quantiles from the null distribution (only affects individual cluster results)
#' @return A list containing the individual cluster outputs from label.transmissions, 
#' the overall matrix of transmission rates, an array of the null transmission rate samples from the permutation test, 
#' the percentiles of each label transmission rate, and the raw p values for
#' @export
#' 
label.transmissions.wrapper <- function(output, df, labelname, 
                                        permute_test = TRUE, nPermute = 10000, percentiles = c(0.025, 0.975)){
  if(!(labelname %in% names(df))) stop("That labelname is not in df")
  #make it load stuff separately to allow larger datasets?
  nClust <- length(output)
  for(i in seq_len(nClust)){
    if(class(output[[i]]$MCMCstate)[1] != "phybreak"){
      stop("output[[i]]$MCMCstate must be a phybreak object for all i")
    }
  }
  
  #make sure labels are a factor
  labels <- as.factor(df[labelname][[1]])
  
  #levels of label factors
  label.levels <- levels(labels)
  
  #number of unique labels
  nLabels <- length(label.levels)
  
  #function to extract label info and call label.transmissions
  lt.miniwrap <- function(x, df, labelname, permute_test, nPermute, percentiles){
    #find number of individuals in cluster
    nInd <- length(unique(x$MCMCstate$d$hostnames))
    #find dataframe indices of individuals in this cluster
    df_indices <- sapply(x$MCMCstate$d$hostnames[1:nInd], 
                         FUN = function(host, df){which(df$patient_ID == host)},
                         df = df)
    label.transmission <- label.transmissions(MCMCstate = x$MCMCstate, 
                                              labels = as.factor(df[labelname][[1]])[df_indices], 
                                              infector.posterior.probabilities = phybreak.infector.posts(x$MCMCstate)$post_support,
                                              permute_test = permute_test,
                                              nPermute = nPermute,
                                              percentiles = percentiles)
  }
  
  lt.wrapped <- parallel::mclapply(output, FUN = lt.miniwrap, 
                                   df = df,
                                   labelname = labelname,
                                   permute_test = permute_test,
                                   nPermute = nPermute,
                                   percentiles = percentiles,
                                   mc.preschedule = FALSE,
                                   mc.cores = parallel::detectCores())
  
  #combine label transmission and permutation test results
  all.label.transmissions <- matrix(0, nrow = dim(lt.wrapped[[1]]$label.transmissions)[1], ncol = dim(lt.wrapped[[1]]$label.transmissions)[1])
  all.label.transmission.null <- array(0, dim = dim(lt.wrapped[[1]]$label.transmissions.null))
  for(i in seq_along(lt.wrapped)){
    all.label.transmissions <- all.label.transmissions + lt.wrapped[[i]]$label.transmissions
    all.label.transmission.null <- all.label.transmission.null + lt.wrapped[[i]]$label.transmissions.null
  }
  
  #find percentiles and p values
  equal_permutes <- vector(mode = "list", length = nLabels^2)
  dim(equal_permutes) <- c(nLabels, nLabels)
  percentiles <- matrix(NA, nrow = nLabels, ncol = nLabels)
  pval_uncor <- matrix(1, nrow = nLabels, ncol = nLabels)
  for(i in seq_len(nLabels)){
    for(j in seq_len(nLabels)){
      null.sorted <- sort(all.label.transmission.null[i,j,])
      equal_permutes[[i,j]] <- which(null.sorted == all.label.transmissions[i,j])
      if(length(equal_permutes[[i,j]]) > 0){
        #print(length(equal_permutes[[i,j]]))
        percentiles[i,j] <- (min(equal_permutes[[i,j]]-1)+max(equal_permutes[[i,j]]-1))/(2*(nPermute-1))
        #print(quant)
        pval_uncor[i,j] <- min(percentiles[i,j], 1-percentiles[i,j])*2
      } else{
        if(all(null.sorted <= all.label.transmissions[i,j])){
          percentiles[i,j] <- 1
          pval_uncor[i,j] <- 2/(nPermute-1)
        } else{
          percentiles[i,j] <- min(which(sort(all.label.transmission.null[i,j,]) > all.label.transmissions[i,j])-1)/(nPermute-1)
          pval_uncor[i,j] <- min(percentiles[i,j], 1-percentiles[i,j])*2
          #don't let p value be 0
          pval_uncor[i,j] <- max(pval_uncor[i,j], 2/(nPermute-1))
        }
      }
    }
  }
  #add names to percentiles and p values
  rownames(percentiles) <- rownames(all.label.transmissions)
  colnames(percentiles) <- colnames(all.label.transmissions)
  rownames(pval_uncor) <- rownames(all.label.transmissions)
  colnames(pval_uncor) <- colnames(all.label.transmissions)
  #add names to null samples
  dimnames(all.label.transmission.null) <- list(to = rownames(all.label.transmissions), 
                                                from = colnames(all.label.transmissions), 
                                                sample = 1:nPermute)
  return(list(lt.wrapped = lt.wrapped, 
              all.label.transmissions = all.label.transmissions, 
              all.label.transmission.null = all.label.transmission.null,
              percentiles = percentiles, pval_uncor = pval_uncor))
}

#'@title Pool smaller groups of individuals together
#'@description Function to add cutoffs to grouping (pool smaller groups of individuals together in an "other" group)
#'@param df A data frame containing infection time distributions as in the output from prepare.HIV.data 
#'@param group_name The name of the column in df to consider
#'@param max_groups The maximum number of groups, including the "other" group
#'@return A vector of the group that each individual belongs to with (potentially) the smallest groups pooled into "other"
group.cutoff <- function(df, group_name, max_groups){
  if(!(group_name %in% names(df))) stop("That group_name is not in df")
  group_table <- table(df[,group_name])
  
  #find cutoff value for pooling groups with smaller numbers of individuals
  if(length(group_table) <= max_groups){
    group_cutoff_count <- 0
  } else{
    group_cutoff_count <- sort(group_table[names(group_table) != "unknown"], decreasing = TRUE)[max_groups]
  }
  
  other_group_names <- c(names(which(group_table <= group_cutoff_count)), "unknown")
  group_cutoff <- df[,group_name]
  group_cutoff[which(group_cutoff %in% other_group_names)] <- "other"
  return(group_cutoff)
}