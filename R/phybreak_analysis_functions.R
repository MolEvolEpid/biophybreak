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
      ggplot2::theme(axis.text.x = ggtext::element_markdown(color = treecolors, angle = angle)) +
      ggplot2::guides(fill = ggplot2::guide_legend(order = 1),
                      color = ggplot2::guide_legend(order = 2, override.aes = list(fill = "#FFFFFF")))
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

#function to generate uninformed random transmission tree
#(infectors chosen uniformly from all possible infectors)
rtrans_tree_unif <- function(nInds){
  infectors.rand <- integer(length = nInds)
  infectors.rand[c(1,2)] <- c(0,1)
  if(nInds > 2){
    for(i in 3:nInds){
      infectors.rand[i] <- sample(i-1,1) #pick random individual from 1 to i-1
    }
  }
  return(infectors.rand)
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
#' @param between.clust.bs Whether or not the bootstrapping should include demographic label randomization between clusters as well as within clusters
#' @param label.levels The labels of the demographic groups of interest. 
#' @param label.probs The overall proportions of labels across the clusters of interest
#' @param quantiles Quantiles at which to find quantiles from the null distribution
#' @return A matrix containing the the probabilities for each individual to be infected by a sampled individual of a certain label
#' and a matrix with the number of transmissions between and within individuals with each label
#' @export
#' 
label.transmissions <- function(MCMCstate, labels, infector.posterior.probabilities = phybreak.infector.posts(MCMCstate)$post_support,
                                permute_test = FALSE, nPermute = 10000, 
                                between.clust.bs = FALSE, label.levels = NULL, label.probs = NULL,
                                quantiles = c(0.025, 0.975)){
  #make sure labels are a factor
  labels <- as.factor(labels)
  #return(label.probs)
  #levels of label factors
  if(is.null(label.levels)) label.levels <- levels(labels)
  
  #number of unique labels
  nLabels <- length(label.levels)
  
  #number of individuals
  nInd <- dim(infector.posterior.probabilities)[2]
  
  #number of MCMC samples
  nSample <- length(MCMCstate$s$logLik)
  
  #labels of infectors for each individual
  infector.labels <- matrix(nrow = nInd, ncol = nLabels)
  for(i in seq_len(nInd)){
    for(j in seq_len(nLabels)){
      infector.labels[i,j] <- infector.posterior.probabilities[2:(nInd+1),i] %*% (labels == label.levels[j])
    }
  }
  
  #matrix for mean number of transmissions between pairs of labels
  label.transmissions.mean <- matrix(0, nrow = nLabels, ncol = nLabels)
  for(i in seq_len(nInd)){
    label.transmissions.mean[which(label.levels == labels[i]),] <- 
      label.transmissions.mean[which(label.levels == labels[i]),] + infector.labels[i,]
  }
  #return(label.transmissions.mean)
  rownames(infector.labels) <- MCMCstate$d$hostnames[1:nInd]
  colnames(infector.labels) <- label.levels
  
  rownames(label.transmissions.mean) <- label.levels
  colnames(label.transmissions.mean) <- label.levels
  
  if(permute_test){
    if(between.clust.bs){
      #draw labels of each individual in each cluster from total label distributions
      labels.one_hot <- rmultinom(n = nPermute*nInd, size = 1, prob = label.probs)
      #convert to name
      labels.vector <- label.levels[apply(labels.one_hot, 2, FUN = function(x) which(x == 1))]
      labels.bcbs <- matrix(labels.vector, nrow = nInd, ncol = nPermute)
    }
    #permuted labels (fixed to exact within-cluster distribution)
    labels.perm <- lapply(seq_len(nPermute), FUN = function(x, labels) sample(labels), labels)
    #draw bootstrap samples from MCMC
    samples <- sample(seq_len(nSample), nPermute, replace = TRUE)
    #find infectors
    infectors <- MCMCstate$s$infectors[,samples]
    #uniformly chosen infectors for null distribution
    infectors.null <- sapply(rep(nInd,nPermute), FUN = rtrans_tree_unif)
    #matrix for number of transmissions between pairs of labels (within-cluster bootstrap)
    label.transmissions.boot <- array(0, dim = c(nLabels, nLabels, nPermute))
    #matrix for number of transmissions between pairs of labels
    label.transmissions.null <- array(0, dim = c(nLabels, nLabels, nPermute))
    #matrix for number of transmissions between pairs of labels (inter-cluster bootstrap)
    if(between.clust.bs) label.transmissions.bcbs <- array(0, dim = c(nLabels, nLabels, nPermute))
    #i'm sure there's a much better way to do this
    for(i in seq_len(nPermute)){
      for(j in seq_len(nInd)){
        if(infectors[j,i] != 0){
          label.transmissions.boot[which(label.levels == labels[j]), which(label.levels == labels[infectors[j,i]]),i] <-
            label.transmissions.boot[which(label.levels == labels[j]), which(label.levels == labels[infectors[j,i]]),i] + 1
        }
        if(infectors.null[j,i] != 0){
          label.transmissions.null[which(label.levels == labels.perm[[i]][j]), which(label.levels == labels.perm[[i]][infectors.null[j,i]]),i] <-
            label.transmissions.null[which(label.levels == labels.perm[[i]][j]), which(label.levels == labels.perm[[i]][infectors.null[j,i]]),i] + 1
        }
        if(between.clust.bs){
          if(infectors.null[j,i] != 0){
            label.transmissions.bcbs[which(label.levels == labels.bcbs[j,i]), which(label.levels == labels.bcbs[[i]][infectors.null[j,i]]),i] <-
              label.transmissions.bcbs[which(label.levels == labels.bcbs[j,i]), which(label.levels == labels.bcbs[[i]][infectors.null[j,i]]),i] + 1
          }
        }
      }
      #labels_permute <- sample(labels, replace = FALSE)
      #label.transmissions.null[,,i] <- label.transmissions(MCMCstate = MCMCstate, 
      #                                                     labels = labels_permute,
      #                                                     infector.posterior.probabilities = infector.posterior.probabilities)[[2]]
    }
    
    #find quantiles for each entry of the matrix
    quantiles.boot <- array(dim = c(nLabels, nLabels, length(quantiles)))
    quantiles.null <- array(dim = c(nLabels, nLabels, length(quantiles)))
    if(between.clust.bs) quantiles.bcbs <- array(dim = c(nLabels, nLabels, length(quantiles)))
    for(i in seq_along(quantiles)){
      quantiles.boot[,,i] <- apply(label.transmissions.boot, MARGIN = c(1,2), FUN = quantile, probs = quantiles[i])
      quantiles.null[,,i] <- apply(label.transmissions.null, MARGIN = c(1,2), FUN = quantile, probs = quantiles[i])
      if(between.clust.bs){
        quantiles.bcbs[,,i] <- apply(label.transmissions.bcbs, MARGIN = c(1,2), FUN = quantile, probs = quantiles[i])
      }
    }
    if(between.clust.bs){
      return(list(infector.labels = infector.labels, 
                  label.transmissions.mean = label.transmissions.mean,
                  label.transmissions.boot = label.transmissions.boot,
                  label.transmissions.null = label.transmissions.null, 
                  label.transmissions.bcbs = label.transmissions.bcbs,
                  quantiles.boot = quantiles.boot,
                  quantiles.null = quantiles.null,
                  quantiles.bcbs = quantiles.bcbs))
    }
    return(list(infector.labels = infector.labels, 
                label.transmissions.mean = label.transmissions.mean,
                label.transmissions.boot = label.transmissions.boot,
                label.transmissions.null = label.transmissions.null, 
                quantiles.boot = quantiles.boot,
                quantiles.null = quantiles.null))
  } else{
    return(list(infector.labels = infector.labels, label.transmissions.mean = label.transmissions.mean))
  }
}

#' @title Find transmission rates between individuals with certain labels over multiple clusters
#' @description Function to find the expected number of transmission between and within individuals 
#' that belong to certain groups over multiple clusters
#' @param output A list output as from run.biophybreak. Does not need to include the 'MCMCstate's if filenames and filepath are provided.
#' For large datasets, use of filenames and filepath is recommended instead
#' @param filenames A vector of file names of output from biophybreak
#' @param filepath The path to the biophybreak MCMC output (should include trailing "/")
#' @param df A dataframe as from prepare.HIV.data
#' @param labelname The name of the column in df to use as the label
#' @param permute_test Whether or not to run a permutation test of how the individuals are labeled 
#' to find a null distribution for the label transmission rates
#' @param nPermute Number of permutations to use in the permutation test
#' @param between.clust.bs Whether or not the bootstrapping should include demographic label randomization between clusters as well as within clusters
#' @param quantiles Quantiles at which to find quantiles from the null distribution (only affects individual cluster results)
#' @param max.cores The maximum number of cores to use for parallel processing
#' @return A list containing the following elements:
#' \itemize{
#'    \item \code{lt.wrapped} individual cluster outputs from label.transmissions
#'    \item \code{all.label.transmissions.mean} the overall matrix of transmission rates
#'    \item \code{all.label.transmissions.null.mean} the matrix of the mean transmission rates in the null distribution using only the clustering information
#'    \item \code{all.label.transmissions.bcbs.mean} the matrix of the mean transmission rates in the null distribution using only the overall demographic label proportions
#'    \item \code{all.label.transmissions.quantile} the transmission rates at the provided quantiles
#'    \item \code{all.label.transmissions.null.quantile} the quantile values of transmission rates in the null distribution using only the clustering information
#'    \item \code{all.label.transmissions.bcbs.quantile} the quantile values of transmission rates in the null distribution using only the overall demographic label proportions
#'    \item \code{all.label.transmissions.boot} all bootstrapped values for the actual group transmission rates using clustering and transmission history information
#'    \item \code{all.label.transmissions.null} all bootstrapped values for the null transmission rates using only the within-cluster demographic label information
#'    \item \code{all.label.transmissions.bcbs} all bootstrapped values for the null transmission rates using only the overall demographic label proportions
#'    \item \code{percentiles} the percentile value of \code{all.label.transmissions.mean} with respect to \code{all.label.transmissions.null}
#'    \item \code{pval_uncor} Conversion of \code{percentiles} to a two-tailed p-value
#'    \item \code{percentiles_bc} the percentile value of \code{all.label.transmissions.mean} with respect to \code{all.label.transmissions.bcbs}
#'    \item \code{pval_uncor_bc} Conversion of \code{percentiles_bc} to a two-tailed p-value
#'    \item \code{percentiles_null_bc} the percentile value of \code{all.label.transmissions.null.mean} with respect to \code{all.label.transmissions.bcbs}
#'    \item \code{pval_uncor_null_bc} Conversion of \code{percentiles_null_bc} to a two-tailed p-value
#' }
#' @export
#' 
label.transmissions.wrapper <- function(output, filenames = NULL, filepath = NULL, 
                                        df, 
                                        labelname, 
                                        permute_test = TRUE, nPermute = 10000, 
                                        between.clust.bs  = FALSE,
                                        quantiles = c(0.025, 0.975),
                                        max.cores = parallel::detectCores()){
  
  if(any(sapply(MCMCs_metadata, FUN = function(x) is.null(x$MCMCstate)))){
    if(is.null(filenames) || is.null(filepath)){
      stop("Either 'MCMCstate' must be present in output or filename and filepath must be provided")
    }
  }
  if(!(labelname %in% names(df))) stop("That labelname is not in df")
  
  #find number of clusters
  nClust <- length(output)
  
  #find IDs of included individuals
  pat_ids <- unlist(sapply(output, FUN = function(x) unique(x$inputs$biophybreakdata$sample.hosts)))
  
  #find indices in df of included individuals
  pat_indices <- which(df$patient_ID %in% pat_ids)
  
  #subset version of data frame
  df.sub <- df[pat_indices,]
  
  #make sure labels are a factor
  labels <- as.factor(df.sub[labelname][[1]])
  
  #levels of label factors
  label.levels <- levels(labels)
  
  #number of unique labels
  nLabels <- length(label.levels)
  
  #find totals of each kind of label
  label.totals <- table(labels)
  
  #find probabilities of each label
  label.probs <- label.totals/sum(label.totals)
  
  #find distribution of cluster sizes
  clust.sizes <- sapply(output, FUN = function(x) x$inputs$nInds)
  
  #function to extract label info and call label.transmissions
  #' @title Intermediary wrapper for label.transmissions
  #' @description
  #' Function to aid in parallelization of label.transmissions by extracting label info and calling label.transmissions.
  #' Not intended for use outside of label.transmission.wrapper.
  #' @export
  lt.miniwrap <- function(x = NULL, filename = NULL, filepath = NULL, df, labelname, 
                          permute_test, nPermute, 
                          between.clust.bs, label.levels, label.probs,
                          quantiles){
    #load MCMC output if not provided
    if(is.null(x$MCMCstate)){
      if(!is.null(filename) && !is.null(filepath)){
        load(paste0(filepath, filename))
        x <- output
      } else{
        stop("Either MCMC output or filename and filepath must be provided")
      }
    }
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
                                              between.clust.bs = between.clust.bs, 
                                              label.levels = label.levels, 
                                              label.probs = label.probs,
                                              quantiles = quantiles)
  }
  
  #setup for parallel processing
  cl <- parallel::makePSOCKcluster(min(max.cores, max(length(output), length(filenames))))
  parallel::clusterEvalQ(cl, library(biophybreak))
  
  lt.wrapped <- parallel::clusterMap(cl,
                                     fun = lt.miniwrap, 
                                     filename = filenames,
                                     filepath = filepath,
                                     df = list(df.sub),
                                     labelname = list(labelname),
                                     permute_test = permute_test,
                                     nPermute = nPermute,
                                     between.clust.bs = between.clust.bs, 
                                     label.levels = list(label.levels), 
                                     label.probs = list(label.probs),
                                     quantiles = list(quantiles))
  parallel::stopCluster(cl)
  #return(lt.wrapped)
  #combine label transmission and permutation test results
  all.label.transmissions.mean <- matrix(0, nrow = dim(lt.wrapped[[1]]$label.transmissions.mean)[1], 
                                         ncol = dim(lt.wrapped[[1]]$label.transmissions.mean)[1])
  if(permute_test){
    all.label.transmissions.boot <- array(0, dim = dim(lt.wrapped[[1]]$label.transmissions.boot))
    all.label.transmissions.null <- array(0, dim = dim(lt.wrapped[[1]]$label.transmissions.null))
    if(between.clust.bs) all.label.transmissions.bcbs <- array(0, dim = dim(lt.wrapped[[1]]$label.transmissions.bcbs))
  }
  for(i in seq_along(lt.wrapped)){
    all.label.transmissions.mean <- all.label.transmissions.mean + lt.wrapped[[i]]$label.transmissions.mean
    if(permute_test){
      all.label.transmissions.boot <- all.label.transmissions.boot + lt.wrapped[[i]]$label.transmissions.boot
      all.label.transmissions.null <- all.label.transmissions.null + lt.wrapped[[i]]$label.transmissions.null
      if(between.clust.bs) all.label.transmissions.bcbs <- all.label.transmissions.bcbs + lt.wrapped[[i]]$label.transmissions.bcbs
    }
  }
  
  if(permute_test){
    #find percentiles and p values
    equal_permutes <- vector(mode = "list", length = nLabels^2)
    dim(equal_permutes) <- c(nLabels, nLabels)
    percentiles <- matrix(NA, nrow = nLabels, ncol = nLabels)
    pval_uncor <- matrix(1, nrow = nLabels, ncol = nLabels)
    #find quantiles of bootstrapped transmission rates (using all info)
    all.label.transmissions.quantile <- apply(all.label.transmissions.boot, MARGIN = c(1,2), FUN = quantile, probs = quantiles)
    #find means and quantiles of nulls keeping within-cluster label distributions
    all.label.transmissions.null.mean <- apply(all.label.transmissions.null, MARGIN = c(1,2), FUN = mean)
    all.label.transmissions.null.quantile <- apply(all.label.transmissions.null, MARGIN = c(1,2), FUN = quantile, probs = quantiles)
    if(between.clust.bs){
      equal_permutes_bc <- vector(mode = "list", length = nLabels^2)
      dim(equal_permutes_bc) <- c(nLabels, nLabels)
      percentiles_bc <- matrix(NA, nrow = nLabels, ncol = nLabels)
      pval_uncor_bc <- matrix(1, nrow = nLabels, ncol = nLabels)
      equal_permutes_null_bc <- vector(mode = "list", length = nLabels^2)
      dim(equal_permutes_null_bc) <- c(nLabels, nLabels)
      percentiles_null_bc <- matrix(NA, nrow = nLabels, ncol = nLabels)
      pval_uncor_null_bc <- matrix(1, nrow = nLabels, ncol = nLabels)
      #find means and quantiles of nulls with only overall label distributions
      all.label.transmissions.bcbs.mean <- apply(all.label.transmissions.bcbs, MARGIN = c(1,2), FUN = mean)
      all.label.transmissions.bcbs.quantile <- apply(all.label.transmissions.bcbs, MARGIN = c(1,2), FUN = quantile, probs = quantiles)
    }
    for(i in seq_len(nLabels)){
      for(j in seq_len(nLabels)){
        null.sorted <- sort(all.label.transmissions.null[i,j,])
        equal_permutes[[i,j]] <- which(null.sorted == all.label.transmissions.mean[i,j])
        if(length(equal_permutes[[i,j]]) > 0){
          #print(length(equal_permutes[[i,j]]))
          percentiles[i,j] <- (min(equal_permutes[[i,j]]-1)+max(equal_permutes[[i,j]]-1))/(2*(nPermute-1))
          #print(quant)
          pval_uncor[i,j] <- min(percentiles[i,j], 1-percentiles[i,j])*2
        } else{
          if(all(null.sorted <= all.label.transmissions.mean[i,j])){
            percentiles[i,j] <- 1
            pval_uncor[i,j] <- 2/(nPermute-1)
          } else{
            percentiles[i,j] <- min(which(sort(all.label.transmissions.null[i,j,]) > all.label.transmissions.mean[i,j])-1)/(nPermute-1)
            pval_uncor[i,j] <- min(percentiles[i,j], 1-percentiles[i,j])*2
            #don't let p value be 0
            pval_uncor[i,j] <- max(pval_uncor[i,j], 2/(nPermute-1))
          }
        }
        if(between.clust.bs){ #between cluster measures
          null.sorted_bc <- sort(all.label.transmissions.bcbs[i,j,])
          equal_permutes_bc[[i,j]] <- which(null.sorted_bc == all.label.transmissions.mean[i,j])
          equal_permutes_null_bc[[i,j]] <- which(null.sorted_bc == all.label.transmissions.null.mean[i,j])
          if(length(equal_permutes_bc[[i,j]]) > 0){
            #print(length(equal_permutes[[i,j]]))
            percentiles_bc[i,j] <- (min(equal_permutes_bc[[i,j]]-1)+max(equal_permutes_bc[[i,j]]-1))/(2*(nPermute-1))
            #print(quant)
            pval_uncor_bc[i,j] <- min(percentiles_bc[i,j], 1-percentiles_bc[i,j])*2
          } else{
            if(all(null.sorted_bc <= all.label.transmissions.mean[i,j])){
              percentiles_bc[i,j] <- 1
              pval_uncor_bc[i,j] <- 2/(nPermute-1)
            } else{
              percentiles_bc[i,j] <- min(which(sort(all.label.transmissions.bcbs[i,j,]) > all.label.transmissions.mean[i,j])-1)/(nPermute-1)
              pval_uncor_bc[i,j] <- min(percentiles_bc[i,j], 1-percentiles_bc[i,j])*2
              #don't let p value be 0
              pval_uncor_bc[i,j] <- max(pval_uncor_bc[i,j], 2/(nPermute-1))
            }
          }
          if(length(equal_permutes_null_bc[[i,j]]) > 0){
            #print(length(equal_permutes[[i,j]]))
            percentiles_null_bc[i,j] <- (min(equal_permutes_bc[[i,j]]-1)+max(equal_permutes_bc[[i,j]]-1))/(2*(nPermute-1))
            #print(quant)
            pval_uncor_bc[i,j] <- min(percentiles_bc[i,j], 1-percentiles_bc[i,j])*2
          } else{
            if(all(null.sorted_bc <= all.label.transmissions.null.mean[i,j])){
              percentiles_null_bc[i,j] <- 1
              pval_uncor_null_bc[i,j] <- 2/(nPermute-1)
            } else{
              percentiles_null_bc[i,j] <- min(which(sort(all.label.transmissions.bcbs[i,j,]) > all.label.transmissions.null.mean[i,j])-1)/(nPermute-1)
              pval_uncor_null_bc[i,j] <- min(percentiles_null_bc[i,j], 1-percentiles_null_bc[i,j])*2
              #don't let p value be 0
              pval_uncor_null_bc[i,j] <- max(pval_uncor_bc[i,j], 2/(nPermute-1))
            }
          }
        }
      }
    }
    #add names to percentiles, p values, and null means
    rownames(percentiles) <- rownames(all.label.transmissions.mean)
    colnames(percentiles) <- colnames(all.label.transmissions.mean)
    rownames(pval_uncor) <- rownames(all.label.transmissions.mean)
    colnames(pval_uncor) <- colnames(all.label.transmissions.mean)
    rownames(all.label.transmissions.null.mean) <- rownames(all.label.transmissions.mean)
    colnames(all.label.transmissions.null.mean) <- colnames(all.label.transmissions.mean)
    if(between.clust.bs){
      #actual compared to between cluster bootstrapping
      rownames(percentiles_bc) <- rownames(all.label.transmissions.mean)
      colnames(percentiles_bc) <- colnames(all.label.transmissions.mean)
      rownames(pval_uncor_bc) <- rownames(all.label.transmissions.mean)
      colnames(pval_uncor_bc) <- colnames(all.label.transmissions.mean)
      #clustering info compared to between cluster bootstrapping
      rownames(percentiles_null_bc) <- rownames(all.label.transmissions.mean)
      colnames(percentiles_null_bc) <- colnames(all.label.transmissions.mean)
      rownames(pval_uncor_null_bc) <- rownames(all.label.transmissions.mean)
      colnames(pval_uncor_null_bc) <- colnames(all.label.transmissions.mean)
      #means
      rownames(all.label.transmissions.bcbs.mean) <- rownames(all.label.transmissions.mean)
      colnames(all.label.transmissions.bcbs.mean) <- colnames(all.label.transmissions.mean)
    }
    #add names to bootstrapped and null samples
    dimnames(all.label.transmissions.boot) <- list(to = rownames(all.label.transmissions.mean), 
                                                   from = colnames(all.label.transmissions.mean), 
                                                   sample = 1:nPermute)
    dimnames(all.label.transmissions.null) <- list(to = rownames(all.label.transmissions.mean), 
                                                   from = colnames(all.label.transmissions.mean), 
                                                   sample = 1:nPermute)
    if(between.clust.bs){
      dimnames(all.label.transmissions.bcbs) <- list(to = rownames(all.label.transmissions.mean), 
                                                     from = colnames(all.label.transmissions.mean), 
                                                     sample = 1:nPermute)
    }
    if(between.clust.bs){
      return(list(lt.wrapped = lt.wrapped, 
                  all.label.transmissions.mean = all.label.transmissions.mean, 
                  all.label.transmissions.null.mean = all.label.transmissions.null.mean, 
                  all.label.transmissions.bcbs.mean = all.label.transmissions.bcbs.mean, 
                  all.label.transmissions.quantile = all.label.transmissions.quantile,
                  all.label.transmissions.null.quantile = all.label.transmissions.null.quantile,
                  all.label.transmissions.bcbs.quantile = all.label.transmissions.bcbs.quantile,
                  all.label.transmissions.boot = all.label.transmissions.boot,
                  all.label.transmissions.null = all.label.transmissions.null,
                  all.label.transmissions.bcbs = all.label.transmissions.bcbs,
                  percentiles = percentiles, pval_uncor = pval_uncor,
                  percentiles_bc = percentiles_bc, pval_uncor_bc = pval_uncor_bc,
                  percentiles_null_bc = percentiles_null_bc, pval_uncor_null_bc = pval_uncor_null_bc))
    }
    return(list(lt.wrapped = lt.wrapped, 
                all.label.transmissions.mean = all.label.transmissions.mean, 
                all.label.transmissions.null.mean = all.label.transmissions.null.mean,
                all.label.transmissions.quantile = all.label.transmissions.quantile,
                all.label.transmissions.null.quantile = all.label.transmissions.null.quantile,
                all.label.transmissions.boot = all.label.transmissions.boot,
                all.label.transmissions.null = all.label.transmissions.null,
                percentiles = percentiles, pval_uncor = pval_uncor))
  } else{
    return(list(lt.wrapped = lt.wrapped, 
                all.label.transmissions.mean = all.label.transmissions.mean))
  }
}

#'@title Pool smaller groups of individuals together
#'@description Function to add cutoffs to grouping (pool smaller groups of individuals together in an "other" group)
#'@param df A data frame containing infection time distributions as in the output from prepare.HIV.data 
#'@param group_name The name of the column in df to consider
#'@param max_groups The maximum number of groups, including the "other" group
#'@return A vector of the group that each individual belongs to with (potentially) the smallest groups pooled into "other"
#'@export
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

#
#' @title Plot group transmission rates
#' @description Function to make violin plots of transmission rates 
#' between individuals with certain labels over multiple clusters
#' @param transmission.list A list output as from "label.transmission.wrapper"
#' @param include.null Whether or not to include the null distributions of the transmission rates given the within-cluster label demographics
#' @param bcbs Whether or not to include the null distributions of the transmission rates given only the overall label demographics
#' @param title A string for the title of the plot
#' @export
transmission.rate.violins <- function(transmission.list, include.null = FALSE, bcbs = FALSE, title = NULL){
  rates.long <- reshape2::melt(transmission.list$all.label.transmissions.mean)
  names(rates.long) <- c("to", "from", "value")
  boot.long <- reshape2::melt(transmission.list$all.label.transmissions.boot)
  null.long <- reshape2::melt(transmission.list$all.label.transmissions.null)
  if(bcbs) bcbs.long <- reshape2::melt(transmission.list$all.label.transmissions.bcbs)
  #add column for whether it's actual or null
  boot.long$set <- rep("All Information", dim(boot.long)[1])
  null.long$set <- rep("Clustering Only", dim(null.long)[1])
  if(bcbs) bcbs.long$set <- rep("Overall Only", dim(null.long)[1])
  #concatenate actual and null (and between cluster if present)
  df.long <- rbind(boot.long, null.long)
  if(bcbs) df.long <- rbind(df.long, bcbs.long)
  
  if(include.null){
    ggplot2::ggplot(data = df.long, ggplot2::aes(x = to, fill = from, color = set, y = value)) + 
      ggplot2::geom_violin(adjust = 6, scale = "width") + 
      ggplot2::scale_color_manual(values = c("#000000", "#777777", "#EEEEEE")) + 
      ggplot2::xlab("Recipient") + 
      ggplot2::ylab("Number of Transmissions") +
      ggplot2::labs(fill = "Source", color = "Information") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_bw() +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  } else{
    fromTo <- paste(boot.long$from, boot.long$to, sep = " to ")
    boot.long$fromTo <- fromTo
    fromTo_unique <- unique(fromTo)
    empty_indices <- c()
    for(i in seq_along(fromTo_unique)){
      indices <- which(fromTo == fromTo_unique[i])
      if(sum(boot.long$value[indices]) == 0){
        empty_indices <- c(empty_indices,indices)
      }
    }
    if(is.null(empty_indices)){
      df.fromTo <- boot.long
    } else{
      df.fromTo <- boot.long[-empty_indices,]
    }
    
    #decide whether angle should be changed (based on length of fromTo_unique)
    if(length(fromTo_unique) > 8){
      label_angle <- 90
    } else{
      label_angle <- 0
    }
    
    #plot actual numbers of transmissions, omitting categories with 0 transmissions
    ggplot2::ggplot(data = df.fromTo, aes(x = fromTo, y = value, fill = fromTo)) + 
      ggplot2::geom_violin(adjust = 6, scale = "width") +
      ggplot2::xlab("Transmission Type") +
      ggplot2::ylab("Number of Transmissions") +
      ggplot2::ggtitle(title) +
      ggplot2::theme_bw() + 
      ggplot2::scale_fill_discrete(guide = "none") +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = label_angle)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  }
}