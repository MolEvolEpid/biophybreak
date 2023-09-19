#simulation functions for transmission histories for use with phybreak

#' @title Generate a transmission tree
#' @description  Function to generate a simple random transmission tree
#' @param nInd The number of individuals
#' @param time.span The amount of time between the first and last infection.
#' @param gen.dist.type The type of distribution for the generation function 
#' (either "gamma" or "function", a user supplied function)
#' @param gen.shape Shape of the gamma distribution for the generation function, if a gamma function is used.
#' @param gen.mean Mean of the gamma distribution for the generation function, if a gamma function is used.
#' @param gen.dist Generation distribution, if gen.dist.type = "function" is used.
#' @param infection.ages A number or vector for the amount of time between infection and sampling. Must be positive.
#' @param tr_het Type of transmission heterogeneity. Does not directly change transmission heterogeneity, 
#' but makes some transmission heterogeneities more or less likely. 
#' Can be "uniform_fixed" (no relative weight differences), "uniform_rnd" (relative weights are a uniform random variable), 
#' or "power" (weights follow a power law)
#' @param power.base Exponent base for the weights if tr_het = "power" is used.
#' @param sep The minimum amount of time (in years) between any two infections.
#' @param post.sam.trans.weight How much transmission after sampling is penalized. 
#' 1 means no penalty. 0.1 Would be 10 times less likely.
#' @param seed The pseudorandom number generator seed.
#' @return A transmission history with individual names, donors,infection times, and sampling times.
#' @export
#' 
rtrans_tree <- function(nInd = 10, time.span = (nInd-1)*0.5, gen.dist.type = "gamma", gen.shape = 1, gen.mean = 100, gen.dist = NULL,
                        infection.ages = 1.5, tr_het = "uniform_fixed", power.base = 2, sep = .1, 
                        post.sam.trans.rate = 1,
                        seed = sample(2e9, 1)){
  
  #gen.dist.type: if "gamma", use gen.shape and gen.mean; otherwise specify "user.fn" to use a custom function in gen.dist, 
  #such as stepfun(c(0, 0.5), c(0, 3, 1)) (can be improper prior)
  set.seed(seed)
  
  if(nInd < 2) stop("There must be at least 2 individuals")
  
  #roll infection ages
  t_inf <- numeric(length = nInd)
  t_inf[1] <- time.span
  t_inf[nInd] <- 0
  
  if(nInd > 2){
    for(i in 2:(nInd-1)){
      t_inf[i] <- runif(1, 0, time.span)
      while(min(abs(t_inf[i] - t_inf[-i])) < sep){
        t_inf[i] <- runif(1, 0, time.span)
      }
    }
  }
  
  #t_inf <- runif(nInd - 2, 0, time.span) #infection times are uniform 0 to the number of hosts
  #order times from oldest to most recent
  t_inf <- sort(t_inf, decreasing = TRUE)
  
  #change so that first host is infected at time 0 and positive numbers are after that
  t_inf <- max(t_inf) - t_inf
  
  #find sampling times
  t_sam <- t_inf + infection.ages
  
  #set up post sampling transmission rate stuff
  if(!(length(post.sam.trans.rate) == 1 | length(post.sam.trans.rate == nInd))){
    stop("Length of post.sam.trans.rate must be either 1 or the number of individuals")
  } else if(length(post.sam.trans.rate) == 1){
    post.sam.trans.rate <- rep(post.sam.trans.rate, nInd) #make it the same length as the number of individuals
  }
  
  #expected amount of heterogeneity of number of infectors (relative)
  if(tr_het == "uniform_fixed"){
    rel_tr_weight <- rep(1, nInd)
  } else if(tr_het == "uniform_rnd"){
    rel_tr_weight <- runif(nInd)
  } else if(tr_het == "power"){
    rel_tr_weight <- power.base^(rexp(nInd))
    #print(rel_tr_weight)
  }
  
  gen.rate <- gen.shape/gen.mean
  
  donor <- integer(nInd)
  
  donor[1] <- 0 #infected from outside
  donor[2] <- 1 #individual infected 2nd must have been infected by the first individual
  
  if(nInd >= 3){
    for(i in 3:nInd){
      #weights
      if(gen.dist.type == "gamma"){
        dens <- dgamma(t_inf[i] - t_inf[1:(i-1)], shape = gen.shape, rate = gen.rate)
      } else if(gen.dist.type == "function"){
        dens <- gen.dist(t_inf[i] - t_inf[1:(i-1)])
      } else{
        stop("invalid gen.dist.type")
      }
      
      probs <- dens*rel_tr_weight[1:(i-1)]*post.sam.penalty(t_sam[1:(i-1)], t_inf[i], post.sam.trans.rate[1:(i-1)])
      probs <- probs/sum(probs) #normalize
      donor[i] <- sample((i-1), size = 1, prob = probs) #uniformly choose an individual that was already infected to be donor
    }
  }
  
  #names of hosts
  names <- paste0("P.", 1:nInd)
  
  tt <- data.frame(names, donor, t_inf, t_sam, stringsAsFactors = FALSE)
  
  return(tt)
}

post.sam.penalty <- function(t_sam, t_inf, post.sam.trans.rate){
  n <- length(t_sam)
  rel.rate <- numeric(n)
  for(i in 1:n){
    if(t_inf < t_sam[i]){
      rel.rate[i] <- 1
    } else if(t_inf >= t_sam[i]){
      rel.rate[i] <- post.sam.trans.rate[1]
    }
  }
  return(rel.rate)
}

post.sam.trans.adjust <- function(t_sam, t_inf, post.sam.trans.rate){
  if(t_sam > t_inf){
    rate <- 1
  } else if(t_sam <= t_inf){
    rate <- post.sam.trans.rate
  } 
  return(rate)
}

#' @title Simulate a phyreakdata transmission history
#' @description Function to simulate phylogenic tree and sequences in phybreakdata format from a transmission histoey
#' @param tt A transmission tree data frame. Must contain the elements "names", "donor", "t_inf", and "t_sam" (see below)
#' May be left blank if these are provided separately.
#' @param a A vector for the effective population size (inverse of coalescence rate) for each patient at the time of infection.
#' A value of 0 corresponds to a complete transmission bottleneck (only 1 lineage transmitted), 
#' with larger values indicating a wider transmission bottleneck (potentially multiple sampled lineages transmitted)
#' Typical values for a wide bottleneck may be around 20.
#' @param b A vector for the linear growth rate of the effective population size per generation (usually 1.5 days for HIV). 
#' Faster growth rates indicated faster growth of effective population size, 
#' resulting in lower coalescence rates near sampling (longer branch lengths).
#' Typical values are around 3.
#' @param rhoD A vector for the rates for REVERSE TIME migration from donor to recipient during an ongoing transmission window.
#' If this is non-zero, it indicates that it is possible for a recipient to transmit a lineage back to the individual who infected them.
#' This is typically set to 0. 
#' This value does not do anything for individuals who did not infect anyone else.
#' @param rhoR A vector for the rates of REVERSE TIME migration from recipient to donor during an ongoing transmission window.
#' If this is non-zero, it indicates that the recipient can be infected by the donor at multiple times during an ongoing transmission window.
#' If ongoing transmission is used, typical values are 0.1
#' This value does not do anything for the index case, as their infector is not simulated.
#' @param nSamples A vector for the number of sequences from each patient.
#' Can be left blank if sample.times is provided.
#' If neither nSamples nor sample.times are provided, it will default to one sequence per individual.
#' If nSamples is provided but not sample.times, every sequence will be assumed to come from t_sam for that individual
#' @param sample.times A list of vectors for the sample times for each individual.
#' For example, sample.times = list(c(2.2, 2.3, 2.7), c(2.5, 3.5)) 
#' would mean that the first patient has one sequence from time 2.2, one from 2.3, and one from 2.7
#' and the second patient has one sequence from 2.5 and one from 3.5.
#' @param gen.rate The generation rate in number of generations per year. 
#' Defaults to 365/1.5 (One generation is approximately 1.5 days)
#' @param tr_window A vector for the length of the transmission window a recipient has with their donor, 
#' used when rhoD or rhoR is non-zero, given in years after the recient is infected.
#' Defaults to infinity, indicating that the ongoing transmission windowis open indefinitely.
#' Has no effect if rhoD and rhoR are zero.
#' Has no effect for the index case.
#' @param mut_model The mutation model to be used if seq-gen is used for sequences, 
#' which can be "env" or "pol" for the HIV-1 envelope or polymerase genes
#' @param user_seqs User provided sequences to be used if seq-gen is not used, in order of sampling time
#' @param user_seq_hosts The hosts of the user-provided sequences
#' @param user_seq_names The names of the user-provided sequences
#' @param anc_seq The ancestral sequence for seq-gen to use as the root state
#' @param sim_index The name of the tree for saving files. 
#' @param tree_loc The directory in which to output saved tree files.
#' Defaults to an empty character, resulting in saving in the current directory.
#' @return A phybreakdata object representing the transmission history, sample times, phylogeny, and sequences
#' @export
#' 
sim.coal.phybreak <- function(tt, 
                              a = rep(0, dim(tt)[1]), b = rep(3, dim(tt)[1]), 
                              rhoD = rep(0, dim(tt)[1]), rhoR = rep(0, dim(tt)[1]), 
                              nSamples = rep(1, dim(tt)[1]),
                              sample.times = NULL,
                              gen.rate = 365/1.5,
                              tr_window = rep(Inf, dim(tt)[1]), #can transmit after sampling
                              mut_model = "env",
                              user_seqs = NULL,
                              user_seq_hosts = NULL,
                              user_seq_names = NULL,
                              anc_seq = NULL,
                              sim_index = 1,
                              tree_loc = "", #location to store the tree and sequence files
                              seed = sample(2e9, 1)){
  #find number of individuals
  nInd <- length(tt$names)
  
  #transmission tree (and other parameters) must be ordered from oldest to most recent infection date
  if(is.unsorted(tt$t_inf)){
    tt_old <- tt
    new_order <- order(tt$t_inf)
    tt <- tt[new_order,]
    a <- a[new_order]
    b <- b[new_order]
    rhoR <- rhoR[new_order]
    rhoD <- rhoD[new_order]
    nSamples <- nSamples[new_order]
    sample.times <- sample.times[new_order]
    tr_window <- tr_window[new_order]
    donor_shuff <- tt$donor
    for(i in 1:nInd){
      if(tt$donor[i] != 0){
        tt$donor[i] <- which(tt$names == tt_old$names[donor_shuff[i]])
      }
    }
  }
  
  #call sim.coal.tree.R to make the phylogenetic tree
  treeAndTips <- sim.coal.tree(tt = tt,
                               a = a, b = b, 
                               rhoD = rhoD, rhoR = rhoR,
                               nSamples = nSamples, 
                               sample.times = sample.times,
                               gen.rate = gen.rate,
                               tr_window = tr_window,
                               file_path = tree_loc,
                               tree_name = sim_index,
                               save_tree = FALSE,
                               plot_tree = FALSE,
                               seed = seed)
  tree <- treeAndTips[1]
  sample.names.tips <- simplify2array(strsplit(treeAndTips[2:length(treeAndTips)], split = " "))[1,]
  
  if(is.null(user_seqs)){
    #make sequences with seq-gen if none are provided
    #tree_file <- paste0(tree_loc, "trees/tree_newick_", nInd, "_", sim_index, "_", seed, ".tree")
    #seq_file <- paste0(tree_loc, "trees/seqs_", nInd, "_", sim_index, "_", seed, ".dat")
    #write(tree, file = tree_file) 
    #modify tree file to include ancestral sequence, if provided
    if(!is.null(anc_seq)){
      #find length of ancestral sequences
      seq.length <- length(strsplit(anc_seq, split = "")[[1]])
      #include ancestral sequence and other formatting needed with it
      tree <- paste(paste0("1 ", seq.length), 
                    paste0("Anc ", anc_seq), 
                    1, 
                    tree, 
                    sep = "\n")
      #flag for including ancestral sequence
      anc.seq.flag <- " -k1"
    } else{
      #default to 1000 if no sequence is given
      seq.length <- 1000
      #do not include flag to use ancestral sequence
      anc.seq.flag <- ""
    }
    if(mut_model == "env"){
      #envelope values
      #system(paste0("./seq-gen -m GTR -l 1000 -a .38 -i .3 -s .0083",
      #              " -r 1.71, 2.88, .535, .359, 2.37, .677", #upper triangle of rate matrix
      #              " -f .4627,.1474,.1598,.2302 -z ", seed, #nucleotide frequencies
      #              " < ", tree_file, " > ", seq_file))
      seqs_phyclust <- phyclust::seqgen(opts = paste0("-mGTR -l", seq.length, " -g4 -a0.38 -i0.3 -s0.0083",  anc.seq.flag,
                                                      " -r1.71,2.88,0.535,0.359,2.37,0.677", #upper triangle of rate matrix
                                                      " -f0.4627,0.1474,0.1598,0.2302 -or"),
                                        newick.tree = tree)
    } else if(mut_model == "pol"){
      #pol values
      #system(paste0("./seq-gen -m GTR -l 1000 -g 4 -a .015 -i .255 -s .00241",
      #              " -r 1.09637, 23.97305, 0.46634, 0.00023, 13.74563, 1.00000", #upper triangle of rate matrix
      #              " -f 0.39110, 0.17228,0.21614,0.22048 -z ", seed, #nucleotide frequencies
      #              " < ", tree_file, " > ", seq_file))
      seqs_phyclust <- phyclust::seqgen(opts = paste0("-mGTR -l", seq.length, " -g4 -a0.015 -i0.255 -s0.00241",  anc.seq.flag,
                                                      " -r1.09637,23.97305,0.46634,0.00023,13.74563,1.00000", #upper triangle of rate matrix
                                                      " -f0.39110,0.17228,0.21614,0.22048 -or"),
                                        newick.tree = tree)
      #seqs_phyclust <- phyclust::seqgen(opts = paste0("-mGTR -l1000 -g4 -a0.015 -i0.255 -s0.00241",
      #                                                " -r1.09637,23.97305,0.46634,0.00023,13.74563,1.00000", #upper triangle of rate matrix
      #                                                " -f0.39110,0.17228,0.21614,0.22048 -or"),
      #                                  newick.tree = tree)
    }
    #system(paste0("./seq-gen -m GTR -l 1000 -s .001 -f .389,.165,.228,.219 -z ", seed, " < ", tree_file, " > ", seq_file))
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3511177/
    
    #seqs <- as.matrix(read.table(seq_file))
    #remove first row
    #seqs <- seqs[-1,]
    
    #return(list(seqs, seqs_phyclust))
    #return(seqs_phyclust)
    seqs <- simplify2array(strsplit(seqs_phyclust[-1], split = " "))[2,]
    sequences <- tolower(t(simplify2array(lapply(seqs, FUN = function(x) strsplit(x, split = "")[[1]]))))
    sample.names.seq <- simplify2array(strsplit(seqs_phyclust[-1], split = " "))[1,]
    
    extract_hostname <- function(seq_name){
      seq_name_components <- strsplit(seq_name, split = "_")[[1]]
      nComponents <- length(seq_name_components)
      hostname <- paste(seq_name_components[1:(nComponents-1)], collapse = "_")
      return(hostname)
    }
    hosts <- sapply(sample.names.seq, FUN = extract_hostname)
    
    #sequences <- tolower(t(simplify2array(strsplit(seqs[,2], split = ""))))
    #hosts <- simplify2array(strsplit(seqs[,1], split = "_"))[1,]
    #sample.names.seq <- seqs[,1]
  } else{
    #use provided sequences
    sequences <- user_seqs
    hosts <- user_seq_hosts
    #name samples
    sample.names.seq <- user_seq_names
    #unique.hosts <- unique(hosts)
    #nsamples_total <- 0
    #for(i in unique.hosts){
    #  nsamplesi <- length(which(hosts == i))
    #  sample.names[which(hosts == i)] <- paste(sample.names[which(hosts == i)], nsamples_total+(1:nsamplesi), sep = "_")
    #  nsamples_total <- nsamples_total +nsamplesi
    #}
  }
  
  
  #hosts <- simplify2array(strsplit(seqs[,1], split = "_"))[1,]
  #host_index <- as.integer(simplify2array(strsplit(hosts, split = "[.]"))[2,])
  host_index <- integer(length(hosts))
  for(i in 1:length(hosts)){
    host_index[i] <- which(hosts[i] == tt$names)
  }
  
  if(is.list(sample.times)){
    sample.times.unlist <- unlist(sample.times)
    sample.times.reorder <- numeric(length = length(sample.times.unlist))
    for(i in 1:length(sample.times.unlist)){
      sample.times.reorder[i] <- sample.times.unlist[which(sample.names.tips == sample.names.seq[i])]
    }
    sample.times <- sample.times.reorder
  } else{
    sample.times <- tt$t_sam[host_index]
  }
  host.names <- hosts
  #sim.infection.times <- tt$t_inf
  #sim.infectors <- tt$donor
  
  #order transmission tree by sample time
  tt_old <- tt
  tt <- tt[order(tt$t_sam),]
  donor_shuff <- tt$donor
  for(i in 1:nInd){
    if(tt$donor[i] != 0){
      tt$donor[i] <- which(tt$names == tt_old$names[donor_shuff[i]])
    }
  }
  
  #make ordering transformations (borrowed from phybreakdata)
  allhosts <- unique(host.names)
  allfirsttimes <- rep(FALSE, length(sample.times))
  sapply(allhosts, function(x) allfirsttimes[which(min(sample.times[host.names == x]) == sample.times & (host.names == x))[1]] <<- TRUE)
  outputorderhosts <- order(sample.times[allfirsttimes])
  orderedhosts <- host.names[allfirsttimes][outputorderhosts]
  outputordersamples <- order(!allfirsttimes, match(host.names, orderedhosts), sample.times)
  sequences <- sequences[outputordersamples, ]
  sample.times <- sample.times[outputordersamples]
  sample.names.seq <- sample.names.seq[outputordersamples]
  host.names <- host.names[outputordersamples]
  
  new.order <- integer(length = nInd)
  for(i in 1:nInd){
    new.order[i] <- which(tt$names == orderedhosts[i])
  }
  tt_old <- tt
  tt <- tt[new.order,]
  donor_shuff <- tt$donor
  for(i in 1:nInd){
    if(tt$donor[i] != 0){
      tt$donor[i] <- which(tt$names == tt_old$names[donor_shuff[i]])
    }
  }
  
  sim_data <- phybreakdata(sequences = sequences, 
                           sample.times = sample.times, 
                           sample.names = sample.names.seq, 
                           host.names = host.names,
                           sim.infection.times = tt$t_inf,
                           sim.infectors = tt$donor, 
                           sim.tree = ape::read.tree(text = tree))
  #plotPhyloTrans(sim_data)
  
  return(sim_data)
}

#' @title Simulate a coalescence tree
#' @description Simulation of a coalescence tree given a transmission tree
#' @param tt A transmission tree data frame. Must contain the elements "names", "donor", "t_inf", and "t_sam" (see below)
#' May be left blank if these are provided separately.
#' @param names A vector of the names of the patients
#' @param donor A vector for the infectors for each patient. 0 indicated the index case (the infector was not sampled) 
#' and all other values must be positive integers corresponding to each individual's infector.
#' @param t_inf The infection times of each of the patients in years, usually relative to the infection time of the index patient.
#' @param a A vector for the effective population size (inverse of coalescence rate) for each patient at the time of infection.
#' A value of 0 corresponds to a complete transmission bottleneck (only 1 lineage transmitted), 
#' with larger values indicating a wider transmission bottleneck (potentially multiple sampled lineages transmitted)
#' Typical values for a wide bottleneck may be around 20.
#' @param b A vector for the linear growth rate of the effective population size per generation (usually 1.5 days for HIV). 
#' Faster growth rates indicated faster growth of effective population size, 
#' resulting in lower coalescence rates near sampling (longer branch lengths).
#' Typical values are around 3.
#' @param rhoD A vector for the rates for REVERSE TIME migration from donor to recipient during an ongoing transmission window.
#' If this is non-zero, it indicates that it is possible for a recipient to transmit a lineage back to the individual who infected them.
#' This is typically set to 0. 
#' This value does not do anything for individuals who did not infect anyone else.
#' @param rhoR A vector for the rates of REVERSE TIME migration from recipient to donor during an ongoing transmission window.
#' If this is non-zero, it indicates that the recipient can be infected by the donor at multiple times during an ongoing transmission window.
#' If ongoing transmission is used, typical values are 0.1
#' This value does not do anything for the index case, as their infector is not simulated.
#' @param nSamples A vector for the number of sequences from each patient.
#' Can be left blank if sample.times is provided.
#' If neither nSamples nor sample.times are provided, it will default to one sequence per individual.
#' If nSamples is provided but not sample.times, every sequence will be assumed to come from t_sam for that individual
#' @param sample.times A list of vectors for the sample times for each individual.
#' For example, sample.times = list(c(2.2, 2.3, 2.7), c(2.5, 3.5)) 
#' would mean that the first patient has one sequence from time 2.2, one from 2.3, and one from 2.7
#' and the second patient has one sequence from 2.5 and one from 3.5.
#' @param gen.rate The generation rate in number of generations per year. 
#' Defaults to 365/1.5 (One generation is approximately 1.5 days)
#' @param tr_window A vector for the length of the transmission window a recipient has with their donor, 
#' used when rhoD or rhoR is non-zero, given in years after the recient is infected.
#' Defaults to infinity, indicating that the ongoing transmission windowis open indefinitely.
#' Has no effect if rhoD and rhoR are zero.
#' Has no effect for the index case.
#' @param file_path The directory in which to output saved tree files and plots.
#' Defaults to an empty character, resulting in saving in the current directory.
#' @param tree_name The name of the tree for saving files. 
#' Defaults to 1, resulting in the file "tree1.txt" (if save_tree == TRUE and tree_name_digits is 3)
#' @param tree_name_digits The number of characters in the tree name. 
#' Defaults to the number of characters in the given tree name, so the name with both defaults for name and digits 
#' would be "tree1.txt". If 3 is supplied, the name would be tree001.txt. This does not do anything if tree_name is a character.
#' @param save_tree if TRUE, it will save the newick tree along with the tip labels and states in a text file
#' @param plot_tree If TRUE, this will plot the coalescent tree the pdf file "treeplot(tree_name).pdf"
#' @param plot_colors The colors to use for the tip labels for each individual. 
#' If there are only two individuals, it defaults to red for teh donor and blue for the recipient.
#' Otherwise defaults to the r base graphics colors.
#' @return a coalescent tree in the Newick format (with branch lengths) and tip labels and states.
#' @examples 
#' tree <- sim.coal.tree(names = c("D", "R"), donor = c(0,1), t_inf = c(0,2), t_sam = c(3,4),
#'                       a = rep(20, 2), b = rep(3, 2),
#'                       nSamples = rep(500, 2))
#' tree <- sim.coal.tree(names = c("D", "R"), donor = c(0,1), t_inf = c(0,2), t_sam = c(3,4),
#'                       rhoD = rep(0, 2), rhoR = rep(0.1, 2),
#'                       sample.times = list(c(rep(2.2, 200), rep(2.3, 100), rep(3, 200)), c(rep(2.5, 250), rep(4, 250))),
#'                       tr_window = rep(2, 2))
#' @export
#' 

sim.coal.tree <- function(tt = NULL, 
                          names = NULL, donor = NULL, t_inf = NULL, t_sam = NULL,
                          a = rep(0, dim(tt)[1]), b = rep(3, dim(tt)[1]), 
                          rhoD = rep(0, dim(tt)[1]), rhoR = rep(0, dim(tt)[1]), 
                          nSamples = NULL,
                          sample.times = NULL,
                          gen.rate = 365/1.5,
                          tr_window = rep(Inf, dim(tt)[1]), #can transmit after sampling
                          file_path = "", 
                          tree_name = 1, 
                          tree_name_digits = nchar(as.character(tree_name)), 
                          save_tree = FALSE, 
                          plot_tree = FALSE,
                          plot_colors = NULL,
                          seed = sample(2e9, 1)){
  # multiple coalescent simulations
  
  #make sure there is a valid tree or way to build one
  if(!is.null(names) & !is.null(donor) & !is.null(t_inf) & !is.null(t_sam)){
    nInd <- length(names)
    if(nInd < 2) stop("There must be at least two individuals in the transmission tree")
    else if(nInd != length(donor) | nInd != length(t_inf) | nInd != length(t_sam)){
      stop("Lenth of 'names', 'donor', 't_inf', and 't_sam' must all be the same")
    } else if(!is.null(tt)){
      stop("Conflicting inputs: must only input one of transmission tree or names, donor, t_inf, and t_sam")
    } else{
      tt <- data.frame(names, donor, t_inf, t_sam, stringsAsFactors = FALSE)
    }
  }
  if(is.null(tt)){
    stop("No way to build transmission tree: Either transmission tree or names, donor, t_inf, and t_sam must be supplied")
  }
  
  #make sure t_sam matches sample.times if provided?
  
  set.seed(seed)
  
  #number of individuals in transmission tree
  nInd <- dim(tt)[1]
  
  #make sure inputs have the right number of elements
  if(length(a) != nInd) stop("a must be same length as number of individuals")
  if(length(b) != nInd) stop("b must be same length as number of individuals")
  if(length(rhoD) != nInd) stop("rhoD must be same length as number of individuals")
  if(length(rhoR) != nInd) stop("rhoR must be same length as number of individuals")
  
  if(any(c(a, b, rhoD, rhoR) < 0)) stop("Values of a, b, rhoD, and rhoR must be non-negative")
  
  #check whether (forwards time) transmission from recipient back to donor is ever possible
  rhoD_any <- any(rhoD > 0)
  
  #make sure nSamples and sample.times inputs are consistent
  if(is.null(nSamples) & is.null(sample.times)){
    nSamples <- rep(1, nInd) #default to one sample per individual if neither are provided
  }
  if(is.null(nSamples)){ #sample.times provided but not nSamples
    if(!is.list(sample.times)){
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else if(length(sample.times) != nInd){
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else{
      nSamples <- sapply(sample.times, length) #number of samples for each individual
    }
  }
  if(is.null(sample.times)){ #all sample times from each individual
    sample.times <- list()
    for(i in 1:nInd){
      sample.times[[i]] <- rep(tt$t_sam[i], nSamples[i])
    }
  } 
  if(!is.null(nSamples) & !is.null(sample.times)){ #make sure it is a list and lengths are correct if they are provded
    if(!is.list(sample.times)){
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else if(length(sample.times) != nInd){
      stop("sample.times must be given as a list of vectors with the sample times of each sequence for each patient")
    } else{
      for(i in 1:nInd){
        if(length(sample.times[[i]]) != nSamples[i]){
          stop("Length of sample.times[[i]] must be equal to nSamples[i]")
        }
      }
    }
  }
  
  
  
  #transmission tree (and other parameters) must be ordered from oldest to most recent infection date
  if(is.unsorted(tt$t_inf)){
    tt_old <- tt
    new_order <- order(tt$t_inf)
    tt <- tt[new_order,]
    a <- a[new_order]
    b <- b[new_order]
    rhoR <- rhoR[new_order]
    rhoD <- rhoD[new_order]
    nSamples <- nSamples[new_order]
    sample.times <- sample.times[new_order]
    tr_window <- tr_window[new_order]
    donor_shuff <- tt$donor
    for(i in 1:nInd){
      if(tt$donor[i] != 0){
        tt$donor[i] <- which(tt$names == tt_old$names[donor_shuff[i]])
      }
    }
  }
  
  #find which individuals have at least 1 sample
  nSamples_nonzero <- which(nSamples > 0)
  
  #sample times as vector
  sample.times.unlist <- unlist(sample.times)
  #find all unique sample times for easier checking of values
  sample.times.unique <- unique(unlist(sample.times))
  
  #find times when transmission windows open
  tr_window_times <- tt$t_inf + tr_window
  
  #TODO make sure values are possible
  
  #starting guess for how many iterations may be needed (can be extended in loop)
  iter_guess <- length(sample.times.unique)+ nInd #start guessing no migration events
  
  #vector for times of events
  t <- vector(length = iter_guess)
  #find oldest sampling time
  t[1] <- max(sample.times.unique)
  
  #initial lineages present
  lin_init <- -(1:sum(nSamples))
  #initial locations of each lineage (which patient each lineage is in)
  loc_init <- list()
  names <- character(0)
  leaf_times <- unlist(sample.times)
  for(i in 1:nInd){
    loc_init[[i]] <- rep(-i, nSamples[i]) #negative values imply samples were taken from that host at an earlier time
    loc_init[[i]][which(sample.times[[i]] == t[1])] <- i #positive values imply samples are in that host at this time
    names <- c(names, rep(tt$names[i], nSamples[i]))
  }
  #unlist into vector
  loc_init <- unlist(loc_init)
  
  #lins[[1]] <- data.frame(lin = lin_init, loc = loc_init)
  lins <- data.frame(lin = lin_init, loc = loc_init)
  
  #current number of lineages
  k_all <- nSamples #lineages in each host (including ones that have already been sampled and can't coalesce)
  k <- integer(nInd) #lineages available to coalesce in each host
  for(i in 1:nInd){
    k[i] <- length(which(loc_init == i))
  }
  k_tot <- sum(nSamples) #total number of lineages
  
  #names for leaves on newick tree
  names_newick <- paste(names, -lin_init, sep = "_")
  
  #declare matrix for internal nodes in tree
  nodes <- matrix(nrow = k_tot - 1, ncol = 2)
  
  int_node <- 0 #counter for number of internal nodes
  
  #initialize matrices for RVs
  z <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  mD <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  mR <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  t_new_w <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  t_new_s <- matrix(Inf, nrow = iter_guess, ncol = nInd)
  
  #initialize other things to keep track
  t_event <- numeric(iter_guess) #time for each event (not cumulative)
  event <- integer(iter_guess)
  
  #initialize list for migration events
  migrate <- vector("list", iter_guess)
  
  i <- 1 #index for number of events
  
  #loop as long as there is more than one lineage in longest infected patient or it is during the time second patient is infected
  while(k_all[1] > 1 || (rhoD_any && t[i] - tt$t_inf[2] > -1e-15) || sum(k_all[2:nInd]) > 0){
    #double lengths if necessary
    if(i == iter_guess){
      
      t <- c(t, numeric(iter_guess))
      
      z <- rbind(z, matrix(Inf, nrow = iter_guess, ncol = nInd))
      mD <- rbind(mD, matrix(Inf, nrow = iter_guess, ncol = nInd))
      mR <- rbind(mR, matrix(Inf, nrow = iter_guess, ncol = nInd))
      t_new_w <- rbind(t_new_w, matrix(Inf, nrow = iter_guess, ncol = nInd))
      t_new_s <- rbind(t_new_s, matrix(Inf, nrow = iter_guess, ncol = nInd))
      
      t_event <- c(t_event, numeric(iter_guess))
      event <- c(event, integer(iter_guess))
      
      migrate <- c(migrate, vector("list", iter_guess))
      
      iter_guess <- 2*iter_guess
    }
    #coalescence events
    #for(j in 1:nInd){
    #  if(k[j] > 1){ #also include "& t[i] <= tt$t_sam[j]"?
    #    unif_rv <- runif(1)
    #    z[i,j] <- Fz(unif_rv, k[j], a[j], b[j], (t[i] - tt$t_inf[j])*gen.rate)/gen.rate
    #  }
    #  #else{
    #  #  z[i,j] <- Inf
    #  #}
    #}
    k_atleast2 <- which(k >= 2)
    #print(k_atleast2)
    nk_atleast2 <- length(k_atleast2)
    #print(Fz(runif(nk_atleast2), k[k_atleast2], a[k_atleast2], b[k_atleast2], (t[i] - tt$t_inf)*gen.rate)/gen.rate)
    z[i,k_atleast2] <- Fz(runif(nk_atleast2), k[k_atleast2], a[k_atleast2], b[k_atleast2], (t[i] - tt$t_inf[k_atleast2])*gen.rate)/gen.rate
    
    #(reverse time) migrations events from donor to recipient
    if(rhoD_any){
      for(j in 1:nInd){
        if(rhoD[j] > 0 && k[j] >= 1 &&
           any(t[i] - tt$t_inf[which(tt$donor == j)] > 0 && #test to see if any recipients with donor j are in a contact window with j
               t[i] - tt$t_inf[which(tt$donor == j)] <= tr_window[which(tt$donor == j)])){
          mD[i,j] <- Fm(runif(1), k[j], rhoD[j], a[j], b[j], (t[i] - tt$t_inf[j])*gen.rate)/gen.rate
        } #else{
        #  mD[i,j] <- Inf
        #}
      }
    }
    
    #(reverse time) migrations events from recipient to donor
    t_t_inf_diffs <- t[i] - tt$t_inf
    k_nonzero <- which(k >= 1)
    for(j in k_nonzero){
      if(j != 1 && t_t_inf_diffs[j] >= 0 && t[i] <= tr_window_times[j]){
        if(rhoR[j] > 0){
          mR[i,j] <- min(t_t_inf_diffs[j], Fm(runif(1), k[j], rhoR[j], a[j], b[j], (t[i] - tt$t_inf[j])*gen.rate)/gen.rate)
        } else{
          mR[i,j] <- t_t_inf_diffs[j]
        }
        
      } #else{
      #  mR[i,j] <- Inf
      #}
    }
    
    #time when a new transmission window opens
    #for(j in 1:nInd){
    #  if(t[i] > tt$t_inf[j] + tr_window[j]){   
    #    #after end of transmission window (forwards time)
    #    #(potentially allow transmission after sampling)
    #    t_new_w[i,j] <- t[i] - (tt$t_inf[j] + tr_window[j])
    #  } #else{
    #  #  t_new_w[i,j] <- Inf
    #  #}
    #}
    tr_window_times_t_diffs <- t[i] - tr_window_times
    t_new_w[i,tr_window_times_t_diffs > 0] <- tr_window_times_t_diffs[tr_window_times_t_diffs > 0]
    #t_new_w[i,(t[i] > tr_window_times)]
    
    #time when a new sample is taken #i think there is a simpler way to do this
    for(j in nSamples_nonzero){
      if(any(t[i] > sample.times[[j]])){
        #is this after (forwards time) any samples have been taken from patient j?
        t_sam_diff <- t[i] - sample.times[[j]] #find difference between current time and sample times
        t_new_s[i,j] <- min(t_sam_diff[t_sam_diff > 0]) #minimum of values that are positive
      } #else{
      #  t_new_s[i,j] <- Inf
      #}
      
    }
    
    #which event happenes first within each class of event
    z_i <- which.min(z[i,]) #index of first coalescence event
    z_t <- z[i, z_i] #time of first coalescence event
    mD_i <- which.min(mD[i,])
    mD_t <- mD[i, mD_i]
    mR_t <- min(mR[i,])
    mR_i <- which(mR[i,] == mR_t)
    
    t_new_w_i <- which.min(t_new_w[i,])
    t_new_w_t <- t_new_w[i, t_new_w_i]
    t_new_s_i <- which.min(t_new_s[i,])
    t_new_s_t <- t_new_s[i, t_new_s_i]
    
    #print(c(z_t, mD_t, mR_t, t_new_w_t, t_new_s_t))
    
    #which kind of event happens first
    event[i] <- which.min(c(z_t, mD_t, mR_t, t_new_w_t, t_new_s_t))
    #time to first event
    t_event[i] <- c(z_t, mD_t, mR_t, t_new_w_t, t_new_s_t)[event[i]]
    
    #prevent reverse time migration to recipient if that would mean recipient to be infected before transmission
    if(event[i] == 2 & all(t_event[i] > t[i] - tt$t_inf[which(tt$donor == mD_i)])){
      event[i] <- order(c(z_t, mD_t, mR_t, t_new_t))[2] # take second fastest event if fastest cannot happen
      t_event[i] <- sort(c(z_t, mD_t, mR_t, t_new_t))[2]
      #make sure second event is also possible
      if(t_event[i] == Inf) break
    }
    
    #lins[[i+1]] <- lins[[i]] #copy current lineages to new timestep
    lins_prev <- lins #make copy of old lineages
    
    #first event is coalescence somewhere
    if(event[i] == 1){
      int_node <- int_node + 1
      #coalescence takes place in individual z_i
      possible_lins <- lins_prev$lin[which(lins_prev$loc == z_i)] #lineages in z_i that can coalesce
      #choose two lineages to coalesce
      coal <- sample(possible_lins, 2, replace = FALSE)
      nodes[int_node,] <- coal #add coalescence event to nodes in tree
      #add new internal lineage
      #lins[[i+1]] <- rbind(lins[[i+1]], c(int_node, z_i))
      lins <- rbind(lins, c(int_node, z_i))
      #remove coalesced lineages
      lins <- lins[-which(lins_prev$lin == coal[1] | lins_prev$lin == coal[2]),]
    }
    #first event is reverse time migration from donor to recipient
    else if(event[i] == 2){
      #possible_lins <- lins[[i]]$lin[which(lins[[i]]$loc == mD_i)]
      possible_lins <- lins_prev$lin[which(lins_prev$loc == mD_i)]
      migrate <- possible_lins[sample.int(length(possible_lins), 1)]
      #find potential recipients
      pot_recips <- which(tt$donor == mD_i & t_event[i] < t[i] - tt$t_inf)
      recips <- pot_recips[which(t[i] - tt$t_inf[pot_recips] > 0 & 
                                   t[i] - tt$t_inf[pot_recips] <= tr_window[pot_recips])]
      recip <- recips[sample.int(length(recips), 1)] #where does lineage actually go?
      #move migrating lineage (reverse time) from donor to recipient
      lins$loc[which(lins$lin == migrate)] <- recip #pick which recipient to use
    }
    #first event is reverse time migration from recipient to donor
    else if(event[i] == 3){
      #was this the initial transmission event, where all lineages go to donor?
      if(any(t_event[i] == t[i] - tt$t_inf[mR_i])){
        #move migrating lineages (reverse time) from recipient to donor
        for(j in 1:length(mR_i)){
          lins$loc[which(lins_prev$loc == mR_i[j])] <- tt$donor[mR_i[j]]
        }
      }
      else{
        possible_lins <- lins$lin[which(lins_prev$loc == mR_i)]
        migrate <- possible_lins[sample.int(length(possible_lins), 1)]
        #move migrating lineage (reverse time) from recipient to donor
        lins$loc[which(lins$lin == migrate)] <- tt$donor[mR_i]
      }
    }
    #first event is new transmission window opening/new sample
    if(event[i] == 4 | event[i] == 5 | any(abs(t[i] - sample.times.unique - t_event[i]) <= 1e-15)){
      
      #see if event is new sample
      if(any(abs(t[i] - sample.times.unique - t_event[i]) <= 1e-15)){
        new_sam_index <- which((t[i] - sample.times.unlist - t_event[i]) <= 1e-15)
        new_sam <- lin_init[new_sam_index]
        lins$loc[which(lins$lin %in% new_sam)] <- abs(lins$loc[which(lins$lin %in% new_sam)])
      }
    }
    
    #reduce current time by elapsed time during event
    t[i+1] <- t[i] - t_event[i]
    #t[i+1] <- t[1] - sum(sort(t_event))
    
    #find new k's
    #for(j in 1:nInd){
    #  #k[j] <- length(which(lins[[i+1]]$loc == j))
    #  k[j] <- length(which(lins$loc == j))
    #  #k_all[j] <- length(which(lins[[i+1]]$loc == j | lins[[i+1]]$loc == -j))
    #  k_all[j] <- length(which(lins$loc == j | lins$loc == -j))
    #}
    k <- tabulate(lins$loc, nbins = nInd)
    k_all <- k + tabulate(-lins$loc, nbins = nInd)
    
    i <- i + 1
  }
  
  #if(dim(lins[[i]])[1] != 1){
  if(dim(lins)[1] != 1){
    stop("Tree Not Fully Coalesced")
  }
  
  #time of events in reverse time
  rev_time <- cumsum(t_event)
  
  #times of coalescence events
  #coal_times <- rev_time[which(event == 1)]
  
  coal_times <- t[which(event == 1) + 1]
  
  # declare character vector for building up parts of newick tree
  branch <- character(int_node)
  
  #build newick tree
  for(i in 1:int_node){
    if(all(nodes[i,] < 0)){ #both edges go to leaves
      branch[i] <- paste0("(", names_newick[-nodes[i, 1]], ":", leaf_times[-nodes[i, 1]] - coal_times[i], ",",
                          names_newick[-nodes[i, 2]], ":", leaf_times[-nodes[i, 2]] - coal_times[i], ")")
    }
    else if(nodes[i,1] < 0 & nodes[i,2] > 0){ #first edge goes to leaf, second goes to internal node
      branch[i] <- paste0("(", names_newick[-nodes[i, 1]], ":", leaf_times[-nodes[i, 1]] - coal_times[i], ",",
                          branch[nodes[i,2]], ":", -coal_times[i] + coal_times[nodes[i,2]], ")")
    }
    else if(nodes[i,1] > 0 & nodes[i,2] < 0){ #first edge goes to internal node, second goes to leaf
      branch[i] <- paste0("(", branch[nodes[i,1]], ":", -coal_times[i] + coal_times[nodes[i,1]], ",",
                          names_newick[-nodes[i, 2]], ":", leaf_times[-nodes[i, 2]] - coal_times[i], ")")
    }
    else if(all(nodes[i,] > 0)){ #both edges go to internal nodes
      branch[i] <- paste0("(", branch[nodes[i,1]], ":", -coal_times[i] + coal_times[nodes[i,1]], ",",
                          branch[nodes[i,2]], ":", -coal_times[i] + coal_times[nodes[i,2]], ")")
    }
  }
  
  #add semicolon to finish tree
  tree_newick <- paste0(branch[int_node], ";")
  
  #tip states
  #tip_states <- c(rep(0, nSamples[1]), rep(1, nSamples[2]))
  tip_states <- vector()
  for(i in 1:nInd){
    tip_states <- c(tip_states, rep(i-1, nSamples[i]))
  }
  #make text object with newick tree and tip names/states
  treeAndTips <- c(tree_newick, paste(t(names_newick), t(tip_states)))
  
  #format tree name
  if(!is.character(tree_name)){
    #add left padding zeros if the tree name is not a character
    tree_name <- formatC(tree_name, width = tree_name_digits, format = "d", flag = "0")
  }
  
  if(save_tree == TRUE){
    #add trailing "/" to file_path if it needs it
    #is the string not empty? is the last char not a "/"?
    if(nchar(file_path) != 0 & substr(file_path,nchar(file_path),nchar(file_path)) != "/"){ 
      file_path <- paste0(file_path, "/")
    }
    write(treeAndTips, file = file.path(paste0(file_path,"tree", tree_name, ".txt")))
  }
  
  #plot tree if requested
  if(plot_tree == TRUE){
    #convert to phylo class
    tree <- ape::read.tree(text = tree_newick)
    
    #find colors for tip labels
    if(is.null(plot_colors)){
      if(nInd == 2){
        plot_colors = c("Red", "Blue") #use red and blue if there is just a donor and a recipient
      } else{
        plot_colors = 1:nInd
      }
    }
    
    #find which tips are for recipient and color appropriately #TODO make it work for more than two individuals
    #tip_colors <- rep("red", sum(nSamples[1:2]))
    #tip_colors[grep(tt$names[2], tree$tip.label)] <- "blue"
    
    tip_colors <- vector(length = sum(nSamples))
    for(i in 1:nInd){
      tip_colors[grep(tt$names[i], tree$tip.label)] <- plot_colors[i]
    }
    
    #find infection ages at time of last sample
    i_age <- max(tt$t_sam) - tt$t_inf
    
    #print(coal_times[int_node])
    
    #plot to pdf
    pdf_height <- max(8,sqrt(sum(nSamples)))
    pdf(file = file.path(paste0(file_path,"treeplot", tree_name, ".pdf")), width = 11, height = pdf_height)
    
    ape::plot.phylo(tree, tip.color = tip_colors, cex = min(1, 6/sqrt(sum(nSamples))), #label scaling seems to work okay for pdfs as drawn
                    x.lim = c(min(0, tt$t_inf[1] - coal_times[int_node]), max(sample.times.unique) - coal_times[int_node])) 
    #ape::plot.phylo(tree) #label scaling seems to work okay for pdfs as drawn
    title(main = paste0("Tree ", tree_name,
                        ", aD = ", a[1], ", bD = ", b[1], ", aR = ", a[2], ", bR = ", b[2], ", rhoD = ", rhoD[1], ", rhoR = ", rhoR[2], 
                        ", nD = ", nSamples[1], ", nR = ", nSamples[2]),
          adj = 1)
    axis(1) #turn on x axis
    abline(v = tt$t_inf - coal_times[int_node], col = plot_colors, lty = 2)
    
    #Draw rectangle for transmission window if there are only 2 individuals and it is applicable
    if(nInd == 2 & tr_window[2] != Inf & rhoR[2] > 0){
      rect(tt$t_inf[2] - coal_times[int_node], 0, tt$t_inf[2] + tr_window[2] - coal_times[int_node], sum(nSamples), 
           col = "#0000FF18", border = NA)
    }
    
    legend("topleft", legend = paste0(c(tt$names), " Infected"), lty = 2, col = plot_colors,
           xpd = TRUE, inset = c(0,-.08*(10/pdf_height)))
    
    dev.off()
  }
  
  return(treeAndTips)
}

#write(tree_newick, file = "tree_newick.tree")

#functions for inverse cdf sampling for coalescence and migration 
#coalescence
Fz <- function(u, k, a, b, t){
  Fz <- (1 - (1 - u)^(b/choose(k,2)))*(a + b*t)/b
}
#migration
Fm <- function(u, k, rho, a, b, t){
  Fm <- (1 - (1 - u)^(b/(k*rho)))*(a + b*t)/b
}

#' @title Subsample phybreakdata
#' @description Subsamples sequences from phybreakdata simulation data
#' @param sim_data A phybreakdata object to subsample sequences from
#' @param nseq_sub A number or vector for the number of sequences to keep from each individual
#' @return A subsampled phybreakdata object
#' @export
#' 
subsam.seqs <- function(sim_data, nseq_sub){
  #number of individuals
  nInds <- length(sim_data$sim.infection.times)
  #new number of sequences
  host.names <- unique(sim_data$sample.hosts)
  
  #make sure nseq_sub is the correct length
  if(length(nseq_sub) == 1){
    nseq_sub <- rep(nseq_sub, nInds)
  } else if(length(nseq_sub) != nInds){
    stop("nseq_sub must have length equal to 1 or the number of individuals")
  }
  
  #find indices which sequences to keep
  sub_indices <- c()
  
  for(i in 1:nInds){
    #subsample specific sequences if nsub_seq is given as a list
    if(is.list(nseq_sub)){
      sub_indices <- c(sub_indices, which(sim_data$sample.hosts == host.names[i])[nseq_sub[[i]]])
    } else{
      sub_indices <- c(sub_indices, which(sim_data$sample.hosts == host.names[i])[1:nseq_sub[i]])
    }
  }
  sub_indices <- sort(sub_indices)
  
  sequences <- sim_data$sequences[sub_indices]
  sample.times <- sim_data$sample.times[sub_indices]
  sample.hosts <- sim_data$sample.hosts[sub_indices]
  sample.names <- names(sim_data$sample.times[sub_indices])
  sim.tree <- ape::keep.tip(sim_data$sim.tree, sub_indices)
  
  #reorder stuff to be consistent
  #make ordering transformations (borrowed from phybreakdata)
  allhosts <- unique(sample.hosts)
  allfirsttimes <- rep(FALSE, length(sample.times))
  sapply(allhosts, function(x) allfirsttimes[which(min(sample.times[sample.hosts == x]) == sample.times & (sample.hosts == x))[1]] <<- TRUE)
  outputorderhosts <- order(sample.times[allfirsttimes])
  orderedhosts <- host.names[allfirsttimes][outputorderhosts]
  outputordersamples <- order(!allfirsttimes, match(sample.hosts, orderedhosts), sample.times)
  sequences <- sequences[outputordersamples, ]
  sample.times <- sample.times[outputordersamples]
  sample.names <- sample.names[outputordersamples]
  host.names <- host.names[outputordersamples]
  
  sim_data_sub <- phybreakdata(sequences = sequences, 
                               sample.times = sample.times, 
                               sample.names = sample.names, 
                               host.names = sample.hosts,
                               sim.infection.times = sim_data$sim.infection.times[orderedhosts],
                               sim.infectors = sim_data$sim.infectors[orderedhosts], 
                               sim.tree = sim.tree)
  
  return(sim_data_sub)
}

