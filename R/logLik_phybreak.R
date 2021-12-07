#' Log-likelihood of a phybreak-object.
#' 
#' The likelihood of a \code{phybreak}-object is calculated, with the option to include or exclude parts of the 
#' likelihood for genetic data, phylogenetic tree (within-host model), sampling times and generation times.
#' 
#' The sequence likelihood is calculated by Felsenstein's pruning algorithm, assuming a prior probability of 0.25 
#' for each nucleotide. The within-host likelihood is the likelihood of coalescence times given the within-host model 
#' and slope. The generation interval and sampling interval likelihood are log-densities of the gamma distributions 
#' for these variables.
#' 
#' @param object An object of class \code{phybreak}.
#' @param genetic Whether to include the likelihood of the mutation model.
#' @param withinhost Whether to include the likelihood of within-host (coalescent) model.
#' @param sampling Whether to include the likelihood of the sampling model (sampling intervals).
#' @param generation Whether to include the likelihood of the transmission model (generation intervals).
#' @param ... Some methods for this generic require additional arguments. None are used in this method.
#' @return The log-likelihood as an object of class logLik.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' logLik(MCMCstate)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' logLik(MCMCstate)
#' 
#' tree0 <- get_phylo(MCMCstate)
#' seqdata <- get_data(MCMCstate)$sequences
#' phangorn::pml(tree0, seqdata, rate = 0.75*get_parameters(MCMCstate, "mu") 
#' logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
#'        sampling = FALSE, generation = FALSE) #should give the same result as 'pml'
#' @export
logLik.phybreak <- function(object, genetic = TRUE, withinhost = TRUE, sampling = TRUE, generation = TRUE, 
                            distance = TRUE, ...) {
  res <- 0
  if (genetic) {
    if(object$h$use.pml == FALSE){
      res <- res + with(object, .likseq(matrix(unlist(d$sequences), ncol = d$nsamples), 
                                        attr(d$sequences, "weight"), 
                                        v$nodeparents, v$nodetimes, p$mu, d$nsamples))
    } else if(object$h$use.pml == TRUE){
      res <- res + with(object, phangorn::pml(get_phylo(object), d$sequences, rate = p$mu, 
                                              bf = p$base.freq, Q = p$sub.rates, 
                                              shape = p$het.shape, inv = d$inv.sites))
    }
  }
  if (generation) {
    res <- res + with(object, lik_gentimes(p$gen.shape, p$gen.mean, v$inftimes, v$infectors,
                                           p$gen.nonpar, p$gen.pdf, d$sample.times, p$gen.dens.scale, p$post.sam.trans.rate))
  }
  if (sampling) {
    res <- res + with(object, lik_sampletimes(p$obs, p$sample.shape, p$sample.mean, v$nodetimes, v$inftimes, 
                                              p$sample.nonpar, p$sample.pdf))
  }
  if (withinhost) {
    objectenv <- object
    objectenv$v <- phybreak2environment(objectenv$v)
    res <- res + with(object, lik_coaltimes(objectenv))
  }
  if(distance && !is.null(object$d$distances)) {
    res <- res + with(object, lik_distances(p$dist.model, p$dist.exponent, p$dist.scale, p$dist.mean, 
                                            v$infectors, d$distances))
  }
  attributes(res) <- list(
    nobs = object$p$obs,
    df = 1 + object$h$est.mG + object$h$est.mS + object$h$est.wh.s + object$h$est.wh.e + object$h$est.wh.0 +
      object$h$est.dist.e + object$h$est.dist.s + object$h$est.dist.m,
    genetic = genetic, withinhost = withinhost, sampling = sampling, generation = generation, distance = distance
  )
  class(res) <- "logLik"
  return(res)
}


### calculate the log-likelihood of generation intervals 
lik_gentimes <- function(shapeG, meanG, inftimes, infectors, gen.nonpar, gen.pdf, sample.times, gen.dens.scale, post.sam.trans.rate) {
  if(gen.nonpar == FALSE){
    sum(dgamma(inftimes[infectors > 0] - 
                 inftimes[infectors[infectors > 0]], 
               shape = shapeG, scale = meanG/shapeG, log = TRUE))
  } else{
    #for individual gen densities
    #sum(log(mapply(FUN = function(gen.pdf, inftimes, infectors){gen.pdf(inftimes[infectors > 0] - 
    #                                                                      inftimes[infectors[infectors > 0]])}, 
    #               gen.pdf = gen.pdf, inftimes = inftimes, infectors = infectors)))
    #find penalties for transmission after sampling
    rel.rate <- mapply(FUN = post.sam.trans.adjust, 
                       t_sam = sample.times[infectors[infectors > 0]] - sample.times[1], 
                       t_inf = inftimes[infectors > 0],
                       post.sam.trans.rate = post.sam.trans.rate[infectors > 0])
    sum(log(gen.pdf(inftimes[infectors > 0] - inftimes[infectors[infectors > 0]])*gen.dens.scale*rel.rate))
  }
  
}

#function to for adjustment of transmission rate for individuals who have already been sampled
#post.sam.trans.rate is the relative likelihood of transmission after sampling compared to before
#rates are scaled to sum to 1
#for example, if post.sam.trans.rate is 0.1, then rate = 10/11 before sampling and 1/11 after sampling
post.sam.trans.adjust <- function(t_sam, t_inf, post.sam.trans.rate){
  if(t_sam > t_inf){
    rate <- 1/(1+post.sam.trans.rate)
  } else if(t_sam <= t_inf){
    rate <- post.sam.trans.rate/(1+post.sam.trans.rate) 
  } 
  return(rate)
}

### calculate the log-likelihood of sampling intervals 
lik_sampletimes <- function(obs, shapeS, meanS, nodetimes, inftimes, sample.nonpar, sample.pdf) {
  if(sample.nonpar == FALSE) {
    sum(dgamma(nodetimes[1:obs] - inftimes, shape = shapeS, scale = meanS/shapeS, log = TRUE))
  }
  else{
    sum(log(mapply(FUN = function(sample.pdf, nodetimes, inftimes) {max(sample.pdf(nodetimes - inftimes), 1e-16)}, #TODO right way to handle zeros here
                   sample.pdf = sample.pdf, nodetimes = nodetimes[1:obs], inftimes = inftimes)))
  }
}

### calculate the log-likelihood of distances 
lik_distances <- function(dist.model, dist.exponent, dist.scale, dist.mean, infectors, distances) {
  if(dist.model == "none") return(0)
  distancevector <- distances[cbind(1:length(infectors), infectors)]
  switch(dist.model,
         power = sum(log(
           dist.exponent * sin(pi/dist.exponent) / 
             (dist.scale * pi * (1 + (distancevector/dist.scale)^dist.exponent))
           )),
         exponential = sum(
           log(dist.exponent) - dist.exponent * distancevector
         ),
         poisson = sum(
           -dist.mean + distancevector * log(dist.mean) - lgamma(1 + distancevector)
         )
  )
}

### calculate the log-likelihood of coalescent intervals 
lik_coaltimes <- function(phybreakenv) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  if(phybreakenv$p$wh.model == "linear" && phybreakenv$p$wh.bottleneck == "wide") {
    if(min(phybreakenv$v$inftimes) - min(phybreakenv$v$nodetimes[phybreakenv$v$nodetypes == "c"]) > 
       phybreakenv$p$sample.mean + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) return(-Inf)
  }
  
  remove0nodes <- phybreakenv$v$nodetypes != "0"
  nodetypes <- phybreakenv$v$nodetypes[remove0nodes]
  nodehosts <- phybreakenv$v$nodehosts[remove0nodes]
  nodetimes <- phybreakenv$v$nodetimes[remove0nodes]
  inftimes <- c(min(phybreakenv$v$inftimes) - phybreakenv$p$sample.mean, phybreakenv$v$inftimes)
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodehosts, nodetimes)
  
  coalnodes <- coalnodes[orderednodes]
  orderedhosts <- nodehosts[orderednodes]
  
  bottlenecks <- sapply(0:phybreakenv$p$obs, function(i) sum((orderedhosts == i) * (1 - 2 * coalnodes))) - 1
  dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage[!duplicated(orderedhosts)] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  whtimes <- nodetimes[orderednodes] - inftimes[orderedhosts + 1]
  whtimes[c(!duplicated(orderedhosts)[-1], FALSE)] <- 0
  
  logcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                         linear = -log(phybreakenv$p$wh.level + phybreakenv$p$wh.slope * whtimes[coalnodes]),
                         exponential = 
                           -log(phybreakenv$p$wh.level * 
                                  exp(phybreakenv$p$wh.exponent * 
                                        whtimes[coalnodes])),
                         constant = -log(phybreakenv$p$wh.level) * coalnodes)
  cumcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite=,
                         linear = log(whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope + 
                                        ((whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) == 0)) / phybreakenv$p$wh.slope,
                         exponential  = -1/(phybreakenv$p$wh.level * phybreakenv$p$wh.exponent * 
                                              exp(phybreakenv$p$wh.exponent * whtimes)),
                         constant = whtimes/phybreakenv$p$wh.level)
  coalratediffs <- cumcoalrates - c(0, head(cumcoalrates, -1))
  logcoalescapes <- -coalratediffs * choose(nrlineages, 2)
  
  return(sum(logcoalrates) + sum(logcoalescapes))
}


### calculate the log-likelihood of coalescent intervals in a single host
lik_coaltimes_host <- function(phybreakenv, hostID) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  if(phybreakenv$p$wh.model == "linear" && phybreakenv$p$wh.bottleneck == "wide") {
    if(min(phybreakenv$v$inftimes) - min(phybreakenv$v$nodetimes[phybreakenv$v$nodetypes == "c"]) > 
       phybreakenv$p$sample.mean + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) return(-Inf)
  }
  
  selecthostnodes <- phybreakenv$v$nodehosts == hostID
  nodetypes <- phybreakenv$v$nodetypes[selecthostnodes]
  nodetimes <- phybreakenv$v$nodetimes[selecthostnodes]
  inftime <- phybreakenv$v$inftimes[hostID]
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodetimes)
  
  coalnodes <- coalnodes[orderednodes]

  bottlenecks <- sum(1 - 2 * coalnodes) - 1
  dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage[1] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  whtimes <- nodetimes[orderednodes] - inftime

  logcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                         linear = -log(phybreakenv$p$wh.level + phybreakenv$p$wh.slope * whtimes[coalnodes]),
                         exponential = 
                           -log(phybreakenv$p$wh.level * 
                                  exp(phybreakenv$p$wh.exponent * 
                                        whtimes[coalnodes])),
                         constant = -log(phybreakenv$p$wh.level) * coalnodes)
  cumcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite=,
                         linear = log(whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope + 
                                        ((whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) == 0)) / phybreakenv$p$wh.slope,
                         exponential  = -1/(phybreakenv$p$wh.level * phybreakenv$p$wh.exponent * 
                                              exp(phybreakenv$p$wh.exponent * whtimes)),
                         constant = whtimes/phybreakenv$p$wh.level)
  coalratediffs <- cumcoalrates - c(0, head(cumcoalrates, -1))
  logcoalescapes <- -coalratediffs * choose(nrlineages, 2)
  
  return(sum(logcoalrates) + sum(logcoalescapes))
}


### calculate the log-likelihood of the within-host topology
lik_topology_host <- function(phybreakenv, hostID) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  selecthostnodes <- phybreakenv$v$nodehosts == hostID
  nodetypes <- phybreakenv$v$nodetypes[selecthostnodes]
  nodetimes <- phybreakenv$v$nodetimes[selecthostnodes]
  inftime <- phybreakenv$v$inftimes[hostID]
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodetimes)
  
  coalnodes <- coalnodes[orderednodes]
  
  bottlenecks <- sum(1 - 2 * coalnodes) - 1
  dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage[1] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  logcoalprobabilities <- -log(choose(nrlineages[c(FALSE, head(coalnodes, -1))], 2))
  
  return(sum(logcoalprobabilities))
}
