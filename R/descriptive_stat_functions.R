#'@title Plot densities of predicted infection times
#'@description Function to plot the individual and total densities of predicted infection times on one page (in two panels)
#'@param df A data frame containing infection time distributions as in the output from prepare.HIV.data 
#'with find_infection_age_distributions = TRUE.
#'@param sort_quant The quantile with which to chronologically sort clusters and individuals within clusters. 
#'Also determines the left boundary of the rectangle enclosing each cluster. 
#'For example, if set to 0.01, it uses the first percentile of the predicted infection times for sorting.
#'@param xmin_override If specified, this overrides the default method of determining the lower bound for plotting.
#'@export
infection.distributions.plot <- function(df, 
                                         sort_quant = 0.01,
                                         xmin_override = NULL){
  #find names of all clusters
  cluster_names <- unique(df$cluster_ID)
  
  #find number of clusters
  nClust <- length(cluster_names)
  
  #find 1st and 0.1st percentile infection times
  mins <- sapply(df$cdf, FUN = function(cdf) cdf$x[which(cdf$y >= sort_quant)[1]])
  
  #put mins into data frame
  df$mins <- mins
  
  #split dataframe into clusters
  df_split <- vector(mode = "list", length = nClust)
  for(i in seq_len(nClust)){
    df_split[[i]] <- df[df$cluster_ID == cluster_names[i],]
  }
  #return(df_split)
  #find cluster sizes
  nInds <- sapply(df_split, FUN = function(x) dim(x)[1])
  #return(df_split)
  #find chronological order of clusters
  mins_bycluster <- sapply(df_split, FUN = function(x) min(x$mins))
  #find order of clusters
  cluster_order <- order(mins_bycluster)
  #sort clusters
  df_split_chron <- df_split[cluster_order]
  #sort individuals within cluster by 1st percentile infection time
  for(i in seq_along(df_split_chron)){
    df_split_chron[[i]] <- df_split_chron[[i]][order(df_split_chron[[i]]$mins),]
  }
  
  #find most recent diagnosis times for each cluster
  time_maxes <- sapply(df_split_chron, FUN = function(x) max(sapply(x$cdf, FUN = function(y) max(y$x))))
  
  #y box division
  ymins <- c(0,cumsum(nInds[cluster_order]+1))[seq_along(nInds)]
  ymaxes <- cumsum(nInds[cluster_order]+1)
  
  extrema.df <- data.frame(xmins = mins_bycluster[cluster_order], xmaxes = time_maxes, ymins = ymins, ymaxes = ymaxes)
  
  if(is.null(xmin_override)){
    time_min_all <- min(mins)
  } else{
    time_min_all <- xmin_override
  }
  time_max_all <- max(time_maxes)
  
  #maximum density
  dens_max_all <- max(sapply(df_split_chron, FUN = function(x) max(sapply(x$pdf, FUN = function(y) max(y$y)))))
  
  #find overall average infection times
  dens.fns <- vector(mode = "list", length = sum(nInds))
  informed <- logical(length = sum(nInds))
  #individual index
  k <- 0
  for(i in seq_along(df_split_chron)){
    for(j in seq_len(nInds[cluster_order[i]])){
      k <- k + 1
      informed[k] <- df_split_chron[[i]]$usable_pol[j]
      dens.fns[[k]] <- approxfun(df_split_chron[[i]]$cdf[[j]]$x, rev(df_split_chron[[i]]$infection_age_dists_diag[[j]]$y), yleft = 0, yright = 0)
    }
  }
  #informed <- as.vector(informed)
  times <- seq(time_min_all, time_max_all, length.out = 10000)
  overall_dens_inf <- rep(0, length(times))
  overall_dens_all <- rep(0, length(times))
  for(k in 1:length(dens.fns)){
    if(informed[k] == 1){
      overall_dens_inf <- overall_dens_inf + dens.fns[[k]](times)
    }
    overall_dens_all <- overall_dens_all + dens.fns[[k]](times)
  }
  
  overall_dens.df <- data.frame(time = times, 
                                inf = overall_dens_inf,
                                all = overall_dens_all)
  
  colors <- c("#0000FF", "#000000")
  
  p1 <- ggplot2::ggplot() + 
    ggplot2::scale_alpha_continuous(range = c(0, dens_max_all)) + 
    ggplot2::xlim(time_min_all, time_max_all+1) +
    ggplot2::ylim(0, sum(nInds)+length(nInds)) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Individual Patient Infection Time Densities") +
    ggplot2::theme_classic() + ggplot2::theme(legend.position = "none", 
                                              axis.title.x = ggplot2::element_blank(),
                                              #axis.title.y = ggplot2::element_blank(), 
                                              axis.text.y = ggplot2::element_blank(), 
                                              axis.ticks.y = ggplot2::element_blank(),
                                              plot.margin = ggplot2::margin(l = 5, r = 10)
    )
  y_level <- 0
  
  p1 <- p1 + ggplot2::geom_rect(data = extrema.df, ggplot2::aes(xmin = xmins, xmax = xmaxes, ymin = ymins+0.5, ymax = ymaxes-0.5), 
                                color = "#F2F2F2", fill = "#F2F2F2", size = 0.3) + 
    ggplot2::annotate("text", x = extrema.df$xmaxes+1, y = (extrema.df$ymins+extrema.df$ymaxes)/2, 
                      label = cluster_names[cluster_order], size = 0.9)
  for(i in seq_along(nInds)){
    for(j in seq_len(nInds[cluster_order[i]])){
      #color for line depending on whether distribution is informed or not
      informed <- df_split_chron[[i]]$usable_pol[j]
      if(informed == 1){
        col <- colors[1]
      } else{
        col <- colors[2]
      }
      y_level <- y_level + 1
      #trim pdf dataframe to avoid values being outside plot range
      pdf_df_trim <- as.data.frame(df_split_chron[[i]]$pdf[[j]])
      pdf_df_trim <- pdf_df_trim[pdf_df_trim$x >= time_min_all,]
      p1 <- p1 + 
        ggplot2::geom_line(data = pdf_df_trim, 
                           ggplot2::aes(x = x, alpha = y), y = y_level, size = 0.39, color = col)
    }
    #another space in between clusters
    y_level <- y_level + 1 
  }
  
  #separate plot for overall
  #find how much xmin for total density graph needs to be offset so the x axes match
  xmin_offset <- 0.02*(time_max_all+1 - time_min_all)
  #function to format y ticks
  #format_y_ticks <- function(x) sprintf("%4.1f", x)
  #trim overall dataframe to avoid values outside plot range
  overall_dens.df_trim <- overall_dens.df[overall_dens.df$time >= time_min_all+xmin_offset,]
  #return(overall_dens.df_trim)
  p2 <- ggplot2::ggplot() + ggplot2::geom_line(data = overall_dens.df_trim, ggplot2::aes(x = time, y = inf), color = colors[1]) +
    ggplot2::geom_line(data = overall_dens.df_trim, ggplot2::aes(x = time, y = all), color = colors[1], lty = 2) +
    ggplot2::xlim(time_min_all+xmin_offset, time_max_all+1) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Total Density") + 
    #ggplot2::scale_y_continuous(labels = format_y_ticks)
    ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 10))
  
  #pdf(file = "all_infection_times.pdf", width = 7.5, height = 8.75)
  gridExtra::grid.arrange(p1, p2, nrow = 2, heights = c(1, 0.2))
  #dev.off()
}