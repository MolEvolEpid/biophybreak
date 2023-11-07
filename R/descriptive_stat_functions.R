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

#'@title Plot descriptive/demographics statistics
#'@description Function to plot descriptive and demographic statistics of a dataset like gender and birth location distributions
#'@param df A data frame containing infection time distributions as in the output from prepare.HIV.data 
#'with find_infection_age_distributions = TRUE.
#'@param max_groups The maximum number of birth or suspected transmission locations to be plotting
#'before starting to group less common locations into "other"
#'@param save_plots Whether or not to save the plots as pdfs
#'@param run_set_title The first part of the file names if plots are saved
#'@param width Width of pdf outputs if plots are saved
#'@param height Height of pdf outputs if plots are saved
#'@export
descriptive.plots <- function(df, 
                              max_groups = 10, 
                              save_plots = FALSE, 
                              run_set_title = "1", 
                              width = 9, 
                              height = 5){
  #find names of all clusters
  cluster_names <- unique(df$cluster_ID)
  
  #find number of clusters
  nClust <- length(cluster_names)
  
  #split dataframe into clusters
  df_split <- vector(mode = "list", length = nClust)
  for(i in seq_len(nClust)){
    df_split[[i]] <- df[df$cluster_ID == cluster_names[i],]
  }
  
  #find cluster sizes
  nInds <- sapply(df_split, FUN = function(x) dim(x)[1])
  
  p_cluster_sizes <- ggplot2::ggplot(data = data.frame(nInds = as.factor(nInds)), ggplot2::aes(x = nInds)) + 
    ggplot2::geom_bar(stat = "count") +
    ggplot2::xlab("Number of Individuals in Cluster") +
    ggplot2::ylab("Number of Clusters") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  
  if(save_plots) pdf(file = paste0(run_set_title, "_cluster_sizes.pdf"), width = width, height = height)
  plot(p_cluster_sizes)
  if(save_plots) dev.off()
  
  #find proportion of males in each cluster
  prop_male <- sapply(df_split, FUN = function(x) mean(x$gender == "M"))
  
  #proportion of swedish born people
  prop_swe <- sapply(df_split, FUN = function(x) mean(x$birth_location == "SWE"))
  
  p_gender_dist <- ggplot2::ggplot(data = data.frame(dist = prop_male), ggplot2::aes(x = dist)) + 
    ggplot2::geom_bar(stat = "count") +
    ggplot2::xlab("Proportion Male") +
    ggplot2::ylab("Number of Clusters") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  
  if(save_plots) pdf(file = paste0(run_set_title, "_gender_dist.pdf"), width = width, height = height)
  plot(p_gender_dist)
  if(save_plots) dev.off()
  
  p_gender_dist_clust_size <- ggplot2::ggplot(data = data.frame(dist = prop_male, nInds = nInds), 
                                              ggplot2::aes(x = dist, y = nInds)) +
    ggplot2::geom_jitter(col = "#00000030") +
    ggplot2::xlab("Proportion Male") +
    ggplot2::ylab("Number of Individuals in Cluster") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_log10(breaks = c(2^(0:ceiling(max(nInds)))))
  
  if(save_plots) pdf(file = paste0(run_set_title, "_gender_dist_clust_size.pdf"), width = width, height = height)
  plot(p_gender_dist_clust_size)
  if(save_plots) dev.off()
  
  df$risk_group <- as.factor(df$risk_group)
  
  p_route_dist <- ggplot2::ggplot(df, ggplot2::aes(risk_group)) + 
    ggplot2::geom_bar(stat = "count") +
    ggplot2::xlab("Suspected Transmission Route") +
    ggplot2::ylab("Number of Individuals") +
    ggplot2::theme_bw() + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  
  if(save_plots) pdf(file = paste0(run_set_title, "_route_dist.pdf"), width = width, height = height)
  plot(p_route_dist)
  if(save_plots) dev.off()
  
  #turn "NA"s into "Unknown" (maybe this is not safe if a dataset uses "NA" for "North America"?)
  df$birth_location[df$birth_location == "NA"] <- "unknown"
  df$suspected_infection_location[df$suspected_infection_location == "NA"] <- "unknown"
  
  #maximum number of groups in birth location and transmission location breakdown
  #max_groups <- 10
  #find birth locations and group locations with small numbers of individuals
  birth_table <- table(df$birth_location)
  #find cutoff value for pooling locations with smaller numbers of individuals
  if(length(birth_table) <= max_groups){
    birth_cutoff_count <- 0
  } else{
    birth_cutoff_count <- sort(birth_table[names(birth_table) != "unknown"], decreasing = TRUE)[max_groups]
  }
  other_birth_names <- c(names(which(birth_table <= birth_cutoff_count)), "unknown")
  birth_loc_cutoff <- df$birth_location
  birth_loc_cutoff[which(birth_loc_cutoff %in% other_birth_names)] <- "other"
  
  #find number of each
  birth_table_cutoff <- table(birth_loc_cutoff)
  #put in order from most common to least common
  ordered_levels <- names(birth_table_cutoff)[order(birth_table_cutoff, decreasing = TRUE)]
  #put "other" at the end
  ordered_levels <- ordered_levels[-which(ordered_levels == "other")]
  ordered_levels <- c(ordered_levels, "other")
  #data frame of with ordered factor levels
  df_birth <- data.frame(b = factor(birth_loc_cutoff, 
                                    levels = ordered_levels))
  
  p_birth_locs <- ggplot2::ggplot(data = df_birth, ggplot2::aes(x = b)) + 
    ggplot2::geom_bar(stat = "count") +
    ggplot2::xlab("Birth Location") + 
    ggplot2::ylab("Number of Individuals") + 
    ggplot2::theme_bw() + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  
  if(save_plots) pdf(file = paste0(run_set_title, "_birth_locs.pdf"), width = width, height = height)
  plot(p_birth_locs)
  if(save_plots) dev.off()
  
  inf_table <- table(df$suspected_infection_location)
  #find cutoff value for pooling locations with smaller numbers of individuals
  if(length(inf_table) <= max_groups){
    inf_cutoff_count <- 0
  } else{
    inf_cutoff_count <- sort(inf_table[names(inf_table) != "unknown"], decreasing = TRUE)[max_groups]
  }
  other_inf_names <- c(names(which(inf_table <= inf_cutoff_count)), "unknown")
  inf_loc_cutoff <- df$suspected_infection_location
  inf_loc_cutoff[which(inf_loc_cutoff %in% other_inf_names)] <- "other"
  
  #find number of each
  inf_table_cutoff <- table(inf_loc_cutoff)
  #put in order from most common to least common
  ordered_levels <- names(inf_table_cutoff)[order(inf_table_cutoff, decreasing = TRUE)]
  #put "other" at the end
  ordered_levels <- ordered_levels[-which(ordered_levels == "other")]
  ordered_levels <- c(ordered_levels, "other")
  #data frame of with ordered factor levels
  df_inf <- data.frame(i = factor(inf_loc_cutoff, 
                                  levels = ordered_levels))
  
  p_inf_locs <- ggplot2::ggplot(data = df_inf, ggplot2::aes(x = i)) + 
    ggplot2::geom_bar(stat = "count") +
    ggplot2::xlab("Suspected Infection Location") + 
    ggplot2::ylab("Number of Individuals") + 
    ggplot2::theme_bw() + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.05)), limits = c(0, NA))
  
  if(save_plots) pdf(file = paste0(run_set_title, "_inf_locs.pdf"), width = width, height = height)
  plot(p_inf_locs)
  if(save_plots) dev.off()
  
  return(list(p_cluster_sizes = p_cluster_sizes,
              p_gender_dist = p_gender_dist,
              p_gender_dist_clust_size = p_gender_dist_clust_size,
              p_route_dist = p_route_dist,
              p_birth_locs = p_birth_locs,
              p_inf_locs = p_inf_locs))
}