#'@title Prepare data for/use multiple biomarker model
#'@description Function to put biomarker and other date into a dataframe, calculate the polypmorphism count,
#' and (possibly) run the multiple biomarker model on it
#'@param patient_ID An identifier for each patient
#'@param last_neg_test_date The last negative test date, if available
#'@param first_pos_test_date The first positive test date. It is overriden by the earliest biomarker sample date if that is earlier
#'or if the first positive test date is unavailable.
#'@param ART_start_date The start date of antiretroviral therapy
#'@param BED_dates A list of vectors for the sample dates of the BED values
#'@param BED A list of vectors for values for the BED test
#'@param LAg_dates A list of vectors for the sample dates of the LAg values
#'@param LAg A list of vectors for values for the LAg test
#'@param CD4_dates A list of vectors for the sample dates of the CD4 counts
#'@param CD4 A list of vectors for the CD4+ T-cell counts
#'@param seq_dates A list of vectors for the sample dates of the HIV sequences
#'@param seq A list of HIV sequences (or a list of lists of sequences) for each patient in character or DNAbin format
#'@param seq_names The names of the sequences
#'@param pol_override A list of vectors for values for the polymorphism count if they have already been calculated externally. 
#'These values will replace the pol values calculated by this function.
#'@param seq_length_override A list of vectors for the sequence lengths 
#'(to be used along with pol_override if they have already been calculated externally)
#'@param pol2_dates A list of vectors for the sample dates of the pol2 values
#'@param pol A list of vectors for values for the pol2 value
#'@param VL_dates A list of vectors for the sample dates of the viral load
#'@param VL A list of vectors for values for the viral load
#'@param ART_cutoff_delay The amount of time (in days) after the start of antiretroviral therapy 
#'after which biomarker values should be disregarded
#'@param cluster_ID An identifier for the transmission cluster that the patient belongs to
#'@param gender A vector for the genders of each patient
#'@param age_at_sampling A list of vectors for the age of each patient at the time of each sequence
#'@param birth_location A vector for the birth location of each patient
#'@param suspected_infection_location A vector for the suspected infection location of each patient
#'@param risk_group A vector for the transmission risk group of each patient
#'@param aids_diagnosis_date A vector for the AIDS diagnosis date of each patient if they have been diagnosed with AIDS
#'@param death_date The date of death of each patient if the patient has died
#'@param date_format A character string specifying the format of the dates (to be passed to as.Date())
#'@param find_infection_age_distributions Whether or not to run the multiple biomarker model (MBM) on the data
#'@param prior.type The type of prior used for the time between infection and diagnosis in the multiple biomarker model.
#'1 indicates a gamma distribution with mean 2 years and standard deviation 1.5 years.
#'2 indicates a continuous uniform distribution with minimum 0 and maximum 12 years.
#'3 indicates a distribution for an individual that is HIV positive but has not developed AIDS symptoms, 
#'assuming that AIDS symptoms develop after a length of time according to a gamma distribution with shape 3.349 and rate 0.327
#'4 indicates a user-supplied distribution from the parameter user.prior.pdf
#'@param user.prior.pdf A pdf to use as the prior distribution for the time between infection and diagnosis in the multiple biomarker model.
#'It must be a list with x and y components corresponding to the time before diagnosis and probability density.
#'@param n.adapt Number of adaptation iterations if the MBM is run
#'@param n.burn Number of burn-in iterations if the MBM is run
#'@param n.iter Number of sampling iterations if the MBM is run
#'@param seed The RNG seed
#'@param ... Additional lists of patient data to put into the dataframe
#'@return A dataframe containing the input data, quality controlled versions of the data, possibly infection age distributions
#'@export
prepare.HIV.data <- function(patient_ID,
                             last_neg_test_date = NULL,
                             first_pos_test_date = NULL,
                             ART_start_date = NULL,
                             BED_dates = NULL,
                             BED = NULL,
                             LAg_dates = NULL,
                             LAg = NULL,
                             CD4_dates = NULL,
                             CD4 = NULL,
                             seq_dates,
                             seqs,
                             seq_names = NULL,
                             pol_override = NULL,
                             seq_length_override = NULL,
                             pol2_dates = NULL,
                             pol2 = NULL, 
                             VL_dates = NULL,
                             VL = NULL,
                             ART_cutoff_delay = 3, #in days
                             cluster_ID = NULL,
                             gender = NULL,
                             age_at_sampling = NULL,
                             birth_location = NULL,
                             suspected_infection_location = NULL,
                             risk_group = NULL,
                             aids_diagnosis_date = NULL,
                             death_date = NULL,
                             date_format = "%Y-%m-%d",
                             find_infection_age_distributions = FALSE,
                             prior.type = 1,
                             user.prior.pdf = list(x = c(0, 10), y = c(1/10, 1/10)),
                             n.adapt = 1e4, 
                             n.burn = 1e5, 
                             n.iter = 1e6,
                             seed = sample(2^31-1, 1),
                             ...){
  
  #find total number of individuals
  total_inds <- length(patient_ID)
  
  #check inputs
  
  #convert patient IDs to characters if they are not already
  patient_ID <- as.character(patient_ID)
  
  #convert dates to years (decimal valued, so July 1st 2000 would be approximately 2000.5 (2000.497))
  if(!is.null(last_neg_test_date)){
    last_neg_test_date <- as.double(as.Date(last_neg_test_date, date_format))/365.25 + 1970
  } else{
    last_neg_test_date <- rep(NA, total_inds)
  }
  if(!is.null(first_pos_test_date)){
    first_pos_test_date <- as.double(as.Date(first_pos_test_date, date_format))/365.25 + 1970
  } else{
    first_pos_test_date <- rep(NA, total_inds)
  }
  if(!is.null(ART_start_date)){
    ART_start_date <- as.double(as.Date(ART_start_date, date_format))/365.25 + 1970
  } else{
    ART_start_date <- rep(NA, total_inds)
  }
  if(!is.null(aids_diagnosis_date)){
    aids_diagnosis_date <- as.double(as.Date(aids_diagnosis_date, date_format))/365.25 + 1970
  } else{
    aids_diagnosis_date <- rep(NA, total_inds)
  }
  if(!is.null(death_date)){
    death_date <- as.double(as.Date(death_date, date_format))/365.25 + 1970
  } else{
    death_date <- rep(NA, total_inds)
  }
  
  #BED sample dates and values
  if(!is.null(BED_dates)){
    BED_dates <- lapply(as.list(BED_dates), 
                        FUN = function(x, date_format) as.double(as.Date(x, date_format))/365.25+1970, 
                        date_format = date_format)
  } else{
    BED_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(BED)){
    BED <- as.list(BED)
  } else{
    BED <- as.list(rep(NA, total_inds))
  }
  
  #LAg sample dates and values
  if(!is.null(LAg_dates)){
    LAg_dates <- lapply(as.list(LAg_dates), 
                        FUN = function(x, date_format) as.double(as.Date(x, date_format))/365.25+1970, 
                        date_format = date_format)
  } else{
    LAg_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(LAg)){
    LAg <- as.list(LAg)
  } else{
    LAg <- as.list(rep(NA, total_inds))
  }
  
  #CD4 sample dates and values
  if(!is.null(CD4_dates)){
    CD4_dates <- lapply(as.list(CD4_dates), 
                        FUN = function(x, date_format) as.double(as.Date(x, date_format))/365.25+1970, 
                        date_format = date_format)
  } else{
    CD4_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(CD4)){
    CD4 <- as.list(CD4)
  } else{
    CD4 <- as.list(rep(NA, total_inds))
  }
  
  #sequence sample dates and values
  if(!is.null(seq_dates)){
    seq_dates <- lapply(as.list(seq_dates), 
                        FUN = function(x, date_format) as.double(as.Date(x, date_format))/365.25+1970, 
                        date_format = date_format)
  } else{
    stop("seq_dates must be provided")
  }
  if(!is.null(seqs)){
    seq <- as.list(seqs)
  } else{
    seqs <- as.list(rep("n", total_inds))
    warning("No sequences found. Using blank sequences.")
  }
  
  if(!is.null(seq_names)){
    seq_names <- as.list(seq_names)
  } else{
    n_seqs <- lapply(seqs, FUN = length)
    #print(n_seqs)
    seq_names <- mapply(FUN = function(a,b) if(length(b) > 0) paste(a,b,sep = "_") else character(0), 
                        a = patient_ID, b = lapply(n_seqs, FUN = seq_len))
  }
  
  #pol2 sample dates and values
  if(!is.null(pol2_dates)){
    pol2_dates <- lapply(as.list(pol2_dates), 
                         FUN = function(x, date_format) as.double(as.Date(x, date_format))/365.25+1970, 
                         date_format = date_format)
  } else{
    pol2_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(pol2)){
    pol2 <- as.list(pol2)
  } else{
    pol2 <- as.list(rep(NA, total_inds))
  }
  
  #viral load sample dates and values
  if(!is.null(VL_dates)){
    VL_dates <- lapply(as.list(VL_dates), 
                       FUN = function(x, date_format) as.double(as.Date(x, date_format))/365.25+1970, 
                       date_format = date_format)
  } else{
    VL_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(VL)){
    VL <- as.list(VL)
  } else{
    VL <- as.list(rep(NA, total_inds))
  }
  
  #cluster_ID
  if(!is.null(cluster_ID)){
    cluster_ID <- as.character(cluster_ID)
  } else{
    stop("cluster_ID must be provided")
  }
  
  #gender
  if(!is.null(gender)){
    gender <- as.character(gender)
  } else{
    gender <- as.character(rep(NA, total_inds))
  }
  
  #age at the time of sequence samples
  if(!is.null(age_at_sampling)){
    age_at_sampling <- as.list(age_at_sampling)
  } else{
    age_at_sampling <- as.list(rep(NA, total_inds))
  }
  
  #birth location
  if(!is.null(birth_location)){
    birth_location <- as.character(birth_location)
  } else{
    birth_location <- as.character(rep(NA, total_inds))
  }
  
  #suspected infection location
  if(!is.null(suspected_infection_location)){
    suspected_infection_location <- as.character(suspected_infection_location)
  } else{
    suspected_infection_location <- as.character(rep(NA, total_inds))
  }
  
  #risk group
  if(!is.null(risk_group)){
    risk_group <- as.character(risk_group)
  } else{
    risk_group <- as.character(rep(NA, total_inds))
  }
  
  #use the earliest sampled biomarker as the first positive date if it is blank or after biomarker measurements (add VL dates?)
  first_pos_test_date_adj <- mapply(FUN = min, 
                                    first_pos_test_date, 
                                    ART_start_date,
                                    BED_dates, 
                                    LAg_dates,
                                    CD4_dates,
                                    seq_dates, 
                                    pol2_dates, 
                                    aids_diagnosis_date, 
                                    na.rm = TRUE)
  
  #calculate polymorphism count from sequences if not provided
  if(is.null(pol_override)){
    #collapse inner lists in seqs 
    #TODO check more about what format seqs is in
    seqs <- lapply(seqs, FUN = function(x) if(is.list(x) & class(x) != "phyDat") do.call("rbind", x) else x)
    #function to calculate polymorphism counts
    calculate_pol <- function(sequences){
      if(is.null(dim(sequences))){
        nseq <- 1
      } else{
        nseq <- dim(sequences)[1]
      }
      length_seq <- integer(length = nseq)
      pol <- double(length = nseq)
      for(i in seq_len(nseq)){
        frequencies <- table(as.character(sequences[i,]))
        polymorphic_count <- sum(frequencies[names(frequencies) != "-" & names(frequencies) != "n" & names(frequencies) != "a" & 
                                               names(frequencies) != "c" & names(frequencies) != "g" & names(frequencies) != "t"])
        length_seq[i] <- sum(frequencies[names(frequencies) != "-" & names(frequencies) != "n"]) #length of sequence
        pol[i] <- polymorphic_count/length_seq[i]
      }
      return(list(length = length_seq, pol = pol))
    }
    length_and_pol <- lapply(seqs, FUN = calculate_pol)
    seq_length <- lapply(length_and_pol, FUN = function(x) x$length)
    pol <- lapply(length_and_pol, FUN = function(x) x$pol)
  } else{ #override pol values if provided
    pol <- as.list(pol_override)
    if(!is.null(seq_length_override)){
      if(length(seq_length_override) == total_inds){
        seq_length = as.list(seq_length_override)
      } else{
        stop("If provided, seq_length_override must be the same length as the number of individuals")
      }
    } else{
      seq_length <- as.list(rep(NA, total_inds))
    }
  }
  
  #find number of sequences per individual
  nSeq <- sapply(pol, FUN = length)
  
  #convert ART cutoff delay from days to years
  ART_cutoff_delay_years <- ART_cutoff_delay/365.25
  #remove biomarker values that are taken too late after the start of ART
  
  #function to remove values that are unusable due to being too long after ART start
  remove_late_samples <- function(samples, sample_dates, ART_dates, ART_cutoff_delay_years){
    samples[sample_dates > ART_dates+ART_cutoff_delay_years] <- NA
    return(samples)
  }
  
  BED_qc <- mapply(remove_late_samples,
                   samples = BED, sample_dates = BED_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  LAg_qc <- mapply(remove_late_samples, 
                   samples = LAg, sample_dates = LAg_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  CD4_qc <- mapply(remove_late_samples, 
                   samples = CD4, sample_dates = CD4_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  BED_qc <- mapply(remove_late_samples, 
                   samples = BED, sample_dates = BED_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  pol_qc <- mapply(remove_late_samples, 
                   samples = pol, sample_dates = seq_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  pol2_qc <- mapply(remove_late_samples, 
                    samples = pol2, sample_dates = pol2_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                    SIMPLIFY = FALSE)
  VL_qc <- mapply(remove_late_samples, 
                  samples = VL, sample_dates = VL_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                  SIMPLIFY = FALSE)
  
  #is there at least one usable pol value
  usable_pol <- sapply(pol_qc, FUN = function(x) !all(is.na(x)))
  
  #make sure other things are character vectors
  cluster_ID <- as.character(cluster_ID)
  gender <- as.character(gender)
  birth_location <- as.character(birth_location)
  suspected_infection_location <- as.character(suspected_infection_location)
  risk_group <- as.character(risk_group)
  
  #put everything into a data.frame
  df <- data.frame(patient_ID = patient_ID,
                   last_neg_test_date = last_neg_test_date,
                   first_pos_test_date = first_pos_test_date,
                   first_pos_test_date_adj = first_pos_test_date_adj,
                   ART_start_date = ART_start_date,
                   BED_dates = I(BED_dates),
                   BED = I(BED),
                   BED_qc = I(BED_qc),
                   LAg_dates = I(LAg_dates),
                   LAg = I(LAg),
                   LAg_qc = I(LAg_qc),
                   CD4_dates = I(CD4_dates),
                   CD4 = I(CD4),
                   CD4_qc = I(CD4_qc),
                   seq_dates = I(seq_dates),
                   seqs = I(seqs),
                   seq_names = I(seq_names),
                   seq_length = I(seq_length),
                   nSeq = nSeq,
                   pol = I(pol),
                   pol_qc = I(pol_qc),
                   usable_pol = usable_pol,
                   pol2_dates = I(pol2_dates),
                   pol2 = I(pol2), 
                   pol2_qc = I(pol2_qc),
                   VL_dates = I(VL_dates),
                   VL = I(VL),
                   VL_qc = I(VL_qc),
                   ART_cutoff_delay_years = rep(ART_cutoff_delay_years, length = total_inds), #in days
                   cluster_ID = cluster_ID,
                   gender = gender,
                   age_at_sampling = I(age_at_sampling),
                   birth_location = birth_location,
                   suspected_infection_location = suspected_infection_location,
                   risk_group = risk_group,
                   aids_diagnosis_date = aids_diagnosis_date,
                   death_date = death_date,
                   stringsAsFactors = FALSE)
  
  #add additional user supplied columns
  additional_input <- list(...)
  if("temp_col_name" %in% names(additional_input)) stop("Please use a name other than 'temp_col_name' for additional input columns")
  #return(additional_input)
  for(i in seq_along(additional_input)){
    df$temp_col_name <- I(additional_input[[i]])
    names(df)[names(df) == "temp_col_name"] <- names(additional_input)[i]
  }
  
  if(find_infection_age_distributions){
    #find probability distributions for the time between infection and diagnosis
    infection_age_dists <- run.mbm(df, n.adapt = n.adapt, n.burn = n.burn, n.iter = n.iter, 
                                   prior.type = prior.type, user.prior.pdf = user.prior.pdf,
                                   overall.seed = seed)
    
    #find pdf and cdf in terms of real time
    cdfs <- lapply(infection_age_dists$infection_age_dists_diag, FUN = function(dist){find_cdf(-rev(dist$x), rev(dist$y))})
    real_times <- mapply(FUN = function(diag, dist){rev(diag - dist$x)}, 
                         diag = df$first_pos_test_date_adj, dist = infection_age_dists$infection_age_dists_diag, SIMPLIFY = FALSE)
    #put times and (forward time) cdfs together
    cdf <- vector(mode = "list", length = length(cdfs))
    for(i in seq_along(cdfs)){
      cdf[[i]]$x <- real_times[[i]]
      cdf[[i]]$y <- cdfs[[i]]
    }
    #put times and (forward time) pdfs together
    pdf <- vector(mode = "list", length = length(cdfs))
    for(i in seq_along(cdfs)){
      pdf[[i]]$x <- real_times[[i]]
      pdf[[i]]$y <- rev(infection_age_dists$infection_age_dists_diag[[i]]$y)
    }
    
    #put distributions into dataframe
    #distributions from diagnosis
    df$infection_age_dists_diag <- I(infection_age_dists$infection_age_dists_diag)
    #distributions from first sequence
    df$infection_age_dists_seq <- I(infection_age_dists$infection_age_dists_seq)
    #forwards, real time pdf
    df$pdf <- I(pdf)
    #forwards, real time cdf
    df$cdf <- I(cdf)
  }
  
  return(df)
}

#'@title Wrapper for mbm.predict
#'@description Wrapper function to call mbm.predict in parallel
#'@param df A data frame containing biomarker information as in the output from prepare.HIV.data
#'@param n.adapt Number of iterations for model adaptation
#'@param n.burn Number of burn-in iterations
#'@param n.iter Number of sampling iterations
#'#'@param prior.type The type of prior used for the time between infection and diagnosis in the multiple biomarker model.
#'1 indicates a gamma distribution with mean 2 years and standard deviation 1.5 years.
#'2 indicates a continuous uniform distribution with minimum 0 and maximum 12 years.
#'3 indicates a distribution for an individual that is HIV positive but has not developed AIDS symptoms, 
#'assuming that AIDS symptoms develop after a length of time according to a gamma distribution with shape 3.349 and rate 0.327
#'4 indicates a user-supplied distribution from the parameter user.prior.pdf
#'@param user.prior.pdf A pdf to use as the prior distribution for the time between infection and diagnosis in the multiple biomarker model.
#'It must be a list with x and y components corresponding to the time before diagnosis and probability density.
#'@param overall.seed RNG seed to use to generate other RNG seeds
#'@return Lists for the infection age distributions both in terms of time before diagnosis and time before the first sequence
#'@export
run.mbm <- function(df, n.adapt = 1000, n.burn = 1000, n.iter = 1000, 
                    prior.type = 1, user.prior.pdf = list(x = c(0, 10), y = c(1/10, 1/10)),
                    overall.seed = sample(2^31-1, 1)){
  #find number of individuals
  nInds <- nrow(df)
  
  set.seed(overall.seed)
  seeds <- sample(2^31-1, nInds)
  
  #find number of measurement for each biomarker for each patient
  mBED <- sapply(df$BED_qc, FUN = length)
  mLAg <- sapply(df$LAg_qc, FUN = length)
  mCD4 <- sapply(df$CD4_qc, FUN = length)
  mpol <- sapply(df$pol_qc, FUN = length)
  mpol2 <- sapply(df$pol2_qc, FUN = length)
  
  #find maximum m for each individual (no longer needed)
  #m <- mapply(FUN = max, mBED, mLAg, mCD4, mpol, mpol2)
  
  #calculate delays between infection and sampling
  BED_delays <- mapply(FUN = `-`, df$BED_dates, df$first_pos_test_date_adj, SIMPLIFY = FALSE)
  LAg_delays <- mapply(FUN = `-`, df$LAg_dates, df$first_pos_test_date_adj, SIMPLIFY = FALSE)
  CD4_delays <- mapply(FUN = `-`, df$CD4_dates, df$first_pos_test_date_adj, SIMPLIFY = FALSE)
  pol_delays <- mapply(FUN = `-`, df$seq_dates, df$first_pos_test_date_adj, SIMPLIFY = FALSE)
  pol2_delays <- mapply(FUN = `-`, df$pol2_dates, df$first_pos_test_date_adj, SIMPLIFY = FALSE)
  
  #set delays to 0 if the corresponding biomarker value is NA
  #function to test if biomarker values are NA and set delays if necessary
  find_biomarker_NAs <- function(values, delays){
    #check to make sure they are the same length
    if(length(values) != length(delays)) stop("The number of biomarker samples and sample dates must match")
    #find indices where both biomarker values and biomarker dates are NA
    both_NA <- (is.na(values) & is.na(delays))
    #change delays where both are NA to 0
    delays[both_NA] <- 0
    return(delays)
  }
  
  #apply function
  BED_delays <- mapply(FUN = find_biomarker_NAs, df$BED_qc, BED_delays, SIMPLIFY = FALSE)
  LAg_delays <- mapply(FUN = find_biomarker_NAs, df$LAg_qc, LAg_delays, SIMPLIFY = FALSE)
  CD4_delays <- mapply(FUN = find_biomarker_NAs, df$CD4_qc, CD4_delays, SIMPLIFY = FALSE)
  pol_delays <- mapply(FUN = find_biomarker_NAs, df$pol_qc, pol_delays, SIMPLIFY = FALSE)
  pol2_delays <- mapply(FUN = find_biomarker_NAs, df$pol2_qc, pol2_delays, SIMPLIFY = FALSE)
  
  #calculate time between last negative and first positive test
  last_neg_first_pos_diff <- mapply(FUN = `-`, df$first_pos_test_date_adj, df$last_neg_test_date, SIMPLIFY = FALSE)
  
  #function to put biomarker values into matrices with NAs for missing values
  biomarker_mat <- function(values, m){
    #number of NA values to pad
    n_NA_pad <- m - length(values)
    #put in matrix format for multiple biomarker model
    mat <- matrix(c(values, rep(NA, n_NA_pad)), nrow = 1, ncol = m)
  }
  
  #put all biomarkers into matrices that have the same shape (on a per individual basis)
  t_BED_mats <- mapply(FUN = biomarker_mat, values = BED_delays, m = mBED, SIMPLIFY = FALSE)
  BED_mats <- mapply(FUN = biomarker_mat, values = df$BED_qc, m = mBED, SIMPLIFY = FALSE)
  
  t_LAg_mats <- mapply(FUN = biomarker_mat, values = LAg_delays, m = mLAg, SIMPLIFY = FALSE)
  LAg_mats <- mapply(FUN = biomarker_mat, values = df$LAg_qc, m = mLAg, SIMPLIFY = FALSE)
  
  t_CD4_mats <- mapply(FUN = biomarker_mat, values = CD4_delays, m = mCD4, SIMPLIFY = FALSE)
  CD4_mats <- mapply(FUN = biomarker_mat, values = df$CD4_qc, m = mCD4, SIMPLIFY = FALSE)
  
  t_pol_mats <- mapply(FUN = biomarker_mat, values = pol_delays, m = mpol, SIMPLIFY = FALSE)
  pol_mats <- mapply(FUN = biomarker_mat, values = df$pol_qc, m = mpol, SIMPLIFY = FALSE)
  
  t_pol2_mats <- mapply(FUN = biomarker_mat, values = pol2_delays, m = mpol2, SIMPLIFY = FALSE)
  pol2_mats <- mapply(FUN = biomarker_mat, values = df$pol2_qc, m = mpol2, SIMPLIFY = FALSE)
  
  #load trained MBM parameters
  data("MBM_pars")
  
  #apply mbm to each individual
  infection_age_dists <- parallel::mcmapply(FUN = mbm.predict,
                                            BED = BED_mats, 
                                            LAg = LAg_mats, 
                                            CD4 = CD4_mats, 
                                            pol = pol_mats, 
                                            pol2 = pol2_mats,
                                            prev.neg.time = last_neg_first_pos_diff, 
                                            t.BED.delay = t_BED_mats,
                                            t.LAg.delay = t_LAg_mats,
                                            t.CD4.delay = t_CD4_mats,
                                            t.pol.delay = t_pol_mats,
                                            t.pol2.delay = t_pol2_mats,
                                            mub = rep(list(MBM_pars$mub), nInds), 
                                            Sigmab = rep(list(MBM_pars$Sigmab), nInds), 
                                            sigmae = rep(list(MBM_pars$sigmae), nInds),
                                            n.adapt = n.adapt, n.burn = n.burn, n.iter = n.iter,
                                            prior.type = prior.type,
                                            inf.mean = 2, inf.sd = 1.5, 
                                            max.seroconvert.delay = 2/12,
                                            u1.pdf = list(user.prior.pdf),
                                            seed = seeds,
                                            output.raw = FALSE, 
                                            SIMPLIFY = FALSE, 
                                            mc.cores = parallel::detectCores())
  
  #extract numeric kernel density estimates from output (time before diagnosis)
  infection_age_dists_diag <- lapply(infection_age_dists, FUN = function(x) x$pdf.num$t_pred)
  #extract numeric kernel density estimates from output (time before first sequence)
  infection_age_dists_seq <- lapply(infection_age_dists, FUN = function(x) x$pdf.num.seq$t_pred)
  
  #return distributions from infection times and from sampling times
  return(list(infection_age_dists_diag = infection_age_dists_diag, 
              infection_age_dists_seq = infection_age_dists_seq))
}

#'@title Create list of inputs for run.biophybreak
#'@description Function to create a list of inputs to use with run.biophybreak
#'@param df A dataframe as output from prepare.HIV.data
#'@param wh.model The model for within-host pathogen dynamics (effective pathogen population size = 
#'   N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Names and numbers are allowed.
#'   Options are:
#'   \enumerate{
#'     \item "single": effective size = 0, so coalescence occurs 'just before' transmission in the infector (complete bottleneck)
#'     \item "infinite": effective size = Inf, with complete bottleneck, so coalescence occurs 'just after' transmission in the infectee
#'     \item "linear": effective size at time t after infection = \code{wh.level + wh.slope * t} (complete or wide bottleneck; if complete, \code{wh.level = 0})
#'     \item "exponential": effective size at time t after infection = \code{wh.level * exp(wh.exponent * t)} (wide bottleneck)
#'     \item "constant": effective size = wh.level (wide bottleneck)
#'   }
#'@param wh.bottleneck Whether the bottleneck should be complete or wide, which is only an option if \code{wh.model = "linear"} 
#'   (in that case, \code{"auto"} defaults to \code{"complete"}).
#'@param prior.wh.level.mean Mean of the (gamma) prior distribution of \code{wh.level} 
#'@param gen.type The type of input for the prior for the generation time. Can be 'parametric' (uses gamma distribution), 
#'   'numeric' for use with a user-supplied density, or 'function' if the user supplied density is already a function.
#'@param gen.density Either a numeric distribution with x and y attributes, such as the output of the density() function
#'   (if gen.type is 'numeric'), or a density function (if 'gen.type' is 'function')
#'@param post.sam.trans.rate The rate of transmission after sampling relative to transmission prior to sampling.
#'   "1" means that there is it is not any less likely for an individual to infect another individual after being sampled.
#'   "0" means that it is not possible for an individual to infect another after being sampled.
#'   "0.1" means that an individual is 10 times less likely to infect another after being sampled. 
#'@param est.gen.mean Whether to estimate the mean generation interval or keep it fixed.
#'@param use.pml Whether to use phangorn::pml to calculate sequences likelihoods.
#'@param run.names The names of the clusters to be used to saving files with run.biophybreak
#'@param name.prepend An identifier to add to the beginning of the file names for the output
#'@param overall.seed The seed used to generate the individual run seeds.
#'@export
prepare.biophybreak.data <- function(df, 
                                     wh.model = "linear",
                                     wh.bottleneck = "wide",
                                     prior.wh.level.mean = 0.1,
                                     gen.type = "function",
                                     gen.density = stepfun(c(0, 0.5), c(0, 3, 1)),
                                     post.sam.trans.rate = 1,
                                     est.gen.mean = FALSE,
                                     use.pml = TRUE,
                                     mut.model = "pol",
                                     run.names = "automatic",
                                     name.prepend = "",
                                     overall.seed = sample(2^31-1,1),
                                     ...){
  
  #set the RNG seed
  set.seed(overall.seed)
  
  #names of clusters
  cluster_names <- sort(unique(df$cluster_ID))
  #order numerically instead of alphanumerically if names are all numbers
  if(suppressWarnings(all(!is.na(as.numeric(cluster_names))))){
    cluster_names <- cluster_names[order(as.numeric(cluster_names))]
  }
  #find the number of clusters
  nClust <- length(cluster_names)
  
  #check name of runs
  if(run.names != "automatic"){
    #check to make sure there are the right number of entries
    if(length(run.names) != nCLust){
      stop("If run.names is provided, it must be the same length as the number of clusters")
    }
  }
  
  #find indices of individuals in each cluster
  clust_indiv_indices <- lapply(cluster_names, FUN = function(x, y) which(x == y), y = df$cluster_ID)
  
  #find number of individuals in each cluster
  nInds <- sapply(clust_indiv_indices, FUN = length)
  
  #get data into right format for phybreakdata function
  process_seqs <- function(indices, seqs){
    if(is.list(seqs[indices]) && (class(seqs[[1]]) != "phyDat") && (class(seqs[[1]][[1]]) != "phyDat")){
      do.call("rbind", seqs[indices]) 
    } else if(class(seqs[[1]]) == "phyDat"){
      seqs <- lapply(seqs, as.character)
      do.call("rbind", seqs[indices]) 
    } else if(class(seqs[[1]][[1]]) != "phyDat"){
      #not sure this will work correctly in all situations
      seqs <- rapply(seqs, as.character)
      seqs <- lapply(seqs, FUN = function(x) if(is.list(x)) do.call("rbind", x) else x)
      seqs <- do.call("rbind", seqs[indices])
    } else{
      seqs[indices]
    } 
  }
  sequences <- lapply(clust_indiv_indices, 
                      FUN = process_seqs,
                      seqs = df$seqs)
  #return(sequences)
  #sequence sample times
  seq_sample_times <- lapply(clust_indiv_indices,
                             FUN = function(indices, x) unname(unlist(x[indices])),
                             x = df$seq_dates)
  #sequences names
  seq_names <- unname(lapply(clust_indiv_indices,
                             FUN = function(indices, x) unname(unlist(x[indices])),
                             x = df$seq_names))
  #names of the hosts of each sequences
  seq_host_names <- vector(mode = "list", length = nClust)
  for(i in seq_len(nClust)){
    seq_host_names[[i]] <- unname(unlist(mapply(rep, x = df$patient_ID[clust_indiv_indices[[i]]], df$nSeq[clust_indiv_indices[[i]]])))
  }
  
  #list for data for each cluster
  biophybreakdata <- vector(mode = "list", length = nClust)
  for(i in seq_len(nClust)){
    #create phybreakdata objects
    biophybreakdata[[i]] <- phybreakdata(sequences = sequences[[i]],
                                         sample.times = seq_sample_times[[i]], 
                                         sample.names = seq_names[[i]],
                                         host.names = seq_host_names[[i]])
  }
  
  #find order of individuals in terms of first sequence time
  first_seq_order <- lapply(clust_indiv_indices, FUN = function(indices, x) order(sapply(x[indices], FUN = min)), x = df$seq_dates)
  
  #order infection age distributions
  sample.density <- vector(mode = "list", length = nClust)
  for(i in seq_len(nClust)){
    sample.density[[i]] <- df$infection_age_dists_seq[clust_indiv_indices[[i]][first_seq_order[[i]]]]
  }
  
  #find which clusters have at least 2 individuals
  nInds2p <- which(nInds >= 2)
  
  #roll individual seeds
  seeds <- sample(2^31-1,length(nInds2p))
  
  #put into input list
  inputs <- vector(mode = "list", length = length(nInds2p))
  run_name <- vector(mode = "character", length = length(nInds2p))
  #vector for pml in case values need to be changed
  pml.vec <- rep(use.pml, length(nInds2p))
  
  #make input lists
  for(i in seq_along(nInds2p)){
    if(run.names == "automatic"){
      run_name[i] <- paste0(name.prepend, "_", df$cluster_ID[clust_indiv_indices[[nInds2p[i]]][1]], "_", seeds[i])
    } else{
      run_name[i] <- paste0(name.prepend, "_", run.names[nInds2p[i]], "_", seeds[i])
    }
    #check if sequences are functionally identical
    SNPpatterns <- do.call(rbind, biophybreakdata[[nInds2p[i]]]$sequences)
    nSNPs <- as.integer(
      sum(apply(SNPpatterns, 2, 
                function(x) {
                  max(0, length(unique(x[x < 5])) - 1)
                }
      ) * attr(biophybreakdata[[nInds2p[i]]]$sequences, "weight")
      )
    )
    if(nSNPs == 0){
      pml.vec[i] <- FALSE #don't use pml for seq likelihood if seqs are functionally identical
      message(paste0("Not using pml on ", run_name[i], " because there are no SNPs"))
    } 
    inputs[[i]] <- list(biophybreakdata = biophybreakdata[[nInds2p[i]]], 
                        sample.density = sample.density[[nInds2p[i]]], 
                        run_name = run_name[i],
                        nInds = nInds[nInds2p[i]],
                        seed = seeds[i],
                        wh.model = wh.model,
                        wh.bottleneck = wh.bottleneck,
                        prior.wh.level.mean = prior.wh.level.mean,
                        gen.type = gen.type,
                        gen.density = gen.density,
                        post.sam.trans.rate = post.sam.trans.rate,
                        est.gen.mean = est.gen.mean,
                        use.pml = pml.vec[i],
                        mut.model = mut.model)
  }
  
  return(inputs)
}

#'@title Create phybreak object and run burn-in and sampling
#'@description Wrapper function to facilitate running inference with biophybreak
#'@param inputs A list containing data and parameters to the phybreak function that may change between runs
#'@param burn_iter The number of iterations for burn-in
#'@param iter_max The maximum number of sampling iterations
#'@param iter_block How frequently to check for sufficient effective sample size
#'@param thin Amount of thinning for MCMC samples
#'@param ESS_target Target effective sample size for all sampled parameters. Sampling will stop if this target is reached.
#'@param save Whether the MCMC output should be saved to a file
#'@param outdir The location at which the MCMC output should be saved
#'@return A list containing the MCMC output, the input list, the effective sample sizes over time, 
#'the mixing status, the run name, and the computation time
#'@export
run.biophybreak <- function(inputs, 
                            burn_iter = 10000, 
                            iter_max = 100000, 
                            iter_block = 10000, 
                            thin = 1, 
                            ESS_target = 200,
                            save = TRUE,
                            outdir = "",
                            ...){
  
  #check values in inputs
  if(!("biophybreakdata" %in% names(inputs))) stop("biophybreakdata must be present in inputs")
  if(!("sample.density" %in% names(inputs))) stop("sample.density must be present in inputs")
  if(!("wh.model" %in% names(inputs))) inputs$wh.model <- "linear"
  if(!("wh.bottleneck" %in% names(inputs))) inputs$wh.bottleneck <- "wide"
  if(!("prior.wh.level.mean" %in% names(inputs))) inputs$prior.wh.level.mean <- 0.1
  if(!("gen.type" %in% names(inputs))) inputs$gen.type <- "function"
  if(!("gen.density" %in% names(inputs))) inputs$gen.density <- stepfun(c(0, 0.5), c(0, 3, 1))
  if(!("post.sam.trans.rate" %in% names(inputs))) inputs$post.sam.trans.rate <- 1
  if(!("est.gen.mean" %in% names(inputs))) inputs$est.gen.mean <- FALSE
  if(!("use.pml" %in% names(inputs))) inputs$use.pml <- TRUE
  if(!("mut.model" %in% names(inputs))) inputs$mut.model <- "pol"
  #run name
  if(!("run_name" %in% names(inputs))) stop("run_name must be present in inputs")
  #RNG seed
  if(!("seed" %in% names(inputs))) inputs$seed <- sample(2^31-1,1)
  
  #check to make sure there are at least 2 individuals
  if(length(unique(inputs$biophybreakdata$sample.hosts)) < 2) stop("biophybreakdata must contain at least 2 individuals")
  
  #set RNG seed
  set.seed(inputs$seed)
  MCMCstate <- phybreak(dataset = inputs$biophybreakdata, 
                        sample.nonpar = TRUE, 
                        sample.density = inputs$sample.density, 
                        wh.model = inputs$wh.model, 
                        wh.bottleneck = inputs$wh.bottleneck, 
                        prior.wh.level.mean = inputs$prior.wh.level.mean,
                        est.sample.mean = FALSE, 
                        gen.type = inputs$gen.type, 
                        gen.density = inputs$gen.density,
                        post.sam.trans.rate = inputs$post.sam.trans.rate,
                        est.gen.mean = inputs$est.gen.mean,
                        use.pml = inputs$use.pml, 
                        mut.model = inputs$mut.model,
                        ...)
  
  #initialize time counting
  time <- system.time(0)
  #burn-in
  i <- 0
  while(((i*iter_block) < burn_iter)){
    iters <- min(burn_iter-i*iter_block, iter_block)
    tempseed <- .Random.seed #store RNG state for failsafe
    j <- 0
    repeat{
      j <- j + 1
      #restore RNG state
      .Random.seed <- tempseed
      time <- time + system.time(MCMCstate_temp <- tryCatch(burnin_phybreak(MCMCstate, ncycles = iters), 
                                                            error = function(e) e))
      if(class(MCMCstate_temp)[1] == "phybreak") MCMCstate <- MCMCstate_temp; break
      if(j == 100) stop(MCMCstate)
    }
    i <- i+1
  }
  ESS <- list()
  ESS[[1]] <- 0
  i <- 0
  mixed <- FALSE
  #iter_block <- round(iter/150)
  while(((i*iter_block) < iter_max) & mixed == FALSE){
    iters <- min(iter_max-i*iter_block, iter_block)
    tempseed <- .Random.seed #store RNG state for failsafe
    j <- 0
    repeat{
      j <- j + 1
      #restore RNG state
      .Random.seed <- tempseed
      time <- time + system.time(MCMCstate_temp <- tryCatch(sample_phybreak(MCMCstate, nsample = iters, thin = thin), 
                                                            error = function(e) e))
      #break if successful
      if(class(MCMCstate_temp)[1] == "phybreak") MCMCstate <- MCMCstate_temp; break
      if(j == 100) stop(MCMCstate)
    }
    i <- i+1
    ESS[[i]] <- biophybreak::ESS(MCMCstate)
    mixed_params <- (ESS[[i]] >= ESS_target)
    if(all(mixed_params[!is.na(mixed_params)])){ #are all non-NA ESSs over the target value?
      mixed = TRUE
    }
    output <- list(MCMCstate = MCMCstate, 
                   inputs = inputs, 
                   ESS = ESS,
                   mixed = mixed,
                   run_name = inputs$run_name,
                   time = time)
    if(save) save(output, file = paste0(outdir, "/partial_", inputs$run_name, ".Rdata"))
  }
  output <- list(MCMCstate = MCMCstate, 
                 inputs = inputs, 
                 ESS = ESS,
                 mixed = mixed,
                 run_name = inputs$run_name,
                 time = time)
  if(save){
    #save final file
    save(output, file = paste0(outdir, "/", inputs$run_name, ".Rdata"))
    #remove partial run file
    file.remove(paste0(outdir, "/partial_", inputs$run_name, ".Rdata"))
  } 
  return(output)
}