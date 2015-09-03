

train_propep_hsmm <- function(train_data, aa_group, max_length = 32) {
  train_data <- lapply(train_data, toupper)
  ts <- signalHsmm:::calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  
  propeps <- lapply(train_data, function(i)
    i[attr(i, "propep")[1]:attr(i, "propep")[2]])
  
  t4 <- rep(0, length(aa_group))
  temp <- table(biogram:::degenerate(unlist(propeps), aa_group))
  t4[as.numeric(names(temp))] <- temp
  names(t4) <- 1:length(aa_group)
  
  matures <- lapply(train_data, function(i)
    i[attr(i, "propep")[2]:length(i)])
  t5 <- rep(0, length(aa_group))
  temp <- table(biogram:::degenerate(unlist(matures), aa_group))
  t5[as.numeric(names(temp))] <- temp
  names(t5) <- 1:length(aa_group)
  
  
  overall <- t5 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall_probs_log = log(overall.probs) #for viterbi
  
  lengths <- ts[["lengths"]]
  propepDensity <- signalHsmm:::measure_region(sapply(train_data, function(x){ y=attr(x, "propep"); y[2]-y[1]+1}), max_length)
  params <- apply(lengths, 2, signalHsmm:::measure_region, max_length = max_length)
  params <- cbind(params, propepDensity, rep(1/max_length, max_length))
  
  #setting params for hmm -------
  ngroups <- length(aa_group)
  additional_margin = 10
  pipar <- c(1,0,0,0,0)
  tpmpar <- matrix(c(0, 1, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 1, 0,
                     0, 0, 0, 0, 1,
                     0, 0, 0, 0, 0), 5, byrow = TRUE)
  od <- matrix(c((t1/sum(t1))[1L:ngroups],
                 (t2/sum(t2))[1L:ngroups],
                 (t3/sum(t3))[1L:ngroups],
                 (t4/sum(t4))[1L:ngroups],
                 (t5/sum(t5))[1L:ngroups]), 5, byrow = TRUE)
  
  res <- list(aa_group = aa_group, pipar = pipar, tpmpar = tpmpar, od = od, 
              overall_probs_log = overall_probs_log, params = params)
  class(res) <- "sighsmm_model"
  res
}

