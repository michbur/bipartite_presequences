
signalHsmm_propep_decision <- function(prot, aa_group, pipar, tpmpar, 
                                od, overall_probs_log, params) {
  if (length(prot) == 1) {
    prot <- strsplit(prot, "")[[1]]
    if ("name" %in% names(attributes(prot)))
      attr(prot, "name") <- "undefined_name"
    if (length(prot) == 1)
      stop("Input sequence is too short.")
  }
  if(!is_protein(prot))
    stop("Atypical aminoacids detected, analysis cannot be performed.")
  
  deg_sample <- as.numeric(biogram::degenerate(toupper(prot)[1L:60], aa_group))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  viterbi_res <- duration_viterbi(deg_sample-1, pipar, tpmpar, od, params)
  viterbi_path <- viterbi_res[["path"]]+1
  c_site <- ifelse(any(viterbi_path == 4), 
                   max(which(viterbi_path == 3)), 
                   length(deg_sample))
  propep_end <- ifelse(any(viterbi_path == 5), 
                   max(which(viterbi_path == 4)), 
                   length(deg_sample))
  #get probabilities of signal peptide model
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  #get probabilities of no signal peptide model
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  
  prob.prop.signal <- viterbi_res[["viterbi"]][propep_end, viterbi_path[propep_end]]
  #get probabilities of no signal peptide model
  prob.prop.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:propep_end], 0)
  prob.prop.total <- exp(prob.prop.signal - prob.prop.non)
  
  res <- list(sp_probability = signalHsmm:::rescale(unname(1 - 1/(1 + prob.total))),
              sp_start = 1,
              sp_end = c_site,
              propep_start = c_site+1,
              propep_end = propep_end,
              propep_probability = signalHsmm:::rescale(unname(1 - 1/(1 + prob.prop.total))),
              struc = viterbi_path,
              prot = toupper(prot[1L:70]),
              name = attr(prot, "name"),
              str_approx = 0)
  class(res) <- "hsmm_pred"
  
  #structure approximation - if atypical (normally negative signal peptide)
  while(!all(1L:4 %in% res[["struc"]])) {
    res[["struc"]] <- c(res[["struc"]], which.min(1L:4 %in% res[["struc"]]))
    res[["str_approx"]] <- res[["str_approx"]] + 1
  }
  
  res
}
