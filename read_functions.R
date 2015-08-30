#UniProt query
#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:manual) annotation:(type:propep) AND reviewed:yes


library(signalHsmm)
all_prots <- read_uniprot("Csp_propeptide.txt", c("signal", "propep"))

#good prots - proteins that have single propeptide starting immediately after signal peptide
good_prots <- all_prots[sapply(prots, function(i)
  attr(i, "signal")[2] ==  attr(i, "propep")[1] - 1)]



