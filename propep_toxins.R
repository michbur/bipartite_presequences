#UniProt query
#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:manual) annotation:(type:propep) AND reviewed:yes


library(signalHsmm)
library(seqinr)
library(reshape2)
all_prots <- read_uniprot("sp_propeptide.txt", c("signal", "propep"), kwds = c("toxin", "poison", "venom"))

#good prots - proteins that have single propeptide starting immediately after signal peptide
good_prots <- all_prots[sapply(all_prots, function(i)
  attr(i, "signal")[2] ==  attr(i, "propep")[1] - 1)]

propeps <- lapply(good_prots, function(i)
  i[attr(i, "propep")[1]:attr(i, "propep")[2]])

matures <- lapply(good_prots, function(i)
  i[attr(i, "propep")[2]:length(i)])

calc_comp <- function(x) {
  res <- data.frame(table(factor(x, levels = a()[-1])))
  colnames(res) <- c("aa", "freq")
  res[["freq"]] <- res[["freq"]]/sum(res[["freq"]])
  res
}

comp_propeps <- t(sapply(propeps, function(i)
  calc_comp(i)[["freq"]]))

comp_matures <- t(sapply(matures, function(i)
  calc_comp(i)[["freq"]]))

mcomp <- melt(data.frame(aa = a()[-1], propep = colMeans(comp_propeps), mature = colMeans(comp_matures)))

library(ggplot2)
ggplot(mcomp, aes(x = aa, fill = variable, y = value)) +
  geom_bar(stat = "identity", position = "dodge")
