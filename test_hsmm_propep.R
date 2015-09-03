
# mini test

library(signalHsmm)

all_prots <- read_uniprot("sp_propeptide.txt", c("signal", "propep"), kwds = c("toxin", "poison", "venom"))

# wybór zbioru uczączego
good_prots <- all_prots[sapply(all_prots, function(i)
  attr(i, "signal")[2] ==  attr(i, "propep")[1] - 1)]
propepLength <- sapply(good_prots, function(x){ y=attr(x, "propep"); y[2]-y[1]+1})
train_data <- good_prots[propepLength<=max_length]
test_idx <- sample(1:length(train_data), 100)
test_data <- train_data[test_idx]

max_length = 32
# podstawowe kodowanie z signalHsmm
aa_list <- signalHsmm:::signalHsmm_main_model$aa_group

res <- train_propep_hsmm(train_data, aa_list, max_length)


decisions <- lapply(test_data, function(prot) 
  signalHsmm_propep_decision(prot, aa_group = res[["aa_group"]], 
                           pipar = res[["pipar"]], 
                           tpmpar = res[["tpmpar"]], 
                           od = res[["od"]], 
                           overall_probs_log = res[["overall_probs_log"]], 
                           params = res[["params"]]))

propep_diff <- sapply(test_data, function(x) attr(x, "propep")[2]) - sapply(decisions, function(x) x$proppep_end)
median(abs(propep_diff))
mean(abs(propep_diff))


# nowe kodowanie
aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h", "d", "e"), 
           `3` = c("f", "w", "y", "s", "t", "c", "n", "q"))

res <- train_propep_hsmm(train_data, aa1, max_length)

decisions_aa1 <- lapply(test_data, function(prot) 
  signalHsmm_propep_decision(prot, aa_group = res[["aa_group"]], 
                             pipar = res[["pipar"]], 
                             tpmpar = res[["tpmpar"]], 
                             od = res[["od"]], 
                             overall_probs_log = res[["overall_probs_log"]], 
                             params = res[["params"]]))

propep_diff <- sapply(test_data, function(x) attr(x, "propep")[2]) - sapply(decisions_aa1, function(x) x$proppep_end)
median(abs(propep_diff))
mean(abs(propep_diff))


# nowe kodowanie
aa2 = list(`1` = c("g", "p", "v", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "y", "w"),
           `5` = c("a", "l", "q"))

res <- train_propep_hsmm(train_data, aa2, max_length)

decisions_aa2 <- lapply(test_data, function(prot) 
  signalHsmm_propep_decision(prot, aa_group = res[["aa_group"]], 
                             pipar = res[["pipar"]], 
                             tpmpar = res[["tpmpar"]], 
                             od = res[["od"]], 
                             overall_probs_log = res[["overall_probs_log"]], 
                             params = res[["params"]]))

propep_diff <- sapply(test_data, function(x) attr(x, "propep")[2]) - sapply(decisions_aa2, function(x) x$proppep_end)
median(abs(propep_diff))
mean(abs(propep_diff))

# czy to w ogóle ma sens. Porownajmy p-stwa z sekwencjami bez propeptydu

negative_decision <- lapply(benchmark_dat[1:100], function(prot) 
  signalHsmm_propep_decision(prot, aa_group = res[["aa_group"]], 
                             pipar = res[["pipar"]], 
                             tpmpar = res[["tpmpar"]], 
                             od = res[["od"]], 
                             overall_probs_log = res[["overall_probs_log"]], 
                             params = res[["params"]]))

decisions_comp <- data.frame(rbind(cbind(sapply(decisions_aa2, function(x) x$propep_probability), type="pos"), 
                                   cbind(sapply(negative_decision, function(x) x$propep_probability), type="neg")))

decisions_comp$V1 <- as.numeric(as.character(decisions_comp$V1))

ggplot(decisions_comp, aes(x=V1, fill=type)) +
  geom_density(alpha=.3)

library(ROCR)

pred <- prediction(decisions_comp$V1,decisions_comp$type)

roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)

auc <- performance(pred, measure = "auc")
auc <- unlist(slot(auc, "y.values"))
auc
