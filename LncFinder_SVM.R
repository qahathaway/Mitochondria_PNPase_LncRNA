library(LncFinder)
library(seqinr)

Sequences<- read.fasta(file = "path_to_your_files/Secondary_Structure_Machine_Learning.txt")

SS.seq <- run_RNAfold(Sequences, RNAfold.path = "path_to_your_files/RNAfold", parallel.cores = -1)

##DNA##
my_features <- extract_features(SS.seq, label = "NonCoding",SS.features = TRUE, format = "SS",
                                 frequencies.file = "human",parallel.cores = -1)

EIIP_res <- compute_EIIP(SS.seq, label = "NonCoding", spectrum.percent = 0.25,quantile.probs = seq(0, 1, 0.25))

FickettScore <- compute_FickettScore(Sequences, label = "NonCoding", on.ORF = TRUE,auto.full = TRUE, parallel.cores = -1)
gContent <- compute_GC(Sequences, label = "NonCoding",on.ORF = TRUE,auto.full = TRUE, parallel.cores = -1)
kmer_res <- compute_kmer(Sequences, k = 1:5, step = 1, freq = TRUE,improved.mode = TRUE, on.ORF = TRUE, auto.full = TRUE)


#Write CSV to combine features#
write.csv(t(kmer_res), file = "path_to_your_files/kmer.csv")
write.csv(FickettScore, file = "path_to_your_files/FickettScore.csv")
write.csv(gContent, file = "path_to_your_files/GCcontent.csv")
write.csv(my_features, file = "path_to_your_files/Features.csv")
write.csv(EIIP_res, file = "path_to_your_files/EIIP_res.csv")
write.csv(SS.seq, file = "path_to_your_files/SS_Seq.csv")



###Support Vector Machines (SVM)##
Features_All <- read.csv("path_to_your_files.csv", row.names=1)
Features_All[,2] <- as.numeric(as.character(Features_All[,2]))
Features_All[,1] <- as.factor(as.character(Features_All[,1]))

Feature_cv_tune <- svm_tune(Features_All, folds.num = 10, seed = 100, gamma.range = (2^seq(-5, 0, 1)), cost.range = c(1, 4, 8, 16, 24, 32), parallel.cores = -1)

Feature_cv_res <- svm_cv(Features_All, folds.num = 10, seed = 100, parallel.core = -1, return.model = TRUE)
