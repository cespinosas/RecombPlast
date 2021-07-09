equis <- read.csv(file="Graphs/accAperfromB_5vs50.txt", sep="\t", head=FALSE, colClasses="numeric")
wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative="greater")
a <- "Wilcoxon test, 5pc vs 50pc, for access (pert) to GAP A through plasticity in networks with default B:"
capture.output(a, file="Graphs/stats2.txt", append=TRUE)
capture.output(wtest, file="Graphs/stats2.txt", append=TRUE)
equis <- read.csv(file="Graphs/accAredfromB_5vs50.txt", sep="\t", head=FALSE, colClasses="numeric")
wtest <- wilcox.test(equis$V1, equis$V2, paired=TRUE, alternative="greater")
a <- "Wilcoxon test, 5pc vs 50pc, for access (red) to GAP A through plasticity in networks with default B:"
capture.output(a, file="Graphs/stats2.txt", append=TRUE)
capture.output(wtest, file="Graphs/stats2.txt", append=TRUE)
