#This script plots the scatterpoint graphic accesibility to B through plasticity vs mutation from the data of access.m.

library(ggplot2)

seed = commandArgs(trailingOnly=TRUE)
x <- read.csv(paste0("S", seed[1], "/access_ci_mu.txt"), header = FALSE, sep = " ")
t1 = cor.test(x[,1], x[,2], alternative = "greater", method = "spearman")
letr <- "Congruence between plasticity and mutation"
capture.output(letr, file=paste0("S", seed[1], "/ci_mu_CE.txt"))
capture.output(t1, file=paste0("S", seed[1], "/ci_mu_CE.txt"), append=TRUE)

pdf(paste0("S", seed[1], "/ci_mu_CE.pdf"), width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)

ggplot(data = x, aes(x = x[,1], y = x[,2])) +
    geom_point(alpha = 0.2) + labs(x=expression('Access to '*italic(B)*' through plasticity'), y=expression('Access to '*italic(B)*' through mutation')) + theme_classic() + theme(text=element_text(size=14))
dev.off()
