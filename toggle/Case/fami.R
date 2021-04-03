#This script takes the data from fami.m and plots the scatterpoint graphics that compare the "Parental access to B through plasticity" with the "Offspring's access to B through plasticity", and "Parental access to B through plasticity" with "Offspring's access to B through mutation". Also, it plots a boxplot and the Mann-Whitney test for the parents with 0 offspring with native B and parents with 1 or more offspring with native B.

library(ggplot2)
library(ggExtra)

seed = commandArgs(trailingOnly=TRUE)
x <- read.csv(paste0("S", seed[1], "/pci_pmu_hci_hmu_has.txt"), header = FALSE, sep = " ")




t1 = cor.test(x[,1], x[,3], alternative = "greater", method = "spearman")
letr <- "Association between parental and filial acccess to B through plasticity"
capture.output(letr, file=paste0("S", seed[1], "/parci_offci_CE.txt"))
capture.output(t1, file=paste0("S", seed[1], "/parci_offci_CE.txt"), append=TRUE)

pdf(paste0("S", seed[1], "/parci_offci_CE.pdf"), width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)

p1 <- ggplot(data = x, aes(x = x[,1], y = x[,3])) +
    geom_point(alpha = 0.15) + labs(x=expression('Parental access to '*italic(B)*' through plasticity'), y=expression('Offspring\'s access to '*italic(B)*' through plasticity')) + theme_classic() + theme(text=element_text(size=14))
    ggExtra::ggMarginal(p1, type = "histogram", size=8)
dev.off()






t1 = cor.test(x[,1], x[,4], alternative = "greater", method = "spearman")
letr <- "Association between parental access through plasticity and filial acccess to B through mutation"
capture.output(letr, file=paste0("S", seed[1], "/parci_offmu_CE.txt"))
capture.output(t1, file=paste0("S", seed[1], "/parci_offmu_CE.txt"), append=TRUE)

pdf(paste0("S", seed[1], "/parci_offmu_CE.pdf"), width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)
p1 <- ggplot(data = x, aes(x = x[,1], y = x[,4])) +
    geom_point(alpha = 0.15) + labs(x=expression('Parental access to '*italic(B)*' through plasticity'), y=expression('Offspring\'s access to '*italic(B)*' through mutation')) + theme_classic() + theme(text=element_text(size=14))
    ggExtra::ggMarginal(p1, type = "histogram", size=8)
dev.off()


a <- read.csv(paste0("S", seed[1], "/has_pci.txt"), header = FALSE, sep = "\t")

a$V1 <- factor(a$V1)

pdf(paste0("S", seed[1], "/parci_offas_CE.pdf"), width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)
ggplot(a, aes(x=V1, y=V2)) + geom_point(size=1, shape=16, alpha=0.2, colour="black", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20")+ scale_x_discrete(labels=c(0,expression("">0))) + labs(x=expression('Offspring with native '*italic(B)*' '), y = expression('Parental access to '*italic(B)*' through plasticity')) + stat_summary(geom = "point", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))
dev.off()

s <- read.csv(paste0("S", seed[1], "/has_en_pci0.txt"), header = FALSE, sep = "\t")
c <- read.csv(paste0("S", seed[1], "/has_en_pci1.txt"), header = FALSE, sep = "\t")
wtest <- wilcox.test(s$V1,c$V1, paired=FALSE, alternative="less")
letr <- "Mann-Whitney: Offspring with native B with minimum parental Bn1 less than with higher values of parental Bn1"
capture.output(letr, file=paste0("S", seed[1], "/parci_offas_CE.txt"), append = TRUE)
capture.output(wtest, file=paste0("S", seed[1], "/parci_offas_CE.txt"), append = TRUE)

