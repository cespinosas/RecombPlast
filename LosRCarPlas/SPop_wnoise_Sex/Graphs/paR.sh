a <- "\n\n\n***A founder:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B appears first (two-tailed test), according to a binomial distribution:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "Number of succesful simulations: 871"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
 bt <- binom.test(442, 871, 0.5,alternative="two.sided")
capture.output(bt, file="Graphs/stats.txt", append=TRUE)
a <- "\n\n**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B evolves (5pc) (two-tailed test), according to a binomial distribution:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "Number of succesful simulations: 846"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
 bt <- binom.test(433, 846, 0.5,alternative="two.sided")
capture.output(bt, file="Graphs/stats.txt", append=TRUE)
a <- "\n\n**Probability of obtaining a greater deviation from expected (0.5) frequency of times that B evolves (50pc) (two-tailed test), according to a binomial distribution:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "Number of succesful simulations: 835"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
 bt <- binom.test(426, 835, 0.5,alternative="two.sided")
capture.output(bt, file="Graphs/stats.txt", append=TRUE)
a <- "\n\n\n***AB founder:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "**Probability of obtaining a greater fraction of times that B appears first (two-tailed test), according to a binomial distribution:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "Number of succesful simulations: 1000"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
 bt <- binom.test(1000, 1000, 0.5,alternative="two.sided")
capture.output(bt, file="Graphs/stats.txt", append=TRUE)
a <- "\n\n**Probability of obtaining a greater fraction of times that B evolves (two-tailed test), according to a binomial distribution:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "Number of succesful simulations: 1000"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
 bt <- binom.test(1000, 1000, 0.5,alternative="two.sided")
capture.output(bt, file="Graphs/stats.txt", append=TRUE)
a <- "\n\n**Probability of obtaining a greater fraction of times that B evolves (two-tailed test), according to a binomial distribution:"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
a <- "Number of succesful simulations: 1000"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
 bt <- binom.test(1000, 1000, 0.5,alternative="two.sided")
capture.output(bt, file="Graphs/stats.txt", append=TRUE)
library(ggplot2)
a <- read.csv(file = "Graphs/pacajaprimer.txt", head=FALSE, sep ="\t")
a$V2 <- factor(a$V2)
pdf(file="Graphs/tfirstGA.pdf", width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)
ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour="black", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20")+labs(x=expression('Founder'), y = expression('First genetically assimilated optimum')) + stat_summary(geom = "point", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))
dev.off()
a <- read.csv(file = "Graphs/pacaja5.txt", head=FALSE, sep ="\t")
a$V2 <- factor(a$V2)
pdf(file="Graphs/tuntil5.pdf", width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)
ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour="black", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20")+labs(x=expression('Founder'), y = expression('Time until prevalence (0.05) of assimilated optimum')) + stat_summary(geom = "point", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))
dev.off()
a <- read.csv(file = "Graphs/pacaja50.txt", head=FALSE, sep ="\t")
a$V2 <- factor(a$V2)
pdf(file="Graphs/tuntil50.pdf", width=7.5,height=5.32)
par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)
ggplot(a, aes(x=V2, y=V1)) + geom_point(size=1, shape=16, alpha=0.2, colour="black", position=position_jitter(width=0.3, height=0.001)) + geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20")+labs(x=expression('Founder'), y = expression('Time until prevalence (0.50) of assimilated optimum')) + stat_summary(geom = "point", fun.y=mean, size=3) + theme_classic() + theme(text = element_text(size=14))
dev.off()
equis <- read.csv(file="Graphs/Atp.txt", sep="\t", head=FALSE, colClasses="numeric")
ye <- read.csv(file="Graphs/ABtp.txt", sep="\t", head=FALSE, colClasses="numeric")
wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative="greater")
a <- "Mann-Whitney for time of appearance (one-tailed; time when founder A > time when founder AB):"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
capture.output(wtest, file="Graphs/stats.txt", append=TRUE)
equis <- read.csv(file="Graphs/At5.txt", sep="\t", head=FALSE, colClasses="numeric")
ye <- read.csv(file="Graphs/ABt5.txt", sep="\t", head=FALSE, colClasses="numeric")
wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative="greater")
a <- "Mann-Whitney for time until prevalence (0.05) (one-tailed; time when founder A > time when founder AB):"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
capture.output(wtest, file="Graphs/stats.txt", append=TRUE)
equis <- read.csv(file="Graphs/At50.txt", sep="\t", head=FALSE, colClasses="numeric")
ye <- read.csv(file="Graphs/ABt50.txt", sep="\t", head=FALSE, colClasses="numeric")
wtest <- wilcox.test(equis$V1, ye$V1, paired=FALSE, alternative="greater")
a <- "Mann-Whitney for time until prevalence (0.50) (one-tailed; time when founder A > time when founder AB):"
capture.output(a, file="Graphs/stats.txt", append=TRUE)
capture.output(wtest, file="Graphs/stats.txt", append=TRUE)
