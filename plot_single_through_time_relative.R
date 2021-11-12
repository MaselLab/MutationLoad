library(ggplot2)

setwd("dataforrelativeexponentialmub0.000001chromosomes23popsize20000mud0.000870seed24/")

sb0.0025 <- read.csv("rawdataforNxtimesteps5000popsize20000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.0025.txt", header =  TRUE, sep = ",")

ggplot(data=sb0.0025, aes(x=Nxtimesteps, y=log(Sum.of.wis))) +
  geom_line()+
  ylim(0, 12.5)

ggplot(data=sb0.0025, aes(x=Nxtimesteps, y=Variance.in.log.fitness)) +
  geom_line(color = "blue")
