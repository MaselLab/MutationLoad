library(ggplot2)

#d = 0.5
setwd("C:/Users/walid/OneDrive/Bureau/Work/mutation_load/datafor_absoluteexponentialmub_0.000001chromosomes_23popsize_5000mud_0.000870d0_0.5seed_24/")

sb0.001 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.001.txt", header =  TRUE, sep = "\t")
sb0.005 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.005.txt", header =  TRUE, sep = "\t")
sb0.010 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.010.txt", header =  TRUE, sep = "\t")
sb0.015 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.015.txt", header =  TRUE, sep = "\t")
sb0.020 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.020.txt", header =  TRUE, sep = "\t")

sblength <- c(nrow(sb0.001), nrow(sb0.005), nrow(sb0.010), nrow(sb0.015), nrow(sb0.020))

sbgroups <- rep(c("0.001", "0.005", "0.010", "0.015", "0.020"), sblength)

totaldatas <- rbind(sb0.001, sb0.005, sb0.010, sb0.015, sb0.020)
totaldatas <- cbind(totaldatas, sbgroups)

pdeathrate0.5 <- ggplot(data=totaldatas, aes(x=Time, y=Mean_death_rate, group=sbgroups)) +
  geom_line(aes(color=sbgroups))

ggplot(data=totaldatas, aes(x=Time, y=Mean_birth_rate, group=sbgroups)) +
  geom_line(aes(color=sbgroups))

pdf("popsize_time.pdf")
ppopd0.5
dev.off()

pdf("meandeathrate_time.pdf")
pdeathrate0.5
dev.off()


#d = 0.2
setwd("../")

setwd("dataforabsoluteexponentialmub0.000001chromosomes23popsize5000mud0.000870d00.2seed24/")

sb0.001 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.001.txt", header =  TRUE, sep = "\t")
sb0.005 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.005.txt", header =  TRUE, sep = "\t")
sb0.010 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.010.txt", header =  TRUE, sep = "\t")
sb0.015 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.015.txt", header =  TRUE, sep = "\t")
sb0.020 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.020.txt", header =  TRUE, sep = "\t")

sblength <- c(nrow(sb0.001), nrow(sb0.005), nrow(sb0.010), nrow(sb0.015), nrow(sb0.020))

sbgroups <- rep(c("0.001", "0.005", "0.010", "0.015", "0.020"), sblength)

totaldatas <- rbind(sb0.001, sb0.005, sb0.010, sb0.015, sb0.020)
totaldatas <- cbind(totaldatas, sbgroups)

ppopd0.2 <- ggplot(data=totaldatas, aes(x=Time, y=Pop.size, group=sbgroups)) +
  geom_line(aes(color=sbgroups))

pdeathrate0.2 <- ggplot(data=totaldatas, aes(x=Time, y=Mean_death_rate, group=sbgroups)) +
  geom_line(aes(color=sbgroups))

ggplot(data=totaldatas, aes(x=Time, y=Mean_birth_rate, group=sbgroups)) +
  geom_line(aes(color=sbgroups))

pdf("popsize_time.pdf")
ppopd0.2
dev.off()

pdf("meandeathrate_time.pdf")
pdeathrate0.2
dev.off()

setwd("../")


#########epistasis##########

#d = 0.5
setwd("C:/Users/walid/OneDrive/Bureau/Work/mutation_load/datafor_absoluteexponentialmub_0.000001chromosomes_23popsize_5000mud_0.000870d0_0.5seed_24/epistasis")

sb0.001 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.001.txt", header =  TRUE, sep = "\t")
sb0.005 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.005.txt", header =  TRUE, sep = "\t")
sb0.010 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.010.txt", header =  TRUE, sep = "\t")
sb0.015 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.015.txt", header =  TRUE, sep = "\t")
sb0.020 <- read.csv("rawdataformaxTime8000initialPopSize5000mutrate0.000870chromsize50chromnum23benmutrate0.000001Sb0.020.txt", header =  TRUE, sep = "\t")

sblength <- c(nrow(sb0.001), nrow(sb0.005), nrow(sb0.010), nrow(sb0.015), nrow(sb0.020))

sbgroups <- rep(c("0.001", "0.005", "0.010", "0.015", "0.020"), sblength)

totaldatas <- rbind(sb0.001, sb0.005, sb0.010, sb0.015, sb0.020)
totaldatas <- cbind(totaldatas, sbgroups)

tiff("popsize_time_epistasis.tiff", units="in", width=5, height=5, res=300)
p <- ggplot(data=sb0.001[100:1000,], aes(x=Time, y=Pop.size))
p <- p+geom_line()
p <- p+labs(x="Time", y="Population size")+theme_Publication()
p
dev.off()
