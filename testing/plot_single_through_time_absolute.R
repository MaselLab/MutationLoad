library(ggplot2)

maxTime = 1000
initialPopSize = 1000
mutrate = "2.000000"
chromsize = 10
chromnum = 23
benmutrate = "0.002000"
Sb = "0.010"
bendist = "exponential"
seed = 57
fitness = "absolute"
d_0 = "0.5"

eldir = paste0("datafor_", fitness, "_", bendist, "mub_", benmutrate, "chromosomes_", chromnum, "popsize_", initialPopSize, "mud_", mutrate, "d0_", d_0, "seed_", seed, "/")
setwd(paste0("../", eldir))
filesb1 = paste0("rawdataformaxTime", maxTime, "initialPopSize", initialPopSize, "mutrate", mutrate, "chromsize", chromsize, "chromnum", chromnum, "benmutrate", benmutrate, "Sb", Sb, ".txt")
datasb1 = read.csv(filesb1, header =  TRUE, sep = "\t")


ggplot(data=datasb1, aes(x=Time, y=Pop_size)) +
  geom_line() +
  labs(title = "no modules")

elementsperlb = 5
eldir = paste0("datafor_", fitness, "_",  "elementsperlb_", elementsperlb, bendist, "mub_", benmutrate, "chromosomes_", chromnum, "popsize_", initialPopSize, "mud_", mutrate, "d0_", d_0, "seed_", seed, "/")
setwd(paste0("../", eldir))

filesb1 = paste0("rawdataformaxTime", maxTime, "initialPopSize", initialPopSize, "mutrate", mutrate, "chromsize", chromsize, "chromnum", chromnum, "benmutrate", benmutrate, "Sb", Sb, ".txt")
datasb1 = read.csv(filesb1, header =  TRUE, sep = "\t")


ggplot(data=datasb1, aes(x=Time, y=Pop_size)) +
  geom_line() +
  labs(title = "5 modules")


elementsperlb = 15
eldir = paste0("datafor_", fitness, "_",  "elementsperlb_", elementsperlb, bendist, "mub_", benmutrate, "chromosomes_", chromnum, "popsize_", initialPopSize, "mud_", mutrate, "d0_", d_0, "seed_", seed, "/")
setwd(paste0("../", eldir))

filesb1 = paste0("rawdataformaxTime", maxTime, "initialPopSize", initialPopSize, "mutrate", mutrate, "chromsize", chromsize, "chromnum", chromnum, "benmutrate", benmutrate, "Sb", Sb, ".txt")
datasb1 = read.csv(filesb1, header =  TRUE, sep = "\t")


ggplot(data=datasb1, aes(x=Time, y=Pop_size)) +
  geom_line() +
  labs(title = "15 modules")