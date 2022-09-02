#######set preferred theme########
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

library(grid)
library(ggthemes)
library(scales)
library(ggplot2)
theme_set(theme_Publication())

#Warning, sd and typeofrun are set without checking the simulation data
#Check that the code corresponds to this variables
sb = "0.010"
typeofrun = "single"
slope = "0"
type_fitness = "absolute_"
#elperl = "1"
d0 = "0.2"
r = "0.8"
mud= "2.0"
sdmin = "0.01"
#deldistri = "point_"
bendistri = "exponential_"
mub = "0.0020"
chromosomes = "23"
chromosomesize = "50"
seed = "24"
#sdtosbratio = "0.1"
tskitstatus = "OFF_"

maxTime = "700000"
initialPopsize="20000"
K="25000"
directory = paste0("datafor_", type_fitness, "tskitstatus_", tskitstatus, "bendist_", bendistri, "mub_", mub, "_chromnum_", chromosomes, "_popsize_", initialPopsize, "_mud_", mud, "_d0_", d0, "_seed_", seed, "/")
file1 <- read.csv(paste0(directory, "rawdatafor_maxTime_", maxTime, "_initialPopSize_", initialPopsize, "_mutrate_", mud, "_chromsize_", chromosomesize, "_chromnum_",chromosomes, "_benmutrate_", mub, "_Sb_",sb, ".txt"), header =  TRUE, sep = "\t")
genatend = tail(file1$Time, n = 1)
finalN = tail(file1$Pop_size, n = 1)
printtotable1 <- c(type_fitness, maxTime , as.numeric(K)* (1-as.numeric(d0)), K, mud, chromosomes, chromosomesize , as.numeric(mub)/as.numeric(mud), sb, bendistri, typeofrun, slope, seed, d0, r, sdmin, genatend, finalN)

maxTime = "900000"
initialPopsize="24000"
K="30000"
directory = paste0("datafor_", type_fitness, "tskitstatus_", tskitstatus, "bendist_", bendistri, "mub_", mub, "_chromnum_", chromosomes, "_popsize_", initialPopsize, "_mud_", mud, "_d0_", d0, "_seed_", seed, "/")
file2 <- read.csv(paste0(directory, "rawdatafor_maxTime_", maxTime, "_initialPopSize_", initialPopsize, "_mutrate_", mud, "_chromsize_", chromosomesize, "_chromnum_",chromosomes, "_benmutrate_", mub, "_Sb_",sb, ".txt"), header =  TRUE, sep = "\t")
genatend = tail(file2$Time, n = 1)
finalN = tail(file2$Pop_size, n = 1)
printtotable2 <- c(type_fitness, maxTime , as.numeric(K)* (1-as.numeric(d0)), K, mud, chromosomes, chromosomesize , as.numeric(mub)/as.numeric(mud), sb, bendistri, typeofrun, slope, seed, d0, r, sdmin, genatend, finalN)

maxTime = "1200000"
initialPopsize="28000"
K="35000"
directory = paste0("datafor_", type_fitness, "tskitstatus_", tskitstatus, "bendist_", bendistri, "mub_", mub, "_chromnum_", chromosomes, "_popsize_", initialPopsize, "_mud_", mud, "_d0_", d0, "_seed_", seed, "/")
file3 <- read.csv(paste0(directory, "rawdatafor_maxTime_", maxTime, "_initialPopSize_", initialPopsize, "_mutrate_", mud, "_chromsize_", chromosomesize, "_chromnum_",chromosomes, "_benmutrate_", mub, "_Sb_",sb, ".txt"), header =  TRUE, sep = "\t")
genatend = tail(file3$Time, n = 1)
finalN = tail(file3$Pop_size, n = 1)
printtotable3 <- c(type_fitness, maxTime , as.numeric(K)* (1-as.numeric(d0)), K, mud, chromosomes, chromosomesize , as.numeric(mub)/as.numeric(mud), sb, bendistri, typeofrun, slope, seed, d0, r, sdmin, genatend, finalN)


labelslegend <- c("Census N (without mutation)", "20000", "24000", "28000")
pdflabel <- paste0("_tskitstatus_", tskitstatus , "bendistri_", bendistri)

pdf(paste0("popsize_time", pdflabel, ".pdf"))
p <- ggplot()
p <- p+ geom_line(data = file1[seq(1, nrow(file1), 100),], aes(x = Time/20000, y = Pop_size/20000, color = labelslegend[2]))
p <- p+ geom_line(data = file2[seq(1, nrow(file2), 100),], aes(x = Time/24000, y = Pop_size/24000, color = labelslegend[3]))
p <- p+ geom_line(data = file3[seq(1, nrow(file3), 100),], aes(x = Time/28000, y = Pop_size/28000, color = labelslegend[4]))
p <- p+ labs(x="Time", y="Population size", title = "Population Size v Time, sampled data")
p <- p+ scale_color_discrete(name = labelslegend[1])
p
dev.off()

pdf(paste0("popsize_logtime", pdflabel, ".pdf"))
p <- ggplot()
p <- p+ geom_line(data = file1[seq(1, nrow(file1), 100),], aes(x = log(Time)/20000, y = Pop_size/20000, color = labelslegend[2]))
p <- p+ geom_line(data = file2[seq(1, nrow(file2), 100),], aes(x = log(Time)/24000, y = Pop_size/24000, color = labelslegend[3]))
p <- p+ geom_line(data = file3[seq(1, nrow(file3), 100),], aes(x = log(Time)/28000, y = Pop_size/28000, color = labelslegend[4]))
p <- p+ labs(x="Log(Time)", y="Population size", title = "Population Size v log(Time), sampled data")
p <- p+ scale_color_discrete(name = labelslegend[1])
p
dev.off()

pdf(paste0("popsize_time_smoothed", pdflabel, ".pdf"))
p <- ggplot()
p <- p+ geom_smooth(data = file1, aes(x = Time/20000, y = Pop_size/20000, color = labelslegend[2]))
p <- p+ geom_smooth(data = file2, aes(x = Time/24000, y = Pop_size/24000, color = labelslegend[3]))
p <- p+ geom_smooth(data = file3, aes(x = Time/28000, y = Pop_size/28000, color = labelslegend[4]))
p <- p+ labs(x="Time", y="Population size", title = "Population Size v Time, smoothed line")
p <- p+ scale_color_discrete(name = labelslegend[1])
p
dev.off()

pdf(paste0("popsize_logtime_smoothed", pdflabel, ".pdf"))
p <- ggplot()
p <- p+ geom_smooth(data = file1, aes(x = log(Time)/20000, y = Pop_size/20000, color = labelslegend[2]))
p <- p+ geom_smooth(data = file2, aes(x = log(Time)/24000, y = Pop_size/24000, color = labelslegend[3]))
p <- p+ geom_smooth(data = file3, aes(x = log(Time)/28000, y = Pop_size/28000, color = labelslegend[4]))
p <- p+ labs(x="Log(Time)", y="Population size", title = "Population Size v log(Time), smoothed line")
p <- p+ scale_color_discrete(name = labelslegend[1])
p
dev.off()


datafortable <- data.frame(printtotable1, printtotable2, printtotable3)
datafortable <- t(datafortable)

tablefile = paste0("datafortable", pdflabel, ".txt")
write.table(datafortable, tablefile, append = FALSE, row.names = FALSE, col.names = FALSE)
