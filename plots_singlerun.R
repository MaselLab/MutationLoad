#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}
#######set preferred theme########
theme_Publication <- function(base_size=14, base_family="Cambria") {
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

library(ggplot2)
theme_set(theme_Publication())
#######################main script###########################
print(args[1])
setwd(args[1])

print(args[2])

#main_dir <- "/home/wmawass/MutationLoad"
#dir_list <- list.dirs(main_dir,recursive = FALSE)
#setwd(dir_list[grep("datafor", dir_list)[1]])
#file_list <- list.files(dir_list[grep("datafor", dir_list)[1]], recursive = F)

dat <- read.csv(list.files(getwd(), pattern="rawdatafor"), header =  TRUE, sep = "\t")

tiff(paste("popsize_by_time",as.character(args[2]),".tiff", sep=""), units="in", width=8, height=5, res=300)
p <- ggplot(data=dat[seq(1, nrow(dat), 100),], aes(x=Time, y=Pop_size))
p <- p+geom_line()
p <- p+labs(x="Time", y="Population size")+theme_Publication()
p
dev.off()

tiff(paste("popsize_by_time_log",as.character(args[2]),".tiff", sep=""), units="in", width=8, height=5, res=300)
p <- ggplot(data=dat[seq(1, nrow(dat), 100),], aes(x=log(Time), y=Pop_size))
p <- p+geom_line()
p <- p+labs(x="Log(Time)", y="log(Population size)")+theme_Publication()
p
dev.off()

tiff(paste("popsize_by_time_smoothline",as.character(args[2]),".tiff", sep=""), units="in", width=8, height=5, res=300)
p <- ggplot(data=dat, aes(x=Time, y=Pop_size))
p <- p+geom_smooth()
p <- p+labs(x="Time", y="Population size")+theme_Publication()
p
dev.off()

tiff(paste("popsize_by_time_smoothline_log",as.character(args[2]),".tiff", sep=""), units="in", width=8, height=5, res=300)
p <- ggplot(data=dat, aes(x=log(Time), y=Pop_size))
p <- p+geom_smooth()
p <- p+labs(x="Log(Time)", y="Population size")+theme_Publication()
p
dev.off()

library(R.utils)

dat_subset <- dat[dat$Time > max(dat$Time)-1000,]
scaled_dat_subset <- dat_subset

scaled_dat_subset[,"Time"] <- scale(dat_subset$Time)[,1]
scaled_dat_subset[,"Pop_size"] <- scale(dat_subset$Pop_size)[,1]

fit <- lm(Pop_size~log(Time), data = dat_subset)
sum.fit <- summary(fit)

print(sum.fit$coefficients)

tstats <- (0-sum.fit$coefficients[2,1])/sum.fit$coefficients[2,2]
# Calculates two tailed probability
pval_log<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
print(pval_log)

fit <- lm(Pop_size~Time, data = scaled_dat_subset)
sum.fit <- summary(fit)

print(sum.fit$coefficients)

tstats <- (0-sum.fit$coefficients[2,1])/sum.fit$coefficients[2,2]
# Calculates two tailed probability
pval_scaled <- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
print(pval_scaled)
