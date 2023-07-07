#######set preferred theme########
theme_Publication <- function(base_size=14, base_family="Times New Roman") {
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
library(extrafont)
loadfonts(device = "win")
library(gridExtra)
library(ggpubr)
library(ggridges)
library(plyr)
library(corrplot)
theme_set(theme_Publication())
##################################################
#d=0.5
#pop size
#full simulation
png("changeinpopsize.png", units="in", width=8, height=5, res=300)
p <- ggplot(pop_size_table, aes(x=Time, y=Pop_size))
p <- p+geom_line()+theme_Publication()
p <- p+labs(x="Time", y="Population size")+theme_Publication()
p
dev.off()

#mean death rate or fitness
tiff("popsize_mean_0.5.tiff", units="in", width=5, height=5, res=300)
p1 <- ggplot(data=totaldatas, aes(x=Time, y=Mean_death_rate, group=sbgroups))
p1 <- p1+geom_line(aes(color=sbgroups))
p1 <- p1+labs(x="Time", y="Average death rate/max\nbirth rate")+theme_Publication()
p1
dev.off()

#variance death rate or fitness
tiff("popsize_variance_0.5.tiff", units="in", width=5, height=5, res=300)
p2 <- ggplot(data=totaldatas, aes(x=Time, y=Variance_death_rate, group=sbgroups))
p2 <- p2+geom_line(aes(color=sbgroups))
p2 <- p2+labs(x="Time", y="Variance in death rate/max\nbirth rate")+theme_Publication()
p2
dev.off()

#full simulation, both pop size and fitness
p_arranged <- ggarrange(p1, p2,labels = c("B","C"), ncol = 2, nrow = 1, legend = "none")
tiff("deathrate_0.5.tiff", units="in", width=10, height=7, res=300)
ggarrange(p, p_arranged,labels = c("A","",""),ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
dev.off()

#d=0.2
tiff("popsize_popsize_0.2.tiff", units="in", width=5, height=5, res=300)
p <- ggplot(data=totaldatas, aes(x=Time, y=Pop_size, group=sbgroups))
p <- p+geom_line(aes(color=sbgroups))
p <- p+labs(x="Time", y="Population size")+theme_Publication()
p
dev.off()

#mean death rate or fitness
tiff("popsize_mean_0.2.tiff", units="in", width=5, height=5, res=300)
p1 <- ggplot(data=totaldatas, aes(x=Time, y=Mean_death_rate, group=sbgroups))
p1 <- p1+geom_line(aes(color=sbgroups))
p1 <- p1+labs(x="Time", y="Average death rate/max\nbirth rate")+theme_Publication()
p1
dev.off()

#variance death rate or fitness
tiff("popsize_variance_0.2.tiff", units="in", width=5, height=5, res=300)
p2 <- ggplot(data=totaldatas, aes(x=Time, y=Variance_death_rate, group=sbgroups))
p2 <- p2+geom_line(aes(color=sbgroups))
p2 <- p2+labs(x="Time", y="Variance in death rate/max\nbirth rate")+theme_Publication()
p2
dev.off()

#full simulation, both pop size and fitness
p_arranged <- ggarrange(p1, p2,labels = c("B","C"), ncol = 2, nrow = 1, legend = "none")
tiff("deathrate_0.2.tiff", units="in", width=10, height=7, res=300)
ggarrange(p, p_arranged,labels = c("A","",""),ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
dev.off()
###################################expansive plots###########################
tiff("popsize_time_0.5.tiff", units="in", width=8, height=5, res=300)
p <- ggplot(data=totaldatas, aes(x=Time, y=Pop.size, group=sbgroups))
p <- p+geom_line(aes(color=sbgroups))
p <- p+geom_rect(aes(xmin= -50, xmax = 2000, ymin = 0, ymax = 10500),
                 fill = "transparent", color = "red", size =.5)
p <- p+geom_rect(aes(xmin= 2050, xmax = 7500, ymin = 0, ymax = 10500),
                 fill = "transparent", color = "red", size =.5)
p <- p+geom_rect(aes(xmin= 7550, xmax = 10000, ymin = 0, ymax = 10500),
                 fill = "transparent", color = "red", size =.5)
p <- p+annotate("text", x = ((-50+2000)/2), y=10000, label = "A")
p <- p+annotate("text", x = ((2050+7500)/2), y=10000, label = "B")
p <- p+annotate("text", x = ((7550+10000)/2), y=10000, label = "C")
p <- p+labs(x="Time", y="Population size")+theme_Publication()
p
dev.off()

#first stage: demographic expansion then decline
p_1 <- ggplot(data=totaldatas[which(totaldatas$Time <= 2000),],
              aes(x=Time, y=Pop.size, group=sbgroups))
p_1 <- p_1+geom_line(aes(color=sbgroups))
p_1 <- p_1+labs(x="Time", y="Population size")+theme_Publication()
p_1
#second stage: evolution
p_2 <- ggplot(data=totaldatas[which(totaldatas$Time > 2000 & totaldatas$Time <= 7500),],
              aes(x=Time, y=Pop.size, group=sbgroups))
p_2 <- p_2+geom_line(aes(color=sbgroups))
p_2 <- p_2+labs(x="Time", y="Population size")+theme_Publication()
p_2
#third stage: collapse or stability
p_3 <- ggplot(data=totaldatas[which(totaldatas$Time > 7501),],
              aes(x=Time, y=Pop.size, group=sbgroups))
p_3 <- p_3+geom_line(aes(color=sbgroups))
p_3 <- p_3+labs(x="Time", y="Population size")+theme_Publication()
p_3

p_arranged <- ggarrange(p_1, p_2,p_3,labels = c("A","B","C"),ncol = 3, nrow = 1, legend = "none")
#all plots together
tiff("popsize_0.5_stages.tiff", units="in", width=10, height=6, res=300)
ggarrange(p, p_arranged,ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
dev.off()

#plotting CPU usage

dat_cpu <- data.frame("Final N" = rep(c(16030, 20326, 20206, 24503, 28646),2),
                         "N of gen" = rep(c(500000, 500000, 700000, 900000, 1200000),2),
                         "CPU time" = rep(c(9.40, 12.26, 18.52, 30.03, 46.08),2),
                         "Comp metric" = c(18.09, 29.90, 32.08, 49.06, 66.00, 
                                           11.93, 17.07, 18.29, 26.12, 35.20),
                         "server"= as.factor(c(rep("Fusion", 5), rep("HPC", 5))))

summary(lm(Comp.metric~Final.N, data = subset(dat_cpu, subset = dat_cpu$server=="HPC")))

pdf("computation_metric.pdf")
p <- ggplot(data=dat_cpu, aes(x=Final.N, y=Comp.metric, group=server))
p <- p+ geom_smooth(aes(color=server), method = "lm")
p <- p+ labs(x="Final population size", y="(CPU time (m)/ gen) * Final N", color="Server",
             title = "CPU time scaled by final n \nof gen and final N")
p <- p+theme_Publication() +scale_x_continuous(limits = c(15000, 30000))
p
dev.off()

######
dat_cpuscale <- read.table("~/Work/MutationLoad/MutationLoad/dat_cpu.txt")
colnames(dat_cpuscale) <- c("generations", "N", "scale")

pdf("computation_scale_tskit_off.pdf")
p <- ggplot(data=dat_cpuscale, aes(x=log(as.numeric(N)*as.numeric(generations)), y=scale))
p <- p+ geom_point()
p <- p+ geom_smooth()
p <- p+ labs(x="log(final gen * final N)", y="(CPU time(m)/final gen)/Final N", color="Server",
             title = "Order gen*N")
p <- p+theme_Publication()
p
dev.off()

dat_cpuscale <- read.table("~/Work/MutationLoad/MutationLoad/dat_cpu.txt")
dat_cpuscale[,4] <- rep("OFF", nrow(dat_cpuscale))
colnames(dat_cpuscale) <- c("generations", "N", "scale", "tskit_status")
dat_cpuscale_tskit <- read.table("~/Work/MutationLoad/MutationLoad/dat_cpu_tskit.txt", sep = "\t")
dat_cpuscale_tskit[,4] <- rep("ON", nrow(dat_cpuscale_tskit))
colnames(dat_cpuscale_tskit) <- c("generations", "N", "scale", "tskit_status")

dat_cpu_all <- rbind(dat_cpuscale, dat_cpuscale_tskit)

pdf("computation_scale_tskit_all.pdf")
p <- ggplot(data=dat_cpu_all, aes(x=log(as.numeric(N)*as.numeric(generations)), y=scale, color = tskit_status))
p <- p+ geom_point(aes(color = tskit_status))
p <- p+ geom_smooth()
p <- p+ labs(x="log(final gen * final N)", y="(CPU time(m)/final gen)/Final N", color="tskit_status",
             title = "Order gen*N")
p <- p+theme_Publication()+scale_y_continuous(limits = c(4.0e-08, 3.0e-07))
p
dev.off()


N = 10^4
u = 10^-5
s= seq(from=0, to=0.5, by = 0.0001)
p0=1/N
p = (1 - exp(-4*s*p0*N))/(1 - exp(-4*N*s))

df = data.frame(s, p)

p1 <- ggplot(data= df, aes(x=s, y=p))
p1 <- p1+ geom_line()
p1 <- p1+ labs(x="selection coefficient", y="Fixation probability")
p1 <- p1+theme_Publication()
p1

######################################################title: "Nearly neutral theory"################################
library(ggplot2)
library(reshape2)
# The palette with gray:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Set parameters ----
N = c(1000, 5000, 10000)
s = seq(from = -0.001, to = 0, by = 0.00001)
x = expand.grid(N = N, s = s)
Ns = x$N*x$s
x <- cbind(x, Ns)
Pfix = ifelse(x$Ns != 0, (1-exp(-2*x$Ns/x$N))/(1- exp(-4*x$Ns)), 1/(2*x$N))
Pfix_neutral = 1/(2*x$N)
Pfix_over_Pfixneutral = Pfix/Pfix_neutral
x <- cbind(x, Pfix, Pfix_neutral, Pfix_over_Pfixneutral)

datatoplot <- subset(x, select = c(s, Pfix_over_Pfixneutral, Pfix, N))
tiff("prob_fix.tiff", units="in", width=5, height=5, res=300)
p_Pfix_over_Pfixneutral <- ggplot(datatoplot, aes(x = s, y = Pfix_over_Pfixneutral, color = as.factor(N))) +
  geom_line() +
  labs(x = "s", y = "Probability of Fixation", color = "Population size") +
  scale_colour_manual(values=cbbPalette) +
  theme(text = element_text(size = 20))+theme_Publication()    
p_Pfix_over_Pfixneutral
dev.off()

##########################################################title: flux calculation - analytic########################################
setwd("~/Work/MutationLoad/MutationLoad/Results/")
flux_ud2 <- read.csv("data.csv", header =  FALSE , sep = ",")
colnames(flux_ud2) <- c("delflux","benflux")
flux_ud2$N <- seq(1000,3000,100) 
flux_ud2$delta <- flux_ud2$delflux - flux_ud2$benflux

tiff("flux.tiff", units="in", width=5, height=5, res=300)
flux_ud2_plot <- ggplot(flux_ud2) +
  geom_line(aes(x=N, y=delflux, col="delf"),linewidth = 1,alpha=2) + 
  geom_line(aes(x=N, y=benflux, col="benf"), linewidth = 1,alpha=2)+
  geom_vline(aes(xintercept=1833),linetype="dashed")+
  labs(x = "Population Size", y = "Fitness flux") +
  theme(text = element_text(size = 20))+
  scale_color_manual(values=c("steelblue", "yellow3"), 
                     name="Flux",
                     labels=c("Beneficial", "Deleterious"))+theme_Publication()    
flux_ud2_plot
dev.off()
#change in sensitivity ratio by Ud
flux_ud_2 <- read.csv("data_ud_2.txt", header =  FALSE , sep = ",")
colnames(flux_ud_2) <- c("popsize","delflux","benflux")
flux_ud <- data.frame("Ud" = c(0.5, 1, 1.5, 2, 5, 10),
                         "Ncrit" = c(962,1329,1604,1833,2801,3857),
                         "dfdg" = c(0.836883,0.847712,0.850872,0.852749, 0.855644, 0.855788))

tiff("flux_Ud.tiff", units="in", width=5, height=5, res=300)
flux_ud2_plot <- ggplot(flux_ud, aes(x = Ud)) +
  geom_line(aes(y=dfdg), linewidth = 1) + 
  geom_line(aes(y=Ncrit*0.0001), linewidth = 2, col="red") + 
  scale_y_continuous(
    name = expression(frac("ben. flux sensitivity","del. flux sensitivity") ~ "|" ~ N[crit]),
    sec.axis = sec_axis(~.*10000, name=expression(N[crit]))) + labs(x = expression(U[d])) +
  theme(text = element_text(size = 20))+theme_Publication() 
flux_ud2_plot
dev.off()

tiff("flux_sens_Ud_2.tiff", units="in", width=5, height=5, res=300)
flux_ud2_plot <- ggplot(flux_ud_2, aes(x = popsize)) +
  geom_line(aes(y=delflux, colour = "Deletrious"), linewidth = 1) + 
  geom_vline(aes(xintercept=1833),linetype="dashed")+
  geom_line(aes(y=benflux, colour = "Beneficial"), linewidth = 1) + 
  labs(x="Population size", y="Flux sensitivity")+
  scale_color_manual(values=c("steelblue", "yellow3"), 
                     name="Flux",
                     labels=c("Beneficial", "Deleterious"))+theme_Publication()+
  theme(text = element_text(size = 16))
flux_ud2_plot
dev.off()
#change in sensitivity by Ud
flux_ud <- data.frame("Ud"= seq(0.1, 10.5, 0.5),
                         "ratio"=c(0.755262, 0.84125, 0.849588, 0.852461,
                                   0.853836, 0.854612, 0.855093, 0.855413,
                                   0.855635, 0.855794, 0.855911, 0.856, 0.856067,
                                   0.856119, 0.85616, 0.856192, 0.856217,
                                   0.856236, 0.856251, 0.856263, 0.856272))

tiff("flux_sens_Ud.tiff", units="in", width=6, height=6, res=300)
flux_ud_plot <- ggplot(flux_ud, aes(x = Ud, y=ratio)) +
  geom_line(size = 1) + 
  lims(y = c(0.7, 1))+
  labs(x=expression(U[d]), y=expression(frac("ben. flux sensitivity","del. flux sensitivity") ~ "|" ~ N[crit]))+
  theme_Publication()+
  theme(text = element_text(size = 16))
flux_ud_plot
dev.off()
#change in sensitivity ratio by Ub
flux_ub <- data.frame("Ub" = c(0.0002, 0.001, 0.002, 0.01, 0.02),
                      "Ncrit" = c(5312,2526,1833,866,623),
                      "dfdg" = c(0.856194,0.854896,0.852749,0.830479, 0.801654))

tiff("flux_sens_Ub.tiff", units="in", width=5, height=5, res=300)
flux_ub2_plot <- ggplot(flux_ub, aes(x = Ub)) +
  geom_line(aes(y=dfdg), linewidth = 1) + 
  labs(x = expression(U[b]), y=expression(frac("ben. flux sensitivity","del. flux sensitivity") ~ "|" ~ N[crit])) +
  theme(text = element_text(size = 20))+theme_Publication() + lim(y)+scale_x_continuous(breaks = seq(0, 0.02, 0.002)) 
flux_ub2_plot
dev.off()

#change in sensitivity by Ub/Ud
flux_ub2 <- rbind.data.frame(read.csv("data_ub_ud_0_0001.txt", header =  FALSE , sep = ","),
                             read.csv("data_ub_ud_0_0005.txt", header =  FALSE , sep = ","),
                             read.csv("data_ub_ud_0_001.txt", header =  FALSE , sep = ","),
                             read.csv("data_ub_ud_0_005.txt", header =  FALSE , sep = ","),
                             read.csv("data_ub_ud_0_01.txt", header =  FALSE , sep = ","),
                             read.csv("data_ub_ud_0_05.txt", header =  FALSE , sep = ","))
colnames(flux_ub2) <- c("popsize","delflux","benflux")
flux_ub2$UbUd <- factor(c(rep("0.0001", 11),rep("0.0005",11),rep("0.001",11),rep("0.005",11),rep("0.01",11),rep("0.05",11)),
                        levels = c("0.0001","0.0005","0.001","0.005","0.01","0.05"))
flux_ub2[1:11,]$popsize <- flux_ub2[1:11,]$popsize/5312
flux_ub2[12:22,]$popsize <- flux_ub2[12:22,]$popsize/2526
flux_ub2[23:33,]$popsize <- flux_ub2[23:33,]$popsize/1833
flux_ub2[34:44,]$popsize <- flux_ub2[34:44,]$popsize/866
flux_ub2[45:55,]$popsize <- flux_ub2[45:55,]$popsize/623
flux_ub2[56:66,]$popsize <- flux_ub2[56:66,]$popsize/278

tiff("flux_sens_Ub_Ud.tiff", units="in", width=8, height=8, res=300)
flux_ud2_plot <- ggplot(flux_ub2, aes(x = popsize)) +
  geom_line(aes(y=delflux, group = UbUd, linetype=UbUd), linewidth = 1) + geom_vline(aes(xintercept=rep(1,66)),linetype="dashed")+
  geom_line(aes(y=benflux, group = UbUd, linetype=UbUd), linewidth = 1) + 
  labs(x=expression(frac(N,N[crit])), y="Flux sensitivity")+
  theme_Publication() + theme(text = element_text(size = 16))
flux_ud2_plot
dev.off()

#####################################################title: flux calculation - simulation#############################################
setwd("~/Work/MutationLoad/MutationLoad/Results/HPC")

flux_sim_sb_0.002 <- data.frame("N" = c(900,1000,1100,1200,1300,1400),
                      "delf" = c(0.0001769691, 0.0001618280, 0.0001460149, 0.0001320459,
                                 0.0001217534,0.0001152294),
                      "benf" = c(0.0001023671,0.0001187022,0.0001241141,0.0001358045,
                                 0.0001461675,0.0001539184))
flux_sim_sb_0.004 <- data.frame("N" = c(300,400,500,600,700,800),
                              "delf" = c(0.0005811931,0.0004116249, 0.0003364481, 0.0002824124, 0.0002290653,
                                         0.0002071196),
                              "benf" = c(0.0001437562,0.0002028060,0.0002281810,0.0002825984,0.0003252260,
                                         0.0003465994))
flux_sim_sb_0.008 <- data.frame("N" = c(200,300,400,500,600,700),
                              "delf" = c(0.0008761119, 0.0005951607, 0.0004374721, 0.0003393053,
                                         0.0002842606,0.0002498480),
                              "benf" = c(0.0003435709,0.0005450489,0.0007150669,0.0009169896,
                                         0.0011363158,0.0012798791))

flux_sim_0.016 <- data.frame("N" = c(100,150,200,250,300,350),
                             "delf" = c(0.0017776497, 0.0012567256, 0.0009043413, 0.0007322552,
                                        0.0006115324,0.0005318247),
                             "benf" = c(0.0005900066,0.0011147069,0.0012911635,0.0018375238,
                                        0.0022001508,0.0026153327))

plot(d2$N, d2$expected2)
lines(flux_sim_0.001$N, flux_sim_0.001$benf, col="red")
#####smoothing####
fPow1 <- function(x, a, b) {a * x^b}
fPow2<- function(x, a, b) {a * x + b}
est1 <- coef(nls(delf ~ fPow1(N, a, b),
                 start=c(a=.00001, b=.00001), data=flux_sim_sb_0.004))
nlfit1 <- nls(delf ~ fPow1(N, a, b),
              start=est1, data=flux_sim_sb_0.004)

est2 <- coef(nls(benf ~ fPow2(N, a, b),
                 start=c(a=.00001, b=.00001), data=flux_sim_sb_0.004))
nlfit2 <- nls(benf ~ fPow2(N, a, b),
              start=est2, data=flux_sim_sb_0.004)

est3 <- coef(nls(delflux ~ fPow1(N, a, b),
                 start=c(a=.00001, b=.00001), data=flux_ud2))
nlfit3 <- nls(delflux ~ fPow1(N, a, b),
              start=est3, data=flux_ud2)

est4 <- coef(nls(benflux ~ fPow2(N, a, b),
                 start=c(a=.00001, b=.00001), data=flux_ud2))
nlfit4 <- nls(benflux ~ fPow2(N, a, b),
              start=est4, data=flux_ud2)

d2 <- data.frame(N = seq(300, 800, by = 1))
expected1 <- predict(nlfit1, d2)
expected2 <- predict(nlfit2, d2)
d2 <- data.frame(d2, expected1, expected2)
d1 <- data.frame(N = seq(1000, 6000, by = 1))
expected3 <- predict(nlfit3, d1)
expected4 <- predict(nlfit4, d1)
d1 <- data.frame(d1,expected3, expected4)

#get Ncrit
f1 <- approxfun(d2$expected1 - d2$expected2, d2$N, rule=2)
f1(0)

f2 <- approxfun(d1$expected3 - d1$expected4, d1$N, rule=2)
f2(0)
#plot
tiff("flux_sens_Ud.tiff", units="in", width=8, height=8, res=300)
flux_sim_plot <- ggplot(d2, aes(x=N)) +
  geom_line(aes(y=expected1, color = "delf"), linewidth = 1,alpha=2) + 
  geom_line(aes(y=expected2, color = "benf"), linewidth = 1,alpha=2) +
  geom_line(aes(y=expected3, color = "delf"), linewidth = 1,alpha=0.3) + 
  geom_line(aes(y=expected4, color = "benf"), linewidth = 1,alpha=0.3) +
  geom_vline(xintercept = 1833, linetype="dotted", alpha=0.4)+
  geom_vline(xintercept = 4554, linetype="longdash")+
  labs(x = "Population Size", y = "Fitness flux") +
  theme(text = element_text(size = 20))+
  scale_color_manual(values=c("steelblue", "yellow3"), 
                     name="Flux",
                     labels=c("Beneficial", "Deleterious"))+theme_Publication()
flux_sim_plot
dev.off()

s=0.00948704
Ud=2
Udw=0.00173913
Rw=0.0008695652
#Unlinked
Ne_N=exp(-8 * Ud * s)
#linked+unlined
Ne_N2=exp(-8 * (Ud-Udw)*s) * exp(-Udw/(2*s+Rw))

d2$Ne <- d2$N * Ne_N
d2$Ne <- d2$N * Ne_N2

tiff("flux_sens_Ud_Ne.tiff", units="in", width=8, height=8, res=300)
flux_sim_plot <- ggplot() +
  geom_line(data=d2,aes(x = Ne, y=expected1, color = "delf"), linewidth = 1,alpha=2) + 
  geom_line(data=d2,aes(x = Ne, y=expected2, color = "benf"), linewidth = 1,alpha=2) +
  geom_line(data=d1,aes(x = N, y=expected3, color = "delf"), linewidth = 1,alpha=0.3) + 
  geom_line(data=d1,aes(x = N, y=expected4, color = "benf"), linewidth = 1,alpha=0.3) +
  geom_vline(xintercept = 1833, linetype="dotted", alpha=0.4)+
  geom_vline(xintercept = 3913, linetype="longdash")+
  labs(x = "Effective population Size", y = "Fitness flux") +
  theme(text = element_text(size = 20))+
  scale_color_manual(values=c("steelblue", "yellow3"), 
                     name="Flux",
                     labels=c("Beneficial", "Deleterious"))+theme_Publication()
flux_sim_plot
dev.off()

#derivative calculation at Ncrit

dg_link <- function(N,epsilon) {return((-d2[which(d2$N == N+epsilon),2] + d2[which(d2$N == N-epsilon),2])/2*epsilon)}
df_link <- function(N,epsilon) {return((d2[which(d2$N == N+epsilon),3] - d2[which(d2$N == N-epsilon),3])/2*epsilon)}

dg_nolink <- function(N,epsilon) {return((-d2[which(d2$N == N+epsilon),4] + d2[which(d2$N == N-epsilon),4])/2*epsilon)}
df_nolink <- function(N,epsilon) {return((d2[which(d2$N == N+epsilon),5] - d2[which(d2$N == N-epsilon),5])/2*epsilon)}

#ratio at ncrit no linkage
ratio_nolink <- df_nolink(1833, 10)/dg_nolink(1833, 10)
ratio_nolink
#ratio at ncrit with linkage
ratio_link <- df_link(600, 10)/dg_link(600, 10)
ratio_link
#difference in ratio
ratio_link - ratio_nolink
######simple linkage effect#######
df <- data.frame("ratio" = c(0.85362, 0.9772939), 
                 "Ub" = c(0.002, 0.002))
df$model <- factor(c("No Linkage", "Linkage"),
                levels = c("No Linkage", "Linkage"))

tiff("sens_ratio_simple.tiff", units="in", width=5, height=5, res=300)
flux_sim_plot <- ggplot(df, aes(x = model, y = ratio, group = Ub))+
  geom_point(size=2)+geom_line(size=1)+geom_hline(yintercept = 1, linetype="dashed")+
  annotate("text", x=2, y=0.8, label= expression(U[d]~"= 2"~U[b]~"= 0.002"~s[b]~"= 0.001"))+
  labs(x= "Model", y = expression("Sensitivity ratio at "~N[crit]))+
  lims(y=c(0.8, 1))+theme_Publication()+theme(text = element_text(size = 16))
flux_sim_plot
dev.off()

#############
df <- data.frame("nolinkage" = c(0.845398, 0.857697, 0.878588, 0.916017),
                 "linkage" = c(0.886792, 0.8938948, 1.038609, 1.182704))

df$sb <- factor(c("0.002", "0.004","0.008", "0.016"),
                        levels = c("0.002", "0.004","0.008", "0.016"))

tiff("sens_ratio.tiff", units="in", width=5, height=5, res=300)
flux_sim_plot <- ggplot(df, aes(x = nolinkage, y = linkage, color=sb))+
  geom_point()+geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  annotate("text", x=0.95, y=0.8, label= expression(U[d]~"= 2"~U[b]~"= 0.008"))+
  annotate("text", x=0.97, y=0.96, label= "1:1 ratio",angle=20)+
  labs(x= expression("Sensitivity ratio at "~N[crit]~" with no linkage"),
       y = expression("Sensitivity ratio at "~N[crit]~" with linkage"), color = expression(s[b]))+
  lims(x=c(0.8, 0.98), y=c(0.8, 1.25))+theme_Publication()+theme(text = element_text(size = 16))
flux_sim_plot
dev.off()
#####################################
library(plotly)

M1 <- as.matrix(read_xls("~/Work/MutationLoad/MutationLoad/Results/HPC/mfile.xls"))
colnames(M1) <- NULL

fig <- plot_ly(
  x = seq(0.0001, 0.01, 0.0005), 
  y = seq(0.001, 0.01, 0.0005), 
  z = M1, 
  type = "contour" 
)

fig
