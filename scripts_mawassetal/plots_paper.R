#######set preferred theme########
library(extrafont)
font_import()
loadfonts(device = "win")

theme_Publication <- function(base_size=24, base_family="Times New Roman") {
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
            panel.grid.major = element_blank(),
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
  discrete_scale("fill","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c",
                                       "#662506","#a6cee3","#fb9a99","#984ea3",
                                       "#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#386cb0","#fdb462","#7fc97f",
                                       "#ef3b2c","#662506","#a6cee3","#fb9a99",
                                       "#984ea3","#ffff33")), ...)
  
}

library(ggplot2)
library(scales)
library(gridExtra)
library(ggpubr)
library(ggridges)
library(plyr)
library(corrplot)
library(latex2exp)
library(ggallin)
theme_set(theme_Publication())
setwd("analytical_data/")
######################################Figure 2+Figure 6: flux+derivative example with Ud/Ub=1000##########################################
#analytical data
dat_del <- read.csv("deltable.csv", header = F)
colnames(dat_del) <- c("N", "flux")
dat_ben <- read.csv("bentable.csv", header = F)
colnames(dat_ben) <- c("N", "flux")
dat_del_der <- read.csv("delderivetable.csv", header = F)
colnames(dat_del_der) <- c("N", "der")
dat_ben_der <- read.csv("benderivetable.csv", header = F)
colnames(dat_ben_der) <- c("N", "der")

#sim data
df_sim <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                     "vb" = c(0.0000077289, 0.0000155021, 0.0000213351, 0.0000280516,
                              0.0000341099, 0.0000421500, 0.0000500220),
                     "vd" = c(0.0001508777, 0.0000695224, 0.0000462792, 0.0000347247,
                              0.0000279091, 0.0000237818, 0.0000206981))

#####curve fitting####
# Perform linear regression for vb
linear_model_vb <- lm(vb ~ N, data = df_sim)

# Perform exponential regression for vd
#using an exponential model
#exponential_model_vd <- nls(vd ~ a * exp(b * N), data = df_sim, start = list(a = 1, b = 0))
#using power law function
fPow1 <- function(x, a, b) {a * x^b}
est1 <- coef(nls(vd ~ fPow1(N, a, b),
                 start=c(a=.00001, b=.000001), data=df_sim,control = nls.control(maxiter = 1000)))
power_law_model <- nls(vd ~ fPow1(N, a, b),
                       start=est1, data=df_sim)

# Predict values for the whole range of N
N_range <- seq(1000, 10000, by = 1)  # Adjust the step size as needed
predictions <- data.frame("N" = N_range,
                          "vb_pred" = predict(linear_model_vb, newdata = data.frame(N = N_range)),
                          "vd_pred" = predict(power_law_model, newdata = data.frame(N = N_range)))

#predicting for analytical fluxes
linear_model_vb_ana <- lm(flux ~ N, data = dat_ben)

fPow1 <- function(x, a, b) {a * x^b}
est1 <- coef(nls(flux ~ fPow1(N, a, b),
                 start=c(a=.00001, b=.00001), data=dat_del, control = nls.control(maxiter = 1000)))
power_law_model_ana <- nls(flux ~ fPow1(N, a, b),
                           start=est1, data=dat_del)

#adding environmental deviation, delta=1.5e-05
dat_del$flux_offset <- dat_del$flux+0.000015

est1 <- coef(nls(flux_offset ~ fPow1(N, a, b),
                 start=c(a=.00001, b=.00001), data=dat_del,control = nls.control(maxiter = 1000)))

power_law_model_ana_offset <- nls(flux_offset ~ fPow1(N, a, b),
                                  start=est1, data=dat_del)

predictions_ana <- data.frame("N" = N_range,
                              "vb_pred" = predict(linear_model_vb_ana, newdata = data.frame(N = N_range)),
                              "vd_pred" = predict(power_law_model_ana, newdata = data.frame(N = N_range)),
                              "vd_pred_offset" = predict(power_law_model_ana_offset, newdata = data.frame(N = N_range)))

#get Ncrit
f1 <- approxfun(predictions$vb_pred - predictions$vd_pred, predictions$N, rule=1)
f1(0)

f2 <- approxfun(predictions_ana$vb_pred - predictions_ana$vd_pred, predictions_ana$N, rule=1)
f2(0)

#plotting fluxes with simulation data points and env 
#Figure 2A
tiff("fluxes_analtyical_withsim_withenv.tiff", units="in", width=6, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_del, aes(x=N, y=flux, col="Deleterious"), linewidth =1)
p <- p+geom_line(data = dat_del, aes(x=N, y=flux_offset, col="Deleterious+env. change"),  linewidth =1)
p <- p+geom_line(data = dat_ben, aes(x=N, y=flux, col="Beneficial"), linewidth =1)
p <- p+geom_vline(xintercept = c(1834.058, 2330.3),linewidth=c(1,1), col = c("#ef3b2c","darkorange"),alpha=c(.8, .8), linetype=c("dotted", "dotted"))
#p <- p+geom_point(data = df_sim, aes(x=N, y=vd, col="Deleterious"), size = 3)
#p <- p+geom_point(data = df_sim, aes(x=N, y=vb, col="Beneficial"), size = 3)
p <- p+xlim(c(500, 7000))
p <- p+scale_y_continuous(limits = c(0, 0.00016), 
                          breaks = c(0, 0.00005, 0.0001, 0.00015),
                          labels = c("0", expression("5 ×" ~ 10^-5), expression("1 ×" ~ 10^-4), expression("1.5 ×" ~ 10^-4)))
p <- p+labs(x="Census population size", y="Fitness change\nper generation")+theme_Publication()
p <- p+scale_colour_manual(name = '', 
                           values =c("Deleterious"="#ef3b2c","Deleterious+env. change"="darkorange","Beneficial"="#386cb0"),
                           labels = c('Beneficial','Deleterious+env. change','Deleterious'))
p <- p+guides(color = guide_legend(override.aes = list(size = 4)), 
              linetype=guide_legend(keywidth = 2, keyheight = 1))
p <- p+theme(axis.text.x = element_text(size=18),
             axis.text.y = element_text(size=21),
             axis.title.y = element_text(size=22, ,vjust = -2),
             legend.position = "none")
p
dev.off()

#Figure 6A
tiff("fluxes_analtyical_withsim_withenv_1.tiff", units="in", width=6, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_del, aes(x=N, y=flux, col="Deleterious"),alpha=.5, linewidth =1)
p <- p+geom_line(data = dat_del, aes(x=N, y=flux_offset, col="Deleterious+env. change"),alpha=.5,  linewidth =1)
p <- p+geom_line(data = dat_ben, aes(x=N, y=flux, col="Beneficial"),alpha=.5, linewidth =1)
p <- p+geom_vline(xintercept = c(1834.058, 2330.3),linewidth=c(1,1), col = c("#ef3b2c","darkorange"),alpha=c(.8, .8), linetype=c("dotted", "dotted"))
p <- p+geom_point(data = df_sim, aes(x=N, y=vd, col="Deleterious"), size = 3)
p <- p+geom_point(data = df_sim, aes(x=N, y=vb, col="Beneficial"), size = 3)
p <- p+xlim(c(500, 7000))
p <- p+scale_y_continuous(limits = c(0, 0.00016), 
                          breaks = c(0, 0.00005, 0.0001, 0.00015),
                          labels = c("0", expression("5 ×" ~ 10^-5), expression("1 ×" ~ 10^-4), expression("1.5 ×" ~ 10^-4)))
p <- p+labs(x="Census population size", y="Fitness change\nper generation")+theme_Publication()
p <- p+scale_colour_manual(name = '', 
                           values =c("Deleterious"="#ef3b2c","Deleterious+env. change"="darkorange","Beneficial"="#386cb0"),
                           labels = c('Beneficial','Deleterious+env. change','Deleterious'))
p <- p+guides(color = guide_legend(override.aes = list(size = 4)), 
              linetype=guide_legend(keywidth = 2, keyheight = 1))
p <- p+theme(axis.text.x = element_text(size=18),
             axis.text.y = element_text(size=21),
             axis.title.y = element_text(size=22, ,vjust = -2),
             legend.position = "none")
p
dev.off()

#Figure 2B
tiff("derivative_analtyical_offset.tiff", units="in", width=5.5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_del_der, aes(x=N, y=der, col="Deleterious"), linewidth =1)
p <- p+geom_line(data = dat_ben_der, aes(x=N, y=der, col="Beneficial"), linewidth =1)
p <- p+geom_vline(xintercept = c(1834.058, 2330.3),linewidth=c(1,1), col = c("#ef3b2c","darkorange"),alpha=c(.8, .8), linetype=c("dotted", "dotted"))
p <- p+scale_y_continuous(limits = c(0, 0.00000004),
                          breaks = c(0, 0.00000001, 0.00000002, 0.00000003, 0.00000004),
                          labels = c("0", expression("1 ×" ~ 10^-8), expression("2 ×" ~ 10^-8),expression("3 ×" ~ 10^-8),expression("4 ×" ~ 10^-8)))
p <- p+xlim(c(500, 7000))
p <- p+labs(x="Census population size", y="Fitness flux derivative")+theme_Publication()
p <- p+scale_colour_manual(name = ' ', 
                           values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))
p <- p+guides(color = guide_legend(override.aes = list(size = 5)))
p <- p+theme(axis.text.x = element_text(size=18),
             axis.text.y = element_text(size=21),
             axis.title.y = element_text(size=26),
             legend.position = "none")
p
dev.off()

#Figure 6C
s=0.00948704
Ud=2
L=100
chr=23
Udw=Ud/(L*chr/2)
Udw
Rw=2/(L*chr/2/23)
Rw
#Unlinked
Ne_N=exp(-8 * Ud * s)
Ne_N
#linked+unlined Joseph's equation
Ne_N2=exp(-8 * (Ud-Udw)*s) * exp(-Ud/(2*23))
Ne_N2

Ne_N_sim <- 1834.57/4402.57
Ne_N_sim
Ne_N2 - Ne_N_sim
(Ne_N2 - Ne_N_sim)/Ne_N2*100

predictions$Ne <- predictions$N * Ne_N2
df_sim$Ne <- df_sim$N * Ne_N2

predictions$Ne_ext <- predictions$N * Ne_N_sim
df_sim$Ne_ext <- df_sim$N * Ne_N_sim

#Figure 6C
tiff("fluxes_analtyical_withsim_withinterpolationwithNe_ext.tiff", units="in", width=6, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_del, aes(x=N, y=flux, linetype="With no LD"),col="#ef3b2c", alpha=0.5, linewidth =1)
p <- p+geom_line(data = dat_ben, aes(x=N, y=flux, linetype="With no LD"),col="#386cb0", alpha=0.5, linewidth =1)
p <- p+geom_vline(xintercept = c(1834.058), linewidth = c(1),alpha=c(0.8), linetype=c("dotted"))
p <- p+geom_point(data = df_sim, aes(x=Ne_ext, y=vd, col="Deleterious"), size = 3)
p <- p+geom_point(data = df_sim, aes(x=Ne_ext, y=vb, col="Beneficial"), size = 3)
p <- p+xlim(c(500, 4000))
p <- p+scale_y_continuous(limits = c(0, 0.00016), 
                          breaks = c(0, 0.00005, 0.0001, 0.00015),
                          labels = c("0", expression("5 ×" ~ 10^-5), expression("1 ×" ~ 10^-4), expression("1.5 ×" ~ 10^-4)))
p <- p+labs(x="Effective population size", y="Fitness flux")+theme_Publication()
p <- p+scale_colour_manual(name = ' ', 
                           values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))
p <- p+scale_linetype_manual(name = '',
                             values=c("With LD"=2,"With no LD"=1), labels = c('With LD','With no LD'))
p <- p+guides(linetype=guide_legend(keywidth = 2, keyheight = 1, override.aes=c(col="black")))
p <- p+theme(axis.text.x = element_text(size=18),
             axis.text.y = element_text(size=18),
             legend.position = "none")
p
dev.off()

#derivative cal at ncrit with census pop
dg_sim_n <- function(N, epsilon) {
  index_plus <- which(abs(predictions$N - (N + epsilon)) < 0.0001)
  index_minus <- which(abs(predictions$N - (N - epsilon)) < 0.0001)
  
  return ((predictions$vd_pred[index_plus] - predictions$vd_pred[index_minus]) / (2 * epsilon))
}

df_sim_n <- function(N, epsilon) {
  index_plus <- which(abs(predictions$N - (N + epsilon)) < 0.0001)
  index_minus <- which(abs(predictions$N - (N - epsilon)) < 0.0001)
  
  return ((predictions$vb_pred[index_plus] - predictions$vb_pred[index_minus]) / (2 * epsilon))
}

#derivative cal at ncrit with effective pop
dg_sim_ne <- function(N, epsilon) {
  index_plus <- which(abs(predictions$Ne_ext - (N + epsilon)) < 0.0001)
  index_minus <- which(abs(predictions$Ne_ext - (N - epsilon)) < 0.0001)
  
  return ((predictions$vd_pred[index_plus] - predictions$vd_pred[index_minus]) / (2 * epsilon))
}

df_sim_ne <- function(N, epsilon) {
  index_plus <- which(abs(predictions$Ne_ext - (N + epsilon)) < 0.0001)
  index_minus <- which(abs(predictions$Ne_ext - (N - epsilon)) < 0.0001)
  
  return ((predictions$vb_pred[index_plus] - predictions$vb_pred[index_minus]) / (2 * epsilon))
}

# Set epsilon to 1 for integer data points
epsilon <- 1

# Calculate the change in slope for vb with N using epsilon = 1
change_in_slope_vb <- numeric(length(N_range) - 2)

change_in_slope_vb <- sapply(N_range[2:(length(N_range) - 1)], function(N) df_sim_n(N, epsilon))

change_in_slope_vb_sim <- as.vector(sapply(df_sim$N[-1], function(N) df_sim_n(N, epsilon)))

# Calculate the change in slope for vd with N using epsilon = 1
change_in_slope_vd <- numeric(length(N_range) - 2)

change_in_slope_vd <- sapply(N_range[2:(length(N_range) - 1)], function(N) dg_sim_n(N, epsilon))

change_in_slope_vd_sim <- as.vector(sapply(df_sim$N[-1], function(N) dg_sim_n(N, epsilon)))

# Set epsilon to 0.4167043 for float data points
epsilon <- 0.4167043

# Calculate the change in slope for vb with N using epsilon = 0.4167043
change_in_slope_vb_ne <- numeric(length(N_range) - 2)

change_in_slope_vb_ne <- sapply(predictions$Ne_ext[2:(length(N_range) - 1)], function(N) df_sim_ne(N, epsilon))

change_in_slope_vb_ne_sim <- as.vector(sapply(df_sim$Ne_ext[-1], function(N) df_sim_ne(N, epsilon)))

# Calculate the change in slope for vd with N using epsilon = 0.4167043
change_in_slope_vd_ne <- numeric(length(N_range) - 2)

change_in_slope_vd_ne <- sapply(predictions$Ne_ext[2:(length(N_range) - 1)], function(N) dg_sim_ne(N, epsilon))

change_in_slope_vd_ne_sim <- as.vector(sapply(df_sim$Ne_ext[-1], function(N) dg_sim_ne(N, epsilon)))


#Figure 6B
#plotting fluxes
tiff("derivative_analtyicalwithsim.tiff", units="in", width=6, height=5, res=300)
p1 <- ggplot()
p1 <- p1+geom_line(data = dat_del_der, aes(x=N, y=der, col="Deleterious"),alpha=0.5, linewidth =1)
p1 <- p1+geom_line(data = dat_ben_der, aes(x=N, y=der, col="Beneficial"),alpha=0.5, linewidth =1)
p1 <- p1+geom_vline(xintercept = c(1834.058), linewidth = c(1),alpha=c(0.8), linetype=c("dotted"))
p1 <- p1+geom_point(aes(x=df_sim$N[-1], y=abs(change_in_slope_vd_sim), col="Deleterious"), size=3)
p1 <- p1+geom_point(aes(x=df_sim$N[-1], y=change_in_slope_vb_sim, col="Beneficial"), size=3)
p1 <- p1+scale_y_continuous(limits = c(0, 0.00000004),
                            breaks = c(0, 0.00000001, 0.00000002, 0.00000003, 0.00000004),
                            labels = c("0", expression("1 ×" ~ 10^-8), expression("2 ×" ~ 10^-8),expression("3 ×" ~ 10^-8),expression("4 ×" ~ 10^-8)))
p1 <- p1+xlim(c(500, 7000))
p1 <- p1+labs(x="Effective population size", y="Fitness flux derivative")+theme_Publication()
p1 <- p1+scale_colour_manual(name = ' ', 
                             values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))
p1 <- p1+guides(color = guide_legend(override.aes = list(size = 5)))
p1 <- p1+theme(axis.text.x = element_text(size=18),
               axis.text.y = element_text(size=18),
               legend.position = "none")
p1
dev.off()

#Figure 6D
#plotting fluxes
tiff("derivative_analtyicalwithsim_ext.tiff", units="in", width=6, height=5, res=300)
p1 <- ggplot()
p1 <- p1+geom_line(data = dat_del_der, aes(x=N, y=der, col="Deleterious"),alpha=0.5, linewidth =1)
p1 <- p1+geom_line(data = dat_ben_der, aes(x=N, y=der, col="Beneficial"),alpha=0.5, linewidth =1)
p1 <- p1+geom_vline(xintercept = c(1834.058), linewidth = c(1),alpha=c(0.8), linetype=c("dotted"))
p1 <- p1+geom_point(aes(x=df_sim$Ne_ext[-1], y=abs(change_in_slope_vd_ne_sim), col="Deleterious"), size=3)
p1 <- p1+geom_point(aes(x=df_sim$Ne_ext[-1], y=change_in_slope_vb_ne_sim, col="Beneficial"), size=3)
p1 <- p1+scale_y_continuous(limits = c(0, 0.00000004),
                            breaks = c(0, 0.00000001, 0.00000002, 0.00000003, 0.00000004),
                            labels = c("0", expression("1 ×" ~ 10^-8), expression("2 ×" ~ 10^-8),expression("3 ×" ~ 10^-8),expression("4 ×" ~ 10^-8)))
p1 <- p1+xlim(c(500, 4000))
p1 <- p1+labs(x="Effective population size", y="Fitness flux derivative")+theme_Publication()
p1 <- p1+scale_colour_manual(name = ' ', 
                             values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))
p1 <- p1+guides(color = guide_legend(override.aes = list(size = 5)))
p1 <- p1+theme(axis.text.x = element_text(size=18),
               axis.text.y = element_text(size=18),
               legend.position = "none")
p1
dev.off()

#####################################Figure 3: Env degradation#######################
dat_env_offset_kim_100 <- read.csv("offset_kim_raw_100.csv", header = F)
colnames(dat_env_offset_kim_100) <- c("offset", "del_flux","Ncrit", "bendelratio")
dat_env_offset_kim_100$UbUd <- rep(100,dim(dat_env_offset_kim_100)[1])

dat_env_offset_kim_1000 <- read.csv("offset_kim_raw_1000.csv", header = F)
colnames(dat_env_offset_kim_1000) <- c("offset", "del_flux","Ncrit", "bendelratio")
dat_env_offset_kim_1000$UbUd <- rep(1000,dim(dat_env_offset_kim_1000)[1])

dat_env_offset_kim <- rbind(dat_env_offset_kim_100,dat_env_offset_kim_1000)

dat_env_offset_kim$gen <- 0.1/dat_env_offset_kim$offset
dat_env_offset_kim$ratio <- dat_env_offset_kim$offset/abs((dat_env_offset_kim$del_flux))

custom_breaks <- c(1000.0000, 10000.0000, 100000.0000)  # Custom tick positions
custom_labels <- c(expression(10^3),expression(10^4),expression(10^5))  # Custom tick labels

#plot
tiff("env_offset.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_env_offset_kim, aes(x=0.1/offset, y=bendelratio, group=as.factor(UbUd), col=as.factor(UbUd)), linewidth =1)
p <- p+geom_point(aes(x=0.1/(0.000015),y=1.43395),col="darkorange",shape=42, size=10)
p <- p+scale_x_log10(breaks=custom_breaks,
                     labels=custom_labels)
p <- p+scale_y_log10()+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+geom_hline(yintercept = 1.0, col="red", linetype="dashed")
p <- p+geom_vline(xintercept = c(0.1/2.2e-05), col="#2298e6", linetype="dotdash", linewidth = 1.1, alpha=0.5)
p <- p+geom_vline(xintercept = c(0.1/4.5e-06), col="#9e9e9e", linetype="dotdash", linewidth = 1.1, alpha=0.5)
p <- p+scale_colour_manual(name = expression(U[d]/U[b]),
                           values = c(100, 1000),
                           labels = c("100", "1000"))
p <- p+labs(x="Generations needed for env.\nchange to reduce fitness by 10%",
            y=expression("Drought:Meltdown ratio"))+theme_Publication()
p <- p+guides(colour=guide_legend(keywidth = 1.5, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",axis.title.x=element_text(size = 19.5),
             legend.text = element_text(size=16),axis.text.y = element_text(size = 20))
p
dev.off()

#plot
tiff("env_offset_ratio.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_env_offset_kim, aes(x=0.1/offset, y=ratio, group=as.factor(UbUd), col=as.factor(UbUd)), linewidth =1)
p <- p+geom_point(aes(x=0.1/(0.000015),y=0.000015/0.0000221607),col="darkorange",shape=42, size=8)
p <- p+scale_x_log10(breaks=custom_breaks,
                     labels=custom_labels)
p <- p+geom_abline(intercept = 3, slope = -1, col="black", linetype="dashed")
p <- p+scale_y_log10(breaks=c(0.01,0.1,1, 10, 100),
                     labels=c("0.01","0.1","1","10", "100"))+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+geom_hline(yintercept = 1.0, col="red", linetype="dashed")
p <- p+geom_vline(xintercept = c(0.1/2.2e-05), col="#2298e6", linetype="dotdash", linewidth = 1.1, alpha=0.5)
p <- p+geom_vline(xintercept = c(0.1/4.5e-06), col="#9e9e9e", linetype="dotdash", linewidth = 1.1, alpha=0.5)
p <- p+scale_colour_manual(name = expression(U[d]/U[b]),
                           values = c(100, 1000),
                           labels = c("100", "1000"))
p <- p+labs(x="Generations needed for env.\nchange to reduce fitness by 10%",
            y=expression(atop("Env. change to del. flux ",paste("ratio at ", N[crit]))))+theme_Publication()
p <- p+guides(shape=guide_legend(keywidth = 1, keyheight = 1, ncol =1, nrow=3))
p <- p+theme(legend.position = "none", axis.title.x=element_text(size = 18),
             legend.text = element_text(size = 18),axis.text.y = element_text(size = 18),
             axis.title.y = element_text(size=20,vjust = -2), legend.background = element_blank())
p
dev.off()

#plot
tiff("env_offset_ncrit.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_env_offset_kim, aes(x=(0.1/offset), y=Ncrit, group=as.factor(UbUd), col=as.factor(UbUd)), linewidth =1)
p <- p+geom_point(aes(x=0.1/(0.000015),y=2330.34),col="darkorange",shape=42, size=10)
p <- p+geom_vline(xintercept = c(0.1/2.2e-05), col="#2298e6", linetype="dotdash", linewidth = 1.1, alpha=0.5)
p <- p+geom_vline(xintercept = c(0.1/4.5e-06), col="#9e9e9e", linetype="dotdash", linewidth = 1.1, alpha=0.5)
p <- p+scale_x_log10(breaks=custom_breaks,
                     labels=custom_labels)
p <- p+scale_y_log10(breaks=c(0, 1000, 10000),
                     labels=c("0","1000", "10000"))+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[d]/U[b]),
                           values = c(100, 1000),
                           labels = c("100", "1000"))
p <- p+geom_hline(yintercept = c(623.805, 1834.57), col="red", linetype="dashed")
p <- p+labs(x="Generations needed for env.\nchange to reduce fitness by 10%",
            y=expression(N[crit]))+theme_Publication()
p <- p+guides(shape=guide_legend(keywidth = 1, keyheight = 1, ncol =1, nrow=3))
p <- p+theme(legend.position = "none", axis.title.x=element_text(size = 19.5),
             legend.text = element_text(size = 18),axis.title.y = element_text(size=22,vjust = -2),)
p
dev.off()

###########################################Figure 4: ben:del ratio no env change####################################
#ratio sensitvity to Ud/Ub
dat_ubud_kim <- read.csv("udubtablekim_ratio.csv", header = F)
colnames(dat_ubud_kim) <- c("Ub_Ud", "del_flux","Ncrit", "bendelratio")

custom_breaks <- c(10, 100, 1000, 10000)  # Custom tick positions
custom_labels <- c("10", "100", "1000", "10000")  # Custom tick labels

#without offset
tiff("ubud_analtyical_sim.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ubud_kim, aes(x=Ub_Ud, y=bendelratio), linewidth =1)
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels,
                     name=expression(U[d]~"/"~U[b])) +annotation_logticks(sides = "b", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+ylim(0,1.05)+geom_point(aes(x=1000,y=0.853449),col = "#F05039", shape=42, size=10)
p <- p+labs(y="Drought:Meltdown ratio")+theme_Publication()
p <- p+theme(legend.position.inside = c(0.5, 1), legend.background = element_blank())
p
dev.off()

#ratio sensitvity to mean xb
dat_sb_kim_100 <- read.csv("sbtablekim_100.csv", header = F)
colnames(dat_sb_kim_100) <- c("sb", "del_flux","Ncrit", "bendelratio")
dat_sb_kim_1000 <- read.csv("sbtablekim_1000.csv", header = F)
colnames(dat_sb_kim_1000) <- c("sb", "del_flux","Ncrit", "bendelratio")

custom_breaks <- c(0.001, 0.005, 0.01)  # Custom tick positions
custom_labels <- c("0.001", "0.005", "0.01")  # Custom tick labels

tiff("sb_analtyical.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_sb_kim_100, aes(x=sb, y=bendelratio, linetype="mut. rate ratio = 100"), linewidth =1)
p <- p+geom_line(data = dat_sb_kim_1000, aes(x=sb, y=bendelratio, linetype="mut. rate ratio = 1000"), linewidth =1)
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels)+annotation_logticks(sides = "b", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_linetype_manual(name = '',
                             values=c("mut. rate ratio = 100"=2,"mut. rate ratio = 1000"=1),
                             labels = c("mut. rate ratio = 100", "mut. rate ratio = 1000"))
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+ylim(0,1.05)+geom_point(aes(x=0.001,y=0.853449),col = "#F05039",shape=42, size=10)
p <- p+labs(x=bquote(bar(x)[b]), y="Drought:Meltdown ratio")+theme_Publication()
p <- p+theme(legend.position = "none")
p
dev.off()

#ratio sensitvity to mean xd
dat_sd_kim_100 <- read.csv("sdtablekim_100.csv", header = F)
colnames(dat_sd_kim_100) <- c("sd", "del_flux","Ncrit", "bendelratio")
dat_sd_kim_1000 <- read.csv("sdtablekim_1000.csv", header = F)
colnames(dat_sd_kim_1000) <- c("sd", "del_flux","Ncrit", "bendelratio")

custom_breaks <- c(0.0001, 0.001, 0.01, 0.1)  # Custom tick positions
custom_labels <- c("0.0001", "0.001", "0.01", "0.1")  # Custom tick labels

tiff("sd_analtyical.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_sd_kim_100[-31,], aes(x=sd, y=bendelratio, linetype="mut. rate ratio = 100"), linewidth =1)
p <- p+geom_line(data = dat_sd_kim_1000[-31,], aes(x=sd, y=bendelratio, linetype="mut. rate ratio = 1000"), linewidth =1)
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels)+annotation_logticks(sides = "b", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_linetype_manual(name = '',
                             values=c("mut. rate ratio = 100"=2,"mut. rate ratio = 1000"=1),
                             labels = c("mut. rate ratio = 100", "mut. rate ratio = 1000"))
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+ylim(0,1.5)+geom_point(aes(x=0.0095,y=0.853449),col = "#F05039",shape=42, size=10)
p <- p+labs(x=bquote(-bar(x)[d]), y="Drought:Meltdown ratio")+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 5)))
p <- p+theme(legend.position = "none")
p
dev.off()

tiff("varsd_analtyical.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_sd_kim_100[-31,], aes(x=(sd^2)/0.169, y=bendelratio, linetype="mut. rate ratio = 100"), linewidth =1)
p <- p+geom_line(data = dat_sd_kim_1000[-31,], aes(x=(sd^2)/0.169, y=bendelratio, linetype="mut. rate ratio = 1000"), linewidth =1)
p <- p+scale_x_log10()+annotation_logticks(sides = "b", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_linetype_manual(name = '',
                             values=c("mut. rate ratio = 100"=2,"mut. rate ratio = 1000"=1),
                             labels = c("mut. rate ratio = 100", "mut. rate ratio = 1000"))
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+ylim(0,1.5)+geom_point(aes(x=0.0095,y=0.853449),col = "#F05039",shape=42, size=10)
p <- p+labs(x=bquote(Var(x)[d]), y="Drought:Meltdown ratio")+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 5)))
p <- p+theme(legend.position = "none")
p
dev.off()

###########################################Figure 5: sensitivity to Ub and Ud with env change####################################
#with offset Kim et al. ratio varying Ud fixed Ub
dat_ud_kim_ub0.002_highoffset <- read.csv("offset_udvaryingub0.002_kim.csv", header = F)
colnames(dat_ud_kim_ub0.002_highoffset) <- c("Ud", "del_flux","Ncrit", "bendelratio")

dat_ud_kim_ub0.02_highoffset <- read.csv("offset_udvaryingub0.02_kim.csv", header = F)
colnames(dat_ud_kim_ub0.02_highoffset) <- c("Ud", "del_flux","Ncrit", "bendelratio")

dat_ud_kim_ub0.002_lowoffset <- read.csv("lowoffset_udvaryingub0.002_kim.csv", header = F)
colnames(dat_ud_kim_ub0.002_lowoffset) <- c("Ud", "del_flux","Ncrit", "bendelratio")

dat_ud_kim_ub0.02_lowoffset <- read.csv("lowoffset_udvaryingub0.02_kim.csv", header = F)
colnames(dat_ud_kim_ub0.02_lowoffset) <- c("Ud", "del_flux","Ncrit", "bendelratio")

custom_breaks <- c(1, 3, 10)  # Custom tick positions
custom_labels <- c("1","3","10")  # Custom tick labels

#ben del ratio
tiff("ud_analtyical_offset.tiff", units="in", width=4, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ud_kim_ub0.002_highoffset,
                 aes(x=Ud, y=bendelratio, col="0.002", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_highoffset,
                 aes(x=Ud, y=bendelratio, col="0.02", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.002_lowoffset,
                 aes(x=Ud, y=bendelratio, col="0.002", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_lowoffset,
                 aes(x=Ud, y=bendelratio, col="0.02", linetype="Low offset"), linewidth =1)
p <- p+scale_x_log10(name=expression(U[d]))+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[b]), 
                           values =c("0.002" = "#9e9e9e","0.02"="#2298e6"), 
                           labels = c("0.002","0.02"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('High offset','Low offset'))
p <- p+scale_y_log10(limits=c(0.35,18.5),
                     breaks=custom_breaks,
                     labels=custom_labels)
p <- p+geom_point(aes(x=2,y=0.853449),col = "#F05039",shape=42, size=10)
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+theme_Publication()+labs(y="Drought:Meltdown ratio")
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =2, nrow=1))
p <- p+theme(legend.position = "none",legend.title=element_text(size=16), legend.background = element_blank(),
             legend.text = element_text(size=16), axis.text = element_text(size=20))
p
dev.off()

custom_breaks <- c(300, 1000, 3000, 10000)  # Custom tick positions
custom_labels <- c("300", "1000", "3000", "10000")  # Custom tick labels
#Ncrit
tiff("ud_analtyical_offset_ncrit.tiff", units="in", width=4, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ud_kim_ub0.002_highoffset,
                 aes(x=Ud, y=Ncrit, col="0.002", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_highoffset,
                 aes(x=Ud, y=Ncrit, col="0.02", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.002_lowoffset,
                 aes(x=Ud, y=Ncrit, col="0.002", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_lowoffset,
                 aes(x=Ud, y=Ncrit, col="0.02", linetype="Low offset"), linewidth =1)
p <- p+scale_x_log10(name=expression(U[d]))+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[b]), 
                           values =c("0.002" = "#9e9e9e","0.02"="#2298e6"), 
                           labels = c("0.002","0.02"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('High offset','Low offset'))
p <- p+scale_y_log10(limits=c(120, 21000),
                     breaks=custom_breaks,
                     labels=custom_labels)
p <- p+geom_point(aes(x=2,y=1834.57),col = "#F05039",shape=42, size=10)
p <- p+labs(y=expression(N[crit]))+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",legend.title=element_text(size=14), legend.background = element_blank(),
             legend.text = element_text(size=12), axis.text = element_text(size=20))
p
dev.off()

#
custom_breaks_y <- c(0.01, 0.1, 1)  # Custom tick positions
custom_labels_y <- c(expression(0.01),expression(0.1),expression(1))  # Custom tick labels
#offset del ratio
tiff("ud_analtyical_offset_delratio.tiff", units="in", width=4, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ud_kim_ub0.002_highoffset,
                 aes(x=Ud, y=1.5e-5/abs(del_flux), col="0.002", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_highoffset,
                 aes(x=Ud, y=1.5e-5/abs(del_flux), col="0.02", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.002_lowoffset,
                 aes(x=Ud, y=1e-6/abs(del_flux), col="0.002", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_lowoffset,
                 aes(x=Ud, y=1e-6/abs(del_flux), col="0.02", linetype="Low offset"), linewidth =1)
p <- p+scale_x_log10(name=expression(U[d]))+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[b]), 
                           values =c("0.002" = "#9e9e9e","0.02"="#2298e6"), 
                           labels = c("0.002","0.02"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('High offset','Low offset'))
p <- p+scale_y_log10(limits=c(0.0038,1),
                     breaks=custom_breaks_y,
                     labels=custom_labels_y)
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+labs(y=expression(atop("Env. change to del. flux ",paste("ratio at ", N[crit]))))+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",legend.title=element_text(size=14),
             legend.background = element_blank(), axis.title.y = element_text(size=20, vjust = -1),
             legend.text = element_text(size=12), axis.text = element_text(size=18))
p
dev.off()

#with offset Kim et al. ratio varying Ub fixed Ud
dat_ub_kim_ud2_highoffset <- read.csv("offset_ubvaryingud2_kim.csv", header = F)
colnames(dat_ub_kim_ud2_highoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

dat_ub_kim_ud1_highoffset <- read.csv("offset_ubvaryingud1_kim.csv", header = F)
colnames(dat_ub_kim_ud1_highoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

dat_ub_kim_ud2_lowoffset <- read.csv("lowoffset_ubvaryingud2_kim.csv", header = F)
colnames(dat_ub_kim_ud2_lowoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

dat_ub_kim_ud1_lowoffset <- read.csv("lowoffset_ubvaryingud1_kim.csv", header = F)
colnames(dat_ub_kim_ud1_lowoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

custom_breaks <- c(0.0001, 0.001, 0.01, 0.1)  # Custom tick positions
custom_labels <- c(expression(10^-4),expression(10^-3),expression(10^-2),
                   expression(10^-1))  # Custom tick labels

custom_breaks_y <- c(1, 3, 10)  # Custom tick positions
custom_labels_y <- c("1","3","10")  # Custom tick labels

#ben del ratio
tiff("ub_analtyical_offset.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ub_kim_ud2_highoffset,
                 aes(x=Ub, y=bendelratio, col="2", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_highoffset,
                 aes(x=Ub, y=bendelratio, col="1", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud2_lowoffset,
                 aes(x=Ub, y=bendelratio, col="2", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_lowoffset,
                 aes(x=Ub, y=bendelratio, col="1", linetype="Low offset"), linewidth =1)
p <- p+scale_x_continuous(trans=c("log10","reverse"),name=expression(U[b]),
                          breaks=custom_breaks,
                          labels=custom_labels)+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[d]), 
                           values =c("2" = "#9e9e9e","1"="#2298e6"), 
                           labels = c("1","2"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('Fast env. change','Slow env. change'))
p <- p+scale_y_log10(limits=c(0.35,18.5),
                     breaks=custom_breaks_y,
                     labels=custom_labels_y)
p <- p+geom_point(aes(x=0.002,y=0.853449),col = "#F05039",shape=42, size=10)
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth =.9)
p <- p+labs(y="Drought:Meltdown ratio")+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",legend.title=element_text(size=14), legend.background = element_blank(),
             legend.text = element_text(size=12), axis.text = element_text(size=20))
p
dev.off()

custom_breaks_y <- c(300, 1000, 3000, 10000)  # Custom tick positions
custom_labels_y <- c("300", "1000", "3000", "10000")  # Custom tick labels
#Ncrit
tiff("ub_analtyical_offset_ncrit.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ub_kim_ud2_highoffset,
                 aes(x=Ub, y=Ncrit, col="2", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_highoffset,
                 aes(x=Ub, y=Ncrit, col="1", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud2_lowoffset,
                 aes(x=Ub, y=Ncrit, col="2", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_lowoffset,
                 aes(x=Ub, y=Ncrit, col="1", linetype="Low offset"), linewidth =1)
p <- p+scale_x_continuous(trans=c("log10","reverse"),name=expression(U[b]),
                          breaks=custom_breaks,
                          labels=custom_labels)+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[d]), 
                           values =c("2" = "#9e9e9e","1"="#2298e6"), 
                           labels = c("1","2"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('High offset','Low offset'))
p <- p+scale_y_log10(limits=c(120, 21000),
                     breaks=custom_breaks_y,
                     labels=custom_labels_y)
p <- p+geom_point(aes(x=0.002,y=1834.57),col = "#F05039",shape=42, size=10)
p <- p+labs(y=expression(N[crit]))+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",legend.title=element_text(size=14), legend.background = element_blank(),
             legend.text = element_text(size=12), axis.text = element_text(size=20))
p
dev.off()

custom_breaks_y <- c(0.01, 0.1, 1)  # Custom tick positions
custom_labels_y <- c(expression(0.01),expression(0.1),expression(1))  # Custom tick labels
#offset del ratio
tiff("ub_analtyical_offset_delratio.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ub_kim_ud2_highoffset,
                 aes(x=Ub, y=1.5e-5/abs(del_flux), col="2", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_highoffset,
                 aes(x=Ub, y=1.5e-5/abs(del_flux), col="1", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud2_lowoffset,
                 aes(x=Ub, y=1e-6/abs(del_flux), col="2", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_lowoffset,
                 aes(x=Ub, y=1e-6/abs(del_flux), col="1", linetype="Low offset"), linewidth =1)
p <- p+scale_x_continuous(trans=c("log10","reverse"),name=expression(U[b]),
                          breaks=custom_breaks,
                          labels=custom_labels)+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[d]), 
                           values =c("2" = "#9e9e9e","1"="#2298e6"), 
                           labels = c("1","2"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('High offset','Low offset'))
p <- p+scale_y_log10(limits=c(0.0038, 1),
                     breaks=custom_breaks_y,
                     labels=custom_labels_y)
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth =.9)
p <- p+labs(y=expression(atop("Env. change to del. flux ",paste("ratio at ", N[crit]))))+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",legend.title=element_text(size=14),
             legend.background = element_blank(), axis.title.y = element_text(size=20,vjust = -2),
             legend.text = element_text(size=12), axis.text = element_text(size=18))
p
dev.off()

########################################Figure 7: Reduction in Ne##########################################
#change in Ub wiht Ud=2;sb=0.001
df_red_ub <- data.frame("Ub" = c(0.001, 0.001584893, 0.002511886, 0.003981072, 0.006309573, 0.01),
                        "N_sim_fit" = c(6389.40, 5145.17, 4322.85, 3513.40, 2696.69, 2151.06),
                        "N_sim_sec" = c(6288, 5112, 4060, 3225, 2638, 2100),
                        "N_ana" = c(2528.13, 2043.18, 1650.72, 1332.96, 1075.5, 866.687),
                        "N_e_coal" = c(5212.694265627656, 4325.696819001546, 3604.4934791060996, 
                                       2966.6656237152997, 2302.653578295538, 1865.6689174206622)/2)

custom_breaks <- c(0.001, 0.001584893, 0.002511886, 0.003981072, 0.006309573, 0.01)  # Custom tick positions
custom_labels <- c(expression(10^-3),expression(10^-2.8),expression(10^-2.6),
                   expression(10^-2.4),expression(10^-2.2),expression(10^-2))  # Custom tick labels
Ub <- c(0.001, 0.001584893, 0.002511886, 0.003981072, 0.006309573, 0.01)
Ne_N2 <- rep(0.8227204,6)

#plot
tiff("Ne_reduction_Ub.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_point(data = df_red_ub, aes(x=Ub, y=(N_ana/N_sim_sec), col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), size=3)
p <- p+geom_point(data = df_red_ub, aes(x=Ub, y=(N_e_coal/N_sim_fit), col="Coalescent-Ne", shape="Coalescent-Ne"), size=3)
p <- p+geom_point(aes(x=Ub, y=Ne_N2, col="Del only model", shape="Del only model"), size=3)
p <- p+geom_hline(yintercept = 1.0, col="red", linetype="dashed")
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels)+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+labs(x=bquote(U[b]), y=expression("Factor reduction in "~N[e]))+theme_Publication()
p <- p+scale_shape_manual(name = '', 
                          values =c("Fitness-flux-Ne"=16,"Del only model"=15,"Coalescent-Ne"=17), 
                          labels = c('Coalescent-Ne','Del only model', 'Fitness-flux-Ne'))
p <- p+scale_colour_manual(name = ' ', 
                           values =c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9", "Del only model"="red"),
                           labels = c("Fitness-flux-Ne","Coalescent-Ne", "Del only model"))
p <- p+guides(shape=guide_legend(keywidth = 0.05, keyheight = 0.05, ncol =1, nrow=3))+ylim(c(0,1))
p <- p+theme(legend.position = "none",
             legend.text = element_text(size = 18),
             legend.background = element_blank())
p
dev.off()

#change in Ud wiht Ub=0.02;sb=0.001
df_red_ud_2 <- data.frame("Ud" = c(0.2511886, 0.5011872, 1.0, 1.9952623, 3.9810717),
                          "N_sim_fit" = c(471.78, 695.96, 953.70, 1500, 2488.79),
                          "N_sim_sec" = c(390, 591, 954, 1500, 2396),
                          "N_ana" = c(216.286, 313.904, 445.359, 623.097, 864.754),
                          "N_e_coal" = c(471.0206607437395, 682.0131766892991, 901.7554276638273,
                                         1320.4766134408887, 1847.292478102416)/2)

#plot
tiff("Ne_reduction_Ud_high_Ub.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_point(data = df_red_ud_2, aes(x=Ud, y=(N_ana/N_sim_sec), col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), size=3)
p <- p+geom_point(data = df_red_ud_2, aes(x=Ud, y=(N_e_coal/N_sim_fit), col="Coalescent-Ne", shape="Coalescent-Ne"), size=3)
p <- p+geom_point(aes(x=Ud, y=Ne_N2, col="Del. only", shape="Del. only"), size=3)
p <- p+geom_hline(yintercept = 1.0, col="red", linetype="dashed")
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels)+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+labs(x=bquote(U[d]), y=expression("Factor reduction in "~N[e]))+theme_Publication()
p <- p+scale_shape_manual(name = '', 
                          values =c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17), 
                          labels = c('Del. only','Coalescent-Ne', 'Fitness-flux-Ne'))
p <- p+scale_colour_manual(name = ' ', 
                           values =c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9", "Del. only"="red"),
                           labels = c("Fitness-flux-Ne","Coalescent-Ne", "Del. only"))
p <- p+guides(shape=guide_legend(keywidth = 0.05, keyheight = 0.05, ncol =1, nrow=3))+ylim(c(0,1))
p <- p+theme(legend.position = "none", legend.text = element_text(size = 12),legend.background = element_blank())
p
dev.off()
##################################Figure 8: Change in L#################################
df_sim_L_5 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                         "vb" = c(0.0000052944, 0.0000108097, 0.0000171429, 0.0000230964,
                                  0.0000291456, 0.0000336541, 0.0000375358),
                         "vd" = c(0.0001881372, 0.0000978801, 0.0000685222, 0.0000535122,
                                  0.0000456827, 0.0000404204, 0.0000369434))

df_sim_L_10 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                          "vb" = c(0.0000063669, 0.0000127176, 0.0000193463, 0.0000258238,
                                   0.0000323070, 0.0000385541, 0.0000451127),
                          "vd" = c(0.0001676110, 0.0000807581, 0.0000543916, 0.0000420350,
                                   0.0000347172, 0.0000302558, 0.0000269228))

df_sim_L_15 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                          "vb" = c(0.0000068111, 0.0000140052, 0.0000209353, 0.0000276391,
                                   0.0000330717, 0.0000413129, 0.0000469458),
                          "vd" = c(0.0001564669, 0.0000750004, 0.0000498756, 0.0000371000,
                                   0.0000305829, 0.0000264316, 0.0000234413))

df_sim_L_20 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                          "vb" = c(0.0000066063, 0.0000147613, 0.0000213570, 0.0000288557,
                                   0.0000343266, 0.0000411660, 0.0000481832),
                          "vd" = c(0.0001574305, 0.0000721110, 0.0000470977, 0.0000351276,
                                   0.0000283229, 0.0000241884, 0.0000214138))

df_sim_L_30 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                          "vb" = c(0.0000087184, 0.0000143849, 0.0000208343, 0.0000284751,
                                   0.0000353364, 0.0000426022, 0.0000502448),
                          "vd" = c(0.0001490496, 0.0000686262, 0.0000448703, 0.0000330566,
                                   0.0000265097, 0.0000222820, 0.0000194895))

df_sim_L_40 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000, 7000),
                          "vb" = c(0.0000074434, 0.0000140163, 0.0000222611, 0.0000277963,
                                   0.0000365466, 0.0000421360, 0.0000499012),
                          "vd" = c(0.0001440845, 0.0000681599, 0.0000439414, 0.0000320720,
                                   0.0000252834, 0.0000211924, 0.0000182838))

df_sim_L_50 <- data.frame("N" = c(1000, 2000, 3000, 4000, 5000, 6000),
                          "vb" = c(0.0000075167, 0.0000161510, 0.0000210694, 0.0000289393,
                                   0.0000360204, 0.0000438151),
                          "vd" = c(0.0001494698, 0.0000664589, 0.0000422315, 0.0000312493,
                                   0.0000246996, 0.0000204847))
#####curve fitting####
dg_sim <- function(dat, N, epsilon) {
  index_plus <- which(abs(dat$N - (N + epsilon)) < 0.0001)
  index_minus <- which(abs(dat$N - (N - epsilon)) < 0.0001)
  
  return ((dat$vd_pred[index_plus] - dat$vd_pred[index_minus]) / (2 * epsilon))
}

df_sim <- function(dat, N, epsilon) {
  index_plus <- which(abs(dat$N - (N + epsilon)) < 0.0001)
  index_minus <- which(abs(dat$N - (N - epsilon)) < 0.0001)
  
  return ((dat$vb_pred[index_plus] - dat$vb_pred[index_minus]) / (2 * epsilon))
}

fitting_func <- function(dat, nchr){
  # Perform linear regression for vb
  linear_model_vb <- lm(vb ~ N, data = dat)
  
  # Perform exponential regression for vd
  #using an exponential model
  #exponential_model_vd <- nls(vd ~ a * exp(b * N), data = df_sim, start = list(a = 1, b = 0))
  #using power law function
  fPow1 <- function(x, a, b) {a * x^b}
  est1 <- coef(nls(vd ~ fPow1(N, a, b),
                   start=c(a=.00001, b=.00001), data=dat))
  power_law_model <- nls(vd ~ fPow1(N, a, b),
                         start=est1, data=dat)
  
  # Predict values for the whole range of N
  N_range <- seq(1000, 10000, by = 1)  # Adjust the step size as needed
  predictions <- data.frame("N" = N_range,
                            "vb_pred" = predict(linear_model_vb, newdata = data.frame(N = N_range)),
                            "vd_pred" = predict(power_law_model, newdata = data.frame(N = N_range)))
  
  #get Ncrit
  f1 <- approxfun(predictions$vb_pred - predictions$vd_pred, predictions$N, rule=1)
  print(f1(0))
  
  ####factor reduction in Ne####
  s=0.00948704
  Ud=2
  L=100
  chr=nchr
  Udw=Ud/(L*chr/2)
  Udw
  Rw=2/(L*chr/2/nchr)
  Rw
  #Unlinked
  Ne_N=exp(-8 * Ud * s)
  Ne_N
  #linked+unlined Joseph's equation
  Ne_N2=exp(-8 * (Ud-Udw)*s) * exp(-Ud/(2*nchr))
  print(paste("factor reducation in Ne - BGS ", Ne_N2))
  
  Ne_N_sim <- 1834.57/f1(0)
  print(paste("factor reducation in Ne - fit flux ", Ne_N_sim))
  
  #derivative calculation at Ncrit with census pop
  ratio_sim <- df_sim(predictions, round(f1(0)), 1)/abs(dg_sim(predictions, round(f1(0)), 1))
  
  print(ratio_sim)
}

fitting_func(df_sim_L_5, 5)
fitting_func(df_sim_L_10, 10)
fitting_func(df_sim_L_15, 15)
fitting_func(df_sim_L_20, 20)
fitting_func(df_sim_L_30, 30)
fitting_func(df_sim_L_40, 40)
fitting_func(df_sim_L_50, 50)

#input into data table
change_in_L <- data.frame("nchr"=c(5, 10, 15, 20, 30, 40, 50),
                          "Ncrit"=c(6437.904, 5121.491, 4709.829, 4481.494, 4294.786, 4246.677, 4106.957),
                          "red_NE_BGS"=c(0.703853198490414, 0.777642, 0.8039189, 0.8173884, 0.8310837, 0.838017124554173, 0.842204918723706),
                          "red_Ne_ext"=c(0.284963876284358, 0.3582101, 0.3895194, 0.4093657, 0.4271621, 0.432001291882267, 0.446698076391268),
                          "ratio" = c(1.117898, 1.011884, 0.9634023, 0.9085198,  0.902947, 0.9188318, 0.859119))

#plotting fluxes
tiff("change_L_Ncrit.tiff", units="in", width=5, height=5, res=300)
p1 <- ggplot()
p1 <- p1+geom_point(data = change_in_L, aes(x=nchr, y=Ncrit), size=3)
p1 <- p1+geom_hline(yintercept = 1834.57, col="red", linetype="dashed")
p1 <- p1+labs(x="Number of Chromosomes", y=expression(N[crit]))+theme_Publication()
p1 <- p1+theme(axis.text.x = element_text(size=18),
               axis.text.y = element_text(size=18))
p1
dev.off()

tiff("change_L_ratio.tiff", units="in", width=5, height=5, res=300)
p1 <- ggplot()
p1 <- p1+geom_point(data = change_in_L, aes(x=nchr, y=ratio), size=3)
p1 <- p1+geom_hline(yintercept = 0.85, col="red", linetype="dashed")
p1 <- p1+labs(x="Number of Chromosomes", y=expression("Drought:Meltdown ratio"))+ylim(0, 1.13)+theme_Publication()
p1 <- p1+theme(axis.text.x = element_text(size=18),
               axis.text.y = element_text(size=18))
p1
dev.off()

tiff("change_L_redNe.tiff", units="in", width=5, height=5, res=300)
p1 <- ggplot()
p1 <- p1+geom_point(data = change_in_L, aes(x=nchr, y=100*(1-red_NE_BGS)/(1-red_Ne_ext), shape="Del only model"), size=3)
p1 <- p1+geom_hline(yintercept = 100.0, col="red", linetype="dashed")
p1 <- p1+scale_shape_manual(name = '', 
                            values =c("Fitness-flux-Ne"=16,"Del only model"=15), 
                            labels = c('Del only model','Fitness-flux-Ne'))
p1 <- p1+labs(x="Number of Chromosomes", y=expression(atop(paste("% reduction in fitness-flux ", N[e]),"explained by del only model")))+ylim(0,100)+theme_Publication()
p1 <- p1+guides(shape=guide_legend(ncol =1, nrow=2))
p1 <- p1+theme(legend.position = c(0.83, 1.03), axis.title.y = element_text(size=18),
               legend.text = element_text(size = 12),legend.background = element_blank())
p1
dev.off()
####################################Supplementary Figures####################################################
###############################Figure S1: Comparing DFEs###########################
dat_cdf_kim <- read.csv("kim_comp_CDF.csv", header = F)
colnames(dat_cdf_kim) <- c("s", "density")

custom_breaks <- c(1e-07, 1e-04, 1e-01)  # Custom tick positions
custom_labels <- c(expression(10^-7), expression(10^-4), expression(10^-1))  # Custom tick labels

tiff("dfe_compare_cdf.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_cdf_kim, aes(x=10^s, y=density, col="Kim et al."), linewidth =1)
p <- p+geom_line(data = dat_cdf_boyko, aes(x=10^s, y=density, col="Boyko et al."), linewidth =1)
p <- p +scale_x_log10(breaks=custom_breaks,
                      labels=custom_labels)
p <- p+scale_colour_manual(name = ' ', 
                           values =c("Kim et al."="#000000","Boyko et al."="#CCCCCC"), 
                           labels = c("Boyko et al.","Kim et al."))
p <- p+annotation_logticks()
p <- p+labs(x=bquote(-x[d]), y=expression("Cumulative probability"))+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 5)))
p <- p+theme(legend.position = c(0.5,1.05), legend.background = element_blank())
p
dev.off()

#############################Figure S2: calculating Ncrit for second verison of W(2000)###################################
#Ncrit to Ud/Ub for three models Kim et al. (2017)
dat_ubud_ncrit_compare_kim <- read.csv("Ncrit_compare.csv", header = F)
colnames(dat_ubud_ncrit_compare_kim) <- c("Ub_Ud", "Ncrit_W2","Ncrit_ana")
dat_ubud_ncrit_compare_kim <- dat_ubud_ncrit_compare_kim[-22,]

dat_ubud_ncrit_all_kim <- read.csv("Ncrit_all.csv", header = F)
colnames(dat_ubud_ncrit_all_kim) <- c("Ub_Ud", "Ncrit_ana","Ncrit_W","Ncrit_MD")

dat_ubud_ncrit_compare_kim_udvarying <- read.csv("Ncrit_compare_Udvarying.csv", header = F)
colnames(dat_ubud_ncrit_compare_kim_udvarying) <- c("Ub_Ud", "Ncrit_W2","Ncrit_ana")
dat_ubud_ncrit_compare_kim <- dat_ubud_ncrit_compare_kim[-22,]


custom_breaks <- c(10, 100, 1000, 10000)  # Custom tick positions
custom_labels <- c("10", "100", "1000", "10000")  # Custom tick labels

tiff("ncrit_analtyical_compare.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ubud_ncrit_compare_kim, aes(x=Ub_Ud, y=Ncrit_ana, col="This paper", linetype="Kim"), linewidth =1)
p <- p+geom_line(data = dat_ubud_ncrit_compare_kim, aes(x=Ub_Ud, y=Ncrit_W2, col="Whitlock (2000)V2", linetype="Kim"), linewidth =1)
p <- p+geom_line(data = dat_ubud_ncrit_all_kim, aes(x=Ub_Ud, y=Ncrit_W, col="Whitlock (2000)", linetype="Kim"), linewidth =1)
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels,
                     name=expression(U[d]~"/"~U[b]))
p <- p+geom_point(aes(x=1000,y=1834.57),col="#F05039",shape=42, size=10)
p <- p+scale_colour_manual(name = '', 
                           values =c("This paper"="#000000","Whitlock (2000)V2"="#7fc97f","Whitlock (2000)"="#386cb0"), 
                           labels = c("This paper","Whitlock (2000)V2","Whitlock (2000)"))
#p <- p+annotation_logticks()
p <- p+labs(y=expression(N[crit]))+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 5)))
p <- p+theme(legend.position = "none", legend.background = element_blank())
p
dev.off()

#######################################Figure S3: Ncrit as a function of xb###############################################
#Ncrit sensitivity to xb
dat_sb_kim_100 <- read.csv("sbtablekim_Ncrit_100.csv", header = F)
colnames(dat_sb_kim_100) <- c("sb", "Ncrit")
dat_sb_kim_1000 <- read.csv("sbtablekim_Ncrit_1000.csv", header = F)
colnames(dat_sb_kim_1000) <- c("sb", "Ncrit")
dat_sb_boyko_100 <- read.csv("sbtableboyko_Ncrit_100.csv", header = F)
colnames(dat_sb_boyko_100) <- c("sb", "Ncrit")
dat_sb_boyko_1000 <- read.csv("sbtableboyko_Ncrit_1000.csv", header = F)
colnames(dat_sb_boyko_1000) <- c("sb", "Ncrit")

custom_breaks <- c(0.0001, 0.001, 0.01)  # Custom tick positions
custom_labels <- c("0.0001", "0.001", "0.01")  # Custom tick labels

tiff("sb_analtyical_Ncrit.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_sb_kim_100, aes(x=sb, y=Ncrit, col="Kim et al.", linetype="mut. rate ratio = 100"), linewidth =1)
p <- p+geom_line(data = dat_sb_kim_1000, aes(x=sb, y=Ncrit, col="Kim et al.", linetype="mut. rate ratio = 1000"), linewidth =1)
p <- p+geom_line(data = dat_sb_boyko_100, aes(x=sb, y=Ncrit, col="Boyko et al.", linetype="mut. rate ratio = 100"), linewidth =1)
p <- p+geom_line(data = dat_sb_boyko_1000, aes(x=sb, y=Ncrit, col="Boyko et al.", linetype="mut. rate ratio = 1000"), linewidth =1)
p <- p+scale_x_log10(breaks = custom_breaks,
                     labels = custom_labels)  
p <- p+scale_colour_manual(name = ' ', 
                           values =c("Kim et al."="#000000","Boyko et al."="#CCCCCC"), 
                           labels = c("Boyko et al.","Kim et al."))
p <- p+scale_linetype_manual(name = '',
                             values=c("mut. rate ratio = 100"=2,"mut. rate ratio = 1000"=1),
                             labels = c(expression(U[d]~"/"~U[b]~"= 100"), expression(U[d]~"/"~U[b]~"= 1000")))
p <- p+labs(x=bquote(bar(x)[b]), y=expression(N[crit]))+theme_Publication()
p <- p+guides(color = guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2),
              linetype = guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = c(0.7, 0.8), legend.background = element_blank())
p
dev.off()

#############################Figure S4: calculating dorught:meltdown ratio for changing ud and ub with env change###################################
#with offset Kim et al. ratio varying Ud fixed Ub
dat_ud_kim_ub0.002_highoffset <- read.csv("offset_udvaryingub0.002_kim.csv", header = F)
colnames(dat_ud_kim_ub0.002_highoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_kim_ub0.02_highoffset <- read.csv("offset_udvaryingub0.02_kim.csv", header = F)
colnames(dat_ud_kim_ub0.02_highoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_kim_ub0.002_lowoffset <- read.csv("lowoffset_udvaryingub0.002_kim.csv", header = F)
colnames(dat_ud_kim_ub0.002_lowoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_kim_ub0.02_lowoffset <- read.csv("lowoffset_udvaryingub0.02_kim.csv", header = F)
colnames(dat_ud_kim_ub0.02_lowoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_boyko_ub0.002_highoffset <- read.csv("offset_udvaryingub0.002_boyko.csv", header = F)
colnames(dat_ud_boyko_ub0.002_highoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_boyko_ub0.02_highoffset <- read.csv("offset_udvaryingub0.02_boyko.csv", header = F)
colnames(dat_ud_boyko_ub0.02_highoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_boyko_ub0.002_lowoffset <- read.csv("lowoffset_udvaryingub0.002_boyko.csv", header = F)
colnames(dat_ud_boyko_ub0.002_lowoffset) <- c("Ud", "ncrit", "delratio", "ratio")

dat_ud_boyko_ub0.02_lowoffset <- read.csv("lowoffset_udvaryingub0.02_boyko.csv", header = F)
colnames(dat_ud_boyko_ub0.02_lowoffset) <- c("Ud", "ncrit", "delratio", "ratio")

custom_breaks <- c(1, 3, 10, 30)  # Custom tick positions
custom_labels <- c("1","3","10","30")  # Custom tick labels

#ben del ratio
tiff("ud_analtyical_offset_si.tiff", units="in", width=4, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ud_kim_ub0.002_highoffset,
                 aes(x=Ud, y=ratio, col="0.002", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_highoffset,
                 aes(x=Ud, y=ratio, col="0.02", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.002_lowoffset,
                 aes(x=Ud, y=ratio, col="0.002", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_kim_ub0.02_lowoffset,
                 aes(x=Ud, y=ratio, col="0.02", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ud_boyko_ub0.002_highoffset,
                 aes(x=Ud, y=ratio, col="0.002", linetype="High offset"), linewidth =1)+
  geom_point(data = dat_ud_boyko_ub0.002_highoffset,
             aes(x=Ud, y=ratio), size=1.2, shape=17)
p <- p+geom_line(data = dat_ud_boyko_ub0.02_highoffset,
                 aes(x=Ud, y=ratio, col="0.02", linetype="High offset"), linewidth =1)+
  geom_point(data = dat_ud_boyko_ub0.02_highoffset,
             aes(x=Ud, y=ratio), size=1.2, shape=17)
p <- p+geom_line(data = dat_ud_boyko_ub0.002_lowoffset,
                 aes(x=Ud, y=ratio, col="0.002", linetype="Low offset"), linewidth =1)+
  geom_point(data = dat_ud_boyko_ub0.002_lowoffset,
             aes(x=Ud, y=ratio), size=1.2, shape=17)
p <- p+geom_line(data = dat_ud_boyko_ub0.02_lowoffset,
                 aes(x=Ud, y=ratio, col="0.02", linetype="Low offset"), linewidth =1)+
  geom_point(data = dat_ud_boyko_ub0.02_lowoffset,
             aes(x=Ud, y=ratio), size=1.2, shape=17)
p <- p+scale_x_log10(name=expression(U[d]))+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[b]), 
                           values =c("0.002" = "#9e9e9e","0.02"="#2298e6"), 
                           labels = c("0.002","0.02"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('High offset','Low offset'))
p <- p+scale_y_log10(limits=c(0.35,18.5),
                     breaks=custom_breaks,
                     labels=custom_labels)
p <- p+geom_point(aes(x=2,y=0.853449),col="#F05039",shape=42, size=10)
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=0.9)
p <- p+labs(y="Drought : Meltdown ratio")+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =2, nrow=1))
p <- p+theme(legend.position = "none",legend.title=element_text(size=16), legend.background = element_blank(),
             legend.text = element_text(size=16), axis.text = element_text(size=20))
p
dev.off()

#with offset Kim et al. ratio varying Ub fixed Ud
dat_ub_kim_ud2_highoffset <- read.csv("offset_ubvaryingud2_kim.csv", header = F)
colnames(dat_ub_kim_ud2_highoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_kim_ud1_highoffset <- read.csv("offset_ubvaryingud1_kim.csv", header = F)
colnames(dat_ub_kim_ud1_highoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_kim_ud2_lowoffset <- read.csv("lowoffset_ubvaryingud2_kim.csv", header = F)
colnames(dat_ub_kim_ud2_lowoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_kim_ud1_lowoffset <- read.csv("lowoffset_ubvaryingud1_kim.csv", header = F)
colnames(dat_ub_kim_ud1_lowoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_boyko_ud2_highoffset <- read.csv("offset_ubvaryingud2_boyko.csv", header = F)
colnames(dat_ub_boyko_ud2_highoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_boyko_ud1_highoffset <- read.csv("offset_ubvaryingud1_boyko.csv", header = F)
colnames(dat_ub_boyko_ud1_highoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_boyko_ud2_lowoffset <- read.csv("lowoffset_ubvaryingud2_boyko.csv", header = F)
colnames(dat_ub_boyko_ud2_lowoffset) <- c("Ub", "ncrit", "delratio", "ratio")

dat_ub_boyko_ud1_lowoffset <- read.csv("lowoffset_ubvaryingud1_boyko.csv", header = F)
colnames(dat_ub_boyko_ud1_lowoffset) <- c("Ub", "ncrit", "delratio", "ratio")

custom_breaks <- c(0.0001, 0.001, 0.01, 0.1)  # Custom tick positions
custom_labels <- c(expression(10^-4),expression(10^-3),expression(10^-2),
                   expression(10^-1))  # Custom tick labels

custom_breaks_y <- c(1, 3, 10)  # Custom tick positions
custom_labels_y <- c("1","3","10")  # Custom tick labels

#ben del ratio
tiff("ub_analtyical_offset_si.tiff", units="in", width=5, height=5, res=300)
p <- ggplot()
p <- p+geom_line(data = dat_ub_kim_ud2_highoffset,
                 aes(x=Ub, y=ratio, col="2", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_highoffset,
                 aes(x=Ub, y=ratio, col="1", linetype="High offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud2_lowoffset,
                 aes(x=Ub, y=ratio, col="2", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_kim_ud1_lowoffset,
                 aes(x=Ub, y=ratio, col="1", linetype="Low offset"), linewidth =1)
p <- p+geom_line(data = dat_ub_boyko_ud2_highoffset,
                aes(x=Ub, y=ratio, col="2", linetype="High offset"), linewidth =1)+
geom_point(data = dat_ub_boyko_ud2_highoffset,
          aes(x=Ub, y=ratio), size=1.2, shape=17)
p <- p+geom_line(data = dat_ub_boyko_ud1_highoffset,
                aes(x=Ub, y=ratio, col="1", linetype="High offset"), linewidth =1)+
geom_point(data = dat_ub_boyko_ud1_highoffset,
          aes(x=Ub, y=ratio), size=1.2, shape=17)
p <- p+geom_line(data = dat_ub_boyko_ud2_lowoffset,
                aes(x=Ub, y=ratio, col="2", linetype="Low offset"), linewidth =1)+
geom_point(data = dat_ub_boyko_ud2_lowoffset,
          aes(x=Ub, y=ratio), size=1.2, shape=17)
p <- p+geom_line(data = dat_ub_boyko_ud1_lowoffset,
                aes(x=Ub, y=ratio, col="1", linetype="Low offset"), linewidth =1)+
geom_point(data = dat_ub_boyko_ud1_lowoffset,
          aes(x=Ub, y=ratio), size=1.2, shape=17)
p <- p+scale_x_continuous(trans=c("log10","reverse"),name=expression(U[b]),
                          breaks=custom_breaks,
                          labels=custom_labels)+annotation_logticks(sides = "bl", outside = TRUE)+ coord_cartesian(clip = "off")
p <- p+scale_colour_manual(name = expression(U[d]), 
                           values =c("2" = "#9e9e9e","1"="#2298e6"), 
                           labels = c("1","2"))
p <- p+scale_linetype_manual(name = '',
                             values=c("Low offset"="dashed", "High offset"="solid"),
                             labels = c('Fast env. change','Slow env. change'))
p <- p+scale_y_log10(limits=c(0.35,18.5),
                     breaks=custom_breaks_y,
                     labels=custom_labels_y)
p <- p+geom_point(aes(x=0.002,y=0.853449),col="#F05039",shape=42, size=10)
p <- p+geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth =.9)
p <- p+labs(y="Drought : Meltdown ratio")+theme_Publication()
p <- p+guides(color = guide_legend(override.aes = list(size = 2),ncol =1, nrow=2), 
              linetype=guide_legend(keywidth = 2, keyheight = 1, ncol =1, nrow=2))
p <- p+theme(legend.position = "none",legend.title=element_text(size=14), legend.background = element_blank(),
             legend.text = element_text(size=12), axis.text = element_text(size=20))
p
dev.off()
