#######set preferred theme########
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(ggpubr)
library(ggridges)
library(plyr)
library(corrplot)
library(latex2exp)
library(ggallin)
library(extrafont)
library(patchwork)
library(cowplot)
font_import(prompt = FALSE)
loadfonts(device = "win")
loadfonts(device = "pdf")

theme_Publication <- function(base_size=9, base_family="Times New Roman") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.tag = element_text(family = base_family, face = "italic", size = base_size),
            text = element_text(),
            panel.background = element_rect(fill = "white", colour = NA), 
            plot.background = element_rect(fill = "white", colour = NA),
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

save_procb <- function(plot, filename, col = c("1col", "2col"), height_mm = 100) {
  col <- match.arg(col)
  width_mm <- ifelse(col == "1col", 80, 167)   # Proceedings B column widths
  ggsave(filename, plot = plot,
         width = width_mm/25.4, height = height_mm/25.4,
         units = "in", dpi = 600,
         compression = "lzw", device = "tiff")
}

label_thinspace <- label_number(
  big.mark = "\u2009",  # Unicode thin space
  decimal.mark = "."
)

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

#adding environmental deviation, delta=3e-05
dat_del$flux_offset <- dat_del$flux+1.5e-05

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

f3 <- approxfun(predictions_ana$vb_pred - predictions_ana$vd_pred_offset, predictions_ana$N, rule=1)
f3(0)
#plotting fluxes with simulation data points and env 
# Plot 2A
p1 <- ggplot() +
  geom_line(data = dat_del, aes(x = N, y = flux, col = "Deleterious"), linewidth = 0.7) +
  geom_line(data = dat_del, aes(x = N, y = flux_offset, col = "Deleterious+env. change"), linewidth = 0.7) +
  geom_line(data = dat_ben, aes(x = N, y = flux, col = "Beneficial"), linewidth = 0.7) +
  geom_vline(xintercept = c(3666.41, 4656.43),
             linewidth = 0.7, col = c("#ef3b2c","darkorange"),
             alpha = 0.8, linetype = "dotted") +
  xlim(c(500, 7000)) +
  scale_y_continuous(
    limits = c(0, 0.00031),
    breaks = c(0.00003, 0.00015, 0.0003),
    labels = c( expression("3 ×" ~ 10^-5), expression("1 ×" ~ 10^-4), expression("3 ×" ~ 10^-4))
  ) +
  labs(y = "Fitness change per generation", x = NULL) +   # remove x-label
  scale_colour_manual(name = '', 
                      values = c("Deleterious"="#ef3b2c","Deleterious+env. change"="darkorange","Beneficial"="#386cb0"),
                      labels = c("Beneficial","Deleterious+env. change","Deleterious")) +
  theme_Publication(base_size = 9) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),          # remove x-axis label
        axis.text.x  = element_blank(),          # remove x-axis ticks
        plot.margin = margin(t = .5, r = 1, b = .5, l = 1, unit = "mm"))

# Plot 2B
p2 <- ggplot() +
  geom_line(data = dat_del_der, aes(x = N, y = der, col = "Deleterious"), linewidth = 0.7) +
  geom_line(data = dat_ben_der, aes(x = N, y = der, col = "Beneficial"), linewidth = 0.7) +
  geom_vline(xintercept = c(3666.41, 4656.43),
             linewidth = 0.7, col = c("#ef3b2c","darkorange"),
             alpha = 0.8, linetype = "dotted") +
  scale_y_continuous(
    limits = c(0, 0.00000004),
    breaks = c(0, 0.00000001, 0.00000002, 0.00000003, 0.00000004),
    labels = c("0", expression("1 ×" ~ 10^-8), expression("2 ×" ~ 10^-8),
               expression("3 ×" ~ 10^-8), expression("4 ×" ~ 10^-8))
  ) +
  scale_x_continuous(
    limits = c(500, 7000),
    breaks = c(2000, 4000, 6000),                  # labeled ticks only
    labels = label_number(big.mark = "\u202F")
  ) +
  labs(x = "Census population size", y = "Fitness flux derivative") +
  scale_colour_manual(name = ' ', 
                      values = c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"),
                      labels = c("Beneficial","Deleterious")) +
  theme_Publication(base_size = 9) +
  theme(legend.position = "none",
        plot.margin = margin(t = .5, r = 1, b = .5, l = 1, unit = "mm"))

combined <- (p1 / p2) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "italic", size = 9))

tiff("figure2_combined.tiff",
     width = 80/25.4, height = 140/25.4, units = "in",
     res = 600, compression = "lzw")

print(combined)

dev.off()
##############Figure 6
linewidth=0.7
# Figure 6A
pA <- ggplot() +
  geom_line(data = dat_del, aes(x=N, y=flux, col="Deleterious"), alpha=.5, linewidth =linewidth) +
  geom_line(data = dat_del, aes(x=N, y=flux_offset, col="Del.+env. change"), alpha=.5, linewidth =linewidth) +
  geom_line(data = dat_ben, aes(x=N, y=flux, col="Beneficial"), alpha=.5, linewidth =linewidth) +
  geom_vline(xintercept = c(3666.41, 4656.43), col=c("#ef3b2c","darkorange"), alpha=.8, linetype="dotted") +
  geom_point(data = df_sim, aes(x=N, y=vd, col="Deleterious"), size=1) +
  geom_point(data = df_sim, aes(x=N, y=vb, col="Beneficial"), size=1) +
  xlim(500, 7000) +
  scale_y_continuous(limits = c(0, 0.00031),
                     breaks = c(0.00003, 0.00015, 0.0003),
                     labels = c( expression("3 ×" ~ 10^-5), expression("1 ×" ~ 10^-4), expression("3 ×" ~ 10^-4))) +
  labs(x="Census pop. size", y="Fitness change per gen.") +
  scale_colour_manual(name = '', 
                      values =c("Deleterious"="#ef3b2c","Del.+env. change"="darkorange","Beneficial"="#386cb0"),
                      labels = c('Beneficial','Del.+env. change','Deleterious'))+
  theme_Publication(base_size=9) +
  theme(legend.position="none",
        axis.title.x = element_blank(),          # remove x-axis label
        axis.text.x  = element_blank(),
        axis.title.y = element_text(size=8),
        plot.margin = margin(t = .5, r = 1, b = .5, l = 1, unit = "mm"))

#Figure 6B
s=0.00948704
Ud=2
L=100
chr=23
Udw=Ud/(L*chr/2)
Udw
Rw=2/(L*chr/2/23)
Rw
#Unlinked
Ne_N=exp(-4 * Ud * s)
Ne_N
#linked+unlined Joseph's equation
Ne_N2=exp(-4 * (Ud-Udw)*s) * exp(-Ud/(2*23))
Ne_N2

Ne_N_sim <- 3666.41/4402.57
Ne_N_sim
Ne_N2 - Ne_N_sim
(Ne_N2 - Ne_N_sim)/Ne_N2*100

predictions$Ne <- predictions$N * Ne_N2
df_sim$Ne <- df_sim$N * Ne_N2

predictions$Ne_ext <- predictions$N * Ne_N_sim
df_sim$Ne_ext <- df_sim$N * Ne_N_sim

# Figure 6B
pB <- ggplot() +
  geom_line(data = dat_del, aes(x=N, y=flux, col="Deleterious"), alpha=0.5, linewidth=linewidth) +
  geom_line(data = dat_ben, aes(x=N, y=flux, col="Beneficial"), alpha=0.5, linewidth=linewidth) +
  geom_vline(xintercept=3666.41, linetype="dotted", alpha=0.8, linewidth=linewidth) +
  geom_point(data = df_sim, aes(x=Ne_ext, y=vd, col="Deleterious"), size=1) +
  geom_point(data = df_sim, aes(x=Ne_ext, y=vb, col="Beneficial"), size=1) +
  xlim(500, 7000) +
  scale_y_continuous(limits = c(0, 0.00031),
                     breaks = c(0.00003, 0.00015, 0.0003),
                     labels = c( expression("3 ×" ~ 10^-5), expression("1 ×" ~ 10^-4), expression("3 ×" ~ 10^-4))) +
  labs(x="Effective population size", y="Fitness flux") +
  scale_colour_manual(name = ' ', 
                      values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))+
  theme_Publication(base_size=9) +
  theme(legend.position="none",
        axis.title.x = element_blank(),          # remove x-axis label
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),          # remove x-axis label
        axis.text.y  = element_blank(),
        plot.margin = margin(t = .5, r = .5, b = .5, l = .5, unit = "mm"))


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
epsilon <- 0.8327886

# Calculate the change in slope for vb with N using epsilon = 0.4167043
change_in_slope_vb_ne <- numeric(length(N_range) - 2)

change_in_slope_vb_ne <- sapply(predictions$Ne_ext[2:(length(N_range) - 1)], function(N) df_sim_ne(N, epsilon))

change_in_slope_vb_ne_sim <- as.vector(sapply(df_sim$Ne_ext[-1], function(N) df_sim_ne(N, epsilon)))

# Calculate the change in slope for vd with N using epsilon = 0.4167043
change_in_slope_vd_ne <- numeric(length(N_range) - 2)

change_in_slope_vd_ne <- sapply(predictions$Ne_ext[2:(length(N_range) - 1)], function(N) dg_sim_ne(N, epsilon))

change_in_slope_vd_ne_sim <- as.vector(sapply(df_sim$Ne_ext[-1], function(N) dg_sim_ne(N, epsilon)))


# Figure 6C
pC <- ggplot() +
  geom_line(data = dat_del_der, aes(x=N, y=der, col="Deleterious"), alpha=0.5, linewidth=linewidth) +
  geom_line(data = dat_ben_der, aes(x=N, y=der, col="Beneficial"), alpha=0.5, linewidth=linewidth) +
  geom_vline(xintercept=3666.41, linetype="dotted", alpha=0.8, linewidth=linewidth) +
  geom_point(aes(x=df_sim$N[-1], y=abs(change_in_slope_vd_sim), col="Deleterious"), size=1) +
  geom_point(aes(x=df_sim$N[-1], y=change_in_slope_vb_sim, col="Beneficial"), size=1) +
  scale_y_continuous(limits=c(0, 0.00000004),
                     breaks=c(0, 0.00000001, 0.00000002, 0.00000003, 0.00000004),
                     labels=c("0", expression("1 ×" ~ 10^-8), expression("2 ×" ~ 10^-8),
                              expression("3 ×" ~ 10^-8), expression("4 ×" ~ 10^-8))) +
  scale_x_continuous(
    limits = c(500, 7000),
    breaks = c(2000, 4000, 6000),                  # labeled ticks only
    labels = label_number(big.mark = "\u202F")
  ) +
  scale_colour_manual(name = ' ', 
                      values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))+
  labs(x="Census pop. size", y="Fitness flux derivative") +
  theme_Publication(base_size=9) +
  theme(legend.position="none",axis.title.y = element_text(size=8),
        plot.margin = margin(t = .5, r = 1, b = .5, l = 1, unit = "mm"))

# Figure 6D
pD <- ggplot() +
  geom_line(data = dat_del_der, aes(x=N, y=der, col="Deleterious"), alpha=0.5, linewidth=linewidth) +
  geom_line(data = dat_ben_der, aes(x=N, y=der, col="Beneficial"), alpha=0.5, linewidth=linewidth) +
  geom_vline(xintercept=3666.41, linetype="dotted", alpha=0.8,linewidth=linewidth) +
  geom_point(aes(x=df_sim$Ne_ext[-1], y=abs(change_in_slope_vd_ne_sim), col="Deleterious"), size=1) +
  geom_point(aes(x=df_sim$Ne_ext[-1], y=change_in_slope_vb_ne_sim, col="Beneficial"), size=1) +
  scale_y_continuous(limits=c(0, 0.00000004),
                     breaks=c(0, 0.00000001, 0.00000002, 0.00000003, 0.00000004),
                     labels=c("0", expression("1 ×" ~ 10^-8), expression("2 ×" ~ 10^-8),
                              expression("3 ×" ~ 10^-8), expression("4 ×" ~ 10^-8))) +
  scale_x_continuous(
    limits = c(500, 7000),
    breaks = c(2000, 4000, 6000),                  # labeled ticks only
    labels = label_number(big.mark = "\u202F")
  ) +
  labs(x="Effective pop. size", y="Fitness flux derivative") +
  theme_Publication(base_size=9) +
  scale_colour_manual(name = ' ', 
                      values =c("Deleterious"="#ef3b2c","Beneficial"="#386cb0"), labels = c('Beneficial','Deleterious'))+
  theme(legend.position="none",
        axis.title.y = element_blank(),          # remove x-axis label
        axis.text.y  = element_blank(),
        plot.margin = margin(t = .5, r = 1, b = .5, l = 1, unit = "mm"))

fig6 <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "italic", size = 9))

# Save for 1-column (≈83 mm)
tiff("figure6_combined_1col.tiff",
     width=83/25.4, height=84/25.4, units="in",
     res=600, compression="lzw")

print(fig6)
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

# Reduce line thickness
line_width <- 0.7
# Top plot
p1 <- ggplot(dat_env_offset_kim, aes(x = 0.1/offset, col = as.factor(UbUd))) +
  geom_line(aes(y = bendelratio), linewidth = line_width) +
  geom_point(aes(x = 0.1/(0.000015), y = 1.43349), col = "darkorange", shape = 42, size = 7) +
  geom_hline(yintercept = 1, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(0.1/2.2e-05, 0.1/4.5e-06), col = c("#2298e6","#9e9e9e"),
             linetype = "dotdash", linewidth = line_width, alpha = 0.5) +
  scale_x_log10() +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(labels = label_number(big.mark = "\u202F")) +
  scale_colour_manual(
    name = expression(U[d]/U[b]),
    values = c("100" = "#2298e6", "1000" = "#9e9e9e"),
    labels = c("100", "1000")
  )+
  labs(y = expression(bold("Drought:Meltdown ratio"))) +
  theme_Publication(base_size = 9) +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "mm"),
    legend.position = "none")

# Middle plot
p2 <- ggplot(dat_env_offset_kim, aes(x = 0.1/offset, col = as.factor(UbUd))) +
  geom_line(aes(y = ratio), linewidth = line_width) +
  geom_point(aes(x = 0.1/(0.000015), y = 0.000015/0.0000221823), col = "darkorange", shape = 42, size = 7) +
  geom_abline(intercept = 3, slope = -1, linetype = "dashed") +
  geom_hline(yintercept = 1, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(0.1/2.2e-05, 0.1/4.5e-06), col = c("#2298e6","#9e9e9e"),
             linetype = "dotdash", linewidth = line_width, alpha = 0.5) +
  scale_x_log10() +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(breaks = c(0.01,0.1,1,10,100),
                labels = c("0.01","0.1","1","10","100")) +
  scale_colour_manual(
    name = expression(U[d]/U[b]),
    values = c("100" = "#2298e6", "1000" = "#9e9e9e"),
    labels = c("100", "1000")
  )+
  labs(y = expression(bold("Env. change to del. flux ratio at " * N[crit]))) +
  theme_Publication(base_size = 9) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "mm"),
    legend.position = "none")

# Bottom plot
p3 <- ggplot(dat_env_offset_kim, aes(x = 0.1/offset, col = as.factor(UbUd))) +
  geom_line(aes(y = Ncrit), linewidth = line_width) +
  geom_point(aes(x = 0.1/(0.000015), y = 4656.43), col = "darkorange", shape = 42, size = 7) +
  geom_hline(yintercept = c(1246.54, 3666.41), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(0.1/2.2e-05, 0.1/4.5e-06), col = c("#2298e6","#9e9e9e"),
             linetype = "dotdash", linewidth = line_width, alpha = 0.5) +
  scale_x_log10(labels = label_number(big.mark = "\u202F")) +
  scale_y_log10(labels = label_number(big.mark = "\u202F")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_colour_manual(
    name = expression(U[d]/U[b]),
    values = c("100" = "#2298e6", "1000" = "#9e9e9e"),
    labels = c("100", "1000")
  ) +
  labs(x = "Generations needed for env.\nchange to reduce fitness by 10%", y = expression(bold(N[crit]))) +
  theme_Publication(base_size = 9) +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = "none")

# Combine with patchwork, add panel labels
combined <- (p1 / p2 / p3) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "italic", size = 9))

# Save Proceedings B 1-column TIFF
tiff("figure3_combined.tiff",
     width = 80/25.4, height = 200/25.4, units = "in", # 1 col, 3 rows
     res = 600, compression = "lzw")

print(combined)
dev.off()

###########################################Figure 4: ben:del ratio no env change####################################
#ratio sensitvity to Ud/Ub
dat_ubud_kim <- read.csv("udubtablekim_ratio.csv", header = F)
colnames(dat_ubud_kim) <- c("Ub_Ud", "del_flux","Ncrit", "bendelratio")
#ratio sensitvity to mean sb
dat_sb_kim_100 <- read.csv("sbtablekim_100.csv", header = F)
colnames(dat_sb_kim_100) <- c("sb", "del_flux","Ncrit", "bendelratio")
dat_sb_kim_1000 <- read.csv("sbtablekim_1000.csv", header = F)
colnames(dat_sb_kim_1000) <- c("sb", "del_flux","Ncrit", "bendelratio")
#ratio sensitvity to mean sd
dat_sd_100 <- read.csv("sdtablekim_100.csv", header = F)
colnames(dat_sd_100) <- c("sd", "del_flux","Ncrit", "bendelratio")
dat_sd_100 <- dat_sd_100[-c(19,20,21),]
dat_sd_1000 <- read.csv("sdtablekim_1000.csv", header = F)
colnames(dat_sd_1000) <- c("sd", "del_flux","Ncrit", "bendelratio")
dat_sd_1000 <- dat_sd_1000[-c(19,20,21),]

#Ratio
# Top-left plot
p1 <- ggplot(dat_ubud_kim, aes(x = Ub_Ud, y = bendelratio)) +
  geom_line(linewidth = line_width) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5, col = "red", linewidth = line_width) +
  geom_point(aes(x = 1000, y = 0.853449), col = "#F05039", shape = 42, size = 7) +
  scale_x_log10(labels = label_number(big.mark = "\u202F")) +
  annotation_logticks(sides = "b",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_continuous(breaks= seq(0.2,1, by=.2), limits = c(0,1))+
  labs(x = expression(bold(U[d]/U[b])), y = "Drought:Meltdown ratio") +
  theme_Publication(base_size = 9) +
  theme(plot.margin = margin(1,0,1,1), legend.position = "none")

# Top-right plot
p2 <- ggplot() +
  geom_line(data = dat_sb_kim_100, aes(x = sb, y = bendelratio, linetype = "mut. rate ratio = 100"), linewidth = line_width) +
  geom_line(data = dat_sb_kim_1000, aes(x = sb, y = bendelratio, linetype = "mut. rate ratio = 1000"), linewidth = line_width) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5, col = "red", linewidth = line_width) +
  geom_point(aes(x = 0.001, y = 0.853449), col = "#F05039", shape = 42, size = 7) +
  scale_x_log10(breaks = c(0.001, 0.01), labels = c("0.001","0.01")) +
  annotation_logticks(sides = "b",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_linetype_manual(values = c("mut. rate ratio = 100" = 2, "mut. rate ratio = 1000" = 1)) +
  scale_y_continuous(breaks= seq(.2,1, by=.2), limits = c(0,1))+
  labs(x = bquote(bold(bar(s)[b])), y = "Drought:Meltdown ratio") +
  theme_Publication(base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(1,0,1,0), legend.position = "none")

# Bottom centered plot
p3 <- ggplot() +
  geom_line(data = dat_sd_100, aes(x = sd, y = bendelratio, linetype = "mut. rate ratio = 100"), linewidth = line_width) +
  geom_line(data = dat_sd_1000, aes(x = sd, y = bendelratio, linetype = "mut. rate ratio = 1000"), linewidth = line_width) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5, col = "red", linewidth = line_width) +
  geom_point(aes(x = 0.0095, y = 0.853449), col = "#F05039", shape = 42, size = 7) +
  scale_x_log10(limits=c(0.001, 0.1),breaks = c(0.001, 0.01, 0.1), labels = c("0.001","0.01","0.1")) +
  annotation_logticks(sides = "b",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_linetype_manual(values = c("mut. rate ratio = 100" = 2, "mut. rate ratio = 1000" = 1)) +
  scale_y_continuous(breaks= seq(.2,1, by=.2), limits = c(0,1))+
  labs(x = bquote(bold(-bar(s)[d])), y = "Drought:Meltdown ratio") +
  theme_Publication(base_size = 9) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(1,1,1,0), legend.position = "none")

# Combine rows
combined <- (p1|p2|p3) +
  plot_layout(heights = c(1, 0.9)) +  # bottom row slightly shorter
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "italic", size = 9))

tiff("figure4_combined_2col.tiff",
     width = 167/25.4, height = 80/25.4, units = "in",
     res = 600, compression = "lzw")

print(combined)
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

#with offset Kim et al. ratio varying Ub fixed Ud
dat_ub_kim_ud2_highoffset <- read.csv("offset_ubvaryingud2_kim.csv", header = F)
colnames(dat_ub_kim_ud2_highoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

dat_ub_kim_ud1_highoffset <- read.csv("offset_ubvaryingud1_kim.csv", header = F)
colnames(dat_ub_kim_ud1_highoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

dat_ub_kim_ud2_lowoffset <- read.csv("lowoffset_ubvaryingud2_kim.csv", header = F)
colnames(dat_ub_kim_ud2_lowoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

dat_ub_kim_ud1_lowoffset <- read.csv("lowoffset_ubvaryingud1_kim.csv", header = F)
colnames(dat_ub_kim_ud1_lowoffset) <- c("Ub", "del_flux","Ncrit", "bendelratio")

line_width <- 0.6  # slightly thinner lines for compact figures
text_size <- 9    # base font size for legibility
#ben del ratio
#p1
p1 <- ggplot() +
  geom_line(data = dat_ud_kim_ub0.002_highoffset,
            aes(x=Ud, y=bendelratio, col="0.002", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ud_kim_ub0.02_highoffset,
            aes(x=Ud, y=bendelratio, col="0.02", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ud_kim_ub0.002_lowoffset,
                   aes(x=Ud, y=bendelratio, col="0.002", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_ud_kim_ub0.02_lowoffset,
                 aes(x=Ud, y=bendelratio, col="0.02", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=2,y=0.853532),col = "#F05039",shape=42, size=7)+
  scale_x_log10(name=NULL, breaks = c(0.1, 1, 10), labels=c("0.1", "1", "10")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_y_log10(limits=c(0.35,18.5), breaks=c(1,3,10), labels=c("1","3","10")) +
  scale_colour_manual(values=c("0.002"="#9e9e9e","0.02"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y="Drought:Meltdown ratio") +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))


#p2 #offset del ratio
p2 <- ggplot() +
  geom_line(data = dat_ud_kim_ub0.002_highoffset,
            aes(x=Ud, y=1.5e-5/abs(del_flux-1.5e-5), col="0.002", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ud_kim_ub0.02_highoffset,
            aes(x=Ud, y=1.5e-5/abs(del_flux-1.5e-5), col="0.02", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ud_kim_ub0.002_lowoffset,
            aes(x=Ud, y=1e-6/abs(del_flux-1e-6), col="0.002", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_ud_kim_ub0.02_lowoffset,
            aes(x=Ud, y=1e-6/abs(del_flux-1e-6), col="0.02", linetype="Low offset"), linewidth =line_width)+
  scale_x_log10(name=NULL, breaks = c(0.1, 1, 10), labels=c("0.1", "1", "10")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(0.0038,1), breaks=c(0.01, 0.1, 1), labels=c("0.01","0.1","1")) +
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_colour_manual(values=c("0.002"="#9e9e9e","0.02"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression(bold("Env.: Mut. deg. ratio at "*N[crit])))+
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=8),
        plot.margin=margin(.5,.5,.5,.5))

#p3 #Ncrit
p3 <- ggplot() +
  geom_line(data = dat_ud_kim_ub0.002_highoffset,
            aes(x=Ud, y=Ncrit, col="0.002", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ud_kim_ub0.02_highoffset,
            aes(x=Ud, y=Ncrit, col="0.02", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ud_kim_ub0.002_lowoffset,
            aes(x=Ud, y=Ncrit, col="0.002", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_ud_kim_ub0.02_lowoffset,
            aes(x=Ud, y=Ncrit, col="0.02", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=2,y=3666.41),col = "#F05039",shape=42, size=7)+
  scale_x_log10(name=expression(bold(U[d])), breaks = c(0.1, 1, 10), labels=c("0.1", "1", "10")) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(250, 42000),, labels = label_number(big.mark = "\u202F")) +
  scale_colour_manual(values=c("0.002"="#9e9e9e","0.02"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression(bold(N[crit]))) +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        plot.margin=margin(.5,.5,.5,.5))

#ben del ratio
#p4
p4 <- ggplot() +
  geom_line(data = dat_ub_kim_ud2_highoffset,
            aes(x=Ub, y=bendelratio, col="2", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ub_kim_ud1_highoffset,
            aes(x=Ub, y=bendelratio, col="1", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ub_kim_ud2_lowoffset,
            aes(x=Ub, y=bendelratio, col="2", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_ub_kim_ud1_lowoffset,
            aes(x=Ub, y=bendelratio, col="1", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=0.002,y=0.853532),col = "#F05039",shape=42, size=7)+
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_x_continuous(trans=c("log10","reverse"),
                     name=NULL, breaks = c(0.0001, 0.001, 0.01, 0.1), labels = c("0.0001", "0.001","0.01","0.1") ) +
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_y_log10(limits=c(0.35,18.5), breaks=c(1,3,10), labels=c("1","3","10")) +
  scale_colour_manual(values=c("2"="#9e9e9e","1"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y="Drought:Meltdown ratio") +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))


#p5 #offset del ratio
p5 <- ggplot() +
  geom_line(data = dat_ub_kim_ud2_highoffset,
            aes(x=Ub, y=1.5e-5/abs(del_flux-1.5e-5), col="2", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ub_kim_ud1_highoffset,
            aes(x=Ub, y=1.5e-5/abs(del_flux-1.5e-5), col="1", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ub_kim_ud2_lowoffset,
            aes(x=Ub, y=1e-6/abs(del_flux-1e-6), col="2", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_ub_kim_ud1_lowoffset,
            aes(x=Ub, y=1e-6/abs(del_flux-1e-6), col="1", linetype="Low offset"), linewidth =line_width)+
  scale_x_continuous(trans=c("log10","reverse"),
                     name=NULL, breaks = c(0.0001, 0.001, 0.01, 0.1), labels = c("0.0001", "0.001","0.01","0.1") ) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(0.0038,1), breaks=c(0.01, 0.1, 1), labels=c("0.01","0.1","1")) +
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_colour_manual(values=c("2"="#9e9e9e","1"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression("Env. ch. to mut. deg. ratio at "*N[crit])) +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

#p6 #Ncrit
p6 <- ggplot() +
  geom_line(data = dat_ub_kim_ud2_highoffset,
            aes(x=Ub, y=Ncrit, col="2", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ub_kim_ud1_highoffset,
            aes(x=Ub, y=Ncrit, col="1", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_ub_kim_ud2_lowoffset,
            aes(x=Ub, y=Ncrit, col="2", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_ub_kim_ud1_lowoffset,
            aes(x=Ub, y=Ncrit, col="1", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=0.002,y=3666.41),col = "#F05039",shape=42, size=7)+
  scale_x_continuous(trans=c("log10","reverse"),name=expression(bold(U[b])),
                     breaks = c(0.0001, 0.001, 0.01, 0.1), labels = c("0.0001", "0.001","0.01","0.1") ) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(250, 42000),, labels = label_number(big.mark = "\u202F")) +
  scale_colour_manual(values=c("2"="#9e9e9e","1"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression(N[crit])) +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

# Combine plots in 2 columns x 3 rows
combined <- (p1 | p4) / (p2 | p5) / (p3 | p6) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family="Times New Roman", face="italic", size=9))

# Save as Proceedings B single-column figure (88 mm = 3.46 inches)
tiff("figure5_combined_1col.tiff",
     width = 88/25.4, height = 140/25.4,  # adjust height for readability
     units = "in", res = 600, compression = "lzw")

print(combined)
dev.off()

########################################Figure 7: Reduction in Ne##########################################
#change in Ub wiht Ud=2;sb=0.001
df_red_ub <- data.frame("Ub" = c(0.001, 0.001584893, 0.002511886, 0.003981072, 0.006309573, 0.01),
                        "N_sim_fit" = c(6389.40, 5145.17, 4322.85, 3513.40, 2696.69, 2151.06),
                        "N_sim_sec" = c(6288, 5112, 4060, 3225, 2638, 2100),
                        "N_ana" = c(5052.58, 4083.34, 3298.96, 2663.88, 2149.32, 1731.97),
                        "N_e_coal" = c(5212.694265627656, 4325.696819001546, 3604.4934791060996, 
                                       2966.6656237152997, 2302.653578295538, 1865.6689174206622))

custom_breaks <- c(0.001, 0.005, 0.01)  # Custom tick positions
custom_labels <- c("0.001", "0.005", "0.01")  # Custom tick labels
Ub <- c(0.001, 0.001584893, 0.002511886, 0.003981072, 0.006309573, 0.01)
Ne_N2_ub <- rep(0.8995369,6)

#Fig 7A
p7a <- ggplot() +
  geom_point(data = df_red_ub, aes(x=Ub, y=(N_ana/N_sim_fit), col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), size=1) +
  geom_point(data = df_red_ub, aes(x=Ub, y=(N_e_coal/N_sim_fit), col="Coalescent-Ne", shape="Coalescent-Ne"), size=1) +
  geom_point(aes(x=Ub, y=Ne_N2_ub, col="Del only model", shape="Del only model"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels) +
  labs(x=bquote(bold(U[b])), y=expression(bold("Factor reduction in "~N[e]))) +
  theme_Publication(base_size=9) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del only model"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del only model"="red")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0.2,1,by=0.2),
                     labels=c("0.2", "0.4", "0.6","0.8", "1")) +
  theme(legend.position="none",
        plot.margin=margin(.5,.5,.5,.5))

#change in Ud with Ub=0.02;sb=0.001
df_red_ud_2 <- data.frame("Ud" = c(0.2511886, 0.5011872, 1.0, 1.9952623, 3.9810717),
                          "N_sim_fit" = c(471.78, 695.96, 953.70, 1500, 2488.79),
                          "N_sim_sec" = c(390, 591, 954, 1500, 2396),
                          "N_ana" = c(432.022, 627.141, 889.883, 1245.12, 1728.11),
                          "N_e_coal" = c(471.0206607437395, 682.0131766892991, 901.7554276638273,
                                         1320.4766134408887, 1847.292478102416))
####factor reduction in Ne  based on Del only model
s=0.00948704
Ud=c(0.2511886, 0.5011872, 1.0, 1.9952623, 3.9810717)
L=100
chr=23
Udw=Ud/(L*chr)
Udw
Rw=1/L
Rw
#Unlinked
Ne_N=exp(-4 * Ud * s)
Ne_N
#linked+unlined Joseph's equation
unlinked_term = exp(-4 * (Ud-Udw)*s)
linked_term = exp(-Udw / (2*s + Rw))
Ne_N2_Ud= unlinked_term * linked_term
Ne_N2_Ud

custom_breaks <- c(0.2, 1, 5)  # Custom tick positions
custom_labels <- c("0.2", "1", "5")  # Custom tick labels
#Fig7b
p7b <- ggplot() +
  geom_point(data = df_red_ud_2, aes(x=Ud, y=(N_ana/N_sim_fit), col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), size=1) +
  geom_point(data = df_red_ud_2, aes(x=Ud, y=(N_e_coal/N_sim_fit), col="Coalescent-Ne", shape="Coalescent-Ne"), size=1) +
  geom_point(aes(x=Ud, y=Ne_N2_Ud, col="Del. only", shape="Del. only"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels, limits=c(0.2,5)) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  labs(x=bquote(bold(U[d])), y=expression("Factor reduction in "~N[e])) +
  theme_Publication(base_size=9) +
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del. only"="red")) +
  ylim(c(0,1)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

#change in Sd with Ub=0.02;sb=0.001, Ud=2.0
df_red_sd_2 <- data.frame("Sd" = c(0.05, 0.03, 0.01, 0.005, 0.003),
                          "N_sim_fit" = c(1404.29, 1428.39, 1454.11, 1417.17, 1450),
                          "N_sim_sec" = c(1337, 1452, 1445, 1427, 1450),
                          "N_ana" = c(1094.59, 1140.09, 1241.6, 1306.43, 1352.7),
                          "N_e_coal" = c(1042.8969946062832, 1155.8836408279096, 1277.1174874421363,
                                         1327.5279447736898, 1340.582373437434))

####factor reduction in Ne based on Del only model
Sd=c(0.05, 0.03, 0.01, 0.005, 0.003)
Ud_sd=2
L=100
chr=23
Udw=Ud_sd/(L*chr)
Udw
Rw=1/L
Rw
#Unlinked
Ne_N=exp(-4 * Ud_sd * Sd)
Ne_N
#linked+unlined Joseph's equation
unlinked_term = exp(-4 * (Ud_sd-Udw)*Sd)
linked_term = exp(-Udw / (2*Sd + Rw))
Ne_N2_Sd= unlinked_term * linked_term
Ne_N2_Sd


custom_breaks <- c(0.001, 0.01, 0.05)  # Custom tick positions
custom_labels <- c("0.001", "0.01", "0.05")  # Custom tick labels

p7d <- ggplot() +
  geom_point(data = df_red_sd_2, aes(x=Sd, y=(N_ana/N_sim_fit), col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), size=1) +
  geom_point(data = df_red_sd_2, aes(x=Sd, y=(N_e_coal/N_sim_fit), col="Coalescent-Ne", shape="Coalescent-Ne"), size=1) +
  geom_point(aes(x=Sd, y=Ne_N2_Sd, col="Del. only", shape="Del. only"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  labs(x=bquote(bold(-bar(s)[d])), y=expression("Factor reduction in "~N[e])) +
  theme_Publication(base_size=9) +
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del. only"="red")) +
  ylim(c(0,1)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

#change in Sb with Ub=0.02;sd=0.009, Ud=2.0
df_red_sb_2 <- data.frame("Sb" = c(0.0005, 0.00075, 0.001, 0.0025, 0.005),
                          "N_sim_fit" = c(2939, 1770, 1462, 587, 288),
                          "N_sim_sec" = c(2738, 1799, 1504, 453, 290),
                          "N_ana" = c(2362.45, 1625.63, 1246.54, 533.659, 279.623),
                          "N_e_coal" = c(2487.70697, 1575.4146529, 1269.811607,
                                         534.860132, 264.64819387776504))

####factor reduction in Ne based on Del only model
Sb <- c(0.0005, 0.00075, 0.001, 0.0025, 0.005)
Ne_N2_sb <- rep(0.8995369,5)


custom_breaks <- c(0.0005, 0.005)  # Custom tick positions
custom_labels <- c("0.0005", "0.005")  # Custom tick labels

p7c <- ggplot() +
  geom_point(data = df_red_sb_2, aes(x=Sb, y=(N_ana/N_sim_fit), col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), size=1) +
  geom_point(data = df_red_sb_2, aes(x=Sb, y=(N_e_coal/N_sim_fit), col="Coalescent-Ne", shape="Coalescent-Ne"), size=1) +
  geom_point(aes(x=Sb, y=Ne_N2_sb, col="Del. only", shape="Del. only"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  labs(x=bquote(bold(bar(s)[b])), y=expression(bold("Factor reduction in "~N[e]))) +
  theme_Publication(base_size=9) +
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del. only"="red")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0.2,1,by=0.2),
                     labels=c("0.2", "0.4", "0.6","0.8", "1")) +
  theme(legend.position="none",
        plot.margin=margin(.5,.5,.5,.5))

fig7 <- (p7a | p7b) / (p7c | p7d) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family="Times New Roman", face="italic", size=9))

# Save for single-column (≈83 mm)
tiff("figure7_combined_1col.tiff",
     width=83/25.4, height=100/25.4, units="in",
     res=600, compression="lzw")

print(fig7)
dev.off()

####################################Supplementary Figures####################################################
#############################Figure S1: calculating Ncrit for second verison of W(2000)###################################
#Ncrit to Ud/Ub for three models Kim et al. (2017)
dat_ubud_ncrit_compare_kim <- read.csv("Ncrit_diff_models.csv", header = FALSE)
colnames(dat_ubud_ncrit_compare_kim) <- c("Ub_Ud", "Ncrit_ana", "Ncrit_approx", "Ncrit_full")

tiff("ncrit_analtyical_compare.tiff", units="in", width=3, height=3, res = 600, compression = "lzw")

p <- ggplot(data = dat_ubud_ncrit_compare_kim, aes(x = Ub_Ud))

p <- p + geom_line(aes(y = Ncrit_ana,    color = "This paper"), linewidth = .7)
p <- p + geom_line(aes(y = Ncrit_full,   color = "Whitlock (2000) Full"), linewidth = .7)
p <- p + geom_line(aes(y = Ncrit_approx, color = "Whitlock (2000) Approx"), linewidth = .7)

p <- p + scale_x_log10(
  breaks = c(10, 100, 1000, 10000),
  labels = label_number(big.mark = "\u202F"),
  name = expression(U[d] ~ "/" ~ U[b])
) +
  annotation_logticks(sides = "b",
                      short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside = TRUE) +
  coord_cartesian(clip = "off")
p <- p + geom_point(aes(x = 1000, y = 3666.41), color = "#F05039", shape = 42, size = 5)

p <- p + scale_colour_manual(
  name = '',
  values = c(
    "This paper" = "#000000",
    "Whitlock (2000) Full" = "#7fc97f",    # Green
    "Whitlock (2000) Approx" = "#386cb0"   # Blue
  ),
  breaks = c("This paper", "Whitlock (2000) Full", "Whitlock (2000) Approx") # Enforce order
)
p <- p + scale_y_log10()
p <- p + labs(y = expression(bold(N[crit])))
p <- p + theme_Publication(base_size = 9)
p <- p + theme(legend.position = "none") # Adjust legend position inside plot if desired

# 10. Print and Save
print(p)
dev.off()

#######################################Figure S2: Ncrit as a function of xb###############################################
#Ncrit
# Top-left plot
p1 <- ggplot(dat_ubud_kim, aes(x = Ub_Ud, y = Ncrit)) +
  geom_line(linewidth = line_width) +
  geom_point(aes(x=1000,y=3666.41), col = "#F05039", shape = 42, size = 5) +
  scale_x_log10(breaks = c(10,100,1000,10000), labels = c("10","100","1000","10000")) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  ylim(0,10623.02) +
  labs(x = expression(bold(U[d]/U[b])), y = expression(bold(N[crit]))) +
  theme_Publication(base_size = 9) +
  theme(plot.margin = margin(5,5,5,5), legend.position = "none")

# Top-right plot
p2 <- ggplot() +
  geom_line(data = dat_sb_kim_100, aes(x = sb, y = Ncrit, linetype = "mut. rate ratio = 100"), linewidth = line_width) +
  geom_line(data = dat_sb_kim_1000, aes(x = sb, y = Ncrit, linetype = "mut. rate ratio = 1000"), linewidth = line_width) +
  geom_point(aes(x=0.001,y=3666.41), col = "#F05039", shape = 42, size = 5) +
  scale_x_log10(breaks = c(0.001,0.005,0.01), labels = c("0.001","0.005","0.01")) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_linetype_manual(values = c("mut. rate ratio = 100" = 2, "mut. rate ratio = 1000" = 1)) +
  ylim(0,10623.02) +
  labs(x = bquote(bold(bar(s)[b])), y = expression(N[crit])) +
  theme_Publication(base_size = 9) +
  theme(plot.margin = margin(5,5,5,5),axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(), legend.position = "none")

# Bottom centered plot
p3 <- ggplot() +
  geom_line(data = dat_sd_100, aes(x = sd, y = Ncrit, linetype = "mut. rate ratio = 100"), linewidth = line_width) +
  geom_line(data = dat_sd_1000, aes(x = sd, y = Ncrit, linetype = "mut. rate ratio = 1000"), linewidth = line_width) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5, col = "red", linewidth = line_width) +
  scale_x_log10(limits=c(0.001,0.1), breaks = c(0.001, 0.01, 0.1), labels = c("0.001","0.01","0.1")) +
  geom_point(aes(x=0.009487,y=3666.41), col = "#F05039", shape = 42, size = 5) +
  scale_linetype_manual(values = c("mut. rate ratio = 100" = 2, "mut. rate ratio = 1000" = 1)) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  ylim(0,10623.02) +
  labs(x = bquote(bold(-bar(s)[d])), y = expression(N[crit])) +
  theme_Publication(base_size = 9) +
  theme(plot.margin = margin(5,5,5,5),axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(), legend.position = "none")

# Combine with patchwork
combined <- (p1 | p2 | p3) +   # stack bottom plot
  plot_layout(heights = c(1, 0.7)) + # bottom plot slightly shorter
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "italic", size = 9))

# Save 2-column Proceedings B figure (2-col = 167 mm)
tiff("figure_S2_combined_2col.tiff",
     width = 167/25.4, height = 80/25.4, units = "in",
     res = 600, compression = "lzw")

print(combined)
dev.off()

#############################Figure S3: calculating dorught:meltdown ratio for changing ud and ub with env change###################################
#with offset varying Sd fixed Ud
dat_sd_highoffset_1000 <- read.csv("offset_sdvarying_1000.csv", header = F)
colnames(dat_sd_highoffset_1000) <- c("sd", "del_flux","Ncrit", "bendelratio")

dat_sd_lowoffset_1000 <- read.csv("lowoffset_sdvarying_1000.csv", header = F)
colnames(dat_sd_lowoffset_1000) <- c("sd", "del_flux","Ncrit", "bendelratio")

dat_sd_highoffset_100 <- read.csv("offset_sdvarying_100.csv", header = F)
colnames(dat_sd_highoffset_100) <- c("sd", "del_flux","Ncrit", "bendelratio")

dat_sd_lowoffset_100 <- read.csv("lowoffset_sdvarying_100.csv", header = F)
colnames(dat_sd_lowoffset_100) <- c("sd", "del_flux","Ncrit", "bendelratio")

#with offset varying sb fixed Ud
dat_sb_highoffset_1000 <- read.csv("offset_sbvarying_1000.csv", header = F)
colnames(dat_sb_highoffset_1000) <- c("sb", "del_flux","Ncrit", "bendelratio")

dat_sb_lowoffset_1000 <- read.csv("lowoffset_sbvarying_1000.csv", header = F)
colnames(dat_sb_lowoffset_1000) <- c("sb", "del_flux","Ncrit", "bendelratio")

dat_sb_highoffset_100 <- read.csv("offset_sbvarying_100.csv", header = F)
colnames(dat_sb_highoffset_100) <- c("sb", "del_flux","Ncrit", "bendelratio")

dat_sb_lowoffset_100 <- read.csv("lowoffset_sbvarying_100.csv", header = F)
colnames(dat_sb_lowoffset_100) <- c("sb", "del_flux","Ncrit", "bendelratio")

# Add columns for offset type, Ub, and sd_type
dat_sd_highoffset_1000$offset <- "Faster env. change";dat_sd_highoffset_1000$UdUb <- "1000";
dat_sd_highoffset_100$offset <- "Faster env. change";dat_sd_highoffset_100$UdUb <- "100";
dat_sd_lowoffset_1000$offset <- "Slower env. change";dat_sd_lowoffset_1000$UdUb <- "1000";
dat_sd_lowoffset_100$offset <- "Slower env. change";dat_sd_lowoffset_100$UdUb <- "100";

# Add columns for offset type, Ub, and sd_type
dat_sb_highoffset_1000$offset <- "Faster env. change";dat_sb_highoffset_1000$UdUb <- "1000";
dat_sb_highoffset_100$offset <- "Faster env. change";dat_sb_highoffset_100$UdUb <- "100";
dat_sb_lowoffset_1000$offset <- "Slower env. change";dat_sb_lowoffset_1000$UdUb <- "1000";
dat_sb_lowoffset_100$offset <- "Slower env. change";dat_sb_lowoffset_100$UdUb <- "100";

# Combine datasets
dat_all_sd <- bind_rows(
  dat_sd_highoffset_1000[-c(19,20,21),],
  dat_sd_highoffset_100[-c(19,20,21),],
  dat_sd_lowoffset_1000[-c(19,20,21),],
  dat_sd_lowoffset_100[-c(19,20,21),]
)

dat_all_sb <- bind_rows(
  dat_sb_highoffset_1000,
  dat_sb_highoffset_100,
  dat_sb_lowoffset_1000,
  dat_sb_lowoffset_100
)

#ben del ratio
#p1
p1 <- ggplot() +
  geom_line(data = dat_sd_highoffset_1000[-c(19,20,21),],
            aes(x=sd, y=bendelratio, col="1000", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sd_highoffset_100[-c(19,20,21),],
            aes(x=sd, y=bendelratio, col="100", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sd_lowoffset_1000[-c(19,20,21),],
            aes(x=sd, y=bendelratio, col="1000", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_sd_lowoffset_100[-c(19,20,21),],
            aes(x=sd, y=bendelratio, col="100", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=0.009487042205869916,y=0.853532),col = "#F05039",shape=42, size=7)+
  scale_x_log10(name=NULL, limits=c(0.001,0.1), breaks = c(0.001, 0.01, 0.1), labels=c("0.001", "0.01", "0.1")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_y_log10(limits=c(0.35,18.5), breaks=c(1,3,10), labels=c("1","3","10")) +
  scale_colour_manual(values=c("1000"="#9e9e9e","100"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y="Drought:Meltdown ratio") +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))


#p2 #offset del ratio
p2 <- ggplot() +
  geom_line(data = dat_sd_highoffset_1000[-c(19,20,21),],
            aes(x=sd, y=1.5e-5/abs(del_flux-1.5e-5), col="1000", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sd_highoffset_100[-c(19,20,21),],
            aes(x=sd, y=1.5e-5/abs(del_flux-1.5e-5), col="100", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sd_lowoffset_1000[-c(19,20,21),],
            aes(x=sd, y=1e-6/abs(del_flux-1e-6), col="1000", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_sd_lowoffset_100[-c(19,20,21),],
            aes(x=sd, y=1e-6/abs(del_flux-1e-6), col="100", linetype="Low offset"), linewidth =line_width)+
  scale_x_log10(name=NULL, limits=c(0.001,0.1), breaks = c(0.001, 0.01, 0.1), labels=c("0.001", "0.01", "0.1")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(0.0038,1), breaks=c(0.01, 0.1, 1), labels=c("0.01","0.1","1")) +
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_colour_manual(values=c("1000"="#9e9e9e","100"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression(bold("Env.: Mut. deg. ratio at "*N[crit])))+
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=8),
        plot.margin=margin(.5,.5,.5,.5))

#p3 #Ncrit
p3 <- ggplot() +
  geom_line(data = dat_sd_highoffset_1000[-c(19,20,21),],
            aes(x=sd, y=Ncrit, col="1000", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sd_highoffset_100[-c(19,20,21),],
            aes(x=sd, y=Ncrit, col="100", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sd_lowoffset_1000[-c(19,20,21),],
            aes(x=sd, y=Ncrit, col="1000", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_sd_lowoffset_100[-c(19,20,21),],
            aes(x=sd, y=Ncrit, col="100", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=0.009487042205869916,y=3666.41),col = "#F05039",shape=42, size=7)+
  scale_x_log10(name=expression(bold(-bar(s[d]))),  limits=c(0.001,0.1), breaks = c(0.001, 0.01, 0.1),
                labels=c("0.001", "0.01", "0.1")) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(140, 42000), labels = label_number(big.mark = "\u202F"))  +
  scale_colour_manual(values=c("1000"="#9e9e9e","100"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression(bold(N[crit]))) +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        plot.margin=margin(.5,.5,.5,.5))

#ben del ratio
#p4
p4 <- ggplot() +
  geom_line(data = dat_sb_highoffset_1000,
            aes(x=sb, y=bendelratio, col="1000", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sb_highoffset_100,
            aes(x=sb, y=bendelratio, col="100", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sb_lowoffset_1000,
            aes(x=sb, y=bendelratio, col="1000", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_sb_lowoffset_100,
            aes(x=sb, y=bendelratio, col="100", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=0.001,y=0.853449),col = "#F05039",shape=42, size=7)+
  scale_x_log10(name=NULL, breaks = c(0.001, 0.01), labels = c("0.001", "0.01")) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_y_log10(limits=c(0.35,18.5), breaks=c(1,3,10), labels=c("1","3","10")) +
  scale_colour_manual(values=c("1000"="#9e9e9e","100"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y="Drought:Meltdown ratio") +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))


#p5 #offset del ratio
p5 <- ggplot() +
  geom_line(data = dat_sb_highoffset_1000,
            aes(x=sb, y=1.5e-5/abs(del_flux-1.5e-5), col="1000", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sb_highoffset_100,
            aes(x=sb, y=1.5e-5/abs(del_flux-1.5e-5), col="100", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sb_lowoffset_1000,
            aes(x=sb, y=1e-6/abs(del_flux-1e-6), col="1000", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_sb_lowoffset_100,
            aes(x=sb, y=1e-6/abs(del_flux-1e-6), col="100", linetype="Low offset"), linewidth =line_width)+
  scale_x_log10(breaks = c(0.001, 0.01), labels = c("0.001", "0.01")) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(0.0038,1), breaks=c(0.01, 0.1, 1), labels=c("0.01","0.1","1")) +
  geom_hline(aes(yintercept =1), linetype="dotted", alpha=0.5, col="red", linewidth=line_width)+
  scale_colour_manual(values=c("1000"="#9e9e9e","100"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression("Env. ch. to mut. deg. ratio at "*N[crit])) +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.x=element_blank(),  # remove x-axis for top row
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

#p6 #Ncrit
p6 <- ggplot() +
  geom_line(data = dat_sb_highoffset_1000,
            aes(x=sb, y=Ncrit, col="1000", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sb_highoffset_100,
            aes(x=sb, y=Ncrit, col="100", linetype="High offset"), linewidth=line_width) +
  geom_line(data = dat_sb_lowoffset_1000,
            aes(x=sb, y=Ncrit, col="1000", linetype="Low offset"), linewidth =line_width)+
  geom_line(data = dat_sb_lowoffset_100,
            aes(x=sb, y=Ncrit, col="100", linetype="Low offset"), linewidth =line_width)+
  geom_point(aes(x=0.001,y=3666.41),col = "#F05039",shape=42, size=7)+
  scale_x_log10(name=expression(bold(s[b])), breaks = c(0.001, 0.01), labels = c("0.001", "0.01")) +
  annotation_logticks(sides = "bl",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_y_log10(limits=c(140, 42000), labels = label_number(big.mark = "\u202F")) +
  scale_colour_manual(values=c("1000"="#9e9e9e","100"="#2298e6")) +
  scale_linetype_manual(values=c("High offset"="solid", "Low offset"="dashed")) +
  labs(y=expression(N[crit])) +
  theme_Publication(base_size=text_size) +
  theme(legend.position="none",
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

# Combine plots in 2 columns x 3 rows
combined <- (p1 | p4) / (p2 | p5) / (p3 | p6) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family="Times New Roman", face="italic", size=9))

# Save as Proceedings B single-column figure (88 mm = 3.46 inches)
tiff("figureS3_combined_1col.tiff",
     width = 4, height = 6,  # adjust height for readability
     units = "in", res = 600, compression = "lzw")

print(combined)
dev.off()

##################################Figure S4: Change in L#################################
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
  Udw=Ud/(L*chr)
  Udw
  Rw=1/L
  Rw
  #Unlinked
  Ne_N=exp(-4 * Ud * s)
  Ne_N
  #linked+unlined Joseph's equation
  unlinked_term = exp(-4 * (Ud-Udw)*s)
  linked_term = exp(-Udw / (2*s + Rw))
  Ne_N2= unlinked_term * linked_term
  Ne_N2
  print(paste("factor reducation in Ne - BGS ", Ne_N2))
  
  Ne_N_sim <- 4656.43/f1(0)
  print(paste("factor reducation in Ne - fit flux ", Ne_N_sim))
  
  
  Udw=Ud/(L*chr)
  Udw
  Rw=1/L
  Rw
  #Unlinked
  Ne_N=exp(-4 * Ud * s)
  Ne_N
  #linked+unlined Joseph's equation
  unlinked_term = exp(-4 * (Ud-Udw)*s)
  linked_term = exp(-Udw / (2*s + Rw))
  Ne_N2_Ud= unlinked_term * linked_term
  Ne_N2_Ud
  
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
                          "ratio" = c(1.117898, 1.011884, 0.9634023, 0.9085198,  0.902947, 0.9188318, 0.859119))

#plotting fluxes
p1 <- ggplot(change_in_L, aes(x = nchr, y = Ncrit)) +
  geom_point(size=2) +
  geom_hline(yintercept = 3666.41, col="red", linetype="dotted", linewidth = line_width) +
  labs(x = "Number of Chromosomes", y = expression(bold(N[crit]))) +
  theme_Publication(base_size = 9) +
  theme(plot.margin = margin(5,5,5,5))

p2 <- ggplot(change_in_L, aes(x = nchr, y = ratio)) +
  geom_point(size=2) +
  geom_hline(yintercept = 0.85, col="red", linetype="dotted", linewidth = line_width) +
  ylim(0, 1.13) +
  labs(x = "Number of Chromosomes", y = "Drought:Meltdown ratio") +
  theme_Publication(base_size = 9) +
  theme(plot.margin = margin(5,5,5,5))

# Combine plots in 2 columns x 3 rows
combined <- (p1 / p2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family="Times New Roman", face="italic", size=9))

# Save as Proceedings B single-column figure (88 mm = 3.46 inches)
tiff("figureS4_combined_1col.tiff",
     width = 3, height = 5,  # adjust height for readability
     units = "in", res = 600, compression = "lzw")

print(combined)
dev.off()

###################################Extra plot Ratio vs Ncrit for reviewer########################################
# Varying Ub/Ud (Single curve)
d4_ratio <- dat_ubud_kim %>% 
  select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ub/Ud Ratio")

# Varying Sb (Split by Ratio)
d4_sb_100 <- dat_sb_kim_100 %>% 
  select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Sb (Ud/Ub=100)")

d4_sb_1000 <- dat_sb_kim_1000 %>% 
  select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Sb (Ud/Ub=1000)")

# Combine Static (Excluding Sd)
df_static <- rbind(d4_ratio, d4_sb_100, d4_sb_1000)

# Hex Codes
grey_col <- "#9B9B9B"
blue_col <- "#1F97E5"


# --- Mappings for p_env ---
# Keys must match the 'LineGroup' text exactly
env_colors <- c(
  "Varying Ub (Fixed Ud=2)"     = grey_col,
  "Varying Ub (Fixed Ud=1)"     = blue_col,
  "Varying Ud (Fixed Ub=0.002)" = grey_col,
  "Varying Ud (Fixed Ub=0.02)"  = blue_col
)

env_linetypes <- c(
  "Varying Ub (Fixed Ud=2)"     = "solid",
  "Varying Ub (Fixed Ud=1)"     = "solid",
  "Varying Ud (Fixed Ub=0.002)" = "dashed",
  "Varying Ud (Fixed Ub=0.02)"  = "dashed"
)

# --- Mappings for p_traj ---
traj_colors <- c(
  "Ud/Ub = 1000" = grey_col,
  "Ud/Ub = 100" = blue_col,
  "Ud/Ub = 10" = "black"
)
# --- 2. Process Environmental Data (Figure 5) ---
# Goal: Distinguish Fast/Slow Env AND the fixed parameter (Ub or Ud)

# Group A: Varying Ud (Fixed Ub)
d5_ud_fast_lowUb <- dat_ud_kim_ub0.002_highoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ud (Fixed Ub=0.002)", Speed = "Fast Env Change")

d5_ud_fast_highUb <- dat_ud_kim_ub0.02_highoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ud (Fixed Ub=0.02)", Speed = "Fast Env Change")

d5_ud_slow_lowUb <- dat_ud_kim_ub0.002_lowoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ud (Fixed Ub=0.002)", Speed = "Slow Env Change")

d5_ud_slow_highUb <- dat_ud_kim_ub0.02_lowoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ud (Fixed Ub=0.02)", Speed = "Slow Env Change")

# Group B: Varying Ub (Fixed Ud)
d5_ub_fast_lowUd <- dat_ub_kim_ud1_highoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ub (Fixed Ud=1)", Speed = "Fast Env Change")

d5_ub_fast_highUd <- dat_ub_kim_ud2_highoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ub (Fixed Ud=2)", Speed = "Fast Env Change")

d5_ub_slow_lowUd <- dat_ub_kim_ud1_lowoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ub (Fixed Ud=1)", Speed = "Slow Env Change")

d5_ub_slow_highUd <- dat_ub_kim_ud2_lowoffset %>% select(Ncrit, bendelratio) %>% 
  mutate(LineGroup = "Varying Ub (Fixed Ud=2)", Speed = "Slow Env Change")

# Combine Env
df_env <- rbind(
  d5_ud_fast_lowUb, d5_ud_fast_highUb, d5_ud_slow_lowUb, d5_ud_slow_highUb,
  d5_ub_fast_lowUd, d5_ub_fast_highUd, d5_ub_slow_lowUd, d5_ub_slow_highUd
)

#add Ud/Ub=10 data
dat_env_offset_kim_10 <- read.csv("offset_kim_raw_10.csv", header = F)
colnames(dat_env_offset_kim_10) <- c("offset", "del_flux","Ncrit", "bendelratio")
dat_env_offset_kim_10$UbUd <- rep(10,dim(dat_env_offset_kim_10)[1])

# Combine the 100 and 1000 datasets if not already done
dat_env_offset_kim_10$Group <- "Ud/Ub = 10"
dat_env_offset_kim_100$Group <- "Ud/Ub = 100"
dat_env_offset_kim_1000$Group <- "Ud/Ub = 1000"

df_traj <- rbind(dat_env_offset_kim_10, dat_env_offset_kim_100, dat_env_offset_kim_1000)

# Filter out any weird outliers if necessary, but usually raw data is fine
# We select just what we need
df_traj_plot <- df_traj %>% 
  select(Ncrit, bendelratio, Group, offset) %>%
  arrange(Group, offset) # Ensure ordered by offset for geom_path

# --- Theme Setup (Increased font sizes) ---
theme_S4 <- theme_Publication(base_size = 11) + 
  theme(
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "none"
  )

y_breaks <- c(1, 10^0.5, 10)
y_labels <- c(expression(10^0), expression(10^0.5), expression(10^1))

# --- PANEL A: Ratio vs Ncrit (The "Big Picture") ---
p_a <- ggplot(df_traj_plot, aes(x = Ncrit, y = bendelratio, color = Group)) +
  geom_line(linewidth = 1.2) + 
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.5) +
  scale_y_log10(limits = c(0.5, 10), breaks = y_breaks, labels = y_labels) +
  scale_x_log10(labels = label_number(big.mark = "\u202F")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_color_manual(values = traj_colors) +
  labs(x = expression(bold(N[crit]~from~env.~change~rate)), 
       y = "Drought : Meltdown Ratio") +
  theme_S4

df_ub_only <- df_env %>% filter(grepl("Varying Ub", LineGroup))

# --- PANEL B: Detailed Mutation Rate Scenarios ---
p_b <- ggplot(df_ub_only, aes(x = Ncrit, y = bendelratio, color = LineGroup)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.5) +
  facet_wrap(~Speed) +
  scale_y_log10(limits = c(0.4, 10), breaks = y_breaks, labels = y_labels) +
  scale_x_log10(labels = label_number(big.mark = "\u202F")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_color_manual(values = env_colors) +
  labs(x = expression(bold(N[crit]~from~varying~U[b])), y = NULL) +
  theme_S4 +
  theme()

# 3. PANEL C: Ncrit from varying Ud (Fixed Ub)
df_ud_only <- df_env %>% filter(grepl("Varying Ud", LineGroup))

p_c <- ggplot(df_ud_only, aes(x = Ncrit, y = bendelratio, color = LineGroup, linetype = LineGroup)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", alpha = 0.5) +
  facet_wrap(~Speed) +
  scale_y_log10(limits = c(0.4, 10), breaks = y_breaks, labels = y_labels) +
  scale_x_log10(labels = label_number(big.mark = "\u202F")) +
  annotation_logticks(sides = "bl",short = unit(0.5, "mm"),
                      mid = unit(0.75, "mm"),
                      long = unit(1.25, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_color_manual(values = env_colors) +
  scale_linetype_manual(values = env_linetypes) +
  labs(x = expression(bold(N[crit]~from~varying~U[d])), y = NULL) +
  theme_S4 +
  theme()

# --- Define the centering layout ---
# Spacers act as empty "columns" to push p_a to the middle
top_row <- plot_spacer() + p_a + plot_spacer() + 
  plot_layout(widths = c(0.5, 1, 0.5))

# Combine with bottom row
final_layout <- top_row / (p_b | p_c) + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "italic", size = 11))

# --- Save with a balanced aspect ratio ---
tiff("figureS4_centered_layout.tiff", 
     width = 200/25.4, height = 150/25.4, units = "in", 
     res = 600, compression = "lzw")
print(final_layout)
dev.off()

##############################Figure 7 revision###################################
#change in Ub wiht Ud=2;sb=0.001
df_red_ub <- data.frame("Ub" = c(0.001, 0.005, 0.01),
                        "N_sim_fit" = c(6288, 3005.78, 1934.77),
                        "N_sim_fit_se" = c(32.64, 26.78, 26.05),
                        "N_ana" = c(5052.58, 2395.76, 1731.97),
                        "N_e_coal" = c(5138.39, 2573.871, 1708.0679),
                        "fitfluxratio"=c(0.80377, 0.79761, 0.89655),
                        "coalratio"=c(0.81742, 0.85691, 0.88418),
                        "no_del_Ne"= c(5198.78, 2580.54, 1693.057))

custom_breaks <- c(0.001, 0.005, 0.01)  # Custom tick positions
custom_labels <- c("0.001", "0.005", "0.01")  # Custom tick labels
Ub <- c(0.001, 0.005, 0.01)
Ne_N2_ub <- rep(0.8995369,3)
#Fig 7A
p7a <- ggplot() +
  geom_point(data = df_red_ub, 
             aes(x=Ub, y=fitfluxratio, col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), 
             size=1) +
  geom_point(data = df_red_ub, 
             aes(x=Ub, y=coalratio, col="Coalescent-Ne", shape="Coalescent-Ne"), 
             size=1) +
  geom_point(data = df_red_ub, aes(x=Ub, y=no_del_Ne/N_sim_fit, col="Del only model", shape="Del only model"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels) +
  labs(x=bquote(bold(U[b])), y=expression(bold("Factor reduction in "~N[e]))) +
  theme_Publication(base_size=9) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del only model"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del only model"="red")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0.2,1,by=0.2),
                     labels=c("0.2", "0.4", "0.6","0.8", "1")) +
  theme(legend.position="none",
        plot.margin=margin(.5,.5,.5,.5))

#change in Ud with Ub=0.02;sb=0.001
df_red_ud_2 <- data.frame("Ud" = c(0.50, 1.0, 5.0),
                          "N_sim_fit" = c(2065.87, 2859.92, 9777.55),
                          "N_sim_fit_se" = c(87.16, 120.00, 43.89),
                          "N_ana" = c(1923.34, 2658.01, 5601.41),
                          "N_e_coal" = c(2010, 2657.655, 5822.071),
                          "fitfluxratio"=c(0.94656, 0.94621, 0.57299),
                          "coalratio"=c(0.98921, 0.94608, 0.59556),
                          "no_del_Ne"= c(2007.59, 2676.85, 5942.30))
####factor reduction in Ne  based on Del only model
s=0.00948704
Ud=c(0.50, 1.0, 5.0)
L=100
chr=23
Udw=Ud/(L*chr)
Udw
Rw=1/L
Rw
#Unlinked
Ne_N=exp(-4 * Ud * s)
Ne_N
#linked+unlined Joseph's equation
unlinked_term = exp(-4 * (Ud-Udw)*s)
linked_term = exp(-Udw / (2*s + Rw))
Ne_N2_Ud= unlinked_term * linked_term
Ne_N2_Ud

custom_breaks <- c(0.5, 1, 5)  # Custom tick positions
custom_labels <- c("0.5", "1", "5")  # Custom tick labels
#Fig7b
p7b <- ggplot() +
  geom_point(data = df_red_ud_2, 
             aes(x=Ud, y=fitfluxratio, col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), 
             size=1) +
  geom_point(data = df_red_ud_2, 
             aes(x=Ud, y=coalratio, col="Coalescent-Ne", shape="Coalescent-Ne"), 
             size=1) +
  geom_point(data = df_red_ud_2, aes(x=Ud, y=no_del_Ne/N_sim_fit, col="Del. only", shape="Del. only"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels, limits=c(0.2,5)) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  labs(x=bquote(bold(U[d])), y=expression("Factor reduction in "~N[e])) +
  theme_Publication(base_size=9) +
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del. only"="red")) +
  ylim(c(0,1)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

#change in Sd with Ub=0.02;sb=0.001, Ud=2.0
df_red_sd_2 <- data.frame("Sd" = c(0.005, 0.01, 0.05),
                          "N_sim_fit" = c(4655.56, 4434.29, 4433.55),
                          "N_sim_fit_se" = c(45.54, 160.13, 187.85),
                          "N_ana" = c(3848.67, 3651.69, 3224.47),
                          "N_e_coal" = c(4000.036, 3758.098, 3191.463),
                          "fitfluxratio"=c(0.82742, 0.83689, 0.74242),
                          "coalratio"=c(0.85996, 0.86128, 0.73482),
                          "no_del_Ne"= c(4027.83, 3762.34, 3217.62))

####factor reduction in Ne based on Del only model
Sd=c(0.005, 0.01, 0.05)
Ud_sd=2
L=100
chr=23
Udw=Ud_sd/(L*chr)
Udw
Rw=1/L
Rw
#Unlinked
Ne_N=exp(-4 * Ud_sd * Sd)
Ne_N
#linked+unlined Joseph's equation
unlinked_term = exp(-4 * (Ud_sd-Udw)*Sd)
linked_term = exp(-Udw / (2*Sd + Rw))
Ne_N2_Sd= unlinked_term * linked_term
Ne_N2_Sd

custom_breaks <- c(0.005, 0.01, 0.05)  # Custom tick positions
custom_labels <- c("0.005", "0.01", "0.05")  # Custom tick labels

p7d <- ggplot() +
  geom_point(data = df_red_sd_2, 
             aes(x=Sd, y=fitfluxratio, col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), 
             size=1) +
  geom_point(data = df_red_sd_2, 
             aes(x=Sd, y=coalratio, col="Coalescent-Ne", shape="Coalescent-Ne"), 
             size=1) +
  geom_point(data = df_red_sd_2, aes(x=Sd, y=no_del_Ne/N_sim_fit, col="Del. only", shape="Del. only"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  labs(x=bquote(bold(-bar(s)[d])), y=expression("Factor reduction in "~N[e])) +
  theme_Publication(base_size=9) +
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del. only"="red")) +
  ylim(c(0,1)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),  # remove y-axis for top row
        axis.text.y=element_blank(),
        plot.margin=margin(.5,.5,.5,.5))

#change in Sb with Ub=0.02;sd=0.009, Ud=2.0
df_red_sb_2 <- data.frame("Sb" = c(0.0005, 0.001, 0.005),
                          "N_sim_fit" = c(9127, 4613.18, 933.35),
                          "N_sim_fit_se" = c(93.86, 31.00, 10.24),
                          "N_ana" = c(6947.83, 3666.41, 830.541),
                          "N_e_coal" = c(7179.021, 3866.158, 845.365),
                          "fitfluxratio"=c(0.76197, 0.79509, 0.89083),
                          "coalratio"=c(0.78733, 0.83841, 0.90673),
                          "no_del_Ne"= c(7357.56, 3861.46, 820.71))

####factor reduction in Ne based on Del only model
Sb <- c(0.0005, 0.001, 0.005)
Ne_N2_sb <- rep(0.8995369 ,3)

custom_breaks <- c(0.0005, 0.001, 0.005)  # Custom tick positions
custom_labels <- c("0.0005", "0.001", "0.005")  # Custom tick labels

p7c <- ggplot() +
  geom_point(data = df_red_sb_2, 
             aes(x=Sb, y=fitfluxratio, col="Fitness-flux-Ne", shape="Fitness-flux-Ne"), 
             size=1) +
  geom_point(data = df_red_sb_2, 
             aes(x=Sb, y=coalratio, col="Coalescent-Ne", shape="Coalescent-Ne"), 
             size=1) +
  geom_point(data = df_red_sb_2, aes(x=Sb, y=no_del_Ne/N_sim_fit, col="Del. only", shape="Del. only"), size=1) +
  scale_x_log10(breaks=custom_breaks, labels=custom_labels) +
  annotation_logticks(sides = "b",short = unit(0.25, "mm"),
                      mid = unit(0.5, "mm"),
                      long = unit(1, "mm"),
                      outside=TRUE)+ coord_cartesian(clip="off")+
  labs(x=bquote(bold(bar(s)[b])), y=expression(bold("Factor reduction in "~N[e]))) +
  theme_Publication(base_size=9) +
  scale_shape_manual(name='', values=c("Fitness-flux-Ne"=16,"Del. only"=15,"Coalescent-Ne"=17)) +
  scale_colour_manual(name=' ', values=c("Fitness-flux-Ne"="black","Coalescent-Ne"="#56b4e9","Del. only"="red")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0.2,1,by=0.2),
                     labels=c("0.2", "0.4", "0.6","0.8", "1")) +
  theme(legend.position="none",
        plot.margin=margin(.5,.5,.5,.5))

fig7 <- (p7a | p7b) / (p7c | p7d) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(family="Times New Roman", face="italic", size=9))

# Save for single-column (≈83 mm)
tiff("figure7_combined_1col.tiff",
     width=83/25.4, height=100/25.4, units="in",
     res=600, compression="lzw")

print(fig7)
dev.off()

################ratio sensitvity to CV############################################ 
line_width <- 1

# ==========================================
# 1. HELPER FUNCTION
# ==========================================
read_format_data <- function(filename, condition_name, env_status) {
  df <- read.csv(filename, header = F)
  colnames(df) <- c("cv", "del_flux", "Ncrit", "bendelratio")
  df$Condition <- condition_name
  df$Env <- env_status
  return(df)
}

# ==========================================
# 2. LOAD DATA
# ==========================================
# Base Case: Ratio 1000 (Ud = 2, Ub = 0.002)
df_base_noenv <- read_format_data("noenv_FixedUd_R1000.csv", "Ud=2.0, Ub=0.002 (Ratio 1000)", "No Change")
df_base_env   <- read_format_data("env_FixedUd_R1000.csv", "Ud=2.0, Ub=0.002 (Ratio 1000)", "Change")

# Variant A: Ratio 100 (Ud = 2, Ub = 0.02)
df_varA_noenv <- read_format_data("noenv_FixedUd_R100.csv", "Ud=2.0, Ub=0.02 (Ratio 100)", "No Change")
df_varA_env   <- read_format_data("env_FixedUd_R100.csv", "Ud=2.0, Ub=0.02 (Ratio 100)", "Change")

# Combine all into one dataframe
df_combined <- rbind(df_base_noenv, df_base_env, 
                     df_varA_noenv, df_varA_env)

# Order the factors for a logical legend
df_combined$Condition <- factor(df_combined$Condition, 
                                levels = c("Ud=2.0, Ub=0.002 (Ratio 1000)", 
                                           "Ud=2.0, Ub=0.02 (Ratio 100)"))
df_combined$Env <- factor(df_combined$Env, levels = c("Change", "No Change"))

# ==========================================
# 3. BUILD PLOT
# ==========================================
p_unified <- ggplot(df_combined, aes(x = cv, y = bendelratio, color = Condition, linetype = Env)) +
  geom_line(linewidth = line_width) +
  geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5, col = "red", linewidth = line_width) +
  
  # Red Star
  geom_point(aes(x = 2.432521, y = 0.853449), col = "#F05039", shape = 42, size = 7, inherit.aes = FALSE) +
  geom_point(aes(x = 2.432521, y = 1.43349), col = "darkorange", shape = 42, size = 7, inherit.aes = FALSE) +
  
  scale_x_log10() +
  annotation_logticks(sides = "b", short = unit(0.5, "mm"), mid = unit(0.75, "mm"), long = unit(1.25, "mm"), outside = TRUE) + 
  coord_cartesian(clip = "off") +
  
  # Colors: Grey for Base, Blue for Varying Ub, Orange for Varying Ud
  scale_color_manual(name = "Parameters", 
                     values = c("Ud=2.0, Ub=0.002 (Ratio 1000)" = "#9B9B9B", 
                                "Ud=2.0, Ub=0.02 (Ratio 100)" = "#1F97E5")) +
  scale_linetype_manual(name = "Environment", 
                        values = c("Change" = "solid", "No Change" = "dashed")) +
  
  scale_y_continuous(breaks = seq(0.5, 2, by = 0.5), limits = c(0.48, 2.000157)) + # Adjust limits based on your actual data range
  labs(x = bquote(bold(CV)), y = "Drought:Meltdown ratio") +
  
  theme_Publication() +
  theme(legend.position = "none")

# ==========================================
# 4. EXPORT
# ==========================================
tiff("figureS6_unified.tiff",
     width = 80/25.4, height = 80/25.4, units = "in", # Wider to fit the single legend
     res = 600, compression = "lzw")

print(p_unified)
dev.off()

