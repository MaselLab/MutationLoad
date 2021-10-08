rm(list=ls())

files1 = list.files("Dead100xover/Data/pre-fenwick text files")
files1 <- paste0("Dead100xover/Data/pre-fenwick text files/", files1)
data.for.plot1 = data.frame(param = c(0.00030625, 0.0006125, 0.00091875, 0.001225, 0.00125, 0.001725, 0.00225, 0.0025, 0.00375, 0.005), slope = NA)

files2 = list.files("Dead100xover/Data/post-fenwick text files")
files2 <- paste0("Dead100xover/Data/post-fenwick text files/", files2)
data.for.plot2 = data.frame(param = c(0.00030625, 0.0006125, 0.00091875, 0.001225, 0.00125, 0.001725, 0.00225, 0.0025, 0.00375, 0.005), slope = NA)

for (file.index in 1:length(files1)) {
  curr.data1 <- read.csv(file = files1[file.index], header = TRUE)
  curr.model1 <- lm(Sum.of.death.rates ~ Generation, data = curr.data1)
  #this is appropriate for very linear data, but may not always be appropriate.
  curr.summary1 <- summary(curr.model1) #may be ways to do without this.
  data.for.plot1[file.index, "slope"] <- (curr.summary1$coeff[2] / 1000) / data.for.plot1[file.index, "param"] #check that [2] is actually the slope.
  }

for (file.index in 1:length(files2)) {
  curr.data2 <- read.csv(file = files2[file.index], header = TRUE)
  curr.model2 <- lm(Sum.of.death.rates ~ Generation, data = curr.data2)
  #this is appropriate for very linear data, but may not always be appropriate.
  curr.summary2 <- summary(curr.model2) #may be ways to do without this.
  data.for.plot2[file.index, "slope"] <- (curr.summary2$coeff[2] / 1000) / data.for.plot2[file.index, "param"] #check that [2] is actually the slope.
}

xticks <- seq(0, 0.005, 0.001)
yticks <- seq(1.1, 1.6, 0.05)
plot(x = data.for.plot1$param, y = data.for.plot1$slope, col="red", ylim = c(1.1, 1.6), axes = FALSE)
points(x = data.for.plot2$param, y = data.for.plot2$slope, col="blue")
axis(2, at = yticks, labels = yticks)
axis(1, at = xticks, labels = xticks)
