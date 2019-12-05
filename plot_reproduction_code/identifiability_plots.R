#Code for reproducing graphs involving Frequentist
#estimates of parameters from Gao et. al (2015) data

#Check working directory is the plot_reproduction_code folder
#change if necessary!
getwd()
#remotes::install_github("jumpingrivers/prettyB")

library(readxl)
library(prettyB)
dd = read_xlsx("../getting_started/WarwickModelGaoLOGFitting.xlsx")

#Histograms of k and rd
par(mfrow=c(1,2))
hist_p(log10(dd$K), main=expression("Histogram of k (log scale)"),
       xlab=expression(log[10](k)), ylim=c(0,85))
hist_p(log10(dd$Rd), main=expression("Histogram of" ~ r[d] ~ "(log scale)"),
       xlab=expression(log[10](r[d])), ylim=c(0,85))
par(mfrow=c(1,1))

dev.copy(postscript,"../TL_presentation/k_rd_histograms.eps")
dev.off()


#Plot of rd against k (log scale)
par(mfrow=c(1,2))
plot_p(dd$K, dd$Rd, log = "xy", main=expression(r[d] ~ "against" ~ k),
       xlab=expression(k), ylab=expression(r[d]),
       xaxt="n", yaxt="n",
       ylim=c(1e-15,1e5))
add_x_axis(c(1e-15,1e-8,1e-4,1e-0,1e+05), expression(10^-15,10^-10,10^-5,10^0,10^5))
add_y_axis(c(1e-15,1e-10,1e-05,1e-0,1e+05), expression(10^-15,10^-10,10^-5,10^0,10^5))

dev.copy(postscript,"../TL_presentation/rd_v_k.eps")
dev.off()

#plot of rd * k (log scale)
par(mfrow=c(1,2))
plot_p(dd$K, dd$Rd, log = "xy", main=expression(r[d] ~ "against" ~ k),
       xlab=expression(k), ylab=expression(r[d]),
       xaxt="n", yaxt="n",
       ylim=c(1e-15,1e5))
add_x_axis(c(1e-15,1e-8,1e-4,1e-0,1e+05), expression(10^-15,10^-10,10^-5,10^0,10^5))
add_y_axis(c(1e-15,1e-10,1e-05,1e-0,1e+05), expression(10^-15,10^-10,10^-5,10^0,10^5))
plot_p(dd$Rd*dd$K, log = "y", main=expression(r[d]%*%k),
     xlab="Tumour no.", ylab=expression(r[d]%*%k),
     yaxt = "n")

add_y_axis(c(1e-30,1e-22,1e-14,1e-6,1e2), labels=expression(10^-30,10^-22,10^-14,10^-6,10^2))

dev.copy(postscript,"../TL_presentation/rd_x_k.eps")
dev.off()
