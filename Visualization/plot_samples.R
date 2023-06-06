#!/usr/bin/env Rscript
#source("helpers/kdepairs.default.R")
#source("helpers/kdepairs.R")
source("helpers/kde_ggpairs.R")

library(tools)
library(bayestestR)

cat("Did you delete all the variables you wanted? If no, run 'rm(list=ls())'\n")

# ----------------------------------------------------------------------
# name of the data file
#fname <- "../CAL_KSC_AE_SET_1_RUN_3_IC_ID-22852920/curgen_db_014.txt"

args = commandArgs(trailingOnly=TRUE)
fname <- args[1]
#fres <-"res_tmcmc_FW.tex"

#fname <- "../Summary/DSR00_RUN_1_C33A_controls_Rall_obs_sigma/Calibration/final.txt"
#l = c(expression(r[g]),expression(K))

if(args[2] == "00"){
    l = c(expression(r[g]),expression(K))
}else if(args[2] == "00_treatment"){
    l = c(expression(r[g]),expression(K),expression(u[eff]))
}else if(args[2] == "00_treatment_expanded"){
    l = c(expression(r[g]),expression(K),expression(u[eff]))
}else if(args[2] == "00_treatment_expanded_inv" || args[2] == "00_treatment_expanded_sym"){
    l = c(expression(r[g]),expression(K),expression(u[eff]))
}else if(args[2] == "00_treatment_expanded_sym_no_gen"){
	l = c(expression(r[g]),expression(K),expression(u[eff]),expression(A),expression(B))
}else if(args[2] == "G00"){
    l = c(expression(r[g]),expression(K),expression(theta))
}else if(args[2] == "F00"){
    l = c(expression(r[g]),expression(K),expression(alpha))
}else if(args[2] == "F00_treatment" || args[2] == "F00_treatment_expanded_inv"|| args[2] == "F00_treatment_expanded_sym"){
    l = c(expression(r[g]),expression(K),expression(u[eff]),expression(alpha))
}else if(args[2] == "F01_treatment_expanded_inv" || args[2] == "F01_treatment_expanded"){
    l = c(expression(r[g]),expression(u[eff]),expression(alpha))
}else if(args[2] == "01"){
	l = c(expression(alpha),expression(phi),expression(mu[sr]),expression(mu[rs]),expression(r[1]),expression(r[2]))
}else if(args[2] == "02"){
    l = c(expression(r[s]),expression(r[r]),expression(K),expression(R[0]))
}else if(args[2] == "02_treatment"){
    l = c(expression(r[s]),expression(r[r]),expression(K),expression(alpha),expression(d[s]),expression(d[r]),expression(R[0]))
}else if(args[2] == "03"){
	l = c(expression(r[s]),expression(r[r]),expression(K),expression(epsilon),expression(gamma))
}else if(args[2] == "01_s"){
    l = c(expression(alpha),expression(phi),expression(mu[sr]),expression(mu[rs]),expression(r[1]),expression(r[2]),expression(sigma))
}else if(args[2] == "02_s"){
    l = c(expression(r[s]),expression(r[r]),expression(K),expression(sigma))
}else if(args[2] == "03_s"){
    l = c(expression(r[s]),expression(r[r]),expression(K),expression(epsilon),expression(gamma),expression(sigma))
}else {
    print("Something is wrong...")
}

paste(c("The model parameters are: ", l), collapse=" ")

# names of variables
#l = c(expression(D[u]),expression(s),expression(chi),expression(D[f]),expression(r), expression(sigma ^2))

discrepancy <- 0  # is the last column of the data discrepancy of likelihood?
truelik <- 1  # use likelihood or log-likelihood for coloring

# ----------------------------------------------------------------------
print_tab <- function(text, x, nd, w) {
  cat(text, formatC(x, digits=nd, format='e', width=w), "\n", sep="\t")
}

print_nd <- function(x, nd) {
  cat(format(round(x, digits=nd), nsmall=nd))
}

cat("Reading file", fname,"\n")
data <- read.table(fname)
data <- array(data=unlist(data), dim=dim(data))
nd <- dim(data)[2]  # dimension of the samples

# mean
md <- colMeans(data[, 1:nd-1])
print_tab("means:\t", md, 4, 8)

# standard deviation
m2d <- colMeans(data[, 1:nd-1]^2)
sd <- sqrt(m2d-md^2)
print_tab("stds:\t", sd, 4, 8)

# most probable parameters
best_id <- which.max(data[, nd-1])
best <- data[best_id, 1:nd-1]
print_tab("MPVs:\t", best, 4, 8)

# quantiles
q1 <- c(); for (i in seq(1, nd-1)) q1[i] <- quantile(data[, i], probs=c(0.05))
print_tab("q-0.05:  ", q1, 4, 8)
q2 <- c(); for (i in seq(1, nd-1)) q2[i] <- quantile(data[, i], probs=c(0.95))
print_tab("q-0.95:  ", q2, 4, 8)

# print table
fname <- paste0(file_path_sans_ext(fname), "_tmcmc_params.txt")
sink(fname)
print_tab("means:", md, 4, 8)
print_tab("stds:", sd, 4, 8)
print_tab("MPVs:", best, 4, 8)
print_tab("q+0.05:", q1, 4, 8)
print_tab("q-0.95:", q2, 4, 8)
sink()


dat <- data.frame(data)
colnames(dat) = l
ci_eti <- ci(dat,ci=0.95, method = "HDI")
fname1 <-paste0(file_path_sans_ext(fname), "_CI95.txt")
sink(fname1)
for(i in (3:nd-2))
  print_tab(ci_eti$Parameter[i], c(ci_eti$CI_low[i], ci_eti$CI_high[i]), 4,8)
sink()


# print for latex table
a <- c(); digits <- 3
for (i in seq(1, nd-1)) {
    t <- capture.output(cat(print_nd(md[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, " & ["), collapse = "")
    t <- capture.output(cat(print_nd(q1[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, ", "), collapse = "")
    t <- capture.output(cat(print_nd(q2[i], digits)))
    a <- paste(c(a, t), collapse = "")
    a <- paste(c(a, "] & "), collapse = "")
}

#write.table(a, file=fres,sep = "\t", row.names = T)
a <- capture.output(cat(substr(a, 1, nchar(a)-3)))
a <- paste("latex:", a)
a <- paste(a, "\\\\ \n")
cat(a)



if(discrepancy) data[, nd] <-    -data[, nd]
if(truelik)     data[, nd] <- exp(data[, nd])

#source("helpers/kde_ggpairs.R")
pname <- paste0(file_path_sans_ext(fname), ".png")
png(pname, width=2500, height=2500, units='px', res=300, pointsize=15, type="cairo", family="times")
kdeggpairs(dat, n_1d=20, n_2d=200, labels=l) #n_2d=200
dev.off()
