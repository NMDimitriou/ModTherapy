## Correlate gene measurements with cell counts
using Revise
include("load_data.jl")
include("plotopt.jl")

#Load control data
cdat = load_data("../Data/C33A_control_R*_werr.txt");
#Load treatment data
x0_01dat = load_data("../Data/C33A_x0_01_R*_werr.txt");
x0_1dat  = load_data("../Data/C33A_x0_1_R*_werr.txt");
x1dat    = load_data("../Data/C33A_x1_R*_werr.txt");
x10dat   = load_data("../Data/C33A_x10_R*_werr.txt");
x100dat  = load_data("../Data/C33A_x100_R*_werr.txt");

## Take the average across replicates
import Statistics as stats
num_sets, num_tp, cidx, icded, cded, σc = prep_dat(cdat)
av_cded,δσc = mean_sem_dat(cded,σc,num_tp,num_sets)

num_sets, num_tp, x0_01cidx, icded, x0_01ded, σx0_01 = prep_dat(x0_01dat)
av_x0_01ded,δσx0_01 = mean_sem_dat(x0_01ded,σx0_01,num_tp,num_sets)

num_sets, num_tp, x0_1cidx, icded, x0_1ded, σx0_1 = prep_dat(x0_1dat)
av_x0_1ded,δσx0_1 = mean_sem_dat(x0_1ded,σx0_1,num_tp,num_sets)

num_sets, num_tp, x1cidx, icded, x1ded, σx1 = prep_dat(x1dat)
av_x1ded,δσx1 = mean_sem_dat(x1ded,σx1,num_tp,num_sets)

num_sets, num_tp, x10cidx, icded, x10ded, σx10 = prep_dat(x10dat)
av_x10ded,δσx10 = mean_sem_dat(x10ded,σx10,num_tp,num_sets)

num_sets, num_tp, x100cidx, icded, x100ded, σx100 = prep_dat(x100dat)
av_x100ded,δσx100 = mean_sem_dat(x100ded,σx100,num_tp,num_sets)

## Normalize treatment data with controls
normx0_01  = av_x0_01ded./av_cded
δnormx0_01 = δnorm(av_x0_01ded,av_cded,δσx0_01,δσc)

normx0_1  = av_x0_1ded./av_cded
δnormx0_1 = δnorm(av_x0_1ded,av_cded,δσx0_1,δσc)

normx1    = av_x1ded./av_cded
δnormx1 = δnorm(av_x1ded,av_cded,δσx1,δσc)

normx10   = av_x10ded./av_cded
δnormx10 = δnorm(av_x10ded,av_cded,δσx10,δσc)

normx100  = av_x100ded./av_cded
δnormx100 = δnorm(av_x100ded,av_cded,δσx100,δσc)

## Plot the normalizations
plot()
plot(icded,
     normx0_01,
     yerror = δnormx0_01,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "0.01μM",)
plot!(icded,
     normx0_1,
     yerror = δnormx0_1,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "0.1μM",)
plot!(icded,
      normx1,
      yerror = δnormx1,
      #st = :scatter,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "1μM",)
plot!(icded,
      normx10,
      yerror = δnormx10,
      #st = :scatter,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "10μM",)
plot!(icded,
      normx100,
      yerror = δnormx100,
      #st = :scatter,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "100μM",title=:C33A,
      margins=10mm)
xlabel!("time (hours)")
ylabel!("normalized cell count")
png("normalized_cell_counts_C33A.png")


## Load PCR of oct4
oct4_x0_01dat = load_data("../Data/PCR/C33A_OCT4_x0_01_R1_werr.txt");
oct4_x0_1dat  = load_data("../Data/PCR/C33A_OCT4_x0_1_R1_werr.txt");
oct4_x1dat    = load_data("../Data/PCR/C33A_OCT4_x1_R1_werr.txt");
oct4_x10dat   = load_data("../Data/PCR/C33A_OCT4_x10_R1_werr.txt");
oct4_x100dat  = load_data("../Data/PCR/C33A_OCT4_x100_R1_werr.txt");

# Load PCR of klf4
klf4_x0_01dat = load_data("../Data/PCR/C33A_KLF4_x0_01_R1_werr.txt");
klf4_x0_1dat  = load_data("../Data/PCR/C33A_KLF4_x0_1_R1_werr.txt");
klf4_x1dat    = load_data("../Data/PCR/C33A_KLF4_x1_R1_werr.txt");
klf4_x10dat   = load_data("../Data/PCR/C33A_KLF4_x10_R1_werr.txt");
klf4_x100dat  = load_data("../Data/PCR/C33A_KLF4_x100_R1_werr.txt");

# Load PCR of sox2
sox2_x0_01dat = load_data("../Data/PCR/C33A_SOX2_x0_01_R1_werr.txt");
sox2_x0_1dat  = load_data("../Data/PCR/C33A_SOX2_x0_1_R1_werr.txt");
sox2_x1dat    = load_data("../Data/PCR/C33A_SOX2_x1_R1_werr.txt");
sox2_x10dat   = load_data("../Data/PCR/C33A_SOX2_x10_R1_werr.txt");
sox2_x100dat  = load_data("../Data/PCR/C33A_SOX2_x100_R1_werr.txt");

# Load PCR of NANOG
nanog_x0_01dat = load_data("../Data/PCR/C33A_NANOG_x0_01_R1_werr.txt");
nanog_x0_1dat  = load_data("../Data/PCR/C33A_NANOG_x0_1_R1_werr.txt");
nanog_x1dat    = load_data("../Data/PCR/C33A_NANOG_x1_R1_werr.txt");
nanog_x10dat   = load_data("../Data/PCR/C33A_NANOG_x10_R1_werr.txt");
nanog_x100dat  = load_data("../Data/PCR/C33A_NANOG_x100_R1_werr.txt");

## Calculate correlations between PCR and cell counts
import Statistics as stats

same_tp =  findall(in(icded),oct4_x0_01dat.Column1)
corr_x0_01_oct4 = stats.cor(normx0_01[same_tp],oct4_x0_01dat.Column2)
corr_x0_1_oct4  = stats.cor(normx0_1[same_tp],oct4_x0_1dat.Column2)
corr_x1_oct4  = stats.cor(normx1[same_tp],oct4_x1dat.Column2)
corr_x10_oct4  = stats.cor(normx10[same_tp],oct4_x10dat.Column2)
corr_x100_oct4  = stats.cor(normx100[same_tp],oct4_x100dat.Column2)

corr_x0_01_klf4 = stats.cor(normx0_01[same_tp],klf4_x0_01dat.Column2)
corr_x0_1_klf4  = stats.cor(normx0_1[same_tp],klf4_x0_1dat.Column2)
corr_x1_klf4  = stats.cor(normx1[same_tp],klf4_x1dat.Column2)
corr_x10_klf4  = stats.cor(normx10[same_tp],klf4_x10dat.Column2)
corr_x100_klf4  = stats.cor(normx100[same_tp],klf4_x100dat.Column2)

corr_x0_01_sox2 = stats.cor(normx0_01[same_tp],sox2_x0_01dat.Column2)
corr_x0_1_sox2  = stats.cor(normx0_1[same_tp],sox2_x0_1dat.Column2)
corr_x1_sox2  = stats.cor(normx1[same_tp],sox2_x1dat.Column2)
corr_x10_sox2  = stats.cor(normx10[same_tp],sox2_x10dat.Column2)
corr_x100_sox2  = stats.cor(normx100[same_tp],sox2_x100dat.Column2)

corr_x0_01_nanog = stats.cor(normx0_01[same_tp],nanog_x0_01dat.Column2)
corr_x0_1_nanog  = stats.cor(normx0_1[same_tp],nanog_x0_1dat.Column2)
corr_x1_nanog  = stats.cor(normx1[same_tp],nanog_x1dat.Column2)
corr_x10_nanog  = stats.cor(normx10[same_tp],nanog_x10dat.Column2)
corr_x100_nanog  = stats.cor(normx100[same_tp],nanog_x100dat.Column2)

@printf("%.2f & %.2f & %.2f & %.2f & %.2f \n",0.01, corr_x0_01_oct4,corr_x0_01_klf4,corr_x0_01_sox2,corr_x0_01_nanog)
@printf("%.2f & %.2f & %.2f & %.2f & %.2f \n",0.1 ,corr_x0_1_oct4,corr_x0_1_klf4,corr_x0_1_sox2,corr_x0_1_nanog)
@printf("%.2f & %.2f & %.2f & %.2f & %.2f \n",1,corr_x1_oct4,corr_x1_klf4,corr_x1_sox2,corr_x1_nanog)
@printf("%.2f & %.2f & %.2f & %.2f & %.2f \n",10,corr_x10_oct4,corr_x10_klf4,corr_x10_sox2,corr_x10_nanog)
@printf("%.2f & %.2f & %.2f & %.2f & %.2f \n",100,corr_x100_oct4,corr_x100_klf4,corr_x100_sox2,corr_x100_nanog)


using StatsPlots
corr_count_matrix = hcat(normx0_01[same_tp],
                         normx0_1[same_tp],
                         normx1[same_tp],
                         normx10[same_tp],
                         normx100[same_tp])

corr_oct4_matrix = hcat(oct4_x0_01dat.Column2,
                         oct4_x0_1dat.Column2,
                         oct4_x1dat.Column2,
                         oct4_x10dat.Column2,
                         oct4_x100dat.Column2)

corr_matrix = hcat(corr_count_matrix,corr_oct4_matrix)

## Plots


plot()
cornerplot(corr_matrix,
           compact=true,
           label=["cells 0.01μM",
                  "cells 0.1μM",
                  "cells 1μM",
                  "cells 10μM",
                  "cells 100μM",
                  "OCT4, 0.01μM",
                  "OCT4, 0.1μM",
                  "OCT4, 1μM",
                  "OCT4, 10μM",
                  "OCT4, 100μM"],
            size = (1000, 1000),
            margins=4mm,
            markersize=4,
            yguidefontsize=5,
            legendfontsize=5,
            xguidefontsize=5,
            xtickfontsize=5,
            ytickfontsize=5)
png("corr_C33A_OCT4.png")

plot()
cornerplot(hcat(normx100[same_tp],oct4_x100dat.Column2),
           compact=true,
           label=["C33A cell count, 5FU 100μM",
                  "OCT4 FC, 5FU 100μM"],
            size = (1000, 1000),nbins=3)
png("corr_C33A_OCT4_x100.png")

plot()
cornerplot(hcat(normx10[same_tp],oct4_x10dat.Column2),
           compact=true,
           label=["C33A cell count, 5FU 10μM",
                  "OCT4 FC, 5FU 10μM"],
            size = (1000, 1000),nbins=3)
png("corr_C33A_OCT4_x10.png")

plot()
cornerplot(hcat(normx1[same_tp],oct4_x1dat.Column2),
           compact=true,
           label=["C33A cell count, 5FU 1μM",
                  "OCT4 FC, 5FU 1μM"],
            size = (1000, 1000),nbins=3)
png("corr_C33A_OCT4_x1.png")

plot()
cornerplot(hcat(normx100[same_tp],klf4_x100dat.Column2),
           compact=true,
           label=["C33A cell count, 5FU 100μM",
                  "KLF4 FC, 5FU 100μM"],
            size = (1000, 1000),nbins=3)
png("corr_C33A_KLF4_x100.png")

## Plots
plot()
plot(icded,
              normx0_01,
     yerror = δnormx0_01,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(            oct4_x0_01dat.Column1,
                  oct4_x0_01dat.Column2,
      yerror =    oct4_x0_01dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm)
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 0.01μM")

plot()
plot(icded,
     normx0_1,
     yerror = δnormx0_1,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x0_1dat.Column1,
               oct4_x0_1dat.Column2,
      yerror = oct4_x0_1dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm)
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 0.1μM")

plot()
plot(icded,
     normx1,
     yerror = δnormx1,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x1dat.Column1,
               oct4_x1dat.Column2,
      yerror = oct4_x1dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm)
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 1μM")

plot()
plot(icded,
     normx10,
     yerror = δnormx10,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x10dat.Column1,
               oct4_x10dat.Column2,
      yerror = oct4_x10dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm)
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 10μM")

plot()
plot(icded,
     normx100,
     yerror = δnormx100,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x100dat.Column1,
               oct4_x100dat.Column2,
      yerror = oct4_x100dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm,
      )
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 100μM")
plot!(legendfontsize=15,titlefontsize=20)
png("norm_counts_oct4_100uM.png")

plot()
plot(icded,
     normx10,
     yerror = δnormx10,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x10dat.Column1,
               oct4_x10dat.Column2,
      yerror = oct4_x10dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm,
      )
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 10μM")
plot!(legendfontsize=15,titlefontsize=20)
png("norm_counts_oct4_10uM.png")

using LsqFit
# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
# model(x, p) will accept the full data set as the first argument `x`.
# This means that we need to write our model function so it applies
# the model to the full dataset. We use `@.` to apply the calculations
# across all rows.
#@. model(x, p) = p[1]*exp(-x*p[2])
@. model(x, p) = p[1]/(x+0.1) + p[2]
p0 = [1.8,0.005]
fit100 = curve_fit(model, oct4_x100dat.Column1, oct4_x100dat.Column2, p0)
sigma100 = stderror(fit100)

fit10 = curve_fit(model, oct4_x10dat.Column1, oct4_x10dat.Column2, p0)
sigma10 = stderror(fit10)

using Measurements
params100 = fit100.param .± sigma100
params10 = fit10.param .± sigma10

t=3:2:72
out100 = model(t,params100)
out10 = model(t,params10)

plot()
plot(icded,
     normx100,
     yerror = δnormx100,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x100dat.Column1,
               oct4_x100dat.Column2,
      yerror = oct4_x100dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm,
      )
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 100μM")


plot!(t,out100,ls=:dash,lw=3,label="fitted gene expression model")

plot!(legendfontsize=15,titlefontsize=20)
png("norm_counts_oct4_100uM_fit.png")



plot()
plot(icded,
     normx10,
     yerror = δnormx10,
     #st = :scatter,
     ls = :solid,
     lw = 3,
     marker = (6, 0.9),
     label = "cell count",)
plot!(         oct4_x10dat.Column1,
               oct4_x10dat.Column2,
      yerror = oct4_x10dat.Column3,
      ls = :solid,
      lw = 3,
      marker = (6, 0.9),
      label = "OCT4 fold change",margins=10mm,
      )
xlabel!("time (hours)")
ylabel!("normalized counts")
title!("C33A, 5FU 10μM")


plot!(t,out10,ls=:dash,lw=3,label="fitted gene expression model")

plot!(legendfontsize=15,titlefontsize=20)
png("norm_counts_oct4_10uM_fit.png")
