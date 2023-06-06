using PackageCompiler

create_sysimage(
    [:CSV, :DataFrames, :Printf, :Glob, :DifferentialEquations, :Plots,
     :Distributions, :DiffEqParamEstim, :NLopt],
    sysimage_path="sys_complete.so"
)
