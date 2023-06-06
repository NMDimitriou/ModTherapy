## Optimization routine
include("load_data.jl")
include("plotopt.jl")

## Load the data
#ARGs=["../Data/calibration_C33A_control_*.txt"
#      "../Verification/output_rep*/*.txt"
#      "../Summary/DSR02_treatment_RUN_1_C33A_x0_01_Rall_obs_sigma/Calibration/final_tmcmc_params.txt"
#      "7"
#      "02_treatment"];
edat = load_data(ARGS[1]);
tmc = load_data(ARGS[2]);
dime = parse(Int, ARGS[4]);
tmcbp = CSV.read(ARGS[3],DataFrame;header=0,delim="\t");
modl = ARGS[5];


mvps  = tmcbp.Column1;
bpid  = [startswith(m,"M") for m in mvps];
bpar  = Array(tmcbp[bpid,2:end-2])

tstart   = minimum(Float64.(edat.Column1))
tend     = maximum(Float64.(edat.Column1))
tspan    = (tstart, tend);
δt       = 1.;
obstimes = tstart:δt:tend;

tmct = Array(tmc[1:Int(tend)+1, 1]); #starts from zero
tmcc = Array(tmc[:, end-1]);
tmcc = reshape(tmcc, (Int(tend)+1, 5000));

num_sets = maximum(Array(edat[1:end, 4]));
num_tp = Int(length(edat[:, 1]) / num_sets - 1);
cidx = (edat[1:end, 1] .> 0.0);
ided = Array(edat[cidx, 1]);
ided = ided[1:num_tp];
eded = Array(edat[cidx, 2]);
σ = Array(edat[cidx, 3]);


## Problem setup
include("modprep.jl")

# Choose model
if (modl == "00")
    import .ModPrep00 as mp
    println("Selecting 00")
	IC = edat.Column2[1];
    prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "G00")
    import .ModPrepG00 as mp
    println("Selecting G00")
	IC = edat.Column2[1];
    prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "00_treatment")
    import .ModPrep00_treatment as mp
    println("Selecting 00_treatment")
	IC = edat.Column2[1];
    prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "00_treatment_expanded")
    import .ModPrep00_treatment_expanded as mp
    println("Selecting 00_treatment_expanded")
	IC = edat.Column2[1];
    prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "00_treatment_expanded_inv")
    import .ModPrep00_treatment_expanded_inv as mp
    println("Selecting 00_treatment_expanded_inv")
	IC = edat.Column2[1];
    prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "00_treatment_expanded_sym")
    import .ModPrep00_treatment_expanded_sym as mp
    println("Selecting 00_treatment_expanded_sym")
	IC = edat.Column2[1];
    prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "00_treatment_expanded_sym_no_gen")
	import .ModPrep00_treatment_expanded_sym_no_gen as mp
	println("00_treatment_expanded_sym_no_gen")
	IC = edat.Column2[1];
	prob_ode = ODEProblem(mp.ode, IC, tspan, mp.p);

elseif (modl == "F00")
    import .ModPrepF00 as mp
    println("Selecting F00")
	IC = edat.Column2[1];
	r=0;
    #prob_ode = ODEProblem(mp.solve, IC, tspan, mp.p);

elseif (modl == "F00_treatment")
    import .ModPrepF00_treatment as mp
    println("Selecting F00_treatment")
	IC = edat.Column2[1];
	r=0;

elseif (modl == "F00_treatment_expanded_inv")
    import .ModPrepF00_treatment_expanded_inv as mp
    println("Selecting F00_treatment_expanded_inv")
    IC = edat.Column2[1];
	if(occursin("x100",ARGS[1]))
		r = [3.17, 0.795]
	else
		r = [1.04, 0.811]
	end

elseif (modl == "F00_treatment_expanded_sym")
    import .ModPrepF00_treatment_expanded_sym as mp
    println("Selecting F00_treatment_expanded_sym")
    IC = edat.Column2[1];
	if(occursin("x100",ARGS[1]))
		r = [3.17, 0.795]
	else
		r = [1.04, 0.811]
	end

elseif (modl == "F01_treatment_expanded_inv")
    import .ModPrepF01_treatment_expanded_inv as mp
    println("Selecting F01_treatment_expanded_inv")
    IC = edat.Column2[1];
	if(occursin("x100",ARGS[1]))
		r = [3.17, 0.795]
	else
		r = [1.04, 0.811]
	end

elseif (modl == "F01_treatment_expanded")
    import .ModPrepF01_treatment_expanded as mp
    println("Selecting F01_treatment_expanded")
    IC = edat.Column2[1];
	if(occursin("x100",ARGS[1]))
		r = [3.17, 0.795]
	else
		r = [1.04, 0.811]
	end

elseif (modl == "02")
    import .ModPrep02 as mp
    println("Selecting 02")
	IC=[edat.Column2[1]-mp.p[end],mp.p[end]];
    prob_ode = ODEProblem(mp.ode!, IC, tspan, mp.p);

elseif (modl == "02_treatment")
	import .ModPrep02_treatment as mp
    println("Selecting 02_treatment")
	IC=[edat.Column2[1]-mp.p[end],mp.p[end]];
    prob_ode = ODEProblem(mp.ode!, IC, tspan, mp.p);

elseif (modl == "01")
    import .ModPrep01 as mp
    println("Selecting 01")
	IC=[edat.Column2[1]-mp.p[end],mp.p[end]];
    prob_ode = ODEProblem(mp.ode!, IC, tspan, mp.p);

end


if (mp.isODE == 1)
	sol_ode = solve(prob_ode,Tsit5(), saveat = δt); #Tsit5()

	## Cost function
	function loss_function(sol)
		cout = Array(sol(ided))
		if(ndims(cout) > 1)
			cout = sum(cout';dims=2)'
			cout = repeat(cout, num_sets)'
		else
			cout = repeat(cout', num_sets)'
		end

	    nzid = (σ.>1e-5);
	    SSN = sum((cout[nzid] .- eded[nzid]) .^ 2 ./ σ[nzid] .^ 2)
	    neg_logL = 0.5 * sum(log(2 * π) .+ log.(σ[nzid] .^ 2)) + 0.5 * SSN

	    return (neg_logL)
	end;

	prob_generator = (prob,q) -> remake(prob_ode,
										u0 = IC,
										p = q);

	cost_function = build_loss_objective(prob_ode,
	                                    Tsit5(),
	                                    loss_function,
	                                    saveat=δt,
	                                    prob_generator = prob_generator,
	                                    maxiters=500,
	                                    verbose=false);


	## Optimization
	# NLopt package
	using NLopt
	#opt               = Opt(:LD_MMA, dime)
	#opt.lower_bounds  = mp.lower
	#opt.upper_bounds  = mp.upper
	#opt.xtol_rel      = 1e-10
	#opt.xtol_abs      = 1e-9
	#opt.maxeval       = 10000
	#opt.min_objective = cost_function;
	opt = Opt(:GN_ISRES, dime)
	local_opt = Opt(:LN_NELDERMEAD,dime)
	opt.lower_bounds  = mp.lower
	opt.upper_bounds  = mp.upper
	opt.xtol_rel      = 1e-10
	opt.xtol_abs      = 1e-9
	opt.maxeval       = 50000
	opt.min_objective = cost_function;
	opt.local_optimizer=local_opt
	(minf,minx,ret)   = NLopt.optimize(opt, mp.initial_x)
	numevals = opt.numevals # the number of function evaluations
	println("got $minf at $minx after $numevals iterations (returned $ret)")

	## solve with optimal parameters
	prob = remake(prob_ode, u0 = IC, p = minx)
	sol  = solve(prob, Tsit5(), saveat = δt)

	# solve with optimal parameters from TMCMC
	prob = remake(prob_ode, u0 = IC, p = bpar)
	soltmc  = solve(prob, Tsit5(), saveat = δt)

	# extract solution for error calculation
	soltp = Array(sol(ided));
	if(ndims(soltp)>1)
		soltp = sum(soltp;dims=1)'
	end
	#for plots
	solp  = Array(sol(obstimes))
	if(ndims(solp)>1)
		solp = sum(solp;dims=1)'
	end
	solmc = Array(soltmc(obstimes))
	if(ndims(solmc)>1)
		solmc = sum(solmc;dims=1)'
	end
	leg="ISRES"

else
	sol = [IC];
	sol = mp.do_solve(sol,obstimes,mp.p,r,Int(tend/δt + 1))

	using NLopt

	function loss_function(x::Vector, grad::Vector)
		solc = [IC];
		solc = mp.do_solve(solc,obstimes,x,r,Int(tend/δt + 1))

		ct = solc[ided]
		cout = repeat(ct, num_sets)

		nzid = (σ.>1e-5);
		SSN = sum((cout[nzid] .- eded[nzid]) .^ 2 ./ σ[nzid] .^ 2)
		neg_logL = 0.5 * sum(log(2 * π) .+ log.(σ[nzid] .^ 2)) + 0.5 * SSN

		return (neg_logL)
	end

	#opt = Opt(:G_MLSL, dime)
	opt = Opt(:GN_ISRES, dime)
	local_opt = Opt(:LN_NELDERMEAD,dime)
	opt.lower_bounds  = mp.lower
	opt.upper_bounds  = mp.upper
	opt.xtol_rel      = 1e-10
	opt.xtol_abs      = 1e-9
	opt.maxeval       = 100000
	opt.min_objective = loss_function;
	opt.local_optimizer=local_opt
	(minf,minx,ret)   = NLopt.optimize(opt, mp.initial_x)
	numevals = opt.numevals # the number of function evaluations
	println("got $minf at $minx after $numevals iterations (returned $ret)")
	#res=Optim.optimize(loss_function, mp.lower, mp.upper, mp.initial_x)

	sol	  = [IC];
	sol   = mp.do_solve(sol,obstimes,minx,r,Int(tend/δt + 1));
	soltp = sol[ided];
	solp  = sol;
	sol	  = [IC];
	solmc = mp.do_solve(sol,obstimes,bpar,r,Int(tend/δt + 1));
	leg="ISRES"
end
#println("got $minf at $minx after $numevals iterations (returned $ret)")


## save optimal parameters
io = open("opt_params_mlsl_lds.txt", "w") do io
	for x in minx
		println(io,x)
	end
end

io = open("mlsl_lds_log.txt", "w") do io
	println(io,"got $minf at $minx after $numevals iterations (returned $ret)")
end


## Error estimation
import .Stats00 as st

# estimate RMSE
#soltp = Array(sol(ided))
#soltmcmc = Array(soltmc(ided))
eded = reshape(eded, (num_tp, num_sets))
Ermse,δrmse = st.calcRMSE(solmc[ided],eded)

# estimate CCC
matC = reduce(hcat, (solmc[ided], eded))
MCCC,δCCC = st.calcCCC(matC)

io = open("rmse_ccc.txt", "w") do io
    println(io, "RMSE = ", Ermse, " ± ", δrmse)
    println(io, "CCC = ", MCCC, " ± ", δCCC)
end


## Plot


#solp = sum(solp;dims=1)
#solmc = sum(solmc;dims=1)'
plot()
y = [1.0, 2.0, 3.0, 4.0, 5.0] * 10^5
x=y
h2=plot!(tmct, tmcc, alpha = 0.01, lw = 1, color = :blue, label = "")
h3=plot!(obstimes, solp, lw = 3, color = :red, label = leg)
h1=plot!(
    obstimes,
    solmc,
    lw = 2,
    ls = :dash,
    color = :blue,
    label = "TMCMC",
    legend = :topleft,
    yticks = y,
    yaxis = (formatter = y -> string(round(Int, y / 10^5))),
    margin = 5mm,
    top_margin = 10mm,
)
h4=plot!(
    edat[:, 1],
    edat[:, 2],
    yerror = edat[:, 3],
    st = :scatter,
    marker = (6, 0.9, :black),
    label = "C33A",
)
yaxis!([0, 5e5])
annotate!([(
    abs(min(xlims(h1)[1],xlims(h4)[1])),
    1.05*max(ylims(h1)[2],ylims(h4)[2]),
    Plots.text(L"\times10^{5}", 20, :black, :center),
)])
xlabel!("time (hours)")
ylabel!("cell count")
plot!(legendfontsize=18,titlefontsize=20)

fname = @sprintf("plot_final.png");
png(fname)


avedat = mean(eded;dims=2)
stddat = reshape(σ, (num_tp, num_sets))
semdat = sqrt.(sum(stddat.^2;dims=2)./num_sets^2)
#stddat = std(eded;dims=2)
stdtmc = std(tmcc;dims=2)
tp     = edat.Column1[2:num_tp+1]

plot()
a=plot( avedat,
        #sum(Array(soltmc[ided]);dims=1)',
		Array(solmc[ided]),
        xerror = semdat,
        yerror = stdtmc[ided],
        st=:scatter,
        marker = (7, 0.9, :black),
        label="")
b=plot!(avedat,
        avedat,
        ls=:dash,
        label="",lw=3,
        margin = 15mm,
        yticks = y,
        xticks = x,
        yaxis = (formatter = y -> string(round(Int, y / 10^5))),
        xaxis = (formatter = x -> string(round(Int, x / 10^5)))
        )
xlabel!("experimental observations")
ylabel!("calibrated model")
annotate!([(
        1.1*min(xlims(a)[1],xlims(b)[1]),
        1.05*max(ylims(a)[2],ylims(b)[2]),
    Plots.text(L"\times10^{5}", 20, :black, :center),
)])
annotate!([(
    1.05*max(xlims(a)[2],xlims(b)[2]),
    1.1 *min(ylims(a)[1],ylims(b)[1]),
    Plots.text(L"\times10^{5}", 20, :black, :center),
)])

txt=("CCC = $(round(MCCC, sigdigits=3)) ± $(round(δCCC, sigdigits=3))")
annotate!([(
    0.48*max(xlims(a)[2],xlims(b)[2]),
    0.75*max(ylims(a)[2],ylims(b)[2]),
    Plots.text(txt, 20, :black, :center)
)])

fname = @sprintf("plot_err.png");
png(fname)
