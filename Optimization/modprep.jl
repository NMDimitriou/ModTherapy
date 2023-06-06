using DifferentialEquations
using Distributions
using DiffEqParamEstim

module ModPrep00

    export p, lower, upper, initial_x, ode, isODE

    using Scanf, Printf

    isODE=1;

    lower=Vector{Float64}(undef,2)
    upper=Vector{Float64}(undef,2)

    open("../src/00/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:1
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2])]; #r_γ, K

    initial_x = p;

    function ode(u,p,t)
        (rᵧ, K) = p
        growth = rᵧ*u*(1-u/K)
        return growth
    end;

    precompile(ode, (Float64,))

end

module ModPrepG00

    export p, lower, upper, initial_x, ode, isODE

    using Scanf, Printf

    isODE=1;

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/G00/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K

    initial_x = p;

    function ode(u,p,t)
        (rᵧ, K, θ) = p
        growth = rᵧ*u*(1.0-(u/K)^θ) / θ
        return growth
    end;

    precompile(ode, (Float64,))

end

module ModPrep00_treatment

    export p, lower, upper, initial_x, ode, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/00_treatment/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, u_eff

    initial_x=p;

    function ode(u,p,t)
        (rᵧ, K, u_eff) = p
        growth = rᵧ*u*(1-u/K) - u_eff*u
        return growth
    end;

    precompile(ode, (Float64,))
end

module ModPrep00_treatment_expanded

    export p, lower, upper, initial_x, ode, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/00_treatment/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, u_eff

    initial_x=p;

    function ode(u,p,t)
        (rᵧ, K, u_eff) = p
        growth = rᵧ*u.*(1-u/K) .- u_eff*(3.05./(t .+ 0.1) + 0.8).*u
        return growth
    end;

    precompile(ode, (Float64,))
end #ModPrep00_treatment_expanded


module ModPrep00_treatment_expanded_sym

    export p, lower, upper, initial_x, ode, isODE

    isODE=1;

    using Scanf, Printf
ARGs=["../Data/calibration_C33A_control_*.txt"
      "../Verification/output_rep*/*.txt"
      "../Summary/DSR02_treatment_RUN_1_C33A_x0_01_Rall_obs_sigma/Calibration/final_tmcmc_params.txt"
      "7"
      "02_treatment"];
edat = load_data(ARGS[1]);
tmc = load_data(ARGS[2]);
dime = parse(Int, ARGS[4]);
tmcbp = CSV.read(ARGS[3],DataFrame;header=0,delim="\t");
modl = "00_treatment_expanded_sym"#ARGS[5];


mvps  = tmcbp.Column1;
bpid  = [startswith(m,"M") for m in mvps];
bpar  = Array(tmcbp[bpid,2:end-2])
    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/00_treatment_expanded_sym/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, u_eff

    initial_x=p;
    r = [3.17, 0.795]
    #r = [1.04, 0.811]

    function ode(u,p,t)
        (rᵧ, K, u_eff) = p

        pα     = (abs.(r[1] ./(t .- 73.1)) .+ r[2]);
        growth = rᵧ*u.*(1-u/K) .- u_eff*pα.*u
        return growth
    end;

    precompile(ode, (Float64,))
end #ModPrep00_treatment_expanded_sym

module ModPrep00_treatment_expanded_sym_no_gen

	export p, lower, upper, initial_x, ode, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,5)
    upper=Vector{Float64}(undef,5)

    open("../src/00_treatment_expanded_sym_no_gen/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:4
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
		 .5*(lower[4].+upper[4]),
		 .5*(lower[5].+upper[5])];
	initial_x=p;

    function ode(u,p,t)
        (rᵧ, K, u_eff, A, B) = p

        pα     = (abs.(A ./(t .- 73.1)) .+ B);
        growth = rᵧ*u.*(1-u/K) .- u_eff*pα.*u
        return growth
    end;

    precompile(ode, (Float64,))
end #ModPrep00_treatment_expanded_sym_no_gen

module ModPrep00_treatment_expanded_inv

    export p, lower, upper, initial_x, ode, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/00_treatment_expanded_inv/tmcmc.par") do f
        for i in eachline(f)
            for j in 0:2
                bound=Char(8);
                bound=Printf.@sprintf("B%d",j);
                if occursin(r"^#", i) continue; end

                if occursin(bound, i)
                   lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                   upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
                end
            end
        end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, u_eff

    initial_x=p;
    r = [3.17, 0.795]
    #r = [1.04, 0.811]
    function ode(u,p,t)
        (rᵧ, K, u_eff) = p
        growth = rᵧ*u.*(1-u/K) .- u_eff*(1 ./(r[1] ./(t .+ 0.1) + r[2])).*u
        return growth;
    end
    precompile(ode, (Float64,));
end # module ModPrep00_treatment_expanded_inv


module ModPrepF00

    export p, lower, upper, initial_x, do_solve, isODE

    isODE=0;

    using Scanf, Printf

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/F00/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, α

    initial_x = p;

    using SpecialFunctions

    function MittagLefler(z,a,n)

        Ea = 0.;
        for k in 0:n
            Ea += z^k / gamma(a*k+1.0);
        end

        return Ea;
    end

    function do_solve(u,t,p,r,Nt)

        for i in 2:Nt
            push!(u, 1.0/(1.0/p[2] + (1.0/u[1] - 1.0/p[2])*MittagLefler(-p[1]*t[i]^p[3],p[3],200)));
        end

        return u
    end

end  # module ModPrepF00



module ModPrepF00_treatment

    export p, lower, upper, initial_x, do_solve, isODE

    isODE=0;

    using Scanf, Printf

    lower=Vector{Float64}(undef,4)
    upper=Vector{Float64}(undef,4)

    open("../src/F00_treatment/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:3
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
         .5*(lower[4].+upper[4])]; #r_γ, K, α

    initial_x = p;

    using SpecialFunctions

    function MittagLefler(z,a,n)

        Ea = 0.;
        for k in 0:n
            Ea += z^k / gamma(a*k+1.0);
        end

        return Ea;
    end

    function do_solve(u,t,p,r,Nt)

        for i in 2:Nt
            push!(u, (p[1] - p[3])/(p[1]/p[2] - (p[1]/p[2] - (p[1]-p[3])/u[1])*MittagLefler(-(p[1]-p[3])*t[i]^p[4],p[4],200)));
        end

        return u
    end

end  # module ModPrepF00_treatment


module ModPrepF00_treatment_expanded_inv

    export p, lower, upper, initial_x, do_solve, isODE

    isODE=0;

    using Scanf, Printf

    lower=Vector{Float64}(undef,4)
    upper=Vector{Float64}(undef,4)

    open("../src/F00_treatment_expanded_inv/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:3
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
         .5*(lower[4].+upper[4])]; #r_γ, K, α

    initial_x = p;

    using SpecialFunctions

    function MittagLefler(z,a,n)

        Ea = 0.;
        for k in 0:n
            Ea += z^k / gamma(a*k+1.0);
        end

        return Ea;
    end


    function do_solve(u,t,p,r,Nt)
        #α = Vector{Float64}(undef,Nt);
        α = p[4].*(1 ./(r[1] ./(t .+ 0.1) .+ r[2]));
        α/= maximum(α);

        for i in 2:Nt
            push!(u, (p[1] - p[3])/(p[1]/p[2] - (p[1]/p[2] - (p[1]-p[3])/u[1])*MittagLefler(-(p[1]-p[3])*t[i]^α[i],α[i],200)));
        end

        return u
    end

end  # module ModPrepF00_treatment_expanded_inv


module ModPrepF00_treatment_expanded_sym

    export p, lower, upper, initial_x, do_solve, isODE

    isODE=0;

    using Scanf, Printf

    lower=Vector{Float64}(undef,4)
    upper=Vector{Float64}(undef,4)

    open("../src/F00_treatment_expanded_sym/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:3
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
         .5*(lower[4].+upper[4])]; #r_γ, K, α

    initial_x = p;

    using SpecialFunctions

    function MittagLefler(z,a,n)

        Ea = 0.;
        for k in 0:n
            Ea += z^k / gamma(a*k+1.0);
        end

        return Ea;
    end


    function do_solve(u,t,p,r,Nt)
        #α = Vector{Float64}(undef,Nt);
        α = p[4].*(abs.(r[1] ./(t .- 73.1)) .+ r[2]);
        #α = p[4].*(1 ./(r[1] ./(t .+ 0.1) .+ r[2]));
        α/= maximum(α);

        for i in 2:Nt
            push!(u, (p[1] - p[3])/(p[1]/p[2] - (p[1]/p[2] - (p[1]-p[3])/u[1])*MittagLefler(-(p[1]-p[3])*t[i]^α[i],α[i],200)));
        end

        return u
    end

end  # module ModPrepF00_treatment_expanded_sym


module ModPrepF01_treatment_expanded_inv

    export p, lower, upper, initial_x, do_solve, isODE

    isODE=0;

    using Scanf, Printf

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/F01_treatment_expanded_inv/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, α

    initial_x = p;

    using SpecialFunctions

    function MittagLefler(z,a,n)

        Ea = 0.;
        for k in 0:n
            Ea += z^k / gamma(a*k+1.0);
        end

        return Ea;
    end


    function do_solve(u,t,p,r,Nt)
        #α = Vector{Float64}(undef,Nt);
        α = p[3].*(1 ./(r[1] ./(t .+ 0.1) .+ r[2]));
        α/= maximum(α);

        for i in 2:Nt
            push!(u, u[1]*MittagLefler((p[1]-p[2])*t[i]^α[i],α[i],200));
        end

        return u
    end

end  # module ModPrepF01_treatment_expanded_inv



module ModPrepF01_treatment_expanded

    export p, lower, upper, initial_x, do_solve, isODE

    isODE=0;

    using Scanf, Printf

    lower=Vector{Float64}(undef,3)
    upper=Vector{Float64}(undef,3)

    open("../src/F01_treatment_expanded/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:2
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
                lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
                upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3])]; #r_γ, K, α

    initial_x = p;

    using SpecialFunctions

    function MittagLefler(z,a,n)

        Ea = 0.;
        for k in 0:n
            Ea += z^k / gamma(a*k+1.0);
        end

        return Ea;
    end


    function do_solve(u,t,p,r,Nt)
        #α = Vector{Float64}(undef,Nt);
        α = p[3].*(r[1] ./(t .+ 0.1) .+ r[2]);
        α/= maximum(α);

        for i in 2:Nt
            push!(u, u[1]*MittagLefler((p[1]-p[2])*t[i]^α[i],α[i],200));
        end

        return u
    end

end  # module ModPrepF01_treatment_expanded



module ModPrep01

    export p, lower, upper, initial_x, ode!, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,7)
    upper=Vector{Float64}(undef,7)

    open("../src/01/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:6
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
         .5*(lower[4].+upper[4]),
         .5*(lower[5].+upper[5]),
         .5*(lower[6].+upper[6]),
         .5*(lower[7].+upper[7])]; #r_γ, K, u_eff

    initial_x = p;

    function ode!(du,u,p,r,t)
        (S, R) = u
        (α, ϕ, μrs, μsr, r₁, r₂, NR) = p
        N = S+R
        growthₛ = -α  *S + (2.0-μsr)*α  *S + μrs*α*ϕ*R - r₁*S*R
        growthᵣ = -α*ϕ*R + (2.0-μrs)*α*ϕ*R + μsr*α  *S - r₂*S*R
        @inbounds begin
            du[1] = growthₛ
            du[2] = growthᵣ
        end
        nothing
        #du .= (growthₛ, growthᵣ)
    end;

    precompile(ode!, (Float64,Float64))

    export u0, p, lower, upper, initial_x, ode!
end


module ModPrep02

    export p, lower, upper, initial_x, ode!, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,4)
    upper=Vector{Float64}(undef,4)

    open("../src/02/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:3
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
         .5*(lower[4].+upper[4])]; #r_γ, K, u_eff

    initial_x=p;

    function ode!(du,u,p,t)
        (rₛ,rᵣ,K,NR) = p
        (S,R)=u
        growthₛ = rₛ*S*(1-(S.+R)/K)
        growthᵣ = rᵣ*R*(1-(S.+R)/K)
        @inbounds begin
            du[1] = growthₛ
            du[2] = growthᵣ
        end
        nothing
        #return S.+R
    end;

    precompile(ode!, (Float64,))
end

module ModPrep02_treatment

    export p, lower, upper, initial_x, ode!, isODE

    isODE=1;

    using Scanf, Printf

    lower=Vector{Float64}(undef,7)
    upper=Vector{Float64}(undef,7)

    open("../src/02_treatment/tmcmc.par") do f
    for i in eachline(f)
        for j in 0:6
            bound=Char(8);
            bound=Printf.@sprintf("B%d",j);
            if occursin(r"^#", i) continue; end

            if occursin(bound, i)
               lower[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end-1]
               upper[j+1] = Scanf.@scanf(i, "%s %f %f", String,Float64,Float64)[end]
            end
        end
    end
    end

    p = [.5*(lower[1].+upper[1]),
         .5*(lower[2].+upper[2]),
         .5*(lower[3].+upper[3]),
         .5*(lower[4].+upper[4]),
         .3*(lower[5].+upper[5]),
         .5*(lower[6].+upper[6]),
         .5*(lower[7].+upper[7])]; #r_γ, K, u_eff

    initial_x=p;

    function ode!(du,u,p,t)
        (rₛ,rᵣ,K,α,dₛ,dᵣ,NR) = p
        (S,R)=u
        growthₛ = rₛ*S*(1-(S.+R)./K) - α*S - dₛ*S
        growthᵣ = rᵣ*R*(1-(S.+R)./K) + α*S - dᵣ*R
        @inbounds begin
            du[1] = growthₛ
            du[2] = growthᵣ
        end
        nothing
        #return S.+R
    end;

    precompile(ode!, (Float64,))
end

module Stats00

    import Statistics as stat

    export Ermse, δrmse, MCCC, δCCC

    function calcRMSE(soltp,eded)
        rmse = sqrt.(stat.mean((soltp .- eded) .^ 2, dims = 1))
        Ermse = stat.mean(rmse)
        δrmse = stat.std(rmse)

        return Ermse,δrmse
    end

    function calcCCC(matC)
        ρ = stat.cor(matC, dims = 1)
        EC = stat.mean(matC, dims = 1)
        δC = stat.std(matC, dims = 1)
        CCC =
            2 * ρ .*
            (((EC' .- EC) .^ 2) ./ (δC' .* δC) + δC' ./ δC + (δC' ./ δC) .^ -1) .^ (-1)

        MCCC = stat.mean(CCC[1, 2:end])
        δCCC = stat.std(CCC[1, 2:end])

        return MCCC,δCCC
    end
end
