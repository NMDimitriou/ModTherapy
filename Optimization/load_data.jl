using Printf, CSV, DataFrames, Glob


function load_data(path)

    @show path

    p  = glob(path)
    ed = DataFrame.( CSV.File.( p ; header=false ) );

    for i in 1:length(ed)
        ed[i][!, :replicate] .= i
    end
    ed = reduce(vcat, ed);

    return(ed)
end

function prep_dat(cdat)

    num_sets = maximum(Array(cdat[1:end, 4]));
    num_tp = Int(length(cdat[:, 1]) / num_sets - 1);
    cidx = (cdat[1:end, 1] .> 0.0);
    icded = Array(cdat[cidx, 1]);
    icded = icded[1:num_tp];
    cded = Array(cdat[cidx, 2]);
    σc = Array(cdat[cidx, 3]);

    return num_sets, num_tp, cidx, icded, cded, σc

end

function mean_sem_dat(cded,σ,num_tp,num_sets)

    cded    = reshape(cded, (num_tp, num_sets))
    av_cded = sum(cded;dims=2)./num_sets
    σcded   = reshape(σc, (num_tp, num_sets))
    δσc     = sqrt.(sum(σcded.^2,dims=2))./num_sets

    return av_cded, δσc

end

function δnorm(x,y,δx,δy)

    dnorm = (x./y) .* sqrt.((δx./x).^2 .+ (δy./y).^2)

    return dnorm
end
