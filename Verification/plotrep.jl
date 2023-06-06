using Printf, CSV, DataFrames, Glob, Plots, LaTeXStrings, Measures
default(legend=false, fontfamily="Computer Modern", size=(800,600),
    guidefont = 30, tickfont=20, linewidth=1, markersize=5, grid=false,
    dpi=300)


function plotrep(args) #edat,rdat,sdir)

    @show args

    edat = args[1]
    #rdat = args[2]
    sdir = args[2]

    #edat = "calibration_C33A_control_R1_werr.txt";
    #rdat = "C33A_control_R1_werr.txt";

    #import experimental replicates
    #ed = CSV.read(path   , header=0, DataFrame);
    apath = @sprintf("../Data/%s*.txt",edat)
    path  = glob(apath)
    ed    = DataFrame.( CSV.File.( path; header=false ) );
    for i in 1:length(ed)
        ed[i][!, :replicate] .= i
    end
    edf = reduce(vcat, ed);

    # import reference replicate (contains all the time-points)
    #pathref = @sprintf("../Data/%s.txt",rdat)
    #er = CSV.read(pathref, header=0, DataFrame);

    #id_ed=Array(ed[:,1]);
    #id_er=Array(er[:,1]);
    #com     = in.(id_ed,[id_er]);
    #dat     = edf[com.==1,:];
    #nulldat = edf[com.==0,:];

    #broadcast(x -> in(x, id_ed), id_er)
    #sdpath="output_rep_calibration_C33A_control_R1_werr_vs_00/output_rep_";
    #sdir = "output_rep_calibration_C33A_control_R1_werr_vs_00"
    sfiles = glob("*.txt",sdir)
    sd = DataFrame.( CSV.File.( sfiles; header=false ) );
    for i in 1:length(sd)
        sd[i][!, :replicate] .= i
    end
    sdf = reduce(vcat, sd);

    for t in 1:3
        comtot = in.(edf[:,4],[t])
        plot!(edf[comtot,1],edf[comtot,2],yerror=edf[comtot,3],st=:scatter,marker = (5,0.9,:black))
        #plot!(dat[comtot==0,1],dat[comtot==0,2],yerror=dat[comtot==0,3],st=:scatter,m=:o,color=:yellow)
    end
    for f in 1:5000
        com = in.(sdf[:,3],[f]);
        plot!(sdf[com,1],sdf[com,end-1],alpha=0.01,lw=2,color=:blue) #color=:blue
    end
    xlabel!("time (hours)")
    ylabel!("cell count")

    fname=@sprintf("plot_%s.png",sdir);
    png(fname)
end
plotrep(ARGS)
