using DelimitedFiles
using Plots
using LaTeXStrings; pgfplots()

# Nek5000 Utilities

function makenek(case); run(pipeline(`./mn $case`,devnull)); end
function nek(case,np); run(`neklmpi $case $np`); end

function cpn_bb(case,d1,d2)
    cp("$d1/$case.rea","$d2/$case.rea",force=true);
    cp("$d1/$case.map","$d2/$case.map",force=true);
    cp("$d1/r0.f00001","$d2/r0.f00001",force=true);
    cp("$d1/nek5000","$d2/nek5000",force=true);
end

# General Functions

function mklist(conf,imu)
    println("pdrive: writing file.list ...");
    tgt=conf["data"]*"/"*lpad(imu,2,"0");
    case=conf["case"];
    snaps=filter(x->occursin(Regex("^$(case)0\\.f.*"),x),readdir(tgt));
    writedlm("file.list",relpath(tgt)*"/" .* reverse(snaps));
    cp(tgt*"/"*snaps[end],"r0.f00001",force=true);
end

function gflog(conf)
    ntrain=length(conf["ptrain"]);
    for i=1:ntrain
        dir=conf["data"]*"/"*lpad(i,2,"0");
        log="$dir/logfile";
        if isfile(log)
            cp("$dir/logfile","fom."*lpad(i,2,"0")*".log",force=true);
        end
    end
end

function fom(conf,imu)
    dir=conf["data"]*"/"*lpad(imu,2,"0");
    if !isdir(dir) && !islink(dir)
        println("pdrive: running FOM ...");
        cwd=pwd(); 
        case=conf["case"];
        mkpath(dir);
        makenek(case*"_fom");
        cpn_bb(case,".",dir);
        cd(dir);
        p=conf["ptrain"][imu];
        open("mu","w") do f; write(f,"$p"); end
        nek(case,conf["np"]);
        cd(cwd);
    else
        println("pdrive: FOM run mu=$(imu) found ...");
        cp("$dir/logfile","fom."*lpad(imu,2,"0")*".log",force=true);
    end
end

function post(ptrain,err,res,imus,ifqoi)
    println("pdrive: writing err.pdf ...");

    (nit,ntrain)=size(res);

    for it=1:nit
        (minv,ind)=findmin(res[1:it,:],dims=1)

        plot(ptrain,vec(err[ind]),m=:circle,lab=L"||u-u_g||");
        plot!(ptrain,vec(minv),yscale=:log10,m=:circle,lab=L"\Delta");
        xlabel!(L"\mu");
        ylabel!("Error");
        title!("Error vs. mu for Iteration $it");
        savefig("err."*lpad(it,2,"0")*".pdf");
    end

    if ifqoi
        println("pdrive: found qoi.jl, writing qoi.pdf ...");
        fqoi=fill(0.,ntrain);
        rqoi=fill(0.,ntrain);
        for i=1:ntrain
            flog="fom."*lpad(i,2,"0")*".log";
            fqoi[i]=gfqoi(flog);
        end
        for it=1:nit
            (minv,ind)=findmin(res[1:it,:],dims=1)
            for i=1:ntrain
                rlog="rom."*lpad(imus[ind[i][1]],2,"0")*
                          "."*lpad(i,2,"0")*".log";
                rqoi[i]=grqoi(rlog);
            end
            plot(ptrain,rqoi,m=:circle,lab="ROM");
            plot!(ptrain,fqoi,m=:circle,lab="FOM");
            xlabel!(L"\mu");
            ylabel!("QOI");
            savefig("qoi."*lpad(it,2,"0")*".pdf");
        end
    end
end

function sweep!(err,res,conf,imu,it)
    ptrain = conf["ptrain"];
    ntrain = length(ptrain);
    case   = conf["case"];

    open("mu","w") do f; write(f,"$imu"); end

    makenek(case*"_offline");
    mklist(conf,imu);
    println("pdrive: running offline-mode ...");
    nek(case,conf["np"]);

    makenek(case*"_online");

    println("pdrive: online ptrain sweep (mu*=$imu) ...");

    for i in 1:ntrain
        println("pdrive: iteration $i, mu = $(ptrain[i]) ...");
        open("mu","w") do f; write(f,"$(ptrain[i])"); end
        nek(case,conf["np"]);
        log=String(read("$case.log.1"));
        fname="rom."*lpad(imu,2,"0")*"."*lpad(i,2,"0")*".log";
        cp("$case.log.1",fname,force=true);
        res[it,i]=parse(Float64,match(r"res: *(.*)\n",log).captures[1]);
        err[it,i]=parse(Float64,match(r"err: *(.*)\n",log).captures[1]);
    end
end
