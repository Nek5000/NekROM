println("""
#####################################

# Julia Driver for NekROM
# v0.2.0

# Kento Kaneko
# 2022-05-06

#####################################
""");

include("functions.jl");

println(" ");

print("now: ");
println(now(Dates.UTC));
println("cwd: "*pwd());

if "mor.yaml" in readdir()
    println(" ");
    println("reading mor.yaml ...");
else
    println("did not find mor.yaml, exiting ...");
    exit
end

readconf("mor.yaml");
initvars();
load(ifrom,ifaug,gxyz != Nothing,inus,ifleray);

println(" ");

if ifrom[1]
    (maxe,e0)=(ceig(cu).*dt);
    println("λΔt_cu: $(@sprintf("%.3e , %.3e",maxe,e0))");
end

if ifrom[2]
    (maxe,e0)=(ceig(ct).*dt);
    println("λΔt_ct: $(@sprintf("%.3e , %.3e",maxe,e0))");
end

println(" ");

if ifcopt
    if ifrom[1]
        global ukmin=minimum(uk,dims=2);
        global ukmax=maximum(uk,dims=2);
    end
    if ifrom[2]
        global tkmin=minimum(tk,dims=2);
        global tkmax=maximum(tk,dims=2);
    end
end

ucoef=zeros(Int(nsteps//iostep),nb+1);
tcoef=zeros(Int(nsteps//iostep),nb+1);

kappa0=1.0/5300;
t[1,1]=kappa0/kappa;

for istep=1:nsteps
    ito=min(istep,3);
    global gamma,gxyz
    if istep <= 3; global htfac=[]; global hufac=[]; end
    if ifrom[2]
        f=zeros(nb,1);
        if gamma != 0.0
            f0=(u[:,1]'*tbn)';
            f=f-gamma*f0[2:end];
            f=f+tsa[2:end]*pi*4.0;
        end
        global (rt,et)=setr(at0,bt,ct,f,u,t,et,kappa,alphas,betas,dt,ito)
        if ifcopt
            if istep < 4
                ht=kappa*at+bt*(betas[1,ito]/dt);
                ht=.5*(ht+ht');
                global vt=zeros(nb+1,nb+1);
                global vt1=zeros(nb,nb);

                (d,vt1)=eigen(ht);
                vt[1,1]=1.0;
                vt[2:end,2:end]=vt1;

                global atp=vt1'*at*vt1;
                global btp=vt1'*bt*vt1;
            end

            rtp=vt1'*rt;
            global (t_new,htfac)=step(rtp,atp,btp,kappa,betas,dt,ito,htfac)

            if istep < 4
                tkp=vt'*tk;
                global tkmin=minimum(tkp,dims=2);
                global tkmax=maximum(tkp,dims=2);
            end

            for i=1:nb
                if t_new[i+1] > tkmax[i+1]; t_new[i+1]=tkmax[i+1]; end
                if t_new[i+1] < tkmin[i+1]; t_new[i+1]=tkmin[i+1]; end
            end
            t_new=vt*t_new;
        else
            global (t_new,htfac)=step(rt,at,bt,kappa,betas,dt,ito,htfac)
        end
        t_new[1]=t[1,1];
    end
    if ifrom[1]
        f=zeros(nb,1);
        if gxyz != Nothing
            f=f-buxt*t[:,1]*gxyz[1];
            f=f-buyt*t[:,1]*gxyz[2];
            f=f-buzt*t[:,1]*gxyz[3];
        end
        global (ru,eu)=setr(au0,bu,cu,f,u,u,eu,nu,alphas,betas,dt,ito)
        if ifcopt
            if istep < 4
                hu=nu*au+bu*(betas[1,ito]/dt);
                hu=.5*(hu+hu');

                global vu=zeros(nb+1,nb+1);
                global vu1=zeros(nb,nb);
                (d,vu1)=eigen(hu);

                vu[1,1]=1.0;
                vu[2:end,2:end]=vu1;

                global aup=vu1'*au*vu1;
                global bup=vu1'*bu*vu1;
            end

            rup=vu1'*ru;
            global (u_new,hufac)=step(rup,aup,bup,nu,betas,dt,ito,hufac)

            if istep < 4
                ukp=vu'*uk;
                global ukmin=minimum(ukp,dims=2);
                global ukmax=maximum(ukp,dims=2);
            end

            for i=1:nb
                if u_new[i+1] > ukmax[i+1]; u_new[i+1]=ukmax[i+1]; end
                if u_new[i+1] < ukmin[i+1]; u_new[i+1]=ukmin[i+1]; end
            end
            u_new=vu*u_new;
        else
            global (u_new,hufac)=step(ru,au,bu,nu,betas,dt,ito,hufac)
        end
    end
    global time=time+dt;

    if ifrom[1]
        global u=shift(u,u_new,3);
        if istep > astep
            global ua=ua+u_new;
            global u2a=u2a+u_new*u_new';
            if istep == nsteps
                global ua=ua/(nsteps-astep);
                global u2a=u2a/(nsteps-astep);
            end
        end
    end

    if ifrom[2]
        global t=shift(t,t_new,3);
        if istep > astep
            global ta=ta+t_new;
            global t2a=t2a+t_new*t_new';
            if istep == nsteps
                global ta=ta/(nsteps-astep);
                global t2a=t2a/(nsteps-astep);
            end
        end
    end

    if  mod(istep,iostep) == 0
        s=@sprintf("time: %.5e",time);
        if ifrom[1]
            uene=u[:,1]'*bu0*u[:,1]; uene=uene[1];
            s=s*@sprintf(" , ue: %.5e",uene);
        end
        if ifrom[2]
            tene=t[:,1]'*bt0*t[:,1]; tene=tene[1];
            s=s*@sprintf(" , te: %.5e",tene);
        end
        if inus == 1
            rnus=nus'*t[:,1];
            s=s*@sprintf(" , nu: %.5e",rnus);
        end
        if inus == 2
            ts=tsa'*t[:,1];
            tb=(u[:,1]'*tbn*t[:,1])/(tbd'*u[:,1])
            rnus=1.0/(kappa*(ts-tb));
            s=s*@sprintf(" , nu: %.5e",rnus);
        end
        println(s);
        ucoef[Int(istep//iostep),:]=u[:,1];
        tcoef[Int(istep//iostep),:]=t[:,1];
    end
    if  mod(istep,100) == 0
        for i=1:nb
            println(@sprintf("%i %.8e %.8e romu",i,time,u[i+1,1]));
        end
        for i=1:nb
            println(@sprintf("%i %.8e %.8e romt",i,time,t[i+1,1]));
        end
    end
end

mkpath("out");

if ifrom[1]
    open("out/uf","w") do io; writedlm(io,u[:,1]); end
    open("out/ua","w") do io; writedlm(io,ua); end
    open("out/u2a","w") do io; writedlm(io,reshape(u2a,(nb+1)^2,1)); end
    open("out/bu0","w") do io; writedlm(io,reshape(bu0,(nb+1)^2,1)); end
    open("out/au0","w") do io; writedlm(io,reshape(au0,(nb+1)^2,1)); end
    open("out/wu","w") do io; writedlm(io,wu); end
    open("out/ucoef","w") do io; writedlm(io,ucoef); end
end

if ifrom[2]
    open("out/tf","w") do io; writedlm(io,t[:,1]); end
    open("out/ta","w") do io; writedlm(io,ta); end
    open("out/t2a","w") do io; writedlm(io,reshape(t2a,(nb+1)^2,1)); end
    open("out/bt0","w") do io; writedlm(io,reshape(bt0,(nb+1)^2,1)); end
    open("out/at0","w") do io; writedlm(io,reshape(at0,(nb+1)^2,1)); end
    open("out/wt","w") do io; writedlm(io,wt); end
    open("out/tcoef","w") do io; writedlm(io,tcoef); end
end

if inus == 1
    open("out/nus","w") do io; writedlm(io,reshape(nus,nb+1,1)); end
end

if inus == 2
    open("out/tsa","w") do io; writedlm(io,reshape(tsa,nb+1,1)); end
    open("out/tbn","w") do io; writedlm(io,reshape(tbn,(nb+1)*(nb+1),1)); end
    open("out/tbd","w") do io; writedlm(io,reshape(tbd,nb+1,1)); end
end

println(" ");

print("now: ");
println(now(Dates.UTC));
println("cwd: "*pwd());
