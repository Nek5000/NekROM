#####################################

# Julia Functions for NekROM
# v0.2.0

# Kento Kaneko
# 2022-05-06

#####################################

using Printf, Dates, YAML, DelimitedFiles, LinearAlgebra
BLAS.set_num_threads(2)

function readconf(fname)
    d=YAML.load_file(fname);
    ks=keys(d);
    global ifrom=[false,false];

    global nsteps  = Nothing;
    global iostep  = Nothing;
    global astep   = 0;
    global dt      = Nothing;
    global nb      = Nothing;
    global gxyz    = Nothing;
    global ifaug   = false;
    global inus    = 0;
    global nu      = Nothing;
    global kappa   = Nothing;
    global base    = "./";
    global ips     = Nothing;
    global gamma   = 0.0;
    global ifleray = false;
    global ifcopt  = false;

    if "nsteps" in ks
        nsteps = d["nsteps"];
        println("nsteps:        $nsteps");
        nsteps = Int(nsteps);
    else
        println("missing nsteps ...");
        exit
    end

    if "iostep" in ks
        iostep = d["iostep"];
        println("iostep:        $iostep");
        iostep = Int(iostep);
    else
        println("missing iostep ...");
        exit
    end

    if "astep"    in ks; astep  = d["astep" ]; astep=Int(astep);  end
    if "dt" in ks
        dt = d["dt"    ];
        println(@sprintf("Δt:            %.3e",dt));
    end
    if "type" in ks
        ips = d["type"];
        println("type:           $ips");
    end

    if "gravity" in ks
        gxyz = d["gravity"];
        println("gravity:       $gxyz");
    end

    if "augment" in ks
        ifaug = d["augment"];
        println("augment:       $ifaug");
    end

    if "nusselt" in ks
        inus = d["nusselt"];
        println("nusselt:       $inus");
    end

    if "gamma" in ks
        gamma = d["gamma"];
        println("gamma:        $gamma");
    end

    if "nb" in ks;
        nb = d["nb"];
        println("nb:            $nb");
        if nb < 0; nb = parse(Int,basename(pwd())); end
    end

    if "conductivity" in ks
        kappa=d["conductivity"];
        println("conductivity:  $kappa");
        if kappa < 0.0; kappa = -1.0 / kappa; end
        ifrom[2]=true;
    end

    if "viscosity" in ks;
        nu=d["viscosity"];
        println("viscosity:     $nu");
        if nu < 0.0; nu = -1.0 / nu; end
        ifrom[1]=true;
    end

    if "base" in ks;
        base=d["base"];
        println("base:          $base");
    end

    if "leray" in ks;
        ifleray=d["leray"];
        println("leray:         $ifleray");
    end

    if "copt" in ks;
        ifcopt=d["copt"];
        println("copt:          $ifcopt");
    end

    println(" ");
end

function initvars()
    global time=0.;
    global (alphas,betas)=setcoef();
end

function load(ifrom,ifaug,ifbuoy,inus,ifleray)
    println("loading ROM operators and vectors ...");

    global nb,mb;

    if ifrom[1]
        global bu,bu0,au,au0,cu,cu0,u,ua,u2a,ru,eu,hufac
    end

    if ifrom[2]
        global bt,bt0,at,at0,ct,ct0,t,ta,t2a,rt,et,htfac
    end
    
    mb=readdlm(base*"ops/nb"); mb=Int(mb[1]); println("mb: $mb");
    if nb == 0; nb = mb; end

    if ifaug
#       nb2=Int((nb-1)//2);
#       mb2=Int((mb-1)//2);
        nb3=Int((nb-1)//3);
        mb3=Int((mb-1)//3);

#       index=vcat(collect(1:nb2+1),collect(mb2+2:mb2+2+nb2));
#       index1=vcat(collect(1:nb2),collect(mb2+1:mb2+nb2+1));

        index=vcat(collect(1:nb3+1),collect(mb3+2:mb3+2+nb3),collect(mb3*2+3:mb3*2+2+nb3));
        index1=vcat(collect(1:nb3),collect(mb3+1:mb3+1+nb3),collect(mb3*2+2:mb3*2+1+nb3));
    else
        index=collect(1:nb+1);
        index1=collect(1:nb);
    end

    if ifrom[1]
        bu0_full=readdlm(base*"ops/bu");
        bu0_full=reshape(bu0_full,mb+1,mb+1);
        bu0=bu0_full[index,index];
        bu0=.5*(bu0+bu0');

        au0_full=readdlm(base*"ops/au");
        au0_full=reshape(au0_full,mb+1,mb+1);
        au0=au0_full[index,index];
        au0=.5*(au0+au0');

        cu_full=readdlm(base*"ops/cu");
        cu_full=reshape(cu_full,mb,mb+1,mb+1);
        cu=cu_full[index1,index,index];
        cu=reshape(cu,nb*(nb+1),nb+1);

        u0_full=readdlm(base*"ops/u0");
        if isfile("ops/u0"); u0_full=readdlm("ops/u0"); end
        u0=u0_full[index];

        global wu=zeros(nb+1,nb+1);
        for i=1:nb+1; wu[i,i]=1.0; end
        if ifaug
            if ips == "H10"
                (d,v)=eigen(au0[2:end,2:end]);
            else
                (d,v)=eigen(bu0[2:end,2:end]);
            end
            wu[2:end,2:end]=v;
        end
        wu1=wu[2:end,2:end];

        au0=wu'*au0*wu;
        bu0=wu'*bu0*wu;

        cu=reshape(wu1'*reshape(cu*wu,nb,(nb+1)*(nb+1)),nb,nb+1,nb+1);
        for i=1:nb+1; cu[:,:,i]=cu[:,:,i]*wu; end
        cu=reshape(cu,nb*(nb+1),nb+1);
        if ifleray; cu=cu[:,1:Int(round(nb//2))+1]; end

        au0=.5*(au0+au0');
        bu0=.5*(bu0+bu0');

        au=au0[2:end,2:end];
        bu=bu0[2:end,2:end];

        u0=wu'*u0;

        if ips == "H10"
            sc=diag(au);
        else
            sc=diag(bu);
        end

        u0[2:end]=u0[2:end]./sc;

        ms=Int(readdlm(base*"ops/ns")[1]);
        if ifcopt
            uk_full=readdlm(base*"ops/uk");
            uk_full=reshape(uk_full,mb+1,ms);
            global uk=wu'*uk_full[index,:];
            for i=1:ms
                uk[2:end,i]=uk[2:end,i]./sc;
            end
        end
    end

    u=zeros(nb+1,3);
    ua=zeros(nb+1,1);
    u2a=zeros(nb+1,nb+1);
    ru=zeros(nb,1);
    eu=zeros(nb,3);
    hufac=[];

    if ifrom[2]
        bt0_full=readdlm(base*"ops/bt");
        bt0_full=reshape(bt0_full,mb+1,mb+1);
        bt0=bt0_full[index,index];
        bt0=.5*(bt0+bt0');

        at0_full=readdlm(base*"ops/at");
        at0_full=reshape(at0_full,mb+1,mb+1);
        at0=at0_full[index,index];
        at0=.5*(at0+at0');

        ct_full=readdlm(base*"ops/ct");
        ct_full=reshape(ct_full,mb,mb+1,mb+1);
        ct=ct_full[index1,index,index];
        ct=reshape(ct,nb*(nb+1),nb+1);

        t0_full=readdlm(base*"ops/t0");
        if isfile("ops/t0"); t0_full=readdlm("ops/t0"); end
        t0=t0_full[index];

        global wt=zeros(nb+1,nb+1);
        for i=1:nb+1; wt[i,i]=1.0; end
        if ifaug
            if ips == "H10"
                (d,v)=eigen(at0[2:end,2:end]);
            else
                (d,v)=eigen(bt0[2:end,2:end]);
            end
            wt[2:end,2:end]=v;
        end
        wt1=wt[2:end,2:end];

        at0=wt'*at0*wt;
        bt0=wt'*bt0*wt;

        ct=reshape(wt1'*reshape(ct*wu,nb,(nb+1)*(nb+1)),nb,nb+1,nb+1);
        for i=1:nb+1; ct[:,:,i]=ct[:,:,i]*wt; end
        ct=reshape(ct,nb*(nb+1),nb+1);
        if ifleray; ct=ct[:,1:Int(round(nb//2))+1]; end

        at0=.5*(at0+at0');
        bt0=.5*(bt0+bt0');

        at=at0[2:end,2:end];
        bt=bt0[2:end,2:end];

        t0=wt'*t0;

        if ips == "H10"
            sc=diag(at);
        else
            sc=diag(bt);
        end

        t0[2:end]=t0[2:end]./sc;

        ms=Int(readdlm(base*"ops/ns")[1]);
        if ifcopt
            tk_full=reshape(readdlm(base*"ops/tk"),mb+1,ms);
            global tk=wt'*tk_full[index,:];
            for i=1:ms
                tk[2:end,i]=tk[2:end,i]./sc;
            end
        end
    end

    ta=zeros(nb+1,1);
    t2a=zeros(nb+1,nb+1);
    t=zeros(nb+1,3);
    rt=zeros(nb,1);
    et=zeros(nb,3);
    htfac=[];

    if ifbuoy
        buxt0_full=reshape(readdlm(base*"ops/buxt"),mb+1,mb+1);
        buyt0_full=reshape(readdlm(base*"ops/buyt"),mb+1,mb+1);
        buzt0_full=reshape(readdlm(base*"ops/buzt"),mb+1,mb+1);

        global buxt0=wu'*buxt0_full[index,index]*wt;
        global buyt0=wu'*buyt0_full[index,index]*wt;
        global buzt0=wu'*buzt0_full[index,index]*wt;

        global buxt=buxt0[2:end,:];
        global buyt=buyt0[2:end,:];
        global buzt=buzt0[2:end,:];
    end

    if inus == 1
        nus_full=readdlm(base*"qoi/nus");
        global nus=wt'*nus_full[index];
    end

    if inus == 2
        tsa_full=readdlm(base*"qoi/tsa");
        tsa_full=tsa_full[index];
        global tsa = wt'*tsa_full;

        tbn_full=reshape(readdlm(base*"qoi/tbn"),mb+1,mb+1);
        tbn_full=tbn_full[index,index];
        global tbn = wu'*tbn_full*wt;

        tbd_full=readdlm(base*"qoi/tbd");
        tbd_full=tbd_full[index];
        global tbd = wu'*tbd_full;
    end

    if ifrom[1]; u[:,1]=u0; end
    if ifrom[2]; t[:,1]=t0; end

    println("done loading ...");
end

function ceig(c);
    nb=Int(round(sqrt(size(c,1))));
    mb=size(c,2)-1;
    maxe=0.;
    e0=0.;
    for i=1:mb+1
        cu1=cu[:,i]; cu1=reshape(cu1,nb,nb+1);
        (d,v)=eigen(cu1[:,2:end]);
        if i == 1; e0=maximum(abs.(d)); end
        maxe=max(maxe,maximum(abs.(d)));
    end
    return (maxe,e0)
end

function shift(a,b,n);
    for i=n:-1:2; a[:,i]=a[:,i-1]; end
    a[:,1]=b;
    return a;
end

function setcoef()
    alphas=zeros(3,3);
    betas=zeros(4,3);
    
    alphas[1,1]=  1.0;
    
    alphas[1,2]=  2.0;
    alphas[2,2]= -1.0;
    
    alphas[1,3]=  3.0;
    alphas[2,3]= -3.0;
    alphas[3,3]=  1.0;
    
    betas[1,1]=  1.0;
    betas[2,1]= -1.0;
    
    betas[1,2]=  1.5;
    betas[2,2]= -2.0;
    betas[3,2]=  0.5;
    
    betas[1,3]=  11.0/6;
    betas[2,3]= -3.0;
    betas[3,3]=  1.5;
    betas[4,3]= -1.0/3;

    return (alphas,betas);
end

function step(rhs,a,b,diff,betas,dt,ito,hfac)
    if hfac == [];
        h=b*betas[1,ito]/dt+a*diff;
        h=.5*(h+h');
        hfac=cholesky(h);
    end
    return ([1;hfac\rhs],hfac);
end

function setr(a,b,c,f,u,t,ext,diff,alphas,betas,dt,ito)
    nb=size(u,1)-1;
    rhs=zeros(nb,1);

    ext[:,3]=ext[:,2];
    ext[:,2]=ext[:,1];
    ext[:,1]=ext[:,1]*0;

    if size(c,2) == Int(round(nb//2))+1;
        ext[:,1]=ext[:,1]-reshape(c*u[1:Int(round(nb//2))+1,1],nb,nb+1)*t[:,1];
    else
        ext[:,1]=ext[:,1]-reshape(c*u[:,1],nb,nb+1)*t[:,1];
    end
    ext[:,1]=ext[:,1]-a[2:end,1]*diff;
    ext[:,1]=ext[:,1]+f;

    rhs=rhs+ext*alphas[:,ito];
    rhs=rhs-b*(t[2:end,:]*betas[2:end,ito])/dt;

    return (rhs,ext);
end
