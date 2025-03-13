%#######################################################
%
%# Matlab Driver for Galerkin-based reduced order model 
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-09-05
%
%#######################################################


%% User Parameters
% Problem parmeters:
%   nsteps : int
%     The total number of time steps to integrate the ROM ODEs
%   dt     : float
%     Time step size used for integrating the ROM ODEs, could be larger than the FOM dt
%   iostep : int
%     Every iostep you want to store the ROM coefficients
%   nu     : float
%     Viscosity
%   nb     : int 
%     Size of your ROM
%
% Stabilization methods:
%   ifcopt : bool 
%     Set to true if you want to use constrained ROM (C-ROM)
%   ifleray : bool 
%     Set to true if you want to use Leray-ROM (L-ROM)
%   ifefr : bool 
%     Set to true if you want to use evolve-filter-relax ROM (EFR-ROM)
%   iftr : bool 
%     Set to true if you want to use time-relaxataion ROM (TR-ROM)

clear all; close all;

nsteps = 1.000E+05;%20000; 
dt     = 1.000000E-03;%0.001;
iostep = 250;%500;%10;
nu     = 1./15000;%0.01;
nb     = 20;%40;


ifcopt  = false;
ifleray = false;
ifefr   = false;
iftr    = true;


if ifcopt
    reg_str = 'CROM';
elseif ifleray
reg_str = 'Leray';
elseif ifefr
    reg_str = 'EFR';
elseif iftr;
    reg_str = 'TR';
else
    reg_str = 'GROM';
end;   

if (ifleray)
   radius = 0.01; % Feel free to change it
   if exist('radius', 'var')
      fprintf('The filter radius used in the L-ROM is: %f\n', radius);
   end
elseif (ifefr) || (iftr)
   radius = 0.01; % Feel free to change it
   relax = dt; % Feel free to change it
   if exist('radius', 'var') && exist('relax', 'var')
      fprintf('The filter radius and relaxation used in the ROM are: %f %f\n', radius, relax);
   end
end

%deims= [1 2 3 4 5 6 7 8 9 10]
%deims= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]
algo='sopt';%"gpode";
oversample_factor=2;
deims= [2,4,8,16,32,64,128,256,512,1024,2000];%[2,4,8,16,20 40 80 100 200 400 800 1000]
% number of deim points
for j=1:length(deims)
deim_pts = deims(j);


[au_full, bu_full, cu_full, u0_full, uk_full, mb, ns] = load_full_ops('../ops');

%size(au_full)

[au, a0, bu, cu, c0, c1, c2, c3, u0, uk, ukmin, ukmax] = get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb);

%size(uk_full)
%size(uk)
%exit;

%size(au)
%exit;

% Initialize variables
time   = 0.;
rhs    = zeros(nb,1);
ext    = zeros(nb,3);
hufac  = [];
ucoef  = zeros((nsteps/iostep),nb+1);
if (ifleray) || (ifefr) || (iftr)
   dfHfac = [];
   [dfHfac] = set_df(au, bu, radius, 1, dfHfac);
end

%% Get the grid and POD bases for plotting purposes
cname='../basldc';
%avg_cname='../avgcyl';
bas_snaps = NekSnaps(cname);
[pod_u, pod_v, x_fom, y_fom] = get_grid_and_pod(bas_snaps);

%% Get the non-linear snapshots and calculate the DEIM points
nl_cname = '../csnldc';
nl_snaps = NekSnaps(nl_cname);

% Begin integrate ROM with BDFk/EXTk
u     = zeros(nb+1,3); % vectors for BDF3/EXT3
u(:,1)=u0;
[alphas, betas] = setcoef();

coef_errs_all = [];
coef_errs_single = [];
for deim_pts_indx=1:size(deims,2);
deim_pts = deims(deim_pts_indx);
tensor_coefs = [];
deim_coefs = [];
for k = 1:size(uk,2)
 u = uk(:,k);
 utmp = u;
 %size(uk)
 %size(c3)
 %exit
 c_coef_tensor = conv_deim_fixed(u(:,1), pod_u, pod_v, nl_snaps,deim_pts,k,algo,oversample_factor)-c1+c2*utmp(:,1)+c3*utmp(:,1); 
 c_coef_deim =  (reshape(c0*utmp(:,1),nb,nb+1)*u(:,1));
 tensor_coefs = [tensor_coefs, c_coef_tensor];
 deim_coefs = [deim_coefs, c_coef_deim];
end
coef_err_vec = [];
for coef = 1:size(tensor_coefs,1);    
    coef_err_vec = [coef_err_vec; norm(tensor_coefs(coef,:) - deim_coefs(coef,:))/norm(tensor_coefs(coef,:))];
end;
coef_errs_single = [coef_errs_single, coef_err_vec];
coef_errs_all = [coef_errs_all, norm(tensor_coefs - deim_coefs)/norm(tensor_coefs)];
end;

%hold on;
coef_errs_all
coef_errs_single
%for coef = 1:size(tensor_coefs,1);
%   semilogx(deims, coef_errs_single(:,coef));
%end;

semilogx(deims, coef_errs_all);
xlabel("Number of hyperreduction points");
ylabel("Coefficient two-norm relative error");
title(sprintf('LDC with %s',algo));
pause(60);
exit;

for istep=1:nsteps
   istep
   ito=min(istep,3);
   if istep<= 3
      hufac = [];
   end

   % Compute the right-handed side of the fully discretized system associated to BDFk/EXTk
   rhs = zeros(nb,1);

   ext(:,3)=ext(:,2);
   ext(:,2)=ext(:,1);
   ext(:,1)=ext(:,1)*0;

   if ifleray
      % Compute filtered u
      utmp = [1;(dfHfac\(dfHfac'\u(2:end,1)))];
   else 
      utmp = u;
   end

   %% Compare the NL evaluation results.

   % ROM convection tensor version
   if deim_pts == 0;
    c_coef = (reshape(c0*utmp(:,1),nb,nb+1)*u(:,1));
%  c_coef' %ext(:,1)=ext(:,1)-reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);

   % Pseudo ROM version
   % This should be identical to the above.
%  c_coef = (conv_fom(u(:,1), pod_u, pod_v, bas_snaps))
   %norm(pod_u(:,1:nb+1)*c_coef)
   else;
   % DEIM version
   % Should be close, but not identical to, the above
    c_coef = conv_deim_fixed(u(:,1), pod_u, pod_v, nl_snaps,deim_pts,istep)-c1+c2*utmp(:,1)+c3*utmp(:,1);
   end;



   %norm(pod_u(:,1:nb+1)*c_coef)
   ext(:,1)=ext(:,1)-c_coef;

   %% End result comparison

   ext(:,1)=ext(:,1)-nu*a0;

   if iftr
      % Compute filtered u
      utmp = [1;(dfHfac\(dfHfac'\u(2:end,1)))];
      ext(:,1) = ext(:,1)-relax*(u(2:end,1)-utmp(2:end));
   end

   rhs=rhs+ext*alphas(:,ito);
   rhs=rhs-bu*(u(2:end,:)*betas(2:end,ito))/dt;

   % Solve the linear system to get next step solution
   if (ifcopt)
%     if any(u(2:end,1) > ukmax(2:end)) || any(u(2:end,1) < ukmin(2:end))
         %fprintf('in constrained %d \n',istep);
         options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
         [x,fval,exitflag,output] = fmincon(@(x)rom_residual(x,au,bu,nu,betas,dt,ito,rhs),u(2:end,1),[],[],[],[],ukmin(2:end),ukmax(2:end),[],options);
         u_new = [1,x'];
%     else
%  nsteps = 80000;%1.25000E+05;%20000; 
dt     = 1.000000E-03;%0.001;
iostep = 40;%250;%500;%10;
      % if constraints are satisfied, do normal solve
%        if isempty(hufac)
%           h=bu*betas(1,ito)/dt+au*nu;
%           hfac=chol(h);
%        end
%        u_new = [1,(hfac\(hfac'\rhs))'];
%     end
   else
      if isempty(hufac)
         h=bu*betas(1,ito)/dt+au*nu;
         hfac=chol(h);
      end
      u_new = [1,(hfac\(hfac'\rhs))'];
   end

   % EFR-ROM additional step
   if (ifefr)
      utmp = u_new;
      utmp = [1,(dfHfac\(dfHfac'\u_new(1,2:end)'))'];
      u_new = (1-relax)*u_new + relax*utmp;
   end
        
   u = shift(u,u_new,3);
   time = time+dt;

   if any(isnan(u(:,1)));
      break;
   end;

   bool_plot = false;
   if (mod(istep,iostep) == 0);
      ucoef(istep/iostep,:)=u(:,1);
      u(:,1)

      if bool_plot;      
      u_proj = pod_u(:,1:nb+1)*u(1:end,1);
      v_proj = pod_v(:,1:nb+1)*u(1:end,1);
     
      u_abs = sqrt(u_proj.^2 + v_proj.^2);
      %norm(u_abs)
      hold off;
      patch_plot(x_fom,y_fom, reshape(u_abs,size(x_fom)), [], 'PlotType', 'surface');
      pause(0.01);
     end;
   end

end

%size(ucoef)
%ucoef
%exit;

%plot(ucoef(:,2)); hold on;
%plot(flip(uk(2,:)));
%pause(60);

if deim_pts > 0;
    casedir= sprintf('nb%d_results_deim_pts%d_%s',nb,deim_pts,reg_str)
else
    casedir= sprintf('nb%d_results_%s',nb,reg_str)
end;
mkdir(casedir);

% Dump out ucoef in casedir
fileID = fopen(casedir+"/ucoef",'w');
fprintf(fileID,"%24.15e\n",ucoef);
fclose(fileID);
end

%#####################################
%
%# Auxiliary functions
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-09-05
%
%#####################################

function [a0_full, b0_full, c0_full, u0_full, uk_full, mb, ms] =  load_full_ops(path)
% load_full_ops function loads the full ROM operators stored in the specified path
% specified by input variable path.
%
% Output:
% a0_full : Full stiffness matrix of size mb+1 x mb+1 (The +1 comes from the zeroth mode)
% b-1_full : Full   mass    matrix of size mb+1 x mb+1
% cu_full : Full advection tensor of size mb x mb+1 x mb+1
% u0_full : Vector of size mb+1, which contains the ROM coefficients of the
%           projection of initial conditons onto mb-reduced space
% uk_full : Matrix of size mb+1 x ns. Each column contains the ROM coefficients
%           of the projection of the snapshot onto mb-reduced space
% mb : total number of modes
% ns : nubmer of snapshots used to create those operators
% 
% This function can take times to load operators with nb >= 300.

   fprintf('Loading ROM operators and vectors... \n');
   fprintf('Currently only support velocity... \n');

   mb=dlmread(fullfile(path,"nb"));

   % load stiffness matrix
   a0_full = dlmread(fullfile(path,"au"));
   a0_full = reshape(a0_full,mb+1,mb+1);

   % load mass matrix
   b0_full = dlmread(fullfile(path,"bu"));
   b0_full = reshape(b0_full,mb+1,mb+1);

   % load advection tensor
   c0_full = dlmread(fullfile(path,"cu"));
   c0_full = reshape(c0_full,mb,mb+1,mb+1);

   u0_full = dlmread(fullfile(path,"u0"));

   ms = dlmread(fullfile(path,"ns"));
   uk_full = dlmread(fullfile(path,"uk"));
   uk_full = reshape(uk_full,mb+1,ms);

   fprintf("done loading ... \n");
   
end

function [a, a0, b, c, c0, c1, c2, c3, u0, uk, ukmin, ukmax] =  get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb)
   index  = [1:nb+1];
   index1 = [1:nb];
   index2 = [2:nb+1];
   a      = au_full(index2,index2);
   a0     = au_full(index2,1);
   b      = bu_full(index2,index2);

   cutmp  = cu_full(index1,index,index);
   c0      = reshape(cutmp,nb*(nb+1),nb+1);
   c1 = cutmp(:,1,1);                                                       
   c2 = reshape(cutmp(:,1,:),nb,nb+1);                                      
   c3 = reshape(cutmp(:,:,1),nb,nb+1); 

   c      = cu_full(index1,index2,index2);
   c      = reshape(c,nb*(nb),nb);

   u0      = u0_full(index);

   uk = uk_full(index,:);
   ukmin = min(uk,[],2);
   ukmax = max(uk,[],2);
end

function [hfac] = set_df(a,b,dfRadius,dfOrder,hfac)
% Construct the mth order differential filter (df) 
% I + (\delta^2 B^{-1}A)^m
% with radius delta
   if isempty(hfac)
      bfac = chol(b);
      h = (dfRadius^2)*(bfac\(bfac'\a));
      for i=2:dfOrder
         h = h*((dfRadius^2)*(bfac\(bfac'\a)));
      end
      h = h + eye(size(a));
      hfac=chol(h);
   end
end

function [alphas, betas] = setcoef()
% Setup BFDk/EXTk coefficients
   alphas=zeros(3,3);
   betas=zeros(4,3);

   alphas(1,1)=  1.0;

   alphas(1,2)=  2.0;
   alphas(2,2)= -1.0;

   alphas(1,3)=  3.0;
   alphas(2,3)= -3.0;
   alphas(3,3)=  1.0;

   betas(1,1)=  1.0;
   betas(2,1)= -1.0;

   betas(1,2)=  1.5;
   betas(2,2)= -2.0;
   betas(3,2)=  0.5;

   betas(1,3)=  11.0/6;
   betas(2,3)= -3.0;
   betas(3,3)=  1.5;
   betas(4,3)= -1.0/3;
end

function a = shift(a,b,n)
   for i=n:-1:2; 
      a(:,i)=a(:,i-1); 
   end
   a(:,1)=b;
end

function F = rom_residual(x,a,b,diff,betas,dt,ito,rhs)                                                                                                                                                                                            
   h=b*betas(1,ito)/dt+a*diff;
   F1 = h*x-rhs;
   F = norm(F1);
end

function [pod_u, pod_v, x_fom, y_fom] = get_grid_and_pod(snaps)
  x_fom = snaps.flds{1}.x;
  y_fom = snaps.flds{1}.y;

  [nr, ns, nE] = size(x_fom);
  nL = nr*ns*nE;
  [nbasis, nbasis1] = size(snaps.flds)
  pod_u = [];
  pod_v = [];
  for i=1:nbasis;
    %u_snap = snaps.flds{i}.u;
    pod_u = [pod_u, reshape(snaps.flds{i}.u, nL,1)];
    pod_v = [pod_v, reshape(snaps.flds{i}.v,nL,1)];
  end;


  % The bases are not 
  %[pod_u;pod_v]'*[pod_u;pod_v]
  %exit;

  %avg_snaps = NekSnaps(avg_cname);
  %u_avg = reshape(avg_snaps.flds{1}.u,nL,1);
  %v_avg = reshape(avg_snaps.flds{1}.v,nL,1);
  %% Technically not the POD anymore
  %exit;
  %pod_u(:,1) = u_avg(:,1);
  %pod_v(:,1) = v_avg(:,1);
end;

% This needs to do the same thing as
% reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
function [out_coef] = conv_fom(ucoef, pod_u, pod_v, snaps)
    x=snaps.flds{1}.x;
    y=snaps.flds{1}.y;
    nx1 = size(x,1);
    [zi, w] = zwgll(nx1-1);
    d = deriv_mat(zi);
    [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
    lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
    nL = prod(size(x));
    Me = reshape(jac.*(w*w'),nL,1);

    if false; % Normal way of calculating the gradient
        [ux_fom, uy_fom] = lgrad(u_fom, 0);
        [vx_fom, vy_fom] = lgrad(v_fom, 0);
    else
        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = [];
        uy_pods = [];
        vx_pods = [];
        vy_pods = [];
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods = [ux_pods, reshape(ux_pod, nL,1)];
            uy_pods = [uy_pods, reshape(uy_pod, nL,1)];
            vx_pods = [vx_pods, reshape(vx_pod, nL,1)];
            vy_pods = [vy_pods, reshape(vy_pod, nL,1)];
        end;
        ux_fom = reshape(ux_pods*ucoef, size(x));
        uy_fom = reshape(uy_pods*ucoef, size(x));
        vx_fom = reshape(vx_pods*ucoef, size(x));
        vy_fom = reshape(vy_pods*ucoef, size(x));
    end

    u_fom = reshape(Me.*pod_u*ucoef, size(x));
    v_fom = reshape(Me.*pod_v*ucoef, size(x));

    conv_u_fom = reshape(u_fom.*ux_fom + v_fom.*uy_fom, nL,1);
    conv_v_fom = reshape(u_fom.*vx_fom + v_fom.*vy_fom, nL,1);
    
    out_coef = [pod_u(:,2:end); pod_v(:,2:end)]'*[conv_u_fom; conv_v_fom];
end

% This needs to do the same thing as
% reshape(cu*utmp(:,1),nb,nb+1)*u(:,1);
function [out_coef] = conv_deim(ucoef, pod_u, pod_v, nl_snaps_obj)
    persistent proj_mat u_deimu v_deimu u_deimv v_deimv ux_deimu uy_deimu vx_deimv vy_deimv;
    persistent u_deim_stack v_deim_stack ux_deim_stack uy_deim_stack

    if isempty(proj_mat)
        % Stuff to be precomputed

        % Read in the snapshots.     
        [nl_u_snaps, nl_v_snaps] = get_grid_and_pod(nl_snaps_obj);
        nl_snaps = [nl_u_snaps; nl_v_snaps];

        x=nl_snaps_obj.flds{1}.x;
        y=nl_snaps_obj.flds{1}.y;
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        nL = prod(size(x));
        Me = reshape(jac.*(w*w'),nL,1);
        
        %% Calculate the POD of the NL snapshots
        gramian = nl_snaps'*diag([Me;Me])*nl_snaps;
        gramian = 0.5*(gramian + gramian');
        [nl_eigvecs, nl_eigvals] = eig(gramian);
        [nl_eigvals_sorted, sort_inds] = sort(diag(abs(nl_eigvals)),'descend');
        nl_eigvecs = nl_eigvecs(:,sort_inds);
        nl_bas = nl_snaps*nl_eigvecs; % Using all of the modes 
        nl_bas = nl_bas(:,1:500);
        % Maybe 500 snapshots is not enough?
        nl_bas = orth(nl_bas,1e-16); % Use all of the snapshots for now.
        %nl_bas = orth(nl_bas(:,1:500),1e-16);
        %for i = 1:size(nl_bas,2)
        %    nl_bas(:,i) = nl_bas(:,i)/norm(nl_bas(:,i));
        %end

        separate=false

        if separate % Calculate the QDEIM points separately for u and v
            % This probably won't work.
            [P_u, nl_u_inds] = calc_qdeim_proj_mat(nl_u_snaps);
            [P_v, nl_v_inds] = calc_qdeim_proj_mat(nl_v_snaps);
            n_deim_points_u = ceil(size(nl_bas,2)/2);
            n_deim_points_v = floor(size(nl_bas,2)/2);
            nl_u_inds = nl_u_inds(1,1:n_deim_points_u); 
            nl_v_inds = nl_v_inds(1,1:n_deim_points_v); 

            inds = [nl_u_inds, nl_v_inds + size(nl_u_snaps,1)];
        else % Use the same QDEIM points for u and v.
            % Maybe this isn't right for vector quantities.
            [P, inds] = calc_qdeim_proj_mat(nl_snaps);
            % For testing, use all rows
            %inds = [1:size(nl_snaps,1)];
            %divider = size(nl_u_snaps,1);
            %nl_u_inds = inds(inds <= divider);
            %nl_v_inds = inds(inds > divider);
            %nl_u_inds
            %exit
            %inds = [nl_u_inds, nl_v_inds];
        end;
        inds = inds(:,1:500);

       
        %size(nl_snaps) 
        %size(nl_bas)
        %size(inds) 
        % Kento's ROM approach. Calculate the gradients of the POD modes
        
        ux_pods = [];
        uy_pods = [];
        vx_pods = [];
        vy_pods = [];
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods = [ux_pods, reshape(ux_pod, nL,1)];
            uy_pods = [uy_pods, reshape(uy_pod, nL,1)];
            vx_pods = [vx_pods, reshape(vx_pod, nL,1)];
            vy_pods = [vy_pods, reshape(vy_pod, nL,1)];
        end;

      
        % Is this problematic? What does it mean to integrate on the DEIM points?
        u_deim = Me.*pod_u;
        v_deim = Me.*pod_v;
%       u_deim = pod_u;
%       v_deim = pod_v;
    
        if separate
            u_deimu = u_deim(nl_u_inds,:);
            v_deimu = v_deim(nl_u_inds,:);
            u_deimv = u_deim(nl_v_inds,:);
            v_deimv = v_deim(nl_v_inds,:);

            ux_deimu = ux_pods(nl_u_inds,:);
            uy_deimu = uy_pods(nl_u_inds,:);
            vx_deimv = vx_pods(nl_v_inds,:);
            vy_deimv = vy_pods(nl_v_inds,:);
        else
            u_deim_stack = [u_deim; u_deim];
            u_deim_stack = u_deim_stack(inds,:);
            v_deim_stack = [v_deim; v_deim];
            v_deim_stack = v_deim_stack(inds,:);


            ux_deim_stack = [ux_pods; vx_pods];
            ux_deim_stack = ux_deim_stack(inds,:);
            uy_deim_stack = [uy_pods; vy_pods];
            uy_deim_stack = uy_deim_stack(inds,:);
        end;

        proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas(inds,:)); 
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*([Me;Me].*nl_bas)*inv(nl_bas(inds,:)); 
        % For testing
        %proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]';
    end;

    separate=false;

    if separate
        conv_u_deim = ((u_deimu*ucoef).*(ux_deimu*ucoef) + (v_deimu*ucoef).*(uy_deimu*ucoef));
        conv_v_deim = ((u_deimv*ucoef).*(vx_deimv*ucoef) + (v_deimv*ucoef).*(vy_deimv*ucoef)); 
        out_coef = proj_mat*[conv_u_deim; conv_v_deim];
    else
        conv_deim = ((u_deim_stack*ucoef).*(ux_deim_stack*ucoef) + (v_deim_stack*ucoef).*(uy_deim_stack*ucoef));
        out_coef = proj_mat*conv_deim;
    end;
end


function[P,indices] = calc_qdeim_proj_mat(U_nl)
    [A,B,P] = qr(U_nl');
    P = P(:,1:size(U_nl,2));
    indices = [];
    for col=1:size(P,2);
        [val, ind] = max(P(:,col));
        indices = [indices, ind];
    end;
end


function[P,indices] = calc_deim_proj_mat(U_nl)
    [maxval, ind] = max(abs(U_nl(:,1)), [], 1);
    [n,m] = size(U_nl);
    P = zeros(n,m);
    indices = zeros(1,m);
    P(ind,1) = 1;
    indices(1,1) = ind;
    for i=2:m;
    % These two ways of calculating c should give the same result
        %c = (P(:,1:i-1)'*U_nl(:,1:i-1)) \ P(:,1:i-1)'*U_nl(:,i)
        c = U_nl(indices(1:i-1),1:i-1) \ U_nl(indices(1:i-1),i);
    r = U_nl(:,i) - U_nl(:,1:i-1)*c;
    [maxval, ind] = max(abs(r), [], 1);
    indices(1,i) = ind;
    P(ind,i) = 1;
    end;
end


function [out_coef] = conv_deim_fixed(ucoef, pod_u, pod_v, nl_snaps_obj, deim_pts,istep,algo,oversample_factor)
    persistent proj_mat u_deimu v_deimu u_deimv v_deimv ux_deimu uy_deimu vx_deimv vy_deimv;
    persistent u_deim_stack v_deim_stack ux_deim_stack uy_deim_stack
    persistent inds

%   if isempty(proj_mat)
    if istep == 1
        % Stuff to be precomputed

        % Read in the snapshots.     
        [nl_u_snaps, nl_v_snaps] = get_grid_and_pod(nl_snaps_obj);
        nl_snaps = [nl_u_snaps; nl_v_snaps];

        x=nl_snaps_obj.flds{1}.x;
        y=nl_snaps_obj.flds{1}.y;
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        nL = prod(size(x));
        Me = reshape(jac.*(w*w'),nL,1);
        
        %% Calculate the POD of the NL snapshots
        gramian = nl_snaps'*diag([Me;Me])*nl_snaps;
        gramian = 0.5*(gramian + gramian');
        [nl_eigvecs, nl_eigvals] = eig(gramian);
        [nl_eigvals_sorted, sort_inds] = sort(diag(abs(nl_eigvals)),'descend');
        fileID = fopen('nl_eigvals_sorted.txt', 'w');
        fprintf(fileID, '%24.15e\n', nl_eigvals_sorted);
        fclose(fileID);
        nl_eigvecs = nl_eigvecs(:,sort_inds);
        nl_bas = nl_snaps*nl_eigvecs; % Using all of the modes 
        nl_bas = nl_bas(:,1:deim_pts);
        % Maybe 500 snapshots is not enough?
        nl_bas = orth(nl_bas,1e-16); % Use all of the snapshots for now.


        nl_bas_u = nl_bas(1:nL,:);
        nl_bas_v = nl_bas(nL+1:end,:);
        nl_bas_u = orth(nl_bas_u,1e-16); % Use all of the snapshots for now.
        nl_bas_v = orth(nl_bas_v,1e-16); % Use all of the snapshots for now.
        %nl_bas = orth(nl_bas(:,1:500),1e-16);
        %for i = 1:size(nl_bas,2)
        %    nl_bas(:,i) = nl_bas(:,i)/norm(nl_bas(:,i));
        %end

        separate=false

        if separate % Calculate the QDEIM points separately for u and v
            % This probably won't work.
%           [P_u, nl_u_inds] = calc_qdeim_proj_mat(nl_u_snaps);
%           [P_v, nl_v_inds] = calc_qdeim_proj_mat(nl_v_snaps);
            [P_u, nl_u_inds] = calc_qdeim_proj_mat(nl_bas_u);
            [P_v, nl_v_inds] = calc_qdeim_proj_mat(nl_bas_v); % not sure why the separate approach is not working
            % maybe because this will not satisfiyes divergence-free constraints...
            n_deim_points_u = ceil(size(nl_bas,2)/2);
            n_deim_points_v = floor(size(nl_bas,2)/2);
            nl_u_inds = nl_u_inds(1,1:n_deim_points_u); 
            nl_v_inds = nl_v_inds(1,1:n_deim_points_v); 

            inds = [nl_u_inds, nl_v_inds + size(nl_u_snaps,1)];
        else % Use the same QDEIM points for u and v.
            % Maybe this isn't right for vector quantities.
%           [P, inds] = calc_qdeim_proj_mat(nl_snaps);
            if strcmp(algo, "sopt");
            % Can oversample if desired.
               inds = s_opt_generator(nl_bas, deim_pts + 0, [])
               inds = inds';
            elseif strcmp(algo, "gpode");
               inds = gpode(nl_bas, oversample_factor*size(nl_bas,2));
               inds = inds';
            else;
               [P, inds] = calc_qdeim_proj_mat(nl_bas); % Should select points based on the basis
            end;
            % For testing, use all rows
            %inds = [1:size(nl_snaps,1)];
            %divider = size(nl_u_snaps,1);
            %nl_u_inds = inds(inds <= divider);
            %nl_v_inds = inds(inds > divider);
            %nl_u_inds
            %exit
            %inds = [nl_u_inds, nl_v_inds];
        end;
        inds = inds(:,1:deim_pts);

       
        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = [];
        uy_pods = [];
        vx_pods = [];
        vy_pods = [];
        for i = 1:size(pod_u,2);
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
            ux_pods = [ux_pods, reshape(ux_pod, nL,1)];
            uy_pods = [uy_pods, reshape(uy_pod, nL,1)];
            vx_pods = [vx_pods, reshape(vx_pod, nL,1)];
            vy_pods = [vy_pods, reshape(vy_pod, nL,1)];
        end;

      
        % Is this problematic? What does it mean to integrate on the DEIM points?
        u_deim = Me.*pod_u;
        v_deim = Me.*pod_v;
%       u_deim = pod_u;
%       v_deim = pod_v;
    
        if separate
            u_deimu = u_deim(nl_u_inds,:);
            v_deimu = v_deim(nl_v_inds,:);
            u_deimv = u_deim(nl_u_inds,:);
            v_deimv = v_deim(nl_v_inds,:);

            ux_deimu = ux_pods(nl_u_inds,:);
            uy_deimu = uy_pods(nl_u_inds,:);
            vx_deimv = vx_pods(nl_v_inds,:);
            vy_deimv = vy_pods(nl_v_inds,:);
%           proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*[nl_bas_u; nl_bas_v]*inv([nl_bas_u(nl_u_inds,:); nl_bas_v(nl_v_inds,:)]); 
        else
            u_deim_stack = [u_deim; u_deim];
            u_deim_stack = u_deim_stack(inds,:);
            v_deim_stack = [v_deim; v_deim];
            v_deim_stack = v_deim_stack(inds,:);


            ux_deim_stack = [ux_pods; vx_pods];
            ux_deim_stack = ux_deim_stack(inds,:);
            uy_deim_stack = [uy_pods; vy_pods];
            uy_deim_stack = uy_deim_stack(inds,:);
%           proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas(inds,:)); 
        end;

        if size(inds,2) == size(nl_bas,2)
            proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas(inds,:));
        else
            proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*pinv(nl_bas(inds,:)); 
        end;
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*nl_bas*inv(nl_bas);
%       proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]'*([Me;Me].*nl_bas)*inv(nl_bas(inds,:)); 
        % For testing
        %proj_mat = [pod_u(:,2:end); pod_v(:,2:end)]';
    end;

    separate=false;

    if separate
        conv_u_deim = ((u_deimu(:,2:end)*ucoef(2:end)).*(ux_deimu(:,2:end)*ucoef(2:end)) + (v_deimu(:,2:end)*ucoef(2:end)).*(uy_deimu(:,2:end)*ucoef(2:end)));
        conv_v_deim = ((u_deimv(:,2:end)*ucoef(2:end)).*(vx_deimv(:,2:end)*ucoef(2:end)) + (v_deimv(:,2:end)*ucoef(2:end)).*(vy_deimv(:,2:end)*ucoef(2:end))); 
        out_coef = proj_mat*[conv_u_deim; conv_v_deim];
    else
        size(u_deim_stack)
        size(v_deim_stack)
        size(ucoef)
        conv_deim = ((u_deim_stack(:,2:end)*ucoef(2:end)).*(ux_deim_stack(:,2:end)*ucoef(2:end)) + (v_deim_stack(:,2:end)*ucoef(2:end)).*(uy_deim_stack(:,2:end)*ucoef(2:end)));
        out_coef = proj_mat*conv_deim;
    end;
end

