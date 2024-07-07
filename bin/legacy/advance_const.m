clear all; close all;
% advance state using cten, amat, bmat, c1mat, c2mat, avec, cvec

t = dlmread('ops/cten');
n = nthroot(length(t),3);
c0 = reshape(t,n,n,n);

t = dlmread('ops/amat');
n = sqrt(length(t));
a0 = reshape(t,n,n);

t = dlmread('ops/bmat');
n = sqrt(length(t));
b0 = reshape(t,n,n);

sample_min = dlmread('sample_min');
sample_max = dlmread('sample_max');

% actual boundary for coefficient
coef_min = 0.9*sample_min;
coef_max = 0.9*sample_max;

% sample_min(1) is positive
coef_min(1) = 1.1*sample_min(1);

nb = n - 1;

a = a0(2:n,2:n); % left-handed side A
b = b0(2:n,2:n); % left-handed side B

a1 = reshape(a0(1,2:n),[nb,1]);
c1 = reshape(c0(2:n,2:n,1),nb,nb);
c2 = reshape(c0(2:n,2:n,1),nb,nb);
c3 = reshape(c0(2:n,1,1),[nb,1]);

%dt = 1e-5;
dt = 5e-5;
re = 5e2;

nstep = 1e+7;
nstep = 40000;
iostep = 20000;

u0 = dlmread('ic');
e0 = [1;zeros(nb,1)];

u = zeros(n,3);
convec = zeros(nb,3);
u(1,1:3) = 1;
u(:,1) = u0;

% BDFk/EXTk coefficients (k=1:3)

beta = zeros(3,4);
alpha = zeros(3,3);

beta(1,:) = [1,-1,0,0];
beta(2,:) = [3/2,-2,1/2,0];
beta(3,:) = [11/6,-3,3/2,-1/3];

alpha(1,:) = [1,0,0];
alpha(2,:) = [2,-1,0];
alpha(3,:) = [3,-3,1];

par = [1 0.1 0.01];

% time step
for istep = 1:nstep
    count = min(istep,3);
    helm = (b * beta(count,1) / dt + a / re);

    t = zeros(nb,n);

    for i = 1:n
        t = t + reshape(c0(2:n,i,:),[nb,n]) * u(i,1);
    end

    convec(:,3) = convec(:,2);
    convec(:,2) = convec(:,1);
    convec(:,1) = t * u(:,1);

    rhs = -b * u(2:n,:) * beta(count,2:4)' / dt;

    rhs = rhs - convec * alpha(count,:)'; % advection contributions
    rhs = rhs - a0(2:n,1:n) * e0 / re; 

    u(:,3) = u(:,2);
    u(:,2) = u(:,1);
    invhelm = inv(helm);

% constrained optimization
    for j = 1:length(par);
    
%        If B = identity, need line search to determine alpha
%        B = eye(nb);

        B = helm;
        f = f_eval(helm,invhelm,u(2:n,1),rhs,par(j),sample_min,sample_max,b);
        gradf = gradf_eval(helm,u(2:n,1),rhs,par(j),sample_min,sample_max,b); 

%        invv = inv(B);
%        invv = eye(nb);
    
       for k =1:500;
    
          xo = u(2:n,1);
          s = -B\gradf;
%          s = invv*(-gradf);
          u(2:n,1) = u(2:n,1) + s;

    % check so that the solution does not exceed boundary
    
          if ((u(2,1)-sample_max(1)) >= 1e-8) 
             u(2,1) = coef_max(1);
          elseif ((sample_min(1)-u(2,1)) >= 1e-8)
             u(2,1) = coef_min(1);
          end
          for i = 3:n
             if ((u(i,1)-sample_max(i-1)) >= 1e-8) 
                u(i,1) = coef_max(i-1);
             elseif ((sample_min(i-1)-u(i,1) >= 1e-8))
                u(i,1) = coef_min(i-1); 
             end
          end
          go = gradf; gradf = gradf_eval(helm,u(2:n,1),rhs,par(j),sample_min,sample_max,b); 
          y = gradf-go;

          B = B + (y*y')/(y'*s) - B*(s*s')*B/(s'*B*s);

%        try to use inverse formula
%          invv = (eye(nb)-(s*y'./s'*y))*invv*(eye(nb)-(y*s'./s'*y)) + (s*s')./(s'*s);
    
          fo = f; f = f_eval(helm,invhelm,u(2:n,1),rhs,par(j),sample_min,sample_max,b);
    
          ns(k) = norm(s); nk(k) = abs(f-fo); kk(k) = k;
          ngf(k) = sqrt(gradf'*invhelm*gradf);

          [istep k ns(k) ngf(k) nk(k) abs(nk(k)/ns(k))]
          iter(istep) = k;

          if (ngf(k) < 1e-8) % (nk(k) < 1e-8)
             break
          end
    
       end % BFGS_step
    end % par_step

    if (mod(istep,iostep) == 0)
        m = (istep/iostep);
        str = num2str(m,'./matout/%04.f');
        fname = strcat(str,'.out');
        fid = fopen(fname,'w');
        fprintf(fid,'%d\n',u(2:n,1));
        fclose(fid);
   
        plot([2:n],u(2:n,1),'ko','MarkerFaceColor','k')
        hold on
        plot([2:n],sample_min,'r-',[2:n],sample_max,'b-')
        hold on
        plot([2:n],coef_min,'r-o',[2:n],coef_max,'b-o','linewidth',1.4)
        hold off
        strr = sprintf('Coefficients solved with constraint, istep: %d',istep)
        title(strr)
        legend({'coef a','sample\_min','sample\_max','coef\_min','coef\_max'})
        xlabel('modes')
        drawnow
    end
end % time_step
