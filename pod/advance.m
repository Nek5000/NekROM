% advance state using cten, amat, bmat, c1mat, c2mat, avec, cvec

t = dlmread('cten');
n = nthroot(length(t),3);
c0 = reshape(t,n,n,n);

t = dlmread('amat');
n = sqrt(length(t));
a0 = reshape(t,n,n);

t = dlmread('bmat');
n = sqrt(length(t));
b0 = reshape(t,n,n);

nb = n - 1;

a = diag(a0(2:n,2:n)); % left-handed side A
b = b0(2:n,2:n); % left-handed side B

a1 = reshape(a0(1,2:n),[nb,1]);
%c = c0(1:nb,1:nb,1:nb)
%c1 = reshape(c0(2:n,1,1:n),[nb,n]);
%c2 = reshape(c0(2:n,1:n,1),[nb,n]);
%c3 = reshape(c0(1,1,2:n),[nb,1]);

dt = 1e-5;
nsteps = 10;
iostep = 1;
%nsteps = 2e5;
%iostep = 1e4;

%nsteps = 2000;
%iostep = 100;

re = 1e3;

u0 = dlmread('ic');
e0 = [1;zeros(nb,1)];

u = zeros(n,3);
convec = zeros(nb,3);
u(1,1:3) = 1;

% BDF3/EXT3 coefficients

beta = zeros(3,4);
alpha = zeros(3,3);

beta(1,:) = [1,-1,0,0];
beta(2,:) = [3/2,-2,1/2,0];
beta(3,:) = [11/6,-3,3/2,-1/3]

alpha(1,:) = [1,0,0];
alpha(2,:) = [2,-1,0];
alpha(3,:) = [3,-3,1];

% do BDF3 
for istep = 1:(nsteps/iostep)
    count = min(istep,3) ;
    helm = (b * beta(count,1) / dt + a / re);

    t = zeros(nb,n);
    for i = 1:n
        t = t + reshape(c0(2:n,:,i),[nb,n]) * u(i,1);
    end
    convec(:,2) = convec(:,1);
    convec(:,3) = convec(:,2);
    convec(:,1) = t * u(:,1);

    rhs = b0(2:n,2:n) * u(2:n,:) * beta(count,2:4)' / dt;

%   rhs = rhs - convec * alpha(count,:)'; % advection contributions
    rhs = rhs - a0(2:n,1:n) * e0 / re; 
    tmp = helm \ rhs;
    u(:,2) = u(:,1)
    u(:,3) = u(:,2) 

    if (mod(istep,iostep) == 0)
        m = (istep/iostep)
        str = num2str(m,'%03.f')
        fname = strcat(str,'.out')
        fid = fopen(fname,'w');
        fprintf(fid,'%d\n',u(2:n,1));
        fclose(fid);
    end
end
