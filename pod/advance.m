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

nb = n - 1

a = diag(a0(1:nb,1:nb))
b = b0(1:nb,1:nb)
c = c0(1:nb,1:nb,1:nb)

dt = 1e-5;
nsteps = 2e5;
iostep = 1e4;

%nsteps = 2000;
%iostep = 100;

re = 1e3;
beta0 = 1;

%u0 = zeros(n,1);
%
%u = u0
%
%for istep = 1:nsteps
%    helm = (b * beta0 / dt + a0 / re);
%
%    t = zeros(n,1);
%    for i = 1:n
%        t = t + reshape(c0(:,i,:),[n,n]) * u(i);
%    end
%
%    rhs = b * u / dt - a1 / re;
%    rhs = rhs - t * u - c1 * u - c2 * u; % advection contributions
%   rhs = rhs - c1 * u - c2 * u; % advection contributions
%    rhs = rhs - c3; % not in Patera 2017
%    u = helm \ rhs;
%    if (mod(istep,iostep) == 0)
%        fname = strcat(num2str(istep/iostep),'.out')
%        fid = fopen(fname,'w');
%        fprintf(fid,'%d\n',u);
%        fclose(fid);
%        u
%    end
%end
