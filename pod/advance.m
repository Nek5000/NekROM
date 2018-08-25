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

aa = a0(2:n,2:n);
a = diag(a0(2:n,2:n)); % left-handed side A
b = b0(2:n,2:n); % left-handed side B

%c = c0(1:nb,1:nb,1:nb)
a1 = reshape(a0(1,2:n),[nb,1]);
c1 = reshape(c0(2:n,1,2:n),[nb,nb]);
c2 = reshape(c0(2:n,2:n,1),[nb,nb]);
c3 = reshape(c0(2:n,1,1),[nb,1]);

dt = 1e-5;
nsteps = 2e5;
iostep = 1e4;

%nsteps = 2000;
%iostep = 100;

re = 1e3;
beta0 = 1;

u0 = zeros(n,1);
e0 = [1;zeros(nb,1)]

u(1:n,1:(nsteps/iostep)) = 0
u(1,1:(nsteps/iostep)) = 1

% do BDF1 to get u1, u2, u3
for istep = 1:4
    helm = (b * beta0 / dt + a / re);

    t = zeros(nb,n);
    for i = 1:n
        t = t + reshape(c0(2:n,:,i),[nb,n]) * u(i,istep);
    end

    rhs = b0(2:n,1:n) * u(:,istep) / dt - a1 / re;
    rhs = rhs - t * u(:,istep); % advection contributions
    rhs = rhs - a0(2:n,1:n) * e0 / re - b0(2:n,1:n) * e0 / dt;
    tmp = helm \ rhs;
    u(2:n,istep) = tmp;
end

%    U = [u0,u1,u2,u3]
%    if (mod(istep,iostep) == 0)
%        fname = strcat(num2str(istep/iostep),'.out')
%        fid = fopen(fname,'w');
%        fprintf(fid,'%d\n',u);
%        fclose(fid);
%        u
%    end
