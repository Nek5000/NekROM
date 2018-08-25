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

a = diag(a0(2:n,2:n)) % left-handed side A
b = b0(2:n,2:n) % left-handed side B

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
beta0 = 1;

u0 = zeros(n,1);
e0 = [1;zeros(nb,1)]

u(1:n,1:(nsteps/iostep)) = 0
u(1,1:(nsteps/iostep)) = 1

% do BDF1 to get u1, u2, u3
for istep = 1:3
    helm = (b * beta0 / dt + a / re);

    t = zeros(nb,n);
    for i = 1:n
        t = t + reshape(c0(2:n,:,i),[nb,n]) * u(i,istep);
    end

    rhs = b0(2:n,1:n) * u(:,istep) / dt - a1 / re;
    rhs = rhs - t * u(:,istep); % advection contributions
    rhs = rhs - a0(2:n,1:n) * e0 / re - b0(2:n,1:n) * e0 / dt;
    tmp = helm \ rhs;
    u(2:n,istep+1) = tmp;

    if (mod(istep,iostep) == 0)
        fname = strcat(num2str(istep/iostep),'.out')
        fid = fopen(fname,'w');
        fprintf(fid,'%d\n',u(1:n,istep));
        fclose(fid);
        u(1:n,istep)
    end

end

% BDF3/EXT3 coefficients
beta = [11/6,-3,3/2,-1/3]
alpha = [3,-3,1]

% do BDF3 
for istep = 4:(nsteps/iostep)
    helm = (b * beta(4) / dt + a / re);

    t1 = zeros(nb,n);
    t2 = zeros(nb,n);
    t3 = zeros(nb,n);
    for i = 1:n
        t3 = t3 + reshape(c0(2:n,:,i),[nb,n]) * u(i,istep-3);
        t2 = t2 + reshape(c0(2:n,:,i),[nb,n]) * u(i,istep-2);
        t1 = t1 + reshape(c0(2:n,:,i),[nb,n]) * u(i,istep-1);
    end


    rhs = -beta(1) * b0(2:n,1:n) * u(:,istep-3) / dt ;
    rhs = rhs - beta(2) * b0(2:n,1:n) * u(:,istep-2) / dt ;
    rhs = rhs - beta(3) * b0(2:n,1:n) * u(:,istep-1) / dt ;
    rhs = rhs - a1 / re;

    rhs = rhs - alpha(1) * t3 * u(:,istep-3); % advection contributions
    rhs = rhs - alpha(2) * t2 * u(:,istep-2); % advection contributions
    rhs = rhs - alpha(3) * t1 * u(:,istep-1); % advection contributions
    rhs = rhs - a0(2:n,1:n) * e0 / re - beta(4) * b0(2:n,1:n) * e0 / dt;
    tmp = helm \ rhs;
    u(2:n,istep) = tmp;

    if (mod(istep,iostep) == 0)
        fname = strcat(num2str(istep/iostep),'.out')
        fid = fopen(fname,'w');
        fprintf(fid,'%d\n',u(1:n,istep));
        fclose(fid);
        u(1:n,istep)
    end
end
