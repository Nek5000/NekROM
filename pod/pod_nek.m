close all; format longe; format compact;
load fort.88;
n=sqrt(size(fort,1));
A=reshape(fort,n,n);

[V,D]=eig(A); D=diag(D);

i=1:n;

semilogy(i,D,'ro-')


N=20;
D(1:N)'
V=V(:,1:N); V=reshape(V,n*N,1);

fid = fopen('evectors.dat','w');
fprintf(fid,'%22.13e\n',V);
fclose(fid);
 



