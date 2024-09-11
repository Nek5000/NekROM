#!/usr/bin/env julia

using DelimitedFiles
using Plots; pgfplotsx();
using LaTeXStrings;

#Plots.scalefontsizes(1.15);

f0=readdlm("fom0/nus.dat");
f4=readdlm("fom4/nus.dat");

r0=readdlm("rom0/nus.dat");
r20=readdlm("rom20/nus.dat");

sp=1;

plot(f0[:,2],f0[:,3],lab=L"\mathrm{FOM},\ \epsilon=1.6");
plot!(f4[:,2],f4[:,3],lab=L"\mathrm{FOM},\ \epsilon=2.6");
plot!(r0[:,2],r0[:,3],lab=L"\mathrm{ROM},\ \epsilon=1.6");
plot!(r20[:,2],r20[:,3],lab=L"\mathrm{ROM},\ \epsilon=2.6",legend=:topright);
xlims!((100,200));
xlabel!("Time");
ylabel!("Nu");
ylims!((1.7,2.4));

savefig("nus.pdf");

#Plots.scalefontsizes(1.0/1.15);
