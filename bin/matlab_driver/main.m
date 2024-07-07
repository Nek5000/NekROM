%#####################################
%
%# Projection Based ROM driver in Matlab for NekROM
%# v0.0.0
%
%# Currently only support Galerkin ROM
%
%# Ping-Hsuan Tsai
%# 2024-07-04
%
%#####################################

clear all; close all;

% Parameters users need to specify
Re = 204; % Reynolds number
nb  = 20; % number of modes
dt      = 0.001; % Simultaion time step size
nsteps  = 20000; % Simulation total time steps
iostep  = 10; % Every #steps to store ROM quantities

rom = grom('../ops/'); % Create an instance of the grom class and load NekROM operators

rom = rom.get_N_dim_ops(nb); % Extract nb-dimensional operators and vectors
rom = rom.initialize_vars(Re); % Initialize variables

timestepper = BDFk_EXTk(nsteps, dt, iostep, nb);
timestepper = timestepper.setup();

ucoef = zeros(rom.nb+1, (timestepper.nsteps/timestepper.iostep));
% Solving ROM using BDFk/EXTk time stepping scheme
timestepper.u(:, 1) = rom.u0;
for istep=1:timestepper.nsteps
    timestepper.ito=min(istep, 3);

    timestepper = timestepper.setrhs(rom, timestepper.u); % Compute the RHS in BDFk/EXTk scheme
    [u_new] = timestepper.advance(rom); % Solve for solution at the next time step
    timestepper.u = timestepper.shift(timestepper.u, u_new, 3);

    rom = timestepper.collect_statistics(rom, u_new); % Compute mean coefficient and mean squared coefficient

    if (mod(istep, timestepper.iostep) == 0)
        ucoef(:, istep/timestepper.iostep) = timestepper.u(:, 1);
     end

end

% Plot the first five modes behavior in time

% Setup time stamp for ROM
% TODO: Make it cleaner
t_rom = linspace(timestepper.dt, timestepper.T_final, timestepper.nsteps/timestepper.iostep);
for i=2:min(rom.nb+1, 6)
    figure(1)
    plot(t_rom, ucoef(i, :), 'r'); hold on
    legend(rom.str(), 'FontSize', 14);
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(['$u_', num2str(i), '$'], 'Interpreter', 'latex', 'FontSize', 14);
    title(['Mode ', num2str(i), ' behavior of ', num2str(rom.nb), '-dimensional ', rom.str()], 'Interpreter', 'latex', 'FontSize', 14);
    saveas(gcf, sprintf('u%d.png', i))
    close(1)
end

dump_rom_statistics(rom, timestepper.nsteps);

figure(1)
plot(rom.uas(1:rom.nb+1), 'k-o'); hold on
plot(rom.ua/timestepper.nsteps, 'r-x'); hold off
legend(rom.str(), 'FOM', 'FontSize', 14);
xlabel('Mode $i$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$u_i$', 'Interpreter', 'latex', 'FontSize', 14);
title(['Averaged coefficients', ' of ', num2str(rom.nb), '-dimensional ', rom.str()], 'Interpreter', 'latex', 'FontSize', 14);
saveas(gcf, sprintf('ua_compare_N%d.png', rom.nb))