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
Re = 200; % Reynolds number
nb  = 20; % number of modes
dt      = 0.001; % Simultaion time step size
nsteps  = 20000; % Simulation total time steps
iostep  = 10; % Every #steps to store ROM quantities

rom = grom('./ops/'); % Create an instance of the grom class and load NekROM operators

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

expected_array = [
    1.000000000000000e+00;
    6.744372524690426e+00;
    1.414966897111988e-01;
    9.927149826323546e-02;
    1.781345058571626e-01;
    4.878748658895765e-01;
    -1.155798970939317e-01;
    9.103068810481321e-02;
    4.967081912919089e-03;
    9.461146705725121e-02;
    -1.515864220108622e-03;
    4.201420332769764e-02;
    7.336015956620934e-03;
    2.969786491032016e-02;
    -5.977978596123081e-03;
    1.320037512282558e-02;
    -5.615288285845466e-04;
    1.630039119415496e-02;
    -2.730687649546992e-03;
    5.574378056426805e-03;
    6.480755825354902e-03
];
tolerance = 1e-10;
assert(all(abs(ucoef(:, end) - expected_array) < tolerance), 'Arrays are not equal within the specified tolerance');

% Display ucoef and expected_array
disp('ucoef:');
disp(ucoef(:, end));
disp('expected_array:');
disp(expected_array);