% ---------------------------------------------------------------
% This code generates plots for a network of 2000 Duffing oscillator subsystems 
% with a ring topology.
% ---------------------------------------------------------------
clc
clear
close all

%% System setup
beta = -2;
delta = 0.5;
alpha = 0.01;
eta = 0.1;
k = 1;

A = [0 1 0; -beta -delta -alpha];

B = [k 0;0 k];

D = [0 0 0;eta 0  0];  % matrix D_ij

N = 2000; % Number of subsystems

m = 2; % Rows of A_i, D_i
n = 3; % Number of columns of A_i and D_i

% Initial state
x = 1 * (-4 + 8 * rand(2 * N, 1));

%% Build A_network (according to topology) and B_network

A_network = zeros(N * m, N * n);
B_network = [];

for i = 1:N
    B_network = blkdiag(B_network, B);
end

for i = 1:N
    A_network((i - 1) * m + 1:i * m, (i - 1) * n + 1:i * n) = A;
    if i < N
        A_network(i * m + 1:(i + 1) * m, (i - 1) * n + 1:i * n) = D;
    end
end

% Wrap around
A_network(1:m, (N - 1) * n + 1:N * n) = D;

%% Symbolic controller
syms x_1 x_2

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

 S = 1;

Control_input(1,1)  = ...
 S * (    -2.0346*x_1^3 - 13.4461*x_1^2*x_2 - 5.1204*x_1*x_2^2 - 10.3586*x_2^3 + 1.2659*x_1^2 - 21.134*x_1*x_2 + 10.1984*x_2^2 - 1168.2521*x_1 + 144.0447*x_2);

Control_input(2,1)  = ...
 S * (   24.9364*x_1^3 - 6.0576*x_1^2*x_2 + 19.0191*x_1*x_2^2 - 5.2088*x_2^3 + 36.8267*x_1^2 - 24.4002*x_1*x_2 + 2.5832*x_2^2 + 133.4182*x_1 - 778.8311*x_2);


M_x = [x_1; x_2; x_1^3];

% === Fast function conversion ===

controller_fast = matlabFunction(Control_input, 'Vars', {[x_1; x_2]});
monomial_fast = matlabFunction(M_x, 'Vars', {[x_1; x_2]});

%% Simulation setup
tau = 1e-3;
tspan = 2;  % Simulation time
T = tspan / tau;

States = x.';

for i = 1:T
    u = zeros(2*N, 1);
    mono = zeros(3 * N, 1);

    for j = 1:N
        idx = (j - 1) * 2 + 1;
        x_local = x(idx:idx + 1);

     
        u(2*(j-1)+1:2*j) = controller_fast(x_local);

        mono((j - 1) * 3 + 1:3 * j) = monomial_fast(x_local);
    end

    x = x + tau * (A_network * mono + B_network * u);

    States = [States; x.'];
end


%% ===================  Plot 2D states of network ==========================
num_plots = 120;
random_subs = randperm(N, num_plots); % Random subsystems

figure('Color', 'w')
hold on

% === Plot Trajectories ===
for i = 1:num_plots
    random_subsystems = random_subs(i);
    plot(States(:, 2 * random_subsystems - 1), States(:, 2 * random_subsystems), ...
        'LineWidth', 2, 'LineStyle', '-');
end

% === Define Unsafe/initial Regions ===
% Regions of interest 

X0_bound =  [-4 4;-4 4]; % initial set
X1_bound1 = [-10 -6;-10 -5]; % unsafe set 1
X1_bound2 = [6 10;5 10]; % unsafe set 2

% === Function to plot a filled box ===
plotRegion = @(box, color) fill( ...
    [box(1,1) box(1,2) box(1,2) box(1,1)], ...
    [box(2,1) box(2,1) box(2,2) box(2,2)], ...
    color,'EdgeColor', 'none');


% The following functions plot the level sets corresponding to the initial and unsafe sets.

% Regions of interest 
x1_range = linspace(-10, 10, 100);
x2_range = linspace(-10, 10, 100);
[x1, x2] = meshgrid(x1_range, x2_range);

quad_base =   8.4503 .* x1.^2 - 2.0223 .* x1 .* x2 + 5.6554 .* x2.^2;

V1 = 2000 * (quad_base -    281.3336);
V2 = 2000 * (quad_base -  291.3207);


contour(x1, x2, V1, [0 0], 'LineWidth', 1.5, 'LineColor', 'b')

contour(x1, x2, V2, [0 0], 'LineWidth', 1.5, 'LineColor', 'r', 'LineStyle', '--')


% === Plot regions ===
h0 = plotRegion(X0_bound, 'b');    % Initial set: blue
h1 = plotRegion(X1_bound1, 'r');   % Unsafe set 1: red
h2 = plotRegion(X1_bound2, 'r');   % Unsafe set 2: red


regionAlphas = [0.2, 0.5, 0.5];  % [initial, unsafe1, unsafe2]

set(h0, 'FaceAlpha', regionAlphas(1));
set(h1, 'FaceAlpha', regionAlphas(2));
set(h2, 'FaceAlpha', regionAlphas(3));

grid on
axis equal;
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;
xlim([-10 10])
ylim([-10 10])