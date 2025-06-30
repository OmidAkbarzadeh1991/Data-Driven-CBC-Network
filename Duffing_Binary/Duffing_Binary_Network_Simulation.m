% ---------------------------------------------------------------
% This code generates plots for a network of 1023 Duffing oscillator subsystems 
% with a binary topology.
% ---------------------------------------------------------------
clc
clear
close all

%% System setup
BETA = -2;
DELTA = 0.5;
ALPHA = 0.01;
ETA = 0.05;
k_C = 1;

A = [0 1 0; -BETA -DELTA -ALPHA];

B = [k_C 0;0 k_C];

D = [0 0 0; ETA 0 0]; % matrix D_ij

num_states = 2; % Number of states per subsystem

M = 3; % Number of columns of A_i and D_i

N = 1023; % Number of subsystems

x = 1 * (-4 +  8 * rand(2 * N, 1)); % Initial state

%% Build A_network (according to topology) and B_network

B_network = [];

for i = 1 : N

    B_network = blkdiag(B_network, B);

end

A_network = zeros(num_states*N, M*N);
for j = 1 : N
    for i = 1 : N
        if i == j
            A_network(num_states * i - 1 : num_states * i, M * j - 2 : M * j) = A;
        elseif i == 2 * j
            if i + 1 > N
                A_network(num_states * i - 1 : num_states * i, M * j - 2 : M * j) = D;
            else
                A_network(num_states * i - 1 : num_states * i + num_states, M * j - 2 : M * j) = [D; D];
            end
        end
    end
end
%% Symbolic controller
syms x_1 x_2

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

S = 1;

Control_input(1,1)  = ...
 S * (   -1.7313*x_1^3 - 14.8234*x_1^2*x_2 - 4.6827*x_1*x_2^2 - 12.1856*x_2^3  + 1.0959*x_1^2 - 18.9204*x_1*x_2 + 2.9079*x_2^2 - 640.2218*x_1 + 91.5416 *x_2  );

Control_input(2,1)  = ...
 S * (   24.1086*x_1^3 - 6.3893*x_1^2*x_2 + 19.748*x_1*x_2^2 - 5.6392*x_2^3 + 31.9524*x_1^2 - 13.2464*x_1*x_2 + 1.8457*x_2^2 + 66.1871*x_1 - 474.8328*x_2);


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
    u = zeros(2 * N, 1);
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

%% ===================  Plot 2D states of network =========================

num_plots = 120;
random_subs = randperm(N, num_plots); % Random subsystems

figure('Color', 'w')
hold on

% === Plot Trajectories ===
for i = 1:num_plots
    random_subsystems = random_subs(i);
    plot(States(:, 2 * random_subsystems - 1), States(:, 2 * random_subsystems), ...
        'LineWidth', 1.5, 'LineStyle', '-');
end

% === Define Unsafe/initial Regions ===
% Regions of interest 

X0_bound =  [-4 4;-4 4];  % initial set
X1_bound1 = [-10 -6;-10 -5];  % unsafe set 1
X1_bound2 = [6 10;5 10];  % unsafe set 2

% === Function to plot a filled box ===
plotRegion = @(box, color) fill( ...
    [box(1,1) box(1,2) box(1,2) box(1,1)], ...
    [box(2,1) box(2,1) box(2,2) box(2,2)], ...
    color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');


% The following functions plot the level sets corresponding to the initial and unsafe sets.
% === Add Both Contours ===
x1_range = linspace(-10, 10, 100);
x2_range = linspace(-10, 10, 100);
[x1, x2] = meshgrid(x1_range, x2_range);

quad_base = 8.0675 * x1.^2 - 2.0414 * x1 .* x2 + 6.0436 * x2.^2;

V1 = 1023 * (quad_base -  280.9653);
V2 = 1023 * (quad_base -  290.7664);
 

contour(x1, x2, V1, [0 0], 'LineWidth', 2, 'LineColor', 'k')

contour(x1, x2, V2, [0 0], 'LineWidth', 1.5, 'LineColor', 'r', 'LineStyle', '--')


% === Plot regions ===
plotRegion(X0_bound, 'b');    % Intial: blue
plotRegion(X1_bound1, 'r');   % Unsafe set 1: red
plotRegion(X1_bound2, 'r');   % Unsafe set 2: red


axis equal;
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;
xlim([-10 10])
ylim([-10 10])