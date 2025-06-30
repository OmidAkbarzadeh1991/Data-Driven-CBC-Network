% ---------------------------------------------------------------
% This code generates plots for a network of 2000 spacecraft subsystems 
% with a line topology.
% ---------------------------------------------------------------
clc
clear
close all
%% System setup

j1 = 20; j2 = 20; j3 = 30;

num_states = 3;


A = [0 0 0 (j2 - j3)/j1 0 0; 
     0 0 0 0 (j3 - j1)/j2 0; 
     0 0 0 0 0 (j1 - j2)/j3];


B = blkdiag(1/j1, 1/j2, 1/j3);

D = 4 * B;

D_expanded = [D, zeros(3, 3)];  %D_ij

N = 2000;  % Number of subsystems


x = -2 + 4 * rand(num_states * N, 1); % Initial state

m = 3; % Number of rows of A_i 

n = 6; % Number of columns of A_i 


%% Build A_network (according to topology) and B_network 
A_network = zeros(N*m, N*n);
B_network = [];


for i = 1:N
    B_network = blkdiag(B_network, B);
end


for i = 1:N
    A_network((i-1)*m+1:i*m, (i-1)*n+1:i*n) = A;
    if i < N
        A_network(i*m+1:(i+1)*m, (i-1)*n+1:i*n) = D_expanded;
    end
end

%% Define symbolic controller and monomials
syms x_1 x_2 x_3

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).
S = 1;

u_symbolic =  S *[       
          -50.6439*x_1^2 - 114.2867*x_1*x_2 + 100.2839*x_1*x_3 + 121.479*x_2^2 - 244.9116*x_2*x_3 + 140.0651*x_3^2 - 18245.2415*x_1 - 1172.2093*x_2 - 6485.0942*x_3
  
           -73.0867*x_1^2 - 149.1322*x_1*x_2 + 181.8506*x_1*x_3 + 72.7091*x_2^2 + 232.6915*x_2*x_3 - 29.0203*x_3^2 - 770.8096*x_1 - 8453.0782*x_2 - 3753.6522*x_3

          116.7209*x_1^2 + 546.4828*x_1*x_2 - 430.1346*x_1*x_3 - 352.0659*x_2^2 + 137.7694*x_2*x_3 - 144.2429*x_3^2 - 10470.1736*x_1 - 5060.5*x_2 - 19797.8908*x_3];

M_x = [x_1; x_2; x_3;x_2*x_3; x_3*x_1; x_1*x_2];

% === Automatic conversion to function handle ===

controller_fast = matlabFunction(u_symbolic, 'Vars', {[x_1; x_2; x_3]});
monomial_fast = matlabFunction(M_x, 'Vars', {[x_1; x_2; x_3]});

%% Simulation
tau = 1e-3;
tspan = 20; % Simulation time
T = tspan / tau;

States = x.';

for i = 1:T
    u = [];
    mono = [];

    for j = 1:N
        idx = (j-1)*3 + 1;
        x_local = x(idx:idx+2);

        u_each = controller_fast(x_local);
        m_each = monomial_fast(x_local);

        u = 1 * [u; u_each(:)];
        mono = [mono; m_each(:)];
    end

    x = x + tau * (A_network * mono + B_network * u);

    States = [States; x.'];
end

%% ===================  Plot 3D states of network ==========================

num_plots = 120;
figure
hold on;

random_subsystems  = randperm(N, num_plots); % Ensuring unique random indices

for i = 1:num_plots
    random_subs = random_subsystems(i);
    plot3(States(:, 3 * random_subs - 2), States(:, 3 * random_subs - 1), States(:, 3 * random_subs),...
          'LineWidth', 2, 'LineStyle', '-');
end

% The following functions plot the level sets corresponding to the initial and unsafe sets.

f1 = @(x1, x2, x3) 2000*( 4.886.*x1.^2 + 0.54151.*x1.*x2 + 3.6833.*x1.*x3 + 2.2014.*x2.^2 + 1.7819.*x2.*x3 + 3.4979.*x3.^2 - 73.3892);

f2 = @(x1, x2, x3) 2000*( 4.886.*x1.^2 + 0.54151.*x1.*x2 + 3.6833.*x1.*x3 + 2.2014.*x2.^2 + 1.7819.*x2.*x3 + 3.4979.*x3.^2 - 75.1046);

x_range = linspace(-10, 10, 200);
y_range = linspace(-10, 10, 200);
z_range = linspace(-10, 10, 200);
[x1, x2, x3] = meshgrid(x_range, y_range, z_range);

v1 = f1(x1, x2, x3);
v2 = f2(x1, x2, x3);


hold on;
s1 = isosurface(x1, x2, x3, v1, 0);
p1 = patch(s1);
isonormals(x1, x2, x3, v1, p1);
set(p1, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

s2 = isosurface(x1, x2, x3, v2, 0);
p2 = patch(s2);
isonormals(x1, x2, x3, v2, p2);
set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6); 
axis equal;
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('$x_3$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;

% Regions of interest 

X0_bound =  2 * [-1 1;-1 1;-1 1]; % initial set
X1_bound1 = [2.5 5;-5 -3;-5 -4]; % unsafe set 1
X1_bound2 = [2.5 5;4 5;2.5 5]; % unsafe set 2
X1_bound3 = [-5 -4;4 5;2.5 5]; % unsafe set 3
hold on;

% Define a function to calculate vertices from bounds
calcVertices = @(b) [b(1,1), b(2,1), b(3,1); ...
                     b(1,2), b(2,1), b(3,1); ...
                     b(1,2), b(2,2), b(3,1); ...
                     b(1,1), b(2,2), b(3,1); ...
                     b(1,1), b(2,1), b(3,2); ...
                     b(1,2), b(2,1), b(3,2); ...
                     b(1,2), b(2,2), b(3,2); ...
                     b(1,1), b(2,2), b(3,2)];

% Define a function to plot a single face
plotFace = @(v, f, c) fill3(v(f, 1), v(f, 2), v(f, 3), c, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

% Define vertices for each cuboid
v0 = calcVertices(X0_bound);
v1 = calcVertices(X1_bound1);
v2 = calcVertices(X1_bound2);
v3 = calcVertices(X1_bound3);

% Define faces
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];


% Plot each cuboid by faces
for i = 1:size(faces, 1)
    plotFace(v0, faces(i, :), 'blue'); % Initial set in red
    plotFace(v1, faces(i, :), 'red'); % Unsafe set 1 in green
    plotFace(v2, faces(i, :), 'red'); % Unsafe set 2 in blue
    plotFace(v3, faces(i, :), 'red'); % Unsafe set 3 in magenta
end

view(3);
xlim([-5 5]); 
ylim([-5 5]);  
zlim([-5 5]); 