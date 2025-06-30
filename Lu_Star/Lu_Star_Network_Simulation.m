% ---------------------------------------------------------------
% This code generates plots for a network of 2000 Lu subsystems 
% with a star topology.
% ---------------------------------------------------------------
clc
clear
close all

%% Create symbolic controller
syms x_1 x_2 x_3

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

S = 1;

Control_input(1,1)  =  ...
 S * ( 59.8425*x_1^2 + 37.8192*x_1*x_2 + 251.9851*x_1*x_3 + 5.5285*x_2^2 + 74.9351*x_2*x_3 + 252.0663*x_3^2 + 245.721*x_1 - 315.0832*x_2 + 853.58*x_3 );

Control_input(2,1)  = ...
 S * ( -91.8658*x_1^2 - 185.3774*x_1*x_2 - 245.4829*x_1*x_3 - 59.3153*x_2^2 - 377.5074*x_2*x_3 - 34.1371*x_3^2 - 1117.3166*x_1 - 1543.8441*x_2 - 370.2911*x_3);

% Convert symbolic controller to function handle
controller_fast = matlabFunction(Control_input, 'Vars', {[x_1; x_2; x_3]});

% Similarly for monomials
monomials_sym = [x_1; x_2; x_3; x_1*x_3; x_1*x_2];
monomial_fast = matlabFunction(monomials_sym, 'Vars', {[x_1; x_2; x_3]});

%% System setup
a = 36; c = 28; b = 20;

A = [-a, a, 0, 0, 0;
     0, c, 0, -1, 0;
     0, 0, -b, 0, 1];

B = [0 1 0; 0 0 1].';

D = 0.001 * diag([-1 -1 -1]); % matrix D_ij

N = 2000; % Number of subsystems

m = 3; % Number of rows of A_i 
n = 5; % Number of columns of A_i 

x = 1*(-5 + 10 * rand(3 * N, 1)); % initial state

%% Build A_network (according to topology) and B_network

A_network = zeros(N*m, N*n);
B_network = [];
for i = 1:N
    B_network = blkdiag(B_network, B);
end

for i = 1:N
    A_network((i-1)*m + 1:i*m, (i-1)*n + 1:i*n) = A;
    if i < N
        A_network(i*m + 1:(i+1)*m, 1:3) = D; 
    end
end

%% Simulation
tau = 1e-4;
tspan = 2; % Simulation time
T = tspan / tau;

States = x.';

for i = 1:T
    u = [];
    mono = [];
    
    for j = 1:N
        idx = (j-1)*3 + 1;
        x_sub = x(idx:idx+2); % x1,x2,x3 for subsystem j
        
        u_each = controller_fast(x_sub);     
        m_each = monomial_fast(x_sub);       
        
        u = [u; u_each(:)];
        mono = [mono; m_each(:)];
    end
    
    x = x + tau * (A_network * mono + B_network * u);
    States = [States; x.'];
end

%% ===================  Plot 3D states of network ==========================
num_plots = 120;
figure
hold on;

random_subs = randperm(N, num_plots); % Select random subsystems
for i = 1:num_plots
    random_subsystems = random_subs(i);
    plot3(States(:, 3 * random_subsystems - 2), States(:, 3 * random_subsystems - 1), States(:, 3 * random_subsystems),...
          'LineWidth', 2, 'LineStyle', '-');
end

% The following functions plot the level sets corresponding to the initial and unsafe sets.

f1 = @(x1, x2, x3) 2000 * ( 4.1535.*x1.^2 + 6.9638.*x1.*x2 + 4.5311.*x1.*x3 + 5.5811.*x2.^2 + 1.0341.*x2.*x3 + 3.7818.*x3.^2  -   729.9267);

f2 = @(x1, x2, x3) 2000 * ( 4.1535.*x1.^2 + 6.9638.*x1.*x2 + 4.5311.*x1.*x3 + 5.5811.*x2.^2 + 1.0341.*x2.*x3 + 3.7818.*x3.^2  -  749.8527);

x_range = linspace(-20, 20, 400);
y_range = linspace(-20, 20, 400);
z_range = linspace(-20, 20, 400);
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


axis equal;
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('$x_3$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;


% Regions of interest 

X0_bound = [-5 5;-5 5;-5 5]; % initial set

X1_bound1 = [-20 -10; -20 -15 ;6.5 20]; % unsafe set 1

X1_bound2 = [10 20; 5.5 20 ; 6.5 20]; % unsafe set 2

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
% v3 = calcVertices(X1_bound3);

% Define faces
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

% Plot each cuboid by faces
for i = 1:size(faces, 1)
    plotFace(v0, faces(i, :), 'b'); % Initial set in red
    plotFace(v1, faces(i, :), 'r'); % Unsafe set 1 in green
    plotFace(v2, faces(i, :), 'r'); % Unsafe set 2 in blue
    % plotFace(v3, faces(i, :), 'r'); % Unsafe set 3 in magenta
end

% Configure the plot
view(3);
xlim([-20 20]); 
ylim([-20 20]);  
zlim([-20 20]); 
