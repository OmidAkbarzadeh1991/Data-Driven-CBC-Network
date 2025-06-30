% ---------------------------------------------------------------
% This code generates plots for a network of 1000 Lorenz subsystems 
% with a fully interconnected topology.
% ---------------------------------------------------------------
clc
clear
close all
%% System setup
sigma = 10;
rho = 28;
beta = 8/3;

A = [-sigma, sigma, 0, 0, 0;
     rho, -1, 0, -1, 0;
     0, 0, -beta, 0, 1];


B = [1 0 0; 0 1 0; 0 0 1].';

D_coeff = 0.00002;

D = D_coeff * [-1,  0,  0;  
                0, -1, 0;
                0, 0, -1];

N = 1000;  % Number of subsystems

x = 3 * (-1 + 2 * rand(3 * N, 1)); % initial state

D_expanded = [D, zeros(3, 2)];   % matrix D_ij

%% Build A_network (according to topology) and B_network
B_network = [];
for i = 1:N
    B_network = blkdiag(B_network, B);
end

A_network = zeros(3*N, 5*N);
for i = 1:N
    for j = 1:N
        if i == j
            A_network(3*(i-1)+1:3*i, 5*(j-1)+1:5*j) = A;
        else
          A_network(3*(i-1)+1:3*i, 5*(j-1)+1:5*j) = D_expanded;
        end
    end
end

%% Define symbolic controller and monomials
syms x_1 x_2 x_3

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

S = 1;

Control_input(1,1)  =    S * ( 21.8533*x_1^2 - 67.7554*x_1*x_2 - 67.6251*x_1*x_3 + 47.7192*x_2^2 + 96.3149*x_2*x_3 + 55.7623*x_3^2 - 1237.1117*x_1 + 197.717*x_2 + 260.4928*x_3  );
 
Control_input(2,1)  =    S * ( 117.0317*x_1^2 - 138.0305*x_1*x_2 - 164.201*x_1*x_3 + 6.5751*x_2^2 + 12.0967*x_2*x_3 + 12.9936*x_3^2 + 51.4842*x_1 - 318.5331*x_2 - 109.3699*x_3);

Control_input(3,1)  =    S * ( 19.0949*x_1^2 - 59.5318*x_1*x_2 - 105.0946*x_1*x_3 + 33.7661*x_2^2 + 53.6538*x_2*x_3 + 27.7845*x_3^2 + 314.9155*x_1 - 270.904*x_2 - 458.6184*x_3 );



% Convert symbolic controller to function handle
controller_fast = matlabFunction(Control_input, 'Vars', {[x_1; x_2; x_3]});

% Similarly for monomials
monomials_sym = [x_1; x_2; x_3; x_1*x_3; x_1*x_2];
monomial_fast = matlabFunction(monomials_sym, 'Vars', {[x_1; x_2; x_3]});

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

%% ===================== Plot 3D states of network =======================
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

f1 = @(x1, x2, x3)  1000 * (18.7668.*x1.^2 - 5.8139.*x1.*x2 - 7.728.*x1.*x3 + 5.5492.*x2.^2 + 6.0589.*x2.*x3 + 6.3392.*x3.^2  -  478.7175);

f2 = @(x1, x2, x3)  1000 * (18.7668.*x1.^2 - 5.8139.*x1.*x2 - 7.728.*x1.*x3 + 5.5492.*x2.^2 + 6.0589.*x2.*x3 + 6.3392.*x3.^2 - 490.1669);


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

ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
zlabel('$x_3$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');


axis equal;
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)

set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;
view(3);

% Regions of interest 

X0_bound  = 3 * [-1 1;-1 1;-1 1]; % initial set
X1_bound1 = [-20 -4;-20 -15;4 20]; % unsafe set 1
X1_bound2 = [8 20;11 20;4 20]; % unsafe set 2
X1_bound3 = [8 20;11 20;-20 -5]; % unsafe set 3 

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
    plotFace(v0, faces(i, :), 'b'); % Initial set in red
    plotFace(v1, faces(i, :), 'r'); % Unsafe set 1 in green
    plotFace(v2, faces(i, :), 'r'); % Unsafe set 2 in blue
    plotFace(v3, faces(i, :), 'r'); % Unsafe set 3 in magenta
end

xlim([-20 20]); 
ylim([-20 20]);  
zlim([-20 20]); 