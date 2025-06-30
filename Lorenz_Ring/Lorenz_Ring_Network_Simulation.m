% ---------------------------------------------------------------
% This code generates plots for a network of 2000 Lorenz subsystems 
% with a ring topology.
% ---------------------------------------------------------------
clc
clear
close all
%% System setup

sigma = 10;
rho = 28;
beta = 8/3;

N = 2000;  % Number of subsystems


m = 3;    % State dimension per subsystem
n = 5;    % Number of columns of A_i


A = [-sigma, sigma, 0, 0, 0;
     rho, -1, 0, -1, 0;
     0, 0, -beta, 0, 1];

B = [1 0 0; 0 1 0; 0 0 1].';

D = 0.01 * diag([-1, -1, -1]); 

D_expanded = [D, zeros(3, 2)];  % matrix D_ij

%% Initial state
x = -3 + 6 * rand(3 * N, 1);

%% Build A_network (according to topology) and B_network
A_network = zeros(m*N, n*N);
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

% Wrap-around coupling: last -> first
A_network(1:m, (N-1)*n+1:N*n) = D_expanded;

%% Symbolic Controller Definition
syms x_1 x_2 x_3

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

S = 1;

Control_input(1,1)  = ...
   S * (  -10.0712*x_1^2 + 23.9253*x_1*x_2 + 8.0498*x_1*x_3 - 25.8062*x_2^2  - 6.5798*x_2*x_3 + 6*x_3^2 - 2734.5265*x_1 + 244.7904*x_2 +  249.8385*x_3);

Control_input(2,1)  = ...
   S * ( -32.134*x_1^2 + 135.5226*x_1*x_2 + 35.9192*x_1*x_3 - 13.1452*x_2^2 - 8.983*x_2*x_3 - 42.2782*x_3^2 + 512.586*x_1 - 536.8531*x_2 - 132.6197*x_3);

Control_input(3,1)  = ...
   S * ( -24.4745*x_1^2 + 1.1856*x_1*x_2 - 46.9177*x_1*x_3 - 7.6293*x_2^2  + 36.0765*x_2*x_3 + 8.4757*x_3^2 - 42.8253*x_1 - 53.8814*x_2 - 583.1774*x_3 );


M_x = [x_1; x_2; x_3; x_1*x_3; x_1*x_2];

% === Automatic Conversion ===
controller_fast = matlabFunction(Control_input, 'Vars', {[x_1; x_2; x_3]});
monomial_fast = matlabFunction(M_x, 'Vars', {[x_1; x_2; x_3]});

%% Simulation
tau = 1e-4;
tspan = 2; % Simulation time
T = tspan / tau;

States = x.';

for i = 1:T
    u = zeros(3*N, 1);
    mono = zeros(5*N, 1);

    for j = 1:N
        idx_x = (j-1)*3 + 1;
        x_local = x(idx_x:idx_x+2);

        u_local = controller_fast(x_local);
        m_local = monomial_fast(x_local);

       u(3*(j-1)+1 : 3*j) = u_local(:);
       
        mono(5*(j-1)+1:5*j) = m_local(:);
    end

    x = x + tau * (A_network * mono + B_network * u);

     States = [States; x.'];
end


%% ===================  Plot 3D states of network ==========================
num_plots = 120;
figure
hold on;
random_subs = randperm(N, num_plots); 
for i = 1:num_plots
    random_subsystems = random_subs(i);
    plot3(States(:, 3 * random_subsystems - 2), States(:, 3 * random_subsystems - 1), States(:, 3 * random_subsystems),...
          'LineWidth', 2, 'LineStyle', '-');
end


% The following functions plot the level sets corresponding to the initial and unsafe sets.

f1 = @(x1, x2, x3)  2000 * ( 26.7769.*x1.^2 - 5.821.*x1.*x2  - 3.9914.*x1.*x3 + 5.1485.*x2.^2 + 1.923.*x2.*x3 + 6.1772.*x3.^2  -    501.4458);

f2 = @(x1, x2, x3)  2000 * ( 26.7769.*x1.^2 - 5.821.*x1.*x2  - 3.9914.*x1.*x3 + 5.1485.*x2.^2 + 1.923.*x2.*x3 + 6.1772.*x3.^2   -    514.5160);


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

X0_bound  = [-3 3;-3 3;-3 3]; % initial set
X1_bound1 = [-20 -10;-20 -5;5 20]; % unsafe set 1
X1_bound2 = [3.5 20;15 20;5 20]; % unsafe set 2
X1_bound3 = [3.5 20;15 20;-20 -5]; % unsafe set 3

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