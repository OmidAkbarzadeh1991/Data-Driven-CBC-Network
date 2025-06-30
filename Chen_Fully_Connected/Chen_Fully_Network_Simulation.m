% ---------------------------------------------------------------
% This code generates plots for a network of 2000 Chen subsystems 
% with a fully connected topology.
% ---------------------------------------------------------------
clc
clear
close all
%% System setup
teta1 = 35;
teta2 = 28;
teta3 = 3;

N = 1000;    % Number of subsystems

m = 3;      % States per subsystem
n = 5;      % Number of columns of A_i

A = [-teta1, teta1, 0, 0, 0;
      teta2-teta1, teta2, 0, -1, 0;
      0, 0, -teta3, 0, 1];

B = [0 1 0; 0 0 1].';   % 3x2 input matrix

D_coeff = 0.00005;

D = D_coeff * [-1, 0, 0; 
               0, -1, 0;
               0, 0, -1];

%% Initial state
x = 1 * (-2.5 + 5 * rand(3 * N, 1));  

%% Build A_network (according to topology) and B_network

A_network = zeros(m*N, n*N);
B_network = [];

for i = 1:N
    B_network = blkdiag(B_network, B);
end

for i = 1:N
    for j = 1:N
        if i == j
            A_network(3*(i-1)+1:3*i, 5*(j-1)+1:5*j) = A;
        else
            A_network(3*(i-1)+1:3*i, 5*(j-1)+1:5*j) = [D, zeros(3,2)]; % matrix D_ij
        end
    end
end

%% Symbolic controller and monomials
syms x_1 x_2 x_3

% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

S = 1;

Control_input(1,1)  =  ...
 S * ( 131.852*x_1^2 + 174.5084*x_1*x_2 + 788.7689*x_1*x_3 + 57.7206*x_2^2  + 495.3951*x_2*x_3 + 874.5197*x_3^2 - 2112.1811*x_1 - 56.1417*x_2 + 1193.6544*x_3);

Control_input(2,1)  = ...
 S * ( -86.624*x_1^2 - 104.6167*x_1*x_2 - 467.8595*x_1*x_3 - 29.3959*x_2^2 - 236.3971*x_2*x_3 - 350.9533*x_3^2 - 612.0299*x_1 - 519.3022*x_2 - 1292.8814*x_3);

% Convert symbolic controller to function handle
controller_fast = matlabFunction(Control_input, 'Vars', {[x_1; x_2; x_3]});

% Similarly for monomials
monomials_sym = [x_1; x_2; x_3; x_1*x_3; x_1*x_2];
monomial_fast = matlabFunction(monomials_sym, 'Vars', {[x_1; x_2; x_3]});
M_x = [x_1; x_2; x_3; x_1*x_3; x_1*x_2];


%% Simulation
tau = 1e-4;
tspan = 2;
T = tspan / tau;

States = x.';

for i = 1:T
    u = zeros(2*N, 1);
    mono = zeros(5*N, 1);

    for j = 1:N
        idx = (j-1)*3 + 1;
        x_local = x(idx:idx+2);

        u_local = controller_fast(x_local);
        m_local = monomial_fast(x_local);

        u(2*(j-1)+1:2*j) = u_local(:);
        mono(5*(j-1)+1:5*j) = m_local(:);
    end

    x = x + tau * (A_network * mono + B_network * u);

     States = [States; x.'];
end

%% ===================  Plot 3D states of network ==========================
num_plots = 120;
figure
hold on;

random_subs = randperm(N, num_plots); % Ensuring unique random indices

for i = 1:num_plots
    random_subsystems = random_subs(i);
    plot3(States(:, 3 * random_subsystems - 2), States(:, 3 * random_subsystems - 1), States(:, 3 * random_subsystems),...
          'LineWidth', 2, 'LineStyle', '-');
end


% The following functions plot the level sets corresponding to the initial and unsafe sets.

f1 = @(x1, x2, x3) 1000 * ( 18.9041.*x1.^2 + 11.2502.*x1.*x2 + 17.1341.*x1.*x3 +  2.5725.*x2.^2 + 10.1225.*x2.*x3 + 12.6188.*x3.^2 -    514.5168);

f2 = @(x1, x2, x3)  1000 * ( 18.9041.*x1.^2 + 11.2502.*x1.*x2 + 17.1341.*x1.*x3 +  2.5725.*x2.^2 + 10.1225.*x2.*x3 + 12.6188.*x3.^2  -  527.5925);


x_range = linspace(-20, 20, 200);
y_range = linspace(-20, 20, 200);
z_range = linspace(-20, 20, 200);
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

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;
axis equal;

% Regions of interest 
X0_bound  = [-2.5 2.5;-2.5 2.5;-2.5 2.5]; % initial set
X1_bound1 = [-20 -4;-20 -5;-20 -4]; % unsafe set 1
X1_bound2 = [3.5 20;5 20;4 20]; % unsafe set 2

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


% Define faces
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];


% Plot each cuboid by faces
for i = 1:size(faces, 1)
    plotFace(v0, faces(i, :), 'blue'); % Initial set in red
    plotFace(v1, faces(i, :), 'red'); % Unsafe set 1 in green
    plotFace(v2, faces(i, :), 'red'); % Unsafe set 2 in blue
end

view(3);
xlim([-20 20]); 
ylim([-20 20]);  
zlim([-20 20]); 
