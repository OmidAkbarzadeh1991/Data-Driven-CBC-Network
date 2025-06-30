% ---------------------------------------------------------------
% This code generates plots for a heterogeneous network of 900 subsystems 
% with a line topology.
% ---------------------------------------------------------------
clc
clear
close all

%% System setup

A1 = [0 1 0;0 0 1];

A2 = [0 1.2 0;0 0 1.2];

A3 = [0 1.4 0;0 0 1.4];

D_coeff_1 = -0.001;

D_1 = D_coeff_1 * [1 0;0 1];

D_coeff_2 = -0.002;

D_2 = D_coeff_2 * [1 0;0 1];

D_coeff_3 = -0.005;

D3 = D_coeff_3 * [1 0;0 1]; 

D_coeff_4 = -0.03;

D_4 = D_coeff_4 * [1 0;0 1]; 

D_coeff_5 = -0.04;

D5 = D_coeff_5 * [1 0;0 1]; 


B = [1 0;0 1];


N = 900; % Number of subsystems

m = 2; % Rows of A_i, D_i
n = 3; % Number of columns of A_i

% Initial state
x = 1 * (-4 + 8 * rand(2 * N, 1));

%% Build A_network (according to topology) and B_network

A_network = zeros(N*m, N*n);
B_network = [];

for i = 1:N
    B_network = blkdiag(B_network, B);
end

for i = 1:N
    if i <= 300
        Ai = A1;  Di = D_1;
    elseif i == 301
        Ai = A2;  Di = D_4;
    elseif i <= 600
        Ai = A2;  Di = D_2;
    elseif i == 601
        Ai = A3;  Di = D5;
    else
        Ai = A3;  Di = D3;
    end

    rows = (i-1)*m + (1:m);
    cols = (i-1)*n + (1:n);

    A_network(rows, cols) = Ai;

    if i < N
        next_rows = i*m + (1:m);
        % only fill the first two monomial columns
        A_network(next_rows, cols(1:2)) = Di;
        % leave the 3rd column as zeros
    end
end

%% Symbolic controller (five of them)
% To observe the open-loop behavior of the system, set the control input to zero (i.e., multiply it by zero (S = 0)).

S = 1;

syms x_1 x_2

% subsystems 1–300
Control_input1 = S * [ ...
  -1.4426*x_1^3 + 0.42323*x_1^2*x_2 - 1.4503*x_1*x_2^2 + 0.35573*x_2^3 ...
  - 2.4993*x_1^2 + 10.9552*x_1*x_2 - 7.6288*x_2^2 - 418.1091*x_1 + 101.5867*x_2; ...
   0.24933*x_1^3 - 1.2156*x_1^2*x_2 + 0.31124*x_1*x_2^2 - 1.2642*x_2^3 ...
  - 10.6783*x_1^2 + 10.5302*x_1*x_2 - 1.9067*x_2^2 + 96.6408*x_1 - 373.059*x_2 ...
];

% subsystem 301
Control_input2 = S * [ ...
  -1.5104*x_1^3 + 0.46391*x_1^2*x_2 - 1.5242*x_1*x_2^2 + 0.37598*x_2^3 ...
  - 2.8769*x_1^2 + 11.8001*x_1*x_2 - 7.1423*x_2^2 - 435.8647*x_1 + 106.0094*x_2; ...
   0.25352*x_1^3 - 1.2152*x_1^2*x_2 + 0.33132*x_1*x_2^2 - 1.2846*x_2^3 ...
  - 12.1112*x_1^2 + 10.5356*x_1*x_2 - 1.8667*x_2^2 + 105.4726*x_1 - 379.4204*x_2 ...
];

% subsystems 302–600
Control_input3 = S * [ ...
  -1.4857*x_1^3 + 0.44405*x_1^2*x_2 - 1.4975*x_1*x_2^2 + 0.36862*x_2^3 ...
  - 2.7181*x_1^2 + 11.2202*x_1*x_2 - 6.9768*x_2^2 - 429.8445*x_1 + 104.4676*x_2; ...
   0.25782*x_1^3 - 1.2020*x_1^2*x_2 + 0.32555*x_1*x_2^2 - 1.2615*x_2^3 ...
  - 11.5691*x_1^2 + 10.2378*x_1*x_2 - 1.8179*x_2^2 + 102.4881*x_1 - 372.4698*x_2 ...
];

% subsystem 601
Control_input4 = S * [ ...
  -1.5642*x_1^3 + 0.49245*x_1^2*x_2 - 1.5838*x_1*x_2^2 + 0.39541*x_2^3 ...
  - 3.1965*x_1^2 + 12.2478*x_1*x_2 - 6.5475*x_2^2 - 450.3345*x_1 + 110.2416*x_2; ...
   0.26560*x_1^3 - 1.1997*x_1^2*x_2 + 0.34929*x_1*x_2^2 - 1.2842*x_2^3 ...
  - 13.2676*x_1^2 + 10.3946*x_1*x_2 - 1.8017*x_2^2 + 113.5013*x_1 - 379.7069*x_2 ...
];

% subsystems 602–900
Control_input5 = S * [ ...
  -1.5358*x_1^3 + 0.46925*x_1^2*x_2 - 1.5521*x_1*x_2^2 + 0.38557*x_2^3 ...
  - 2.9906*x_1^2 + 11.5616*x_1*x_2 - 6.3527*x_2^2 - 443.4257*x_1 + 108.1686*x_2; ...
   0.26899*x_1^3 - 1.1885*x_1^2*x_2 + 0.34254*x_1*x_2^2 - 1.2606*x_2^3 ...
  - 12.5955*x_1^2 + 10.0035*x_1*x_2 - 1.7356*x_2^2 + 109.7795*x_1 - 372.5294*x_2 ...
];

% fast‐function conversion
controller_fast1 = matlabFunction(Control_input1, 'Vars', {[x_1; x_2]});
controller_fast2 = matlabFunction(Control_input2, 'Vars', {[x_1; x_2]});
controller_fast3 = matlabFunction(Control_input3, 'Vars', {[x_1; x_2]});
controller_fast4 = matlabFunction(Control_input4, 'Vars', {[x_1; x_2]});
controller_fast5 = matlabFunction(Control_input5, 'Vars', {[x_1; x_2]});

monomial_fast = matlabFunction([x_1; x_2; x_1^2], 'Vars', {[x_1; x_2]});

%% Simulation
tau    = 1e-3;
tspan  = 2;         
T      = tspan/tau; 

States = x.';

for k = 1:T
    u    = zeros(2*N,1);
    mono = zeros(3*N,1);

    for j = 1:N
        % extract local state
        idx     = (j-1)*2 + (1:2);
        x_local = x(idx);

        % dispatch the right symbolic controller
        if     j <= 300
            u_local = controller_fast1(x_local);
        elseif j == 301
            u_local = controller_fast2(x_local);
        elseif j <= 600      % 302–600
            u_local = controller_fast3(x_local);
        elseif j == 601
            u_local = controller_fast4(x_local);
        else                 % 602–900
            u_local = controller_fast5(x_local);
        end

        % write into the big vectors
        u(idx) = u_local;
        mono((j-1)*3 + (1:3)) = monomial_fast(x_local);
    end

    % network update
    x = x + tau*(A_network*mono + B_network*u);

    States = [States; x.'];
end

%% ===================== Plot 2D states of network =======================
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
% === Add Both Contours ===

x1_range = linspace(-10, 10, 100);
x2_range = linspace(-10, 10, 100);
[x1, x2] = meshgrid(x1_range, x2_range);

quad_base1 =  300 * (3.0597*x1.^2 - 1.4599*x1.*x2 + 2.7439*x2.^2);

quad_base2 =  299 * (3.1479*x1.^2 - 1.5244*x1.*x2 + 2.7371*x2.^2);

quad_base3 =  299 * (3.2505*x1.^2 - 1.6046*x1.*x2 + 2.7344*x2.^2);

quad_base4 = 3.1924*x1.^2 - 1.5561*x1.*x2 + 2.786*x2.^2;

quad_base5 = 3.3016*x1.^2 - 1.6452*x1.*x2 + 2.7846*x2.^2;

quad_base = quad_base1 + quad_base2 + quad_base3 + quad_base4 + quad_base5;

GAMMA = 300 * (121.1384) + 299 * ( 123.1940) + 299 * (125.6914) + 125.1777 + 127.8835;

BETA  = 300 * (123.3770) + 299 * ( 125.3563) + 299 * (127.7447) +  127.3447 + 129.9200;

V1 =  (quad_base -    GAMMA);
V2 =  (quad_base -  BETA);


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