% ---------------------------------------------------------------
% This code implements a subnetwork of 299 subsystems (subsystem 302 to subsystem 600) 
% with a line topology.
% ---------------------------------------------------------------
clc
clear
close all
tic
ops.solver = 'mosek'; % Specify SOS solver
rng(1)  % Set random seed for reproducibility

%% ======================= Regions of interest ============================

X0_bound =  [-4 4;-4 4]; % initial set
X1_bound1 = [-10 -6;-10 -5]; % unsafe set 1
X1_bound2 = [6 10;5 10]; % unsafe set 2
X_bound =   [-10 10;-10 10]; % space state

%======================= System parameters ================================

T = 20;  % Number of collected data

tau = 0.01; % Sampleing time

interconnected_subs = 1;

n_states = 2; % Number of states

initial = 5 * [-0.1; -0.1];  % Initial state vector [x_1; x_2]

%===================== Subsystem matrices (A_i, B_i, D_i) =================

A = [0 1.2 0;0 0 1.2];

B = [1 0;0 1];

D_coeff = -0.002;  

D = D_coeff * [1 0;0 1]; % D_ij


Matrix_D = repmat(D,1 ,interconnected_subs); % Block matrix D_i   

states = initial.';

statesdot = [];

N0T = [];

%=========================== ODE ==========================================

Vector_W = [];

for i = 1:T

    x0 = initial;

      % Original vector of system monomials: 
    M  = [x0(1),x0(2), x0(1)^2]; 

     % Dictionary containing all combinations of system monomials up to degree 2.
    M1 = [x0(1), x0(2), x0(1) * x0(2), x0(1)^2, x0(2)^2];

       % Random input
       u1(i, :) = 100 * (-1 + 2 * rand); 
       u2(i, :) = 100 * (-1 + 2 * rand);

       u(i,:) = [u1(i, :), u2(i, :)];
    
     tempArray = [];

    for k = 1 : interconnected_subs

      % The result must be robust to all state values of other subsystems affecting this subsystem; 
      % hence, random sampling within a safe range is sufficient.
      w = 0.1 * [-5 + 10 * rand, -5 + 10 * rand];

      tempArray = vertcat(tempArray, w.');

    end

    Vector_W = vertcat(Vector_W, tempArray.');

    odeSystem = @(t, x) A * M.' + B * u(i,:).' + Matrix_D * Vector_W(i,:).'; 
    
    
      tspan = [(i-1) * tau  i * tau];

    options = odeset('RelTol',1e-3,'AbsTol',1e-6);

    [t, X] = ode45(odeSystem, tspan, x0, options);
    
    initial = X(end, :).';
    
    states = [states; initial.'];
    
    statesdot = [statesdot; (A*M.' + B * u(i,:).' + Matrix_D *Vector_W(i,:).').'];
    
    N0T = [N0T; M1];
    
end

%===================== Trajectories in equation (9)  ======================

U0T = u.'; % External input trajectory

W0T = Vector_W.';  % Internal input trajectory

X0T = states(1:end-1, :).'; % States trajectory

% State derivatives trajectory influenced by noise, which satisfies the bound specified in Equation (11).
X1T = statesdot.'+ (-0.3 + 0.6 * rand(n_states,T));

%Dictionary containing all combinations of system monomials up to degree 2 trajectory
N0T = N0T.';

clear tau initial  states  statesdot  x0  x01 u t X i odeSystem  tspan



%% ================== SOS Program: Condition (16d & 15) ===================

pvar x_1 x_2  %Polynomial vars 

x = [x_1 x_2].'; % State vector

%Transformation matrix in equation (1b)
Upsilon_x  = [1 0;0 1; x_2 0; x_1 0; 0 x_2];

% ============= Initialization the sum of squares program =================
prog = sosprogram(x);

% ==================== Define Lagransian Monomials ========================
Mon_Lagransian = monomials(x,[0:2]);


% ========================== Computing b_x ================================
b_x_1 = (x_1 - X_bound(1, 1))*(X_bound(1, 2) - x_1);

b_x_2 = (x_2 - X_bound(2, 1))*(X_bound(2, 2) - x_2); 

b_x = [b_x_1; b_x_2];


%========================= Lagrangian multipliers =========================

[prog,L3_1] = sospolyvar(prog,Mon_Lagransian,'wscoeff');

[prog,L3_2] = sospolyvar(prog,Mon_Lagransian,'wscoeff');

Lambda_x = [L3_1; L3_2];


PI  = monomials(x,0); % Constant PI

mu  = monomials(x,0); % Constant mu

[prog,PI]  = sospolyvar(prog,PI);

[prog,mu]  = sospolyvar(prog,mu);

invP = monomials(x, 0); %Positive definite matrix P_i^(-1)

[prog, invP] = sospolymatrixvar(prog, invP, [length(x), length(x)], 'symmetric');

H_x = monomials(x, 0 : 2); %In condition (15)

[prog, H_x] = sospolymatrixvar(prog, H_x, [T, length(x)]);


epsilon = 0.99;  % Decay rate fixed a priori

bar_varkappa = n_states * 0.09 % Noise bound in (11) 


% Condition (20d) guarantees the fulfillment of Condition (16d).
G = - epsilon * invP - ((X1T- Matrix_D * W0T) * H_x) - (H_x.' * (X1T- Matrix_D * W0T).') - mu * (T * bar_varkappa * eye(length(x))) - PI * (eye(length(x)));

Condition_3_creation = [G  H_x.'; H_x  mu * eye(T)] - Lambda_x.' * b_x * eye(T+ length(x));

% Condition for positiveness of Lagransian multiplier
for ii = 1:length(Lambda_x)
    prog = sosineq(prog, Lambda_x(ii));
end

% Positiveness of condition (15)
prog = soseq(prog, N0T * H_x - Upsilon_x * invP);

% Positiveness of matrix P^(-1)
prog = sosmatrixineq(prog, invP - 10^(-6) * eye(length(x)), 'quadraticMineq');

% Positiveness of condition (20d)
prog = sosmatrixineq(prog, Condition_3_creation,'quadraticMineq');

% Conditions for positiveness of coefficients mu and PI
prog = sosineq(prog, PI - 10^(-6));

prog = sosineq(prog, mu - 10^(-6));

%============ Call SOS solver (requires Mosek installed) =================

prog = sossolve(prog,ops);

%========================== Cheking the results ===========================

PI = sosgetsol(prog,PI);

mu = sosgetsol(prog,mu);

invP = sosgetsol(prog, invP);

P = (double(invP))^(-1);

H_x = sosgetsol(prog, H_x);


msg = 'Control input & Barrier certificate';
border = repmat('-', 1, length(msg) + 4);
disp(border);
disp(['* ', msg, ' *']);
disp(border);
Control_input = U0T * H_x * P * x

Barrier_certificate = x.' * P * x

PI = double(PI) 

mu = double(mu)

disp('eigs of P:');

disp(double(eig(P)));

clear prog


%% ===================== SOS Program: Condition (16a)======================
prog = sosprogram(x);


% ==================== Define Lagransian Monomials ========================
Mon_Lagransian2 = monomials(x,[0:2]);


% ========================== Computing b_x ================================
B_x_1 = (x_1 - X_bound(1, 1))*(X_bound(1, 2) - x_1);

B_x_2 = (x_2 - X_bound(2, 1))*(X_bound(2, 2) - x_2); 

B_x = [B_x_1; B_x_2];


%========================= Lagrangian multipliers =========================

[prog,LL3_1] = sospolyvar(prog,Mon_Lagransian2,'wscoeff');

[prog,LL3_2] = sospolyvar(prog,Mon_Lagransian2,'wscoeff');

Lambdaa_x = [LL3_1; LL3_2];


Mon_phi  = monomials(x,[0]);

[prog,Mon_phi]  = sospolyvar(prog,Mon_phi);

% Condition (20a) implies condition (16a)
Condition_1_creation = (x.' *  P * x) - (Mon_phi * (x(1)^2 + x(2)^2)) - Lambdaa_x.' * B_x;

% Condition for positiveness of Lagransian multiplier

for ii = 1:length(Lambdaa_x)
    prog = sosineq(prog, Lambdaa_x(ii));
end

% Positiveness Condition (20a)
prog = sosineq(prog, Condition_1_creation);

prog = sosineq(prog, Mon_phi);

prog = sossolve(prog,ops);

Phi = sosgetsol(prog,Mon_phi);

disp('phi:');

disp(double(Phi));

rho = (1/(PI))*(norm(Matrix_D))^2  % Interaction gain term

delta_i_j = (1/(Phi * PI))*(norm(Matrix_D))^2 % In (22)

clear x prog 

%% ================== SOS Program: Condition (16b & 16c)===================

pvar x1 x2

var_x = [x1 x2].';

% ============= Initialization the sum of squares program =================

prog = sosprogram(var_x);

% ========================= Computing b_0 b_u =============================

%******* Initial set
b0_1 = (x1 - X0_bound(1, 1))*(X0_bound(1, 2) - x1);
b0_2 = (x2 - X0_bound(2, 1))*(X0_bound(2, 2) - x2);


b0 = [b0_1; b0_2];

%******* Unsafe set 1
b_u_1 = (x1 - X1_bound1(1, 1))*(X1_bound1(1, 2) - x1);
b_u_2 = (x2 - X1_bound1(2, 1))*(X1_bound1(2, 2) - x2);


b_u = [b_u_1; b_u_2];

%******* Unsafe set 2
b_u_2_1 = (x1 - X1_bound2(1, 1))*(X1_bound2(1, 2) - x1);
b_u_2_2 = (x2 - X1_bound2(2, 1))*(X1_bound2(2, 2) - x2);

b_uu = [b_u_2_1; b_u_2_2];

% ==================== Define Lagransian Monomials ========================

Mon_Lagransian1 = monomials(var_x,[0 : 2]); % Lagransian multiplier

% ========================= Lagransian multiplier =========================

[prog,L0_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L0_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');


L0 = [L0_1; L0_2];

[prog,L1_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L1_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');


L1 = [L1_1; L1_2];

[prog,L2_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L2_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');


L2 = [L2_1; L2_2];

% ========================= Sub-barrier: B_i(x_i) =========================

Barrier = var_x.' * P * var_x;  % Sub-barrier function

gamma  = monomials(var_x,0);  % Initial set level-set: gamma_i

beta  = monomials(var_x,0); % Unsafe set level-set: beta_i


[prog,gamma]  = sospolyvar(prog,gamma);

[prog,beta]  = sospolyvar(prog,beta);

% Conditions for positiveness of coefficients(enforce β > γ and non-negativity)
prog = sosineq(prog, beta - gamma);

prog = sosineq(prog, gamma);

prog = sosineq(prog, beta);

% ========================= Constraints ===================================

% Condition 20c which implies 16c must be satisfied individually for each unsafe region separately.

prog = sosineq(prog, -Barrier - L0.'*b0 + gamma); % Condition 20b implies 16b

prog = sosineq(prog, Barrier - L1.'*b_u - beta); % Condition 20c for first unsafe region

prog = sosineq(prog, Barrier - L2.'*b_uu - beta); % Condition 20c for second unsafe region


% Conditions for positiveness of Lagransian multiplier

for jj = 1:length(L0)
    prog = sosineq(prog,L0(jj));
end


for kk = 1:length(L1)
    prog = sosineq(prog,L1(kk));
end

for zz = 1:length(L2)
    prog = sosineq(prog,L2(zz));
end

% ============== Call SOS solver (requires Mosek installed) ==============

prog = sossolve(prog,ops);

% ========================= Cheking the results =========================

SOLV0 = sosgetsol(prog,-Barrier - L0.'*b0 + gamma);
SOLV1 = sosgetsol(prog, Barrier - L1.'*b_u - beta);
SOLV2 = sosgetsol(prog, Barrier - L2.'*b_uu - beta);


[a0,m0] = findsos(SOLV0,ops);
[a1,m1] = findsos(SOLV1,ops);
[a2,m2] = findsos(SOLV2,ops);

% 
a = [length(a0),length(a1), length(a2)];
m = [length(m0),length(m1), length(m2)];

% ========================= Results =========================
a

gamma = double(sosgetsol(prog, gamma))

beta = double(sosgetsol(prog, beta))

toc
allvars = whos;
memused = sum([allvars.bytes])
