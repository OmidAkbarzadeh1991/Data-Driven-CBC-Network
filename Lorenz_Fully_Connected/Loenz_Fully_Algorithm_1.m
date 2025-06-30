% ---------------------------------------------------------------
% This code implements a network of 1000 Lorenz subsystems 
% with a fully interconnected topology.
% ---------------------------------------------------------------
clc
clear
close all
tic

ops.solver = 'mosek'; % Specify SOS solver

rng(1)  % Set random seed for reproducibility

%% ======================= Regions of interest ============================

X0_bound  = 3 * [-1 1;-1 1;-1 1]; % initial set
X1_bound1 = [-20 -4;-20 -15;4 20]; % unsafe set 1
X1_bound2 = [8 20;11 20;4 20]; % unsafe set 2
X1_bound3 = [8 20;11 20;-20 -5]; % unsafe set 3
X_bound =   [-20 20;-20 20;-20 20]; % Space state


%======================= System parameters ================================

D_coeff = 0.00002; % Interconnection coefficient

T = 15;     % Number of collected data

tau = 0.01; % Sampleing time

Subsystems = 1000; % Number of subsystems

initial = 1 * [2; 2; 2]; % Initial state vector [x_1; x_2; x_3]

n_states = 3; % Number of states

%=============== Subsystem matrices (A_i, B_i, D_i) (Lorenz) ==============

sigma = 10; % Lorenz parameter σ

r = 28;  % Lorenz parameter ρ

b = 8/3; % Lorenz parameter β

A = [-sigma, sigma, 0, 0, 0;
     r, -1, 0, -1, 0;
     0, 0, -b, 0, 1];

B = [1 0 0; 0 1 0; 0 0 1].';


D = -D_coeff*[1, 0, 0;   % matrix D_ij
              0, 1, 0;
              0, 0, 1];

arrayD = repmat(D,1 ,Subsystems-1); % Block matrix D_i

states = initial.';
statesdot = [];
N0T = [];

%=========================== ODE ==========================================

% We note that the ODE solver in the subsequent procedure is solely intended for the sake of data collection, 
% which means that if the data is given, this block of code can be commented out, and instead, the given data with the appropriate name should be utilized.

expandedArray = [];

for i = 1:T
    x0 = initial;

    % Original vector of system monomials: 
    M = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2)]; 

   % Dictionary containing all combinations of system monomials up to degree 2.
    M_c = [x0(1), x0(2), x0(3), x0(1)*x0(3), x0(1)*x0(2), x0(2)*x0(3), x0(2)^2, x0(1)^2 , x0(3)^2]; 

       % Random input
         u1(i, :) = 100 * (-1 + 2 * rand); 
         u2(i, :) = 100 * (-1 + 2 * rand);
         u3(i, :) = 100 * (-1 + 2 * rand);

    u(i,:) = [u1(i, :), u2(i, :), u3(i, :)]; 
    
     tempArray = [];
    for k = 1 : Subsystems-1
      
      % The result must be robust to all state values of other subsystems affecting this subsystem; 
      % hence, random sampling within a safe range is sufficient.
      w =  0.1 * [-10 + 20 * rand, -10 + 20 * rand, -10 + 20 * rand];

      tempArray = vertcat(tempArray, w.');

    end
    expandedArray = vertcat(expandedArray, tempArray.');

    odeSystem = @(t, x) A * M.' + B * u(i,:).' + arrayD * expandedArray(i,:).';      

    tspan = [(i-1) * tau  i * tau];
    %options = odeset('RelTol',1e-3,'AbsTol',1e-6);

    [t, X] = ode45(odeSystem, tspan, x0);


    initial = X(end, :).';

    states = [states; initial.'];

    statesdot = [statesdot; (A*M.' + B*u(i, :).' + arrayD *expandedArray(i,:).').'];
    
    N0T = [N0T; M_c];

end

%===================== Trajectories in equation (9)  ======================

U0T = u.'; % External input trajectory

W0T = expandedArray.'; % Internal input trajectory

X0T = states(1:end-1, :).'; % States trajectory

% State derivatives trajectory influenced by noise, which satisfies the bound specified in Equation (11).
X1T = statesdot.'+ (-0.1 + 0.2 * rand(n_states ,T)); 

%Dictionary containing all combinations of system monomials up to degree 2 trajectory
N0T = N0T.'; 

clear tau initial states  statesdot  x0  x01 u t X i odeSystem tspan


%% ================== SOS Program: Condition (16d & 15) ===================

pvar x_1 x_2 x_3    %Polynomial vars 

x = [x_1 x_2 x_3].';   % State vector

%Transformation matrix Upsilon_x in equation (1b)

Upsilon_x = [1 0 0;0 1 0;0 0 1; x_3 0 0; x_2 0 0; 0 x_3 0; 0 x_2 0 ; x_1 0 0 ; 0 0 x_3];

% ========== Initialization the sum of squares (SOS) program ==============

prog = sosprogram(x);


% ==================== Define Lagransian Monomials ========================

Mon_Lagransian = monomials(x,[0:2]);


% ========================== Computing b_x ================================

b_x_1 = (x_1 - X_bound(1, 1))*(X_bound(1, 2) - x_1);

b_x_2 = (x_2 - X_bound(2, 1))*(X_bound(2, 2) - x_2); 

b_x_3 = (x_3 - X_bound(3, 1))*(X_bound(3, 2) - x_3);

b_x = [b_x_1; b_x_2; b_x_3];

%========================= Lagrangian multipliers =========================

[prog,lambda_x_1] = sospolyvar(prog,Mon_Lagransian,'wscoeff');

[prog,lambda_x_2] = sospolyvar(prog,Mon_Lagransian,'wscoeff');

[prog,lambda_x_3] = sospolyvar(prog,Mon_Lagransian,'wscoeff');

lambda_x = [lambda_x_1; lambda_x_2; lambda_x_3];

PI  = monomials(x,0); % Constant PI

mu  = monomials(x,0); % Constant mu

[prog,PI]  = sospolyvar(prog,PI);
[prog,mu]  = sospolyvar(prog,mu);

invP = monomials(x, 0);  %Positive definite matrix P_i^(-1)

[prog, invP] = sospolymatrixvar(prog, invP, [length(x), length(x)], 'symmetric');

H_x = monomials(x, 0 : 1);   %In condition (15)

[prog, H_x] = sospolymatrixvar(prog, H_x, [T, length(x)]);


epsilon = 0.99; % Decay rate fixed a priori

bar_varkappa = n_states * 0.01 % In equation (11)

% Condition (20d) guarantees the fulfillment of Condition (16d).
G = - epsilon * invP - ((X1T- arrayD * W0T) * H_x) - (H_x.' * (X1T- arrayD * W0T).') - mu * (T * bar_varkappa * eye(length(x))) - PI * (eye(length(x)));


Condition_3_creation = [G  H_x.'; H_x  mu * eye(T)] - lambda_x.' * b_x * eye(T + length(x));

% Condition for positiveness of Lagransian multiplier

for ii = 1:length(lambda_x)
    prog = sosineq(prog, lambda_x(ii));
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

msg = 'Local control input & Sub-barrier certificate';
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

B_x_3 = (x_3 - X_bound(3, 1))*(X_bound(3, 2) - x_3);

B_x = [B_x_1; B_x_2; B_x_3];

%========================= Lagrangian multipliers =========================

[prog,lambdaa_x_1] = sospolyvar(prog,Mon_Lagransian2,'wscoeff');

[prog,lambdaa_x_2] = sospolyvar(prog,Mon_Lagransian2,'wscoeff');

[prog,lambdaa_x_3] = sospolyvar(prog,Mon_Lagransian2,'wscoeff');

lambdaa_x = [lambdaa_x_1; lambdaa_x_2; lambdaa_x_3];


Mon_phi  = monomials(x,[0]);  % Constant phi

[prog,Mon_phi]  = sospolyvar(prog,Mon_phi);

% Condition (20a) implies (16a)
Condition_1_creation = (x.' *  P * x) - (Mon_phi * (x(1)^2 + x(2)^2 + x(3)^2)) - lambdaa_x.' * B_x; 

% Condition for positiveness of Lagransian multiplier
for ii = 1:length(lambdaa_x)
    prog = sosineq(prog, lambdaa_x(ii));
end

% Positiveness of condition (20a)
prog = sosineq(prog, Condition_1_creation);

% Positiveness phi
prog = sosineq(prog, Mon_phi);

prog = sossolve(prog,ops);

Phi = sosgetsol(prog,Mon_phi);

disp('phi:');
disp(double(Phi));

clear x prog 
%% ================== SOS Program: Condition (16b & 16c)===================

pvar x1 x2 x3

var_x = [x1 x2 x3].';


% ============= Initialization the sum of squares program =================

prog = sosprogram(var_x);


% ========================= Computing b_0 b_u =============================

%********** Initial set
b0_1 = (x1 - X0_bound(1, 1))*(X0_bound(1, 2) - x1);
b0_2 = (x2 - X0_bound(2, 1))*(X0_bound(2, 2) - x2);
b0_3 = (x3 - X0_bound(3, 1))*(X0_bound(3, 2) - x3);

b0 = [b0_1; b0_2 ; b0_3];

%********** Unsafe set 1
b_u_1 = (x1 - X1_bound1(1, 1))*(X1_bound1(1, 2) - x1);
b_u_2 = (x2 - X1_bound1(2, 1))*(X1_bound1(2, 2) - x2);
b_u_3 = (x3 - X1_bound1(3, 1))*(X1_bound1(3, 2) - x3);

b_u = [b_u_1; b_u_2 ; b_u_3];

%********** Unsafe set 2
b_u_2_1 = (x1 - X1_bound2(1, 1))*(X1_bound2(1, 2) - x1);
b_u_2_2 = (x2 - X1_bound2(2, 1))*(X1_bound2(2, 2) - x2);
b_u_2_3 = (x3 - X1_bound2(3, 1))*(X1_bound2(3, 2) - x3);

b_uu = [b_u_2_1; b_u_2_2 ; b_u_2_3];

%********* Unsafe set 3
b_u_3_1 = (x1 - X1_bound3(1, 1))*(X1_bound3(1, 2) - x1);
b_u_3_2 = (x2 - X1_bound3(2, 1))*(X1_bound3(2, 2) - x2);
b_u_3_3 = (x3 - X1_bound3(3, 1))*(X1_bound3(3, 2) - x3);

b_uuu = [b_u_3_1; b_u_3_2 ; b_u_3_3];

% ==================== Define Lagransian Monomials ========================


Mon_Lagransian1 = monomials(var_x,[0 : 2]); 

% ========================= Lagransian multiplier =======================

[prog,L0_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L0_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L0_3] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');

Lambda0 = [L0_1; L0_2 ; L0_3];

[prog,L1_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L1_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L1_3] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');

Lambda_u_1 = [L1_1; L1_2 ; L1_3];

[prog,L2_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L2_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L2_3] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');

Lambda_u_2 = [L2_1; L2_2 ; L2_3];


[prog,L4_1] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L4_2] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');
[prog,L4_3] = sospolyvar(prog,Mon_Lagransian1,'wscoeff');

Lambda_u_3 = [L4_1; L4_2 ; L4_3];
 
% ========================= Sub-barrier: B_i(x_i) =========================

Barrier = var_x.' * P * var_x;  % Sub-barrier function

gamma  = monomials(var_x,0); % Initial set level-set: gamma_i

beta  = monomials(var_x,0);  % Unsafe set level-set:  beta_i


[prog,gamma]  = sospolyvar(prog,gamma);

[prog,beta]  = sospolyvar(prog,beta);

% Conditions for positiveness of coefficients(enforce β_i > γ_i and non-negativity)

prog = sosineq(prog, beta - gamma);

prog = sosineq(prog, gamma);

prog = sosineq(prog, beta);

% ========================= Constraints =================================

% Condition 20c which implies 16c must be satisfied individually for each unsafe region separately.

prog = sosineq(prog, -Barrier - Lambda0.'*b0 + gamma); % Condition 20b implies 16b

prog = sosineq(prog, Barrier - Lambda_u_1.'*b_u - beta); % Condition 20c for first unsafe region

prog = sosineq(prog, Barrier - Lambda_u_2.'*b_uu - beta); %Condition 20c for second unsafe region

prog = sosineq(prog, Barrier - Lambda_u_3.'*b_uuu - beta); % Condition 20c for third unsafe region

% Conditions for positiveness of Lagransian multiplier

for jj = 1:length(Lambda0)
    prog = sosineq(prog,Lambda0(jj));
end

for kk = 1:length(Lambda_u_1)
    prog = sosineq(prog,Lambda_u_1(kk));
end

for zz = 1:length(Lambda_u_2)
    prog = sosineq(prog,Lambda_u_2(zz));
end

for cc = 1:length(Lambda_u_3)
    prog = sosineq(prog,Lambda_u_3(cc));
end

% ============== Call SOS solver (requires Mosek installed) ===============

prog = sossolve(prog,ops);

% ========================= Cheking the results ===========================

SOLV0 = sosgetsol(prog,-Barrier - Lambda0.'*b0 + gamma);
SOLV1 = sosgetsol(prog, Barrier - Lambda_u_1.'*b_u - beta);
SOLV2 = sosgetsol(prog, Barrier - Lambda_u_2.'*b_uu - beta);
SOLV3 = sosgetsol(prog, Barrier - Lambda_u_3.'*b_uuu - beta);
 

[a0,m0] = findsos(SOLV0,ops);
[a1,m1] = findsos(SOLV1,ops);
[a2,m2] = findsos(SOLV2,ops);
[a3,m3] = findsos(SOLV3,ops);
% 
a = [length(a0),length(a1),length(a2) ,length(a3)];
m = [length(m0),length(m1),length(m2) ,length(m3)];

% =============================== Results =================================
a

msg = 'Level sets';
border = repmat('-', 1, length(msg) + 4);
disp(border);
disp(['* ', msg, ' *']);
disp(border);

gamma = double(sosgetsol(prog, gamma))

beta = double(sosgetsol(prog, beta))

gamma_network = Subsystems * gamma

beta_network = Subsystems *  beta

if gamma_network < beta_network

     msg3='Compostional condition (23a) is satisfied.';
end

border = repmat('-', 1, length(msg3) + 4);
disp(border);
disp(['* ', msg3, ' *']);
disp(border);toc
allvars = whos;
memused = sum([allvars.bytes])

%% ============================= Compositionality ===========================

Delta = (1/PI) * (1/Phi) * norm(arrayD)^2; 

%According to topology

Comp = (ones(Subsystems) - eye(Subsystems)) * Delta;

disp('Cheking compostional condition:');

epsil = repmat(epsilon, 1, Subsystems);


Compostion = ones(1, Subsystems) * (Comp - diag(epsil)); % Compositionality condition in equation (23b)

Comp_numeric = double(Compostion);


% Check compositional condition

if isempty(find(Comp_numeric > 0))
    disp('epsilon of network:');
    disp(max(Comp_numeric));  % Displays the maximum value (which will be <= 0)
    msg2 = 'Compostional condition (23b) is satisfied.';
else
    disp(Comp_numeric(1));
    msg2 = 'Compostional condition is NOT satisfied.';
end


border = repmat('-', 1, length(msg2) + 4);
disp(border);
disp(['* ', msg2, ' *']);
disp(border);


rho = (1/PI) * norm(arrayD)^2 % Interaction gain term


