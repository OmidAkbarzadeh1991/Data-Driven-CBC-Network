% -------------------------------------------------------------------------
% This code implements the compositionality condition for a heterogeneous network
% of 900 subsystems arranged in a line topology.
% -------------------------------------------------------------------------
clc
clear
close all
%============================ parameters ==================================
tic
number_subsystems = 900; % Number of subsystems

epsilon_local = 0.99; % Local decay rate epsilon_i

% Interaction gains rho_i reported in Table II
rho = [1.15e-6, 1.0e-3, 4.62e-6, 1.0e-3, 2.89e-5];

% Corresponding phi_i values used in the verification of Condition (16a)
phi = [0.0783, 0.0785, 0.0789, 0.0790, 0.0794];

%% ========================== Compositionality  ===========================

% For the line topology, Delta_{i,i-1} = rho_i / phi_{i-1}.
aVec = zeros(number_subsystems,1);
aVec(2:300)   = rho(1)/phi(1);
aVec(301)     = rho(2)/phi(1);
aVec(302)     = rho(3)/phi(2);
aVec(303:600) = rho(3)/phi(3);
aVec(601)     = rho(4)/phi(3);
aVec(602)     = rho(5)/phi(4);
aVec(603:end) = rho(5)/phi(5);

% Assemble -hat{epsilon} + Delta
D = diag(-epsilon_local*ones(number_subsystems,1));
L = diag(aVec(2:end), -1);
compose_mat = D + L;

% Check compositional condition (24b)
Composition = ones(1, number_subsystems) * compose_mat;

if all(Composition < 0)

    disp('Compositional condition (24b) is satisfied.');

    max_varpi = max(double(Composition));
    epsilon_upper = -max_varpi;

    disp('Maximum varpi_i:');
    disp(max_varpi);

    disp('Admissible upper bound on network epsilon:');
    disp(epsilon_upper);

    epsilon_network = 0.97;

    if epsilon_network < epsilon_upper
        disp('Selected network epsilon:');
        disp(epsilon_network);
    else
        error('Selected network epsilon does not satisfy Theorem 3.');
    end

else

    disp('Compositional condition (24b) is NOT satisfied.');

end

gamma_network = 300 * (121.1384) + 299 * ( 123.1940) + 299 * (125.6914) + 125.1777 + 127.8835

beta_network  = 300 * (123.3770) + 299 * ( 125.3563) + 299 * (127.7447) +  127.3447 + 129.9200

if gamma_network < beta_network 
    msg3='Compostional condition (23a) is satisfied.';
end

border = repmat('-', 1, length(msg3) + 4);
disp(border);
disp(['* ', msg3, ' *']);
disp(border);
toc