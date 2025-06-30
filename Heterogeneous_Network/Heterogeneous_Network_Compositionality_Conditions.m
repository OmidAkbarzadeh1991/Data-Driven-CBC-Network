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

epsilon = 0.99; % Decay rate fixed a priori

delta_i_j_1 = 1.4757e-05; % subsystems: 1 - 300 
delta_i_j_2 =  0.013201; % subsystems:  301
delta_i_j_3 = 0.00036672;  % subsystems: 302 - 600
delta_i_j_4 =  0.02336;  % subsystems:  601
delta_i_j_5 = 0.0002347;  % subsystems: 602 - 900

%% ========================== Compositionality  ===========================
aVec = zeros(number_subsystems,1);
aVec(2:300)   = 1.4757e-05;   % subsystems: 1 - 300 
aVec(301)     = 5.8888e-05;   % subsystems:  301
aVec(302:600) = 0.00036672;   % subsystems: 302 - 600
aVec(601)     = 0.0001325;    % subsystems:  601
aVec(602:end) = 0.0002347;    % subsystems: 602 - 900

% 2) assemble the matrix
D = diag(-epsilon*ones(number_subsystems,1));    % main diagonal
L = diag(aVec(2:end), -1);                       % sub-diagonal
compose_mat = D + L;

compose_mat(1,end) = 0;

% Check compositional condition
disp('Checking compostional condition:');

Compostion=ones(1, number_subsystems)*(compose_mat); % Compositionality condition in equation (23b)

if isempty(find(Compostion>0))

    disp('epsilon of network:');
    disp(max(double(Compostion)));
    msg2='Compostional condition (23b) is satisfied.';

else

    disp(double(Compostion(1)));
    msg2='Compostional condition is NOT satisfied.';

end

border = repmat('-', 1, length(msg2) + 4);
disp(border);
disp(['* ', msg2, ' *']);
disp(border);

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