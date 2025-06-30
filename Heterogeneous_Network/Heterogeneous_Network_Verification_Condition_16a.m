% ---------------------------------------------------------------
% This code generates a plot indicating the satisfaction of Condition (16a) 
% for five different subsystems including subsystems (301 and 601) in a heterogeneous network of 900 subsystems 
% with a line topology.
% ---------------------------------------------------------------
clc; clear; close all

% grid
x1 = linspace(-10,10,200);
x2 = linspace(-10,10,200);
[X1,X2] = meshgrid(x1,x2);

% surfaces
Z{1} = 3.0597*X1.^2 - 1.4599*X1.*X2 + 2.7439*X2.^2 - 0.0783*(X1.^2+X2.^2);
Z{2} = 3.1479*X1.^2 - 1.5244*X1.*X2 + 2.7371*X2.^2 - 0.0785*(X1.^2+X2.^2);
Z{3} = 3.2505*X1.^2 - 1.6046*X1.*X2 + 2.7344*X2.^2 - 0.0789*(X1.^2+X2.^2);
Z{4} = 3.1924*X1.^2 - 1.5561*X1.*X2 + 2.7860*X2.^2 - 0.0790*(X1.^2+X2.^2);
Z{5} = 3.3016*X1.^2 - 1.6452*X1.*X2 + 2.7846*X2.^2 - 0.0794*(X1.^2+X2.^2);

cols = lines(5);

figure('Color','w','Renderer','opengl'); hold on
for k = 1:5
    mesh(X1, X2, Z{k}, 'EdgeColor', cols(k,:), 'LineWidth', 0.7);
end
hold off
pbaspect([1 1 1])  % equal plot-box aspect ratio
view(45,30)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',24)
zlabel( ...
    '\boldmath $x_i^\top\big[\Upsilon_i(x_i)^{\dagger}\mathcal N_i^{0,\mathcal{T}}\mathcal H_i(x_i)\big]^{-1} x_i - \phi_i\|x_i\|^2$', ...
    'Interpreter','latex', ...
    'FontSize',15);
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;

