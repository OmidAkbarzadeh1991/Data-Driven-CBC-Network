% ---------------------------------------------------------------
% This code generates a plot indicating the satisfaction of Condition (16a) 
% for Duffing oscillator subsystems.
% ---------------------------------------------------------------
clc
clear all
close all

% Define the range for x1 and x2
x1 = linspace(-10, 10, 1000);
x2 = linspace(-10, 10, 1000);

% Create a meshgrid
[X1, X2] = meshgrid(x1, x2);

% Define the function
Z =    8.4503 * X1.^2 - 2.0223* X1 .* x2 + 5.6554* X2.^2 - 0.1054 * (X1.^2 + X2.^2);

% Plot the surface
figure;
surf(X1, X2, Z);
pbaspect([1 1 1])  % equal plot-box aspect ratio
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 28)
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

zlabel( ...
    '\boldmath $x_i^\top\big[\Upsilon_i(x_i)^{\dagger}\mathcal N_i^{0,\mathcal{T}}\mathcal H_i(x_i)\big]^{-1} x_i - \phi_i\|x_i\|^2$', ...
    'Interpreter','latex', ...
    'FontSize',15);

shading interp;
colormap("default");
axis tight;
set(gcf, 'Renderer', 'opengl');
set(gcf, 'Color', 'w');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 1, 'LineWidth', 1);
box on;
xlim([-10 10]);
ylim([-10 10]);