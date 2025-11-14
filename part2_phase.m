%% ======================== Constants & Wave Properties ======================== %%
mu0 = 4 * pi * 1e-7;
eps0 = 8.854e-12;
freq = 15e9;
omega = 2 * pi * freq;
k0 = omega * sqrt(mu0 * eps0);
lambda_0 = 2 * pi / k0;
eta = sqrt(mu0 / eps0);
E0 = 1;
n_0 = 1.5;
a = 1;

%% ======================== Phase Tracing Calculations ======================== %%           % X range
theta_deg = 150;   % Incident angles
theta_rad = deg2rad(theta_deg);      % Convert to radians
x_0 = 2;
n_x0 = sqrt(x_0/a + n_0^2);          % Index at x_0
beta = n_x0 * sin(theta_rad);
disp(beta);

v_0 = sqrt(n_0^2 + x_0/a - beta.^2);
x_t = a * (beta^2 - n_0^2);
v_t = sqrt(n_0^2 + x_t/a - beta.^2);
disp(v_t);
v = linspace(v_t, 200, 1000); 
% for psi from v_0 to v_t we will use v = linspace(v_0, v_t, 1000);
x_rays = a .* (v.^2 - n_0^2 + beta.^2);


%% ======================== Phase Plotting ======================== %%
% Compute psi(v)
psi = -a * ((2/3) * v.^3 + 2 * beta^2 * v - (2/3) * v_0^3 - 2 * beta^2 * v_0);

% Plot
figure;
plot(x_rays, psi, 'LineWidth', 2);
xlabel('x');
ylabel('\psi(x)');
title('\psi(x)');
grid on;
