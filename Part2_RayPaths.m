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

%% ======================== Ray Tracing Parameters ======================== %%
v = linspace(-5, 5, 1000);            % X range
theta_deg = linspace(-90, 90, 25);   % Incident angles
theta_rad = deg2rad(theta_deg);      % Convert to radians
x_0 = 2;
n_x0 = sqrt(x_0/a + n_0^2);          % Index at x_0
beta_all = n_x0 * sin(theta_rad);    % Each value of beta corresponds to an angle

%% ======================== Ray Plotting ======================== %%
figure; 
hold on;

for i = 1:length(beta_all)
    beta = beta_all(i);
 
    v_0 = sqrt(n_0^2 + x_0/a - beta.^2);
    x_rays = a .* (v.^2 - n_0^2 + beta.^2);
    z_rays = 2 .* a .* beta .* (v + v_0);
    
    plot(x_rays, z_rays, 'DisplayName', sprintf('\\theta = %.1f°', theta_deg(i)));

    % ======= Plot the black dot at x_t ======= 
    x_t = a * (beta^2 - n_0^2);
    v_t = sqrt(n_0^2 + x_t/a - beta.^2); 
    z_t = 2 .* a .* beta .* (v_t + v_0);
    plot(x_t, z_t, 'ko', 'MarkerFaceColor', 'k');  % black filled circle
end

xlabel('x', 'FontSize', 12);
ylabel('z', 'FontSize', 12);
title('Ray path for Different Angles', 'FontSize', 14);
legend show;
grid on;