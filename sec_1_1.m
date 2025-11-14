%% ======================== Constants & Wave Properties ======================== %%
mu0 = 4 * pi * 1e-7;           % Permeability of free space (H/m)
eps0 = 8.854e-12;              % Permittivity of free space (F/m)           
freq = 15e9;                   % Frequency (Hz)
omega = 2 * pi * freq;         % Angular frequency (rad/s)
k0 = omega * sqrt(mu0 * eps0); % Wave number in free space (1/m)
lambda_0 = (2 * pi) / k0;      % Wavelength in free space (m)
eta = sqrt(mu0 / eps0);        % Impedance of free space (Ohms)
E0 = 1;                        % Incident field magnitude (V/m)
phi_inc = 3 * pi / 2;          % Incident angle (radians)

%% ======================== Simulation Parameters ======================== %%
eps_r_vec = [1.01, 1.1, 2]; % Relative permittivity values
a_1 = [0.1 * lambda_0, 1.5 * lambda_0, 3.5 * lambda_0, 5 * lambda_0]; % Radius of cylinders

% Loop over permittivity and cylinder radius
for eps_r = eps_r_vec
    for a_plot = a_1
        % Update wave properties for cylinder
        k1 = k0 * sqrt(eps_r);                        % Wave number in cylinder 1 (1/m)
        k2 = k1;                                      % Wave number in cylinder 2 (1/m)
        lambda_1 = (2 * pi) / k1;                     % Wavelength in cylinder 1 (m)
        lambda_2 = (2 * pi) / k2;                     % Wavelength in cylinder 2 (m)
        d_plot = 3 * a_plot;                          % Distance between the two circles (m)
        
        % Number of squares in the object corresponding to the size of the object
        if (a_plot / lambda_0) > 2
            N = (a_plot / lambda_0) * 1000;      
        else
            N = 1600;
        end

        delta = sqrt(((2 * a_plot + d_plot)^2) / N);  % Side length of each element in grid
        eta1 = sqrt(mu0 ./ (eps0 .* eps_r));          % Impedance of cylinder 1 (Ohms)
        eta2 = eta1;                                  % Impedance of cylinder 2 (Ohms)
        
        far = 10 * lambda_0; % Far-field distance
        x = -(a_plot + far) : delta : (a_plot + far);
        z = -(a_plot + far) : delta : (a_plot + far + d_plot);
        [Z, X] = meshgrid(z, x);

        % Call the MoM calculation function
        U_s_MoM = MoM_calc(x, z, Z, a_plot, d_plot, delta, k0, k1, k2, E0, phi_inc, lambda_0, eps_r);

    end 
end

%% ======================== Function Definitions ======================== %%

%% Function: Main Calculation %%
function U_s_MoM = MoM_calc(x, z, Z, a_plot, d_plot, delta, k0, k1, k2, E0, phi_inc, lambda_0, eps_r)
    % ----- Grids & Masks Parameters ----- %

    % Center coordinates of two circular dielectric cylinders
    x0 = 0; z0 = 0;         % Center of first cylinder
    x1 = 0; z1 = d_plot;    % Center of second cylinder (offset in z-direction)
    
    % Define computational grid in x and z directions
    x_obj = -(2 * a_plot + d_plot) / 2 : delta : (2 * a_plot + d_plot) / 2;
    z_obj = -a_plot : delta : (a_plot + d_plot);
    [z_grid, x_grid] = meshgrid(z_obj, x_obj);  % 2D grid for the simulation domain
    
    % Generate binary masks representing the circular objects
    % mask0 and mask1 correspond to each cylinder; tot_mask is the union
    masks = compute_masks(a_plot, x0, z0, x1, z1, x_grid, z_grid);
    mask_N = compute_mask_N(masks); % Number of points in each mask
    mask0 = masks(:,:,1);
    mask1 = masks(:,:,2);
    tot_mask = compute_tot_mask(mask0, mask1);

    % ----- Compute & Plot Circular Grids ----- %
    R = compute_circle_points(tot_mask, x_grid, z_grid); % Extract the coordinates (R) of points inside the combined mask
    
    % Plot the circle grid
    if (a_plot == 3.5 * lambda_0) && (eps_r == 1.1)
        plot_grid(R, x_grid, z_grid, tot_mask);
    end

    % ----- Compute & Plot the Scattered Field (MoM) ----- %

    % Compute the incident field at object points
    b_vec = compute_b_vec(E0, k0, phi_inc, R);

    % Build the object matrix
    O_mat = compute_O_mat(k0, k1, k2, mask_N);

    % Compute the Z matrix 
    Z_mat = compute_Z_mat(R, k0, delta);

    % Compute the a vector 
    a_vec = compute_a_vec(Z_mat, O_mat, b_vec);

    % Map a_vec back onto the full 2D grid
    A_grid_mat = compute_A_grid_mat(tot_mask, z_grid, a_vec);

    % Build the object grid matrix
    O_grid_mat = compute_O_grid_mat(k0, k1, k2, mask0, mask1, z_grid);
    
    % Compute the scattered field using MoM
    U_s_MoM = compute_U_s(x_obj, z_obj, x_grid, z_grid, x, z, Z, A_grid_mat, O_grid_mat, k0);
    
    % Plot scattered field
    plot_U_MoM(x, z, U_s_MoM, a_plot, eps_r, lambda_0);

    % Compute and plot the RCS 
    compute_RCS(E0, x_obj, z_obj, x_grid, z_grid, A_grid_mat, O_grid_mat, k0, lambda_0, a_plot, eps_r, 'using MoM calculation')
end

%% Function: Compute Masking %%
function masks = compute_masks(a, x0, z0, x1, z1, x_grid, z_grid)
    % Create binary masks for two circular regions of radius 'a'
    % mask0 is centered at (x0, z0)
    mask0 = sqrt((x_grid - x0).^2 + (z_grid - z0).^2) <= a;
    
    % mask1 is centered at (x1, z1)
    mask1 = sqrt((x_grid - x1).^2 + (z_grid - z1).^2) <= a;
    
    % Stack both masks along the 3rd dimension into one 3D array
    masks = cat(3, mask0, mask1);
end

%% Function: Combine Masks %%
function tot_mask = compute_tot_mask(mask0, mask1)
    % Combine two masks into one using logical OR
    % This gives a total region covered by either mask
    tot_mask = logical(mask0 + mask1);
end 

%% Function: Count Masked Points %%
function mask_N = compute_mask_N(masks)
    mask_N = [sum(masks(:,:,1), 'all'), sum(masks(:,:,2), 'all')];
end

%% Function: Compute Circular Grid Points %%
function R = compute_circle_points(tot_mask, x_grid, z_grid)
    % Extract the (z, x) coordinates of grid points inside the total mask
    R = [z_grid(tot_mask), x_grid(tot_mask)]; 
end

%% Function: Compute Incident Field Vector inside V %%
function b = compute_b_vec(E0, k0, phi_inc, R)
    % Compute the incident field vector 'b' at each point R inside the region V
    b = E0 * exp(-1j * k0 * (R(:,2) * sin(phi_inc) + R(:,1) * cos(phi_inc)));
end

%% Function: Compute Z Matrix %%
function Z = compute_Z_mat(R, k0, delta)
    N = length(R);           % Number of points inside the object region
    Z = zeros(N, N);         % Initialize Z matrix (NxN)

    for m = 1:N
        for n = 1:N
            if m ~= n
                % Off-diagonal terms: approximate interaction between different points
                Z(m, n) = delta^2 * compute_g(R(m), R(n), k0);  % G(r_m, r_n)
                
                % Replace with 0 if any NaN values are detected
                if isnan(real(Z(m, n))) || isnan(imag(Z(m, n)))
                    Z(m, n) = 0;
                end
            else
                % Diagonal terms: compute self-interaction using double integral
                integrand = @(z, x) (-1i / 4) * besselh(0, 2, k0 * sqrt(z.^2 + x.^2));
                z_mm = integral2(integrand, -delta/2, delta/2, -delta/2, delta/2);

                if isnan(real(z_mm)) || isnan(imag(z_mm))
                    z_mm = 0;
                end

                Z(m, m) = z_mm;
            end
        end
    end
end

%% Function: Green's Function Calculation %%
function g = compute_g(r, r_prime, k)
    distance = norm(r - r_prime);           % Distance between points
    H0 = besselh(0, 2, k * distance);       % Hankel function of the second kind
    g = -1j / 4 * H0;                       % 2D scalar Green's function
end

%% Function: Compute Object Function O %%
function O = compute_O_mat(k0, k1, k2, mask_N)
    % Create a diagonal matrix representing the contrast in permittivity
    O_1 = (k1^2 - k0^2) * ones(1, mask_N(1));  % Region 1
    O_2 = (k2^2 - k0^2) * ones(1, mask_N(2));  % Region 2
    O_vec = [O_1, O_2];                        % Combine into one vector

    O = diag(O_vec);                           % Diagonal contrast matrix
end

%% Function: Compute Vector a %%
function a = compute_a_vec(Z_mat, O, b)
    I = eye(size(Z_mat));           % Identity matrix
    M = I - Z_mat * O;              % MoM matrix (I - Z * O)

    % Check for singularity before solving
    if det(M) == 0
        error('Matrix (I - Z * O) is singular and cannot be inverted.');
    end

    a = M \ b;                      % Solve linear system for unknowns
end

%% Function: Compute A_grid_mat %%
function A = compute_A_grid_mat(tot_mask, z_grid, a_vec)
    A = zeros(size(z_grid));       % Initialize full matrix with zeros
    A(tot_mask) = a_vec;           % Assign computed values only inside the object region
end

%% Function: Compute O_grid_mat %%
function O = compute_O_grid_mat(k0, k1, k2, mask0, mask1, z_grid)
    O = zeros(size(z_grid));           % Initialize with zeros (background)
    O(mask0) = k1^2 - k0^2;            % Assign contrast for first object
    O(mask1) = k2^2 - k0^2;            % Assign contrast for second object
end

%% Function: Compute Scattered Field U_s %%
function integral_result = compute_U_s(x_obj, z_obj, x_grid, z_grid, x, z, Z, A_grid, O_grid, k0)
    % Initialize output matrix for scattered field
    integral_result = zeros(size(Z));
    
    % Loop over each observation point (x_j, z_j)
    for z_i = 1:length(z)
        z_j = z(z_i); 
        for x_i = 1:length(x)
            x_j = x(x_i);

            % Compute distance R between observation point and all source points
            R = sqrt((x_j - x_grid).^2 + (z_j - z_grid).^2);
            R(R == 0) = eps;  % Avoid singularity

            % Green’s function (2D Helmholtz) between source and observation point
            G = (-1i / 4) * besselh(0, 2, k0 * R);

            % Compute the double integral: ∫∫ G * O * A dx dz
            % Using trapezoidal numerical integration
            integral_value = trapz(x_obj, trapz(z_obj, G .* O_grid .* A_grid, 2));

            % Store result in corresponding grid location
            integral_result(x_i, z_i) = integral_value;

            % Optional: debug print if result is invalid
            if isnan(integral_value)
                disp("fail");
            end
        end
    end
end

%% Function: Compute and Plot RCS %%
function RCS = compute_RCS(E0, x_obj, z_obj, x_grid, z_grid, A_grid, O_grid, k0, lambda_0,  a_plot, eps_r, titleText)
    r = 500 * lambda_0;                          % Far-field observation radius
    angles = linspace(0, 2 * pi, 500);           % Observation angles in radians

    % Inline function to compute RCS at a single angle
    compute_rcs = @(phi) calculate_rcs(phi, r, x_grid, z_grid, x_obj, z_obj, O_grid, A_grid, E0, k0);
    
    % Evaluate RCS at all angles using vectorized array function
    RCS = arrayfun(compute_rcs, angles);

    % Create polar plot of RCS magnitude
    figure;
    polarplot(angles, abs(RCS), 'Color', 'k');

    s = a_plot / lambda_0; % Normalize radius by wavelength
    % Add subtitle with radius and eps_r values
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    title(['RCS ', titleText]);
    grid on;

    % Create filename string
    filename = sprintf('RCS_MoM_eps_%.2f_rad_%.1f.png', eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end

function rcs_value = calculate_rcs(angle, r, x_grid, z_grid, x_obj, z_obj, O_grid, A_grid, E0, k0)
    % Convert polar coordinates (angle, r) to Cartesian (x, z)
    [z, x] = pol2cart(angle, r);

    % Compute distance from all points in the object region to observation point
    dist = sqrt((x - x_grid).^2 + (z - z_grid).^2);

    % Green’s function between each source point and the far-field point
    green_func = (-1i / 4) * besselh(0, 2, k0 * dist);

    % Compute the scattered field at angle 'angle' using MoM solution
    U_s = trapz(x_obj, trapz(z_obj, green_func .* O_grid .* A_grid, 2));

    % Compute RCS using far-field formula: σ = 2πr * |U_s / E0|^2
    rcs_value = 2 * pi * r * (abs(U_s / E0)^2);
end

%% Functions: Plot Circular Grids %%
function plot_grid(R, x_grid, z_grid, mask)
    figure;
    hold on;
    
    % Plot points inside circles (blue)
    scatter(R(:,1), R(:,2), 120, 's', 'filled', 'MarkerEdgeColor', 'W', 'MarkerFaceColor', [179/255 208/255 255/255]);
    
    % Plot points outside circles (red)
    scatter(z_grid(~mask), x_grid(~mask), 120, 's', 'filled', 'MarkerEdgeColor', 'W', 'MarkerFaceColor', [255/255 102/255 102/255]);

    axis equal;
    grid on;
    xlabel('z');
    ylabel('x');
    title('Object Grid for Radius = 3.5\lambda');    
    hold off;
end

%% Function: Plot U in MoM %%
function plot_U_MoM(z, x, U, a_plot, eps_r, lambda_0)
    s = a_plot / lambda_0; % Normalize radius by wavelength
    figure; 
    imagesc(z, x, real(U)); 
    colorbar; 
    xlabel('z'); 
    ylabel('x'); 
    title('Real Part of U_{s} in MoM'); 
    set(gca, 'YDir', 'normal'); 
    
    % Add subtitle with radius and eps_r values
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    % Create filename string
    filename = sprintf('U_MoM_eps_%.2f_rad_%.1f.png', eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end