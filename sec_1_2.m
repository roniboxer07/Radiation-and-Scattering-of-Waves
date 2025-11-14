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
        
        % Call the main calculation function
        Born_Rytov_calc(x, z, X, Z, a_plot, d_plot, delta, k0, k1, k2, E0, phi_inc, lambda_0, eps_r);
    end 
end

%% ======================== Function Definitions ======================== %%

%% Function: Main Calculation %%
function Born_Rytov_calc(x, z, X, Z, a_plot, d_plot, delta, k0, k1, k2, E0, phi_inc, lambda_0, eps_r) 
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
    mask0 = masks(:,:,1);
    mask1 = masks(:,:,2);
    tot_mask = compute_tot_mask(mask0, mask1);

    % ----- Compute & Plot Circular Grids ----- %
    R = compute_circle_points(tot_mask, x_grid, z_grid); % Extract the coordinates (R) of points inside the combined mask

    % ----- Compute & Plot the Scattered Field ----- %
    
    % Compute the incident field at scatterer points (R)
    b_vec = compute_b_vec(E0, k0, phi_inc, R);

    % Expand b_vec to the full grid using the total mask
    B_grid_mat = compute_B_grid_mat(tot_mask, z_grid, b_vec);

    % Create object matrix
    O_grid_mat = compute_O_grid_mat(k0, k1, k2, mask0, mask1, z_grid);
    
    % Compute the scattered field using the Born approximation
    U_s_Born =  compute_U_s_Born(x_obj, z_obj,x_grid, z_grid, x, z, Z, B_grid_mat, O_grid_mat,k0);
    
    % Plot the scattered field using the Born approximation
    plot_U(x, z, U_s_Born, a_plot, eps_r, lambda_0 , 'Born approximation');
    
    % Compute and plot the RCS using Born
    compute_RCS_Born(E0, x_obj, z_obj, x_grid, z_grid, B_grid_mat, O_grid_mat, k0, lambda_0,  a_plot, eps_r, 'Born approximation')
    
    % Compute the scattered field using the Rytov approximation
    U_inc = compute_U_inc_vec(E0, k0, phi_inc, X, Z); % Compute the total incident field at the observation points
    U_s_Rytov = compute_U_s_Rytov(U_inc, U_s_Born); % Compute Rytov scattered field approximation using Born and incident field

    % Plot the scattered field using the Rytov approximation
    plot_U(x, z, U_s_Rytov, a_plot, eps_r, lambda_0 , 'Rytov approximation');

    % Compute and plot the RCS using Rytov
    compute_RCS_Rytov(phi_inc, E0, x_obj, z_obj, x_grid, z_grid, B_grid_mat, O_grid_mat, k0, lambda_0,  a_plot, eps_r, 'Rytov approximation')
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

%% Function: Compute Incident Field Vector %%
function U_inc = compute_U_inc_vec(E0, k0, phi_inc, X, Z)
    % Compute the incident field U_inc over the full grid (X, Z)
    U_inc = E0 * exp(-1j * k0 * (X * sin(phi_inc) + Z * cos(phi_inc)));
end

%% Function: Compute B_grid_mat %%
function B = compute_B_grid_mat(tot_mask, z_grid, b_vec)
    % Create a full 2D grid matrix 'B' and assign the values of b_vec to
    B = zeros(size(z_grid));        % Initialize with zeros
    B(tot_mask) = b_vec;            % Assign incident field values inside the mask
end

%% Function: Compute O_grid_mat %%
function O = compute_O_grid_mat(k0, k1, k2, mask0, mask1, z_grid)
    % Compute the objevt matrix 'O'
    % For each region (mask0 and mask1), assign the difference in squared wavenumbers
    O = zeros(size(z_grid));        % Initialize with zeros
    O(mask0) = k1^2 - k0^2;         % Region 0: contrast k1^2 - k0^2
    O(mask1) = k2^2 - k0^2;         % Region 1: contrast k2^2 - k0^2
end

%% Function: Compute U_s_Born %%
function integral_result = compute_U_s_Born(x_obj, z_obj, x_grid, z_grid, x, z, Z, B_grid, O_grid, k0)
    % Initialize output matrix with zeros; same size as Z (observation grid)
    integral_result = zeros(size(Z));
    
    % Loop over all observation points (x_j, z_j) where field is to be computed
    for z_i = 1:length(z)
        z_j = z(z_i);  % z-coordinate of the observation point
        
        for x_i = 1:length(x)
            x_j = x(x_i);  % x-coordinate of the observation point

            % Compute distance R between current observation point and all object points
            R = sqrt((x_j - x_grid).^2 + (z_j - z_grid).^2);
            R(R == 0) = eps;  % Avoid division by zero or singularity at R = 0

            % Green's function
            G = (-1i / 4) * besselh(0, 2, k0 * R);

            % Compute the integral using the Born approximation
            integral_value = trapz(x_obj, trapz(z_obj, G .* O_grid .* B_grid, 2));

            % Store the result in the appropriate (x,z) index
            integral_result(x_i, z_i) = integral_value;

            % Debug - show a message if the result is NaN
            if isnan(integral_value)
                disp("fail");
            end
        end
    end
end

%% Function: Compute U_s_Rytov %%
function U_s_Rytov = compute_U_s_Rytov(U_inc, U_s_Born)
    % This formula transforms the Born field into the Rytov scattered field.
    U_s_Rytov = U_inc .* (exp(U_s_Born ./ U_inc) - 1);
end

%% Function: Compute RCS using Born Approximation %%
function RCS = compute_RCS_Born(E0, x_obj, z_obj, x_grid, z_grid, B_grid, O_grid, k0, lambda_0, a_plot, eps_r, titleText)
    r = 1000 * lambda_0;  % Far-field observation distance
    angles = linspace(0, 2 * pi, 500);  % Observation angles from 0 to 2π

    % Define inline function to compute RCS at a given angle
    compute_rcs = @(phi) calculate_rcs_Born(phi, r, x_grid, z_grid, x_obj, z_obj, O_grid, B_grid, E0, k0);
    
    % Compute RCS for each angle
    RCS = arrayfun(compute_rcs, angles);

    % Plot polar diagram of RCS
    figure;
    polarplot(angles, abs(RCS), 'Color', 'k');

    s = a_plot / lambda_0; % Normalize radius by wavelength
    % Add subtitle with radius and eps_r values
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    title(['RCS using ', titleText]);
    grid on;

    % Extract the first word from titleText
    firstWord = regexp(titleText, '\s', 'split');
    firstWord = firstWord{1};

    % Create filename string
    filename = sprintf('RCS_%s_eps_%.2f_rad_%.1f.png', firstWord, eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end

function rcs_value = calculate_rcs_Born(angle, r, x_grid, z_grid, x_obj, z_obj, O_grid, B_grid, E0, k0)
    % Convert polar coordinates (angle, radius) to Cartesian (x,z)
    [z, x] = pol2cart(angle, r);

    % Calculate distance between observation point and all grid points
    dist = sqrt((x - x_grid).^2 + (z - z_grid).^2);

    % Green's function at observation point
    green_func = (-1i / 4) * besselh(0, 2, k0 * dist);

    % Compute scattered field using Born approximation
    U_s = trapz(x_obj, trapz(z_obj, green_func .* O_grid .* B_grid, 2));

    % Compute RCS value using the far-field expression
    rcs_value = 2 * pi * r * (abs(U_s / E0) .^ 2);
end


%% Function: Compute RCS using Rytov Approximation %%
function RCS = compute_RCS_Rytov(phi_inc, E0, x_obj, z_obj, x_grid, z_grid, B_grid, O_grid, k0, lambda_0, a_plot, eps_r, titleText)
    r = 1000 * lambda_0;  % Far-field observation distance
    angles = linspace(0, 2 * pi, 1000);  % Observation angles (finer resolution)

    % Define inline function to compute RCS at a given angle using Rytov
    compute_rcs = @(phi) calculate_rcs_Rytov(phi_inc, phi, r, x_grid, z_grid, x_obj, z_obj, O_grid, B_grid, E0, k0);
    
    % Compute RCS for each angle
    RCS = arrayfun(compute_rcs, angles);

    % Plot RCS polar diagram
    figure;
    polarplot(angles, abs(RCS), 'Color', 'k');

    s = a_plot / lambda_0; % Normalize radius by wavelength
    % Add subtitle with radius and eps_r values
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    title(['RCS using ', titleText]);
    grid on;

    % Extract the first word from titleText
    firstWord = regexp(titleText, '\s', 'split');
    firstWord = firstWord{1};

    % Create filename string
    filename = sprintf('RCS_%s_eps_%.2f_rad_%.1f.png', firstWord, eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end


function rcs_value = calculate_rcs_Rytov(phi_inc, angle, r, x_grid, z_grid, x_obj, z_obj, O_grid, B_grid, E0, k0)
    % Convert polar coordinates to Cartesian
    [z, x] = pol2cart(angle, r);

    % Distance from each grid point to observation point
    dist = sqrt((x - x_grid).^2 + (z - z_grid).^2);

    % Green's function evaluated at observation point
    green_func = (-1i / 4) * besselh(0, 2, k0 * dist);

    % Scattered field using Born approximation
    U_s_Born = trapz(x_obj, trapz(z_obj, green_func .* O_grid .* B_grid, 2));

    % Compute incident field at far-field point (single complex value)
    U_inc = E0 * exp(-1j * k0 * (x * sin(phi_inc) + z * cos(phi_inc)));

    % Apply Rytov approximation formula to get total scattered field
    U_s = U_inc * (exp(U_s_Born / U_inc) - 1); 

    % Compute RCS value
    rcs_value = 2 * pi * r * (abs(U_s / E0) .^ 2);
end


%% Function: Plot U %%
function plot_U(z, x, U, a_plot, eps_r, lambda_0, titleText)
    s = a_plot / lambda_0; % Normalize radius by wavelength
    figure; 
    imagesc(z, x, real(U)); 
    colorbar; 
    xlabel('z'); 
    ylabel('x'); 
    title(['Real Part of U_{s} using', titleText]); 
    set(gca, 'YDir', 'normal'); 
    
    % Add subtitle with radius and eps_r values
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    % Extract the first word from titleText
    firstWord = regexp(titleText, '\s', 'split');
    firstWord = firstWord{1};

    % Create filename string
    filename = sprintf('U_%s_eps_%.2f_rad_%.1f.png', firstWord, eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end

