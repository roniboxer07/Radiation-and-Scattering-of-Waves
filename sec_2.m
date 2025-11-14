%% ======================== Constants & Wave Properties ======================== %%
mu0 = 4 * pi * 1e-7;           % Permeability of free space (H/m)
eps0 = 8.854e-12;              % Permittivity of free space (F/m)           
freq = 15e9;                   % Frequency (Hz)
omega = 2 * pi * freq;         % Angular frequency (rad/s)
k0 = omega * sqrt(mu0 * eps0); % Wave number in free space (1/m)
lambda_0 = (2 * pi) / k0;      % Wavelength in free space (m)
eta = sqrt(mu0 / eps0);        % Impedance of free space (Ohms)
E0 = 1;                        % Incident field magnitude (V/m)
phi_inc = 3 * pi / 2;              % Incident angle (radians)

%% ======================== Simulation Parameters ======================== %%
eps_r = 1.01; % Relative permittivity value
a_1 = [0.5 * lambda_0, 2 * lambda_0]; % Radius of cylinders

% Loop over cylinder radius
for a_plot = a_1
    % Update wave properties for cylinder
    k1 = k0 * sqrt(eps_r);                        % Wave number in cylinder 1 (1/m)
    k2 = k1;                                      % Wave number in cylinder 2 (1/m)
    lambda_1 = (2 * pi) / k1;                     % Wavelength in cylinder 1 (m)
    lambda_2 = (2 * pi) / k2;                     % Wavelength in cylinder 2 (m)
    
    % Distance between the two circles (m)
    if a_plot == 0.5 * lambda_0
        d_plot = 2 * lambda_0;
    end
    if a_plot == 2 * lambda_0
        d_plot = 8 * lambda_0;
    end
                              
    N = 1600;                                     % Number of squares in the object
    delta = sqrt(((2 * a_plot + d_plot)^2) / N);  % Side length of each element in grid
    eta1 = sqrt(mu0 ./ (eps0 .* eps_r));          % Impedance of cylinder 1 (Ohms)
    eta2 = eta1;                                  % Impedance of cylinder 2 (Ohms)

    far = 10 * lambda_0; % Far-field distance
    x = -(a_plot + far) : delta : (a_plot + far);
    z = -(a_plot + far) : delta : (a_plot + far + d_plot);
    [Z, X] = meshgrid(z, x);
                       
    % Call the main calculation function
    O_obj = O_obj_calc(a_plot, d_plot, delta, k0, k1, k2, E0,  lambda_0, eps_r);

end 

%% ======================== Function Definitions ======================== %%

%% Function: Main Calculation (Object Function via Fourier Domain) %%
function O_fft = O_obj_calc(a_plot, d_plot, delta, k0, k1, k2, E0, lambda_0, eps_r)
    
    % ----- Grids & Masks Parameters ----- %
    
    % Center coordinates of two circular dielectric cylinders
    x0 = 0; z0 = 0;         % Center of first cylinder
    x1 = 0; z1 = d_plot;    % Center of second cylinder (offset in z-direction)
    
    % Define computational grid in x and z directions
    x_obj = -(2 * a_plot + d_plot) / 2 : delta : (2 * a_plot + d_plot) / 2;
    z_obj = -a_plot : delta : (a_plot + d_plot);
    [z_grid, x_grid] = meshgrid(z_obj, x_obj);  % Simulation domain
    
    % Generate binary masks for the two circular objects
    masks = compute_masks(a_plot, x0, z0, x1, z1, x_grid, z_grid);
    mask_N = compute_mask_N(masks);            % Number of points in each object
    mask0 = masks(:,:,1);
    mask1 = masks(:,:,2);
    tot_mask = compute_tot_mask(mask0, mask1); % Combined region

    % Extract coordinates (z,x) of points inside the objects
    R = compute_circle_points(tot_mask, x_grid, z_grid);
    
    % ----- Sampling Parameters (k-space) ----- %
    
    % Number of samples
    if a_plot == 0.5 * lambda_0
        N_val = [10, 20, 25, 50];
    end
    if a_plot == 2 * lambda_0
        N_val = [70, 80, 90, 180];
    end

    for N_samples = N_val
        kz = linspace(-2 * k0, 2 * k0, N_samples);   % kz sampling range
        kx = linspace(-2 * k0, 2 * k0, N_samples);   % kx sampling range
        [KZ, KX] = meshgrid(kz, kx);                % k-space grid
        
        % ----- Scattering Angle Parameters ----- %    
        phi_s = linspace(0, 2*pi, N_samples);       % Scattering angles
        phi_inc = linspace(0, pi, N_samples);       % Incident angles
        [PHI_INC, PHI_S] = meshgrid(phi_inc, phi_s);
        
        % Flatten all angle combinations into 1D list
        Angle_comb = [reshape(PHI_INC, [], 1), reshape(PHI_S, [], 1)];
        
        % Initialize output arrays
        kzx_pair = zeros(size(Angle_comb));         % Stores (kz, kx) values
        O_fourier = zeros(size(Angle_comb, 1), 1);  % Stores Fourier values of object
        
        % ----- Precompute Static Matrices ----- %
        O_mat = compute_O_mat(k0, k1, k2, mask_N);                 % Contrast matrix
        Z_mat = compute_Z_mat(R, k0, delta);                       % Impedance matrix
        O_grid_mat = compute_O_grid_mat(k0, k1, k2, mask0, mask1, z_grid);  % Contrast over grid
        
        % ----- Main Loop: Compute FT of Object Function ----- %
        phi_inc_j = PHI_INC(:);     % Flattened incident angles
        phi_s_j = PHI_S(:);         % Flattened scattering angles
        
        for i = 1:length(phi_inc_j)
            phi_inc_val = phi_inc_j(i);
            phi_s_val = phi_s_j(i);
        
            % Compute incident field vector
            b_vec = compute_b_vec(E0, k0, phi_inc_val, R);
        
            % Solve MoM system: Z * a = b
            a_vec = compute_a_vec(Z_mat, O_mat, b_vec);
        
            % Spread solution back to grid
            A_grid_mat = compute_A_grid_mat(tot_mask, z_grid, a_vec);
        
            % Compute far-field projection for current scattering direction
            u_s_part_two = compute_U_s_part_two(x_obj, z_obj, x_grid, z_grid, A_grid_mat, O_grid_mat, k0, phi_s_val, lambda_0);
        
            % Compute spatial frequency pair (kz, kx)
            kx_s = k0 * (sin(phi_inc_val) - sin(phi_s_val));
            kz_s = k0 * (cos(phi_inc_val) - cos(phi_s_val));
            kzx_pair(i,:) = [kz_s, kx_s];
        
            % Compute Fourier coefficient from scattered field
            O_fourier(i) = compute_O_Fourier(k0, a_plot, u_s_part_two, E0, lambda_0);
        end
        
        % ----- Fill Fourier Domain Grid & Interpolate Missing Values ----- %
        O_fourier_FT = Fill_O_Fourier_missing(k0, kzx_pair, O_fourier, KX, KZ);
        
        % Plot the real part of the Fourier space
        plot_O_r(kx, kz, real(O_fourier_FT), a_plot, eps_r, lambda_0, N_samples);
        
        % ----- Inverse FFT to Reconstruct Object Function in Space ----- %
        O_fft = ifftshift(ifft2(O_fourier_FT));    % Convert back to spatial domain
        
        % Define spatial axis for plotting
        x = linspace(-2*k0, 2*k0, N_samples);
        z = linspace(-2*k0, 2*k0, N_samples);
        plot_O_zx(z, x, abs(O_fft), a_plot, eps_r, lambda_0, N_samples);
    end
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
    
%% Function: Fill Missing Parts of the Object's Fourier Transform %%
function O_fourier_FT = Fill_O_Fourier_missing(k0, kzx_pair, O_fourier, KX, KZ)
    % Apply O(-k) = conj(O(k))
    n = size(kzx_pair, 1);
    for j = 1:n
        kzx_pair(n + j, :) = -kzx_pair(j, :);       % Reflect wave vector
        O_fourier(n + j) = conj(O_fourier(j));      % Use complex conjugate
    end

    % Interpolate the known Fourier values onto a full KX-KZ grid
    F = scatteredInterpolant(kzx_pair(:,1), kzx_pair(:,2), O_fourier, 'natural');
    O_fourier_FT = F(KZ, KX);                       % Interpolated Fourier map

    % Apply a cutoff to remove components beyond 2k0
    O_fourier_FT(KZ.^2 + KX.^2 > (2 * k0)^2) = 0;   % Optional filtering
end

%% Function: Compute Fourier Coefficient from Scattered Field %%
function O = compute_O_Fourier(k0, a_plot, u_s, E0, lambda_0)
    O = u_s .* (sqrt(8 * pi * k0 * a_plot) ./ (E0 * exp(-1i * k0 * (1000 * lambda_0) - 1i * pi * 0.25)));
end

%% Function: Compute Far-Field Scattered Field for Angle φ_s %%
function integral_result = compute_U_s_part_two(x_obj, z_obj, x_grid, z_grid, A_grid, O_grid, k0, phi_s, lambda_0)
    % Far-field observation point (r ≈ 1000 * λ) in direction φ_s
    z_j = (1000 * lambda_0) * cos(phi_s);
    x_j = (1000 * lambda_0) * sin(phi_s);

    % Distance from observation point to every point in the object
    R = sqrt((x_j - x_grid).^2 + (z_j - z_grid).^2);
    R(R == 0) = eps; % Avoid singularity

    % 2D Green’s function 
    G = (-1i / 4) * besselh(0, 2, k0 * R);

    % Perform the double integral over object region
    integral_result = trapz(x_obj, trapz(z_obj, G .* O_grid .* A_grid, 2));
end

%% Function: Plot Real Part of O(kz, kx) %%
function plot_O_r(kz, kx, O, a_plot, eps_r, lambda_0, N_samples)
    s = a_plot / lambda_0;  % Normalize radius by wavelength
    figure;
    imagesc(kz, kx, real(O));
    colorbar;
    xlabel('kz');
    ylabel('kx');
    title(sprintf('Real Part of O(kz, kx) for N_s = %.f', N_samples));
    set(gca, 'YDir', 'normal');
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    % Create filename string
    filename = sprintf('O_kz_kx_real_N_%.f_eps_%.2f_rad_%.1f.png', N_samples, eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end

%% Function: Plot Reconstructed Object O(x, z) %%
function plot_O_zx(z, x, O, a_plot, eps_r, lambda_0, N_samples)
    s = a_plot / lambda_0;  % Normalize radius by wavelength
    figure;
    imagesc(z, x, real(O)); 
    colorbar;
    xlabel('z');
    ylabel('x');
    title(sprintf('Reconstructed Object O(x, z) for N_s = %.f', N_samples));
    set(gca, 'YDir', 'normal');
    subtitle(sprintf('Cylinder Radius = %.1f * \\lambda, \\epsilon_r = %.2f', s, eps_r));

    % Create filename string
    filename = sprintf('Reconstructed_O_N_%.f_eps_%.2f_rad_%.1f.png', N_samples, eps_r, s);
    
    % Save the figure
    saveas(gcf, filename);
end
