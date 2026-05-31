function results = simulateTHzFromDiffusion(params)
% SIMULATETHZFROMDIFFUSION Core physics engine for carrier dynamics and THz emission

    % --- 1. Parameter Unpacking ---
    Lx = params.Lx; Ly = params.Ly;
    Nx = params.Nx; Ny = params.Ny;

    if ~isfield(params, 'Lz'), params.Lz = 3.0e-6; end
    if ~isfield(params, 'Nz'), params.Nz = 201; end
    if ~isfield(params, 'Wdep'), params.Wdep = 0.5e-6; end
    if ~isfield(params, 'Emax'), params.Emax = 8e5; end
    if ~isfield(params, 'usePhotoDember'), params.usePhotoDember = true; end
    if ~isfield(params, 'useBuiltIn'), params.useBuiltIn = true; end

    Lz = params.Lz; Nz = params.Nz;

    D_e = params.D_e; D_h = params.D_h;
    mu_e = params.mu_e; mu_h = params.mu_h;
    tau_r = params.tau_r; alpha = params.alpha;
    w0 = params.w0; I0 = params.I0; tau_p = params.tau_p;

    Wdep = params.Wdep; Emax = params.Emax;
    X0 = params.X0; Y0 = params.Y0; z0 = params.z0;

    theta_rad = 0;
    if isfield(params, 'theta'), theta_rad = params.theta; end

    % --- 2. Spatial Grid Setup ---
    x = linspace(-Lx/2, Lx/2, Nx);
    y = linspace(-Ly/2, Ly/2, Ny);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    z = linspace(0, Lz, Nz).';
    dz = z(2) - z(1);
    [XX, YY] = meshgrid(x, y);

    cos_theta = max(cos(theta_rad), 1e-3);
    Gxy = exp(-2 * ( (XX.^2 / (w0 / cos_theta)^2) + (YY.^2 / w0^2) ));

    % --- 3. Time Grid Setup ---
    dt_fine = 5e-16;
    t_fine  = params.Tmin : dt_fine : params.Tmax;
    Nt_fine = numel(t_fine);

    Nt_out = 800;
    t_out  = linspace(params.Tmin, params.Tmax, Nt_out);
    t_ps   = t_out * 1e12;

    % --- 4. Built-in Depletion Field ---
    E_z = zeros(Nz, 1);
    mask_dep = (z <= Wdep);
    E_z(mask_dep) = -Emax * 0.5 .* (1 + cos(pi*z(mask_dep)/Wdep));

    if isfield(params, 'useBuiltIn') && ~params.useBuiltIn
        E_z(:) = 0;
    end

    % --- 5. 1D Depth Dynamics (always runs) ---
    fprintf('Simulating ultrafast 1D Z-axis dynamics...\n');
    n_1D = zeros(Nz, 1); p_1D = zeros(Nz, 1);
    Jz_1D_history = zeros(Nz, Nt_fine);
    n_surf_history = zeros(1, Nt_fine);

    q = 1.602e-19; hbar_omega = 1.55 * q;
    Ttr = 1 - ((1.0 - 3.7351)/(1.0 + 3.7351))^2;
    G0_1D = Ttr * alpha * I0 / hbar_omega;

    sigma_t = tau_p / (2*sqrt(2*log(2)));
    De_eff = D_e; Dh_eff = D_h;

    if isfield(params, 'usePhotoDember') && ~params.usePhotoDember
        De_eff = 0; Dh_eff = 0;
    end

    for k = 2:Nt_fine
        tt = t_fine(k);

        G = G0_1D * exp(-alpha*z) * exp(-((tt)^2)/(2*sigma_t^2));
        Rpair = min(n_1D, p_1D) / tau_r;

        Fn = zeros(Nz+1, 1); Fp = zeros(Nz+1, 1);
        for i = 2:Nz
            Eface = 0.5*(E_z(i-1) + E_z(i));

            vn = -mu_e * Eface;
            n_up = n_1D(i-1); if vn < 0, n_up = n_1D(i); end
            Fn(i) = vn * n_up - De_eff * (n_1D(i) - n_1D(i-1))/dz;

            vp = mu_h * Eface;
            p_up = p_1D(i-1); if vp < 0, p_up = p_1D(i); end
            Fp(i) = vp * p_up - Dh_eff * (p_1D(i) - p_1D(i-1))/dz;
        end

        divFn = (Fn(2:end) - Fn(1:end-1)) / dz;
        divFp = (Fp(2:end) - Fp(1:end-1)) / dz;

        n_1D = max(0, n_1D - dt_fine*divFn + dt_fine*(G - Rpair));
        p_1D = max(0, p_1D - dt_fine*divFp + dt_fine*(G - Rpair));

        Fn_center = 0.5*(Fn(1:end-1) + Fn(2:end));
        Fp_center = 0.5*(Fp(1:end-1) + Fp(2:end));

        Jz_1D_history(:, k) = q * (Fp_center - Fn_center);
        n_surf_history(k) = n_1D(1);
    end

    % --- 6. 3D Projection and Current Density ---
    fprintf('Projecting dynamics to 2D optical envelope...\n');
    Jz_surf_fine = sum(Jz_1D_history, 1) * dz;

    Jz_surf_out = interp1(t_fine, Jz_surf_fine, t_out);
    n_surf_out  = interp1(t_fine, n_surf_history, t_out);

    [dG_dx, dG_dy] = gradient(Gxy, dx, dy);

    Jx_all = zeros(Ny, Nx, Nt_out);
    Jy_all = zeros(Ny, Nx, Nt_out);
    Jz_all = zeros(Ny, Nx, Nt_out);

    if isfield(params, 'mechanism') && strcmp(params.mechanism, 'Custom')
        % Build real n(x,y,t) and spatial gradients from simulated carriers
        fprintf('Evaluating Custom Jx, Jy, Jz formulas with simulated carriers...\n');

        n_3D = zeros(Ny, Nx, Nt_out);
        dnx_3D = zeros(Ny, Nx, Nt_out);
        dny_3D = zeros(Ny, Nx, Nt_out);

        for k = 1:Nt_out
            n_3D(:,:,k) = n_surf_out(k) * Gxy;
            dnx_3D(:,:,k) = n_surf_out(k) * dG_dx;
            dny_3D(:,:,k) = n_surf_out(k) * dG_dy;
        end

        [YY_3D, XX_3D, TT_3D] = ndgrid(y, x, t_out);

        [Jx_all, Jy_all, Jz_all] = params.CustomJxHandle(n_3D, dnx_3D, dny_3D, XX_3D, YY_3D, TT_3D);
    else
        for k = 1:Nt_out
            Jz_all(:,:,k) = Jz_surf_out(k) * Gxy;
            Jx_all(:,:,k) = q * De_eff * n_surf_out(k) * dG_dx;
            Jy_all(:,:,k) = q * De_eff * n_surf_out(k) * dG_dy;
        end
    end

    % GUI Plotting Variables
    results.t_ps = t_ps; results.x = x; results.y = y;
    results.n_center_time = n_surf_out;

    [~, peak_idx] = max(n_surf_out);
    results.peak_center_time_ps = t_ps(peak_idx);
    results.n_at_peak_center = n_surf_out(peak_idx) * Gxy;

    % --- 7. Far-Field Calculation ---
    fprintf('Computing far-field radiation...\n');
    [Ex, Ey, Ez] = compute_THz_E_from_Jx(x, y, t_out, Jx_all, Jy_all, Jz_all, X0, Y0, z0);

    results.Ex_pt = Ex; results.Ey_pt = Ey; results.Ez_pt = Ez;
    results.E_mag = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
    results.Jx_all = Jx_all; results.Jy_all = Jy_all; results.Jz_all = Jz_all;
    results.X0 = X0; results.Y0 = Y0; results.z0 = z0;
    results.theta = theta_rad;

    % --- 8. 2D Spatial Radiation Calculation ---
    fprintf('Computing 2D spatial radiation for vector map...\n');

    [~, peak_E_idx] = max(results.E_mag);
    t_obs_peak = t_out(peak_E_idx);

    grid_size = 50e-6;
    N_grid = 11;
    x0_vec = linspace(-grid_size, grid_size, N_grid);
    y0_vec = linspace(-grid_size, grid_size, N_grid);
    [X0_grid, Y0_grid] = meshgrid(x0_vec, y0_vec);

    [Ex_2D, Ey_2D] = compute_THz_E_spatial_2D(x, y, t_out, Jx_all, Jy_all, Jz_all, X0_grid, Y0_grid, z0, t_obs_peak);

    results.X0_grid = X0_grid;
    results.Y0_grid = Y0_grid;
    results.Ex_2D = Ex_2D;
    results.Ey_2D = Ey_2D;
    results.t_obs_peak = t_obs_peak;

    fprintf('Simulation complete.\n');
end
