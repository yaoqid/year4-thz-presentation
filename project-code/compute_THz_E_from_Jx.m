function [Ex, Ey, Ez] = compute_THz_E_from_Jx(x, y, t, Jx, Jy, Jz, X0, Y0, z0)
% COMPUTE_THZ_E_FROM_JX Vectorized far-field THz emission calculator
% Uses retarded potentials to compute E-field from 3D current density

    % --- Constants & Grid Setup ---
    c0  = 3e8;
    mu0 = 4*pi*1e-7;
    Nx = numel(x); Ny = numel(y); Nt = numel(t);
    dx = x(2) - x(1); dy = y(2) - y(1); dt = t(2) - t(1);

    % --- 1. Compute Time Derivatives of Current Density ---
    dJx = zeros(size(Jx)); dJy = zeros(size(Jy)); dJz = zeros(size(Jz));
    dJx(:,:,2:Nt-1) = (Jx(:,:,3:Nt) - Jx(:,:,1:Nt-2)) / (2*dt);
    dJy(:,:,2:Nt-1) = (Jy(:,:,3:Nt) - Jy(:,:,1:Nt-2)) / (2*dt);
    dJz(:,:,2:Nt-1) = (Jz(:,:,3:Nt) - Jz(:,:,1:Nt-2)) / (2*dt);

    % --- 2. Compute Geometric Factors ---
    [XXs, YYs] = meshgrid(x, y);
    Rx = X0 - XXs; Ry = Y0 - YYs; Rz = z0;
    R  = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
    
    R(R==0) = eps; % Prevent division by zero
    Rhat_x = Rx ./ R; Rhat_y = Ry ./ R; Rhat_z = Rz ./ R;

    Ex = zeros(1,Nt); Ey = zeros(1,Nt); Ez = zeros(1,Nt);

    % Pre-generate spatial grid indices and matrices to optimize loop
    [grid_y, grid_x] = ndgrid(1:Ny, 1:Nx); 
    R_mat = R;
    Rhat_x_mat = Rhat_x; Rhat_y_mat = Rhat_y; Rhat_z_mat = Rhat_z;
    
    % --- 3. Vectorized Integration Loop ---
    for k = 1:Nt
        % Calculate retarded time and map to discrete time indices
        t_ret = t(k) - R_mat / c0; 
        tau   = (t_ret - t(1))/dt + 1;
        k1    = floor(tau);
        k2    = k1 + 1;
        alpha = tau - k1;
        
        % Boundary check mask for valid retarded times
        valid = (k1 >= 1) & (k2 <= Nt);
        
        % Clamp indices to prevent out-of-bounds errors during interpolation
        k1_safe = max(1, min(Nt, k1));
        k2_safe = max(1, min(Nt, k2));
        
        % Fast linear indexing for 3D matrices
        idx1 = grid_y + (grid_x - 1) * Ny + (k1_safe - 1) * (Ny * Nx);
        idx2 = grid_y + (grid_x - 1) * Ny + (k2_safe - 1) * (Ny * Nx);
        
        % Linear interpolation of current derivatives at retarded time
        dJx_ret = (1-alpha) .* dJx(idx1) + alpha .* dJx(idx2);
        dJy_ret = (1-alpha) .* dJy(idx1) + alpha .* dJy(idx2);
        dJz_ret = (1-alpha) .* dJz(idx1) + alpha .* dJz(idx2);
        
        % Zero out contributions from invalid retarded times
        dJx_ret(~valid) = 0; dJy_ret(~valid) = 0; dJz_ret(~valid) = 0;
        
        % --- 4. Cross-Product and Integration ---
        % Using E_rad proportional to ( Rhat*(Rhat·dJ) - dJ ) / R
        dJdot = Rhat_x_mat .* dJx_ret + Rhat_y_mat .* dJy_ret + Rhat_z_mat .* dJz_ret;
        
        Tx = dJdot .* Rhat_x_mat - dJx_ret;
        Ty = dJdot .* Rhat_y_mat - dJy_ret;
        Tz = dJdot .* Rhat_z_mat - dJz_ret;
        
        Ex(k) = (mu0/(4*pi)) * sum(sum( Tx ./ R_mat )) * dx * dy;
        Ey(k) = (mu0/(4*pi)) * sum(sum( Ty ./ R_mat )) * dx * dy;
        Ez(k) = (mu0/(4*pi)) * sum(sum( Tz ./ R_mat )) * dx * dy;
    end
end