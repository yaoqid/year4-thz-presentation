function [Ex_2D, Ey_2D] = compute_THz_E_spatial_2D(x, y, t, Jx, Jy, Jz, X0_grid, Y0_grid, z0, t_obs)
    % Computes Ex and Ey over a 2D observation grid at a specific snapshot in time.
    c0 = 3e8; 
    mu0 = 4*pi*1e-7;
    dx = x(2)-x(1); 
    dy = y(2)-y(1); 
    dt = t(2)-t(1);
    Nx = numel(x); 
    Ny = numel(y); 
    Nt = numel(t);

    % 1. Compute Time Derivatives (Matches your single-point logic)
    dJx = zeros(size(Jx)); dJy = zeros(size(Jy)); dJz = zeros(size(Jz));
    dJx(:,:,2:Nt-1) = (Jx(:,:,3:Nt) - Jx(:,:,1:Nt-2)) / (2*dt);
    dJy(:,:,2:Nt-1) = (Jy(:,:,3:Nt) - Jy(:,:,1:Nt-2)) / (2*dt);
    dJz(:,:,2:Nt-1) = (Jz(:,:,3:Nt) - Jz(:,:,1:Nt-2)) / (2*dt);

    [XXs, YYs] = meshgrid(x, y);

    Ex_2D = zeros(size(X0_grid));
    Ey_2D = zeros(size(X0_grid));

    % 2. Loop over the coarse observation grid
    for i = 1:size(X0_grid, 1)
        for j = 1:size(X0_grid, 2)
            X0 = X0_grid(i,j);
            Y0 = Y0_grid(i,j);

            % Distance from source plane to this specific observation point
            R = sqrt((X0 - XXs).^2 + (Y0 - YYs).^2 + z0.^2);
            R(R==0) = eps;
            Rhat_x = (X0 - XXs)./R; Rhat_y = (Y0 - YYs)./R; Rhat_z = z0./R;

            % Retarded time mapping
            t_ret = t_obs - R/c0;
            tau = (t_ret - t(1))/dt + 1;
            k1 = floor(tau);
            k2 = k1 + 1;
            alpha = tau - k1;

            valid = (k1 >= 1) & (k2 <= Nt);
            k1_safe = max(1, min(Nt, k1));
            k2_safe = max(1, min(Nt, k2));

            % Interpolate current derivatives at retarded time
            dJx_ret = zeros(Ny, Nx); dJy_ret = zeros(Ny, Nx); dJz_ret = zeros(Ny, Nx);
            
            % Using a nested loop here for the 2D slice interpolation to conserve memory
            for row = 1:Ny
                for col = 1:Nx
                    if valid(row, col)
                        idx1 = k1_safe(row, col); 
                        idx2 = k2_safe(row, col);
                        a = alpha(row, col);
                        dJx_ret(row,col) = (1-a)*dJx(row,col,idx1) + a*dJx(row,col,idx2);
                        dJy_ret(row,col) = (1-a)*dJy(row,col,idx1) + a*dJy(row,col,idx2);
                        dJz_ret(row,col) = (1-a)*dJz(row,col,idx1) + a*dJz(row,col,idx2);
                    end
                end
            end

            % Vector Projection & Integration
            dJdot = Rhat_x .* dJx_ret + Rhat_y .* dJy_ret + Rhat_z .* dJz_ret;
            Tx = dJdot .* Rhat_x - dJx_ret;
            Ty = dJdot .* Rhat_y - dJy_ret;

            Ex_2D(i,j) = (mu0/(4*pi)) * sum(sum(Tx ./ R)) * dx * dy;
            Ey_2D(i,j) = (mu0/(4*pi)) * sum(sum(Ty ./ R)) * dx * dy;
        end
    end
end