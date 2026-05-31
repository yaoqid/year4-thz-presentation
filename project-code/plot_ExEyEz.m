%% plot_ExEyEz.m
% Plots Ex(t), Ey(t), Ez(t) at a user-specified observation point.
% All other parameters match the demo.mlapp GUI defaults.
%
% === CHANGE THESE TO MATCH YOUR GUI SETTINGS ===
X0 = 0.003;      % lateral x position (m)   [GUI default: 3.2 mm]
Y0 = -0.003;            % lateral y position (m)   [GUI default: 0]
z0 = -0.007;      % vertical distance  (m)   [GUI default: 7.73 mm]

% Mechanism: 'BuiltIn', 'PhotoDember', or 'Both'
% GUI default is 'BuiltIn' (radio button selection)
mechanism = 'BuiltIn';
% ================================================

clear params; clc; close all;

%% --- 1. GUI default parameters ---
q  = 1.602e-19;
kB = 1.380649e-23;
T  = 300;

params.Nx  = 601;
params.Ny  = 601;
params.Lx  = 1.5e-3;
params.Ly  = 1.5e-3;

params.mu_e = 0.8;
params.mu_h = 0.1 * params.mu_e;
params.D_e  = params.mu_e * kB * T / q;
params.D_h  = params.mu_h * kB * T / q;

params.tau_r = 1e-10;
params.alpha = 1 / 1.4e-7;
params.w0    = 0.000424 / 2;

FWHM   = 1e-13;
params.tau_p = FWHM / sqrt(2*log(2));

F_Jcm2 = 0.0001;
F_Jm2  = F_Jcm2 * 1e4;
params.I0 = F_Jm2 / (4 * params.tau_p);

params.Tmin = -5e-12;
params.Tmax = 40e-12;

params.Wdep  = 0.5e-6;
params.Emax  = 8e5;
params.Lz    = 3.0e-6;
params.Nz    = 201;

switch mechanism
    case 'BuiltIn'
        params.usePhotoDember = false;
        params.useBuiltIn     = true;
    case 'PhotoDember'
        params.usePhotoDember = true;
        params.useBuiltIn     = false;
    case 'Both'
        params.usePhotoDember = true;
        params.useBuiltIn     = true;
    otherwise
        error('mechanism must be ''BuiltIn'', ''PhotoDember'', or ''Both''');
end
params.theta = 0;

params.X0 = X0;
params.Y0 = Y0;
params.z0 = z0;

%% --- 2. Run simulation ---
theta_deg = atan2d(sqrt(X0^2 + Y0^2), z0);
R_mm = sqrt(X0^2 + Y0^2 + z0^2) * 1e3;
fprintf('Observation point: (%.4f, %.4f, %.5f) m\n', X0, Y0, z0);
fprintf('Angle from normal: %.1f deg,  Distance: %.2f mm\n', theta_deg, R_mm);
fprintf('Running simulation (Nx=%d, Ny=%d) — this may take a few minutes...\n', ...
    params.Nx, params.Ny);

tic;
results = simulateTHzFromDiffusion(params);
elapsed = toc;
fprintf('Simulation complete in %.1f s\n', elapsed);

%% --- 3. Extract E-field components ---
t_ps = results.t_ps;
Ex   = results.Ex_pt;
Ey   = results.Ey_pt;
Ez   = results.Ez_pt;
E_mag = results.E_mag;

[~, idx_peak] = max(E_mag);
fprintf('Peak |E| = %.3e V/m at t = %.2f ps\n', E_mag(idx_peak), t_ps(idx_peak));

%% --- 4. Plot Ex and Ey in separate figures ---
title_base = sprintf('(%.1f mm, %.1f mm, %.1f mm)', ...
    X0*1e3, Y0*1e3, z0*1e3);

% --- Figure 1: Ex(t) ---
fig1 = figure('Units', 'centimeters', 'Position', [3 3 14 9], 'Color', 'w');
ax1  = axes(fig1);
plot(ax1, t_ps, Ex, '-', 'Color', [0.85 0.25 0.2], 'LineWidth', 1.8);
xlabel(ax1, 't (ps)');
ylabel(ax1, 'E_x(t) [V/m]');
title(ax1, ['E_x at ' title_base]);
set(ax1, 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'on');

% --- Figure 2: Ey(t) ---
% fig2 = figure('Units', 'centimeters', 'Position', [3 15 14 9], 'Color', 'w');
% ax2  = axes(fig2);
% plot(ax2, t_ps, Ey, '-', 'Color', [0.2 0.6 0.3], 'LineWidth', 1.8);
% xlabel(ax2, 't (ps)');
% ylabel(ax2, 'E_y(t) [V/m]');
% title(ax2, ['E_y at ' title_base]);
% set(ax2, 'FontSize', 11, 'FontName', 'Times New Roman', ...
%     'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'on');

%% --- 5. Save ---
% exportgraphics(fig1, 'figure_Ex.png', 'Resolution', 600);
% exportgraphics(fig2, 'figure_Ey.png', 'Resolution', 600);
% fprintf('Saved: figure_Ex.png, figure_Ey.png\n');
