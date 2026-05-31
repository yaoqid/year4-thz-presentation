%% plot_5angles_Ex.m
% Plot Ex(t) at 5 observation angles in one figure,
% similar to the style of Fig. 8(a).
% Also prints distance-corrected Epp values for Ex, Ez, and E_det.

clear; clc; close all;

%% --- 0. Choose 5 angles ---
theta_list = [18, 22, 25, 30, 50, 80,];   % change these as you want
X0 = 0.0032;   % fixed lateral x position (m)
Y0 = 0;

mechanism = 'PhotoDember';   % 'BuiltIn', 'PhotoDember', or 'Both'

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

%% --- 2. Run simulation for each angle ---
nAngles = numel(theta_list);
Ez_all = cell(1, nAngles);
t_ps_all = cell(1, nAngles);

% Store Epp results
Epp_Ex_list   = zeros(1, nAngles);
Epp_Ez_list   = zeros(1, nAngles);
Epp_Edet_list = zeros(1, nAngles);
R_list        = zeros(1, nAngles);
z0_list       = zeros(1, nAngles);

for i = 1:nAngles
    theta_deg = theta_list(i);

    if theta_deg == 90
        z0 = 0;
    else
        z0 = X0 / tand(theta_deg);
    end

    params.z0 = z0;

    R = sqrt(X0^2 + Y0^2 + z0^2);
    R_mm = R * 1e3;

    R_list(i) = R;
    z0_list(i) = z0;

    fprintf('Running angle = %.1f deg, z0 = %.4f mm, R = %.4f mm\n', ...
        theta_deg, z0*1e3, R_mm);

    results = simulateTHzFromDiffusion(params);

    t_ps = results.t_ps;
    Ex   = results.Ex_pt;
    Ez   = results.Ez_pt;

    % Rotating transverse detected field in the x-z plane
    if R == 0
        E_det = zeros(size(Ex));
    else
        ux = z0 / R;
        uz = -X0 / R;
        E_det = Ex * ux + Ez * uz;
    end

    % Store waveform for plotting
    t_ps_all{i} = t_ps;
    Ez_all{i}   = Ez;   % use fixed Ex for stacked waveform plot

    % Distance-corrected peak-to-peak values
    Epp_Ex   = (max(Ex)    - min(Ex))    * R;
    Epp_Ez   = (max(Ez)    - min(Ez))    * R;
    Epp_Edet = (max(E_det) - min(E_det)) * R;

    Epp_Ex_list(i)   = Epp_Ex;
    Epp_Ez_list(i)   = Epp_Ez;
    Epp_Edet_list(i) = Epp_Edet;

    fprintf('\n--- Current observation point ---\n');
    fprintf('theta = %.1f deg\n', theta_deg);
    fprintf('X0 = %.4f mm, Y0 = %.4f mm, Z0 = %.4f mm\n', X0*1e3, Y0*1e3, z0*1e3);
    fprintf('R  = %.4f mm\n', R*1e3);
    fprintf('Epp_Ex   * R = %.4f\n', Epp_Ex);
    fprintf('Epp_Ez   * R = %.4f\n', Epp_Ez);
    fprintf('Epp_Edet * R = %.4f\n\n', Epp_Edet);
end

%% --- 3. Print summary table ---
fprintf('\n================ Summary Table ================\n');
fprintf('theta(deg)    z0(mm)      R(mm)      Epp_Ex*R      Epp_Ez*R      Epp_Edet*R\n');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:nAngles
    fprintf('%8.1f   %8.4f   %8.4f   %12.4f   %12.4f   %13.4f\n', ...
        theta_list(i), z0_list(i)*1e3, R_list(i)*1e3, ...
        Epp_Ex_list(i), Epp_Ez_list(i), Epp_Edet_list(i));
end

fprintf('================================================\n');

%% --- 4. Find a good vertical offset ---
maxAmp = 0;
for i = 1:nAngles
    maxAmp = max(maxAmp, max(abs(Ez_all{i})));
end
offset_step = 1.4 * maxAmp;

%% --- 5. Plot all waveforms in one figure ---
fig = figure('Units', 'centimeters', 'Position', [3 3 16 12], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');

for i = 1:nAngles
    t_ps = t_ps_all{i};
    Ez   = Ez_all{i};

    offset = (nAngles - i) * offset_step;
    
    plot(ax, t_ps, Ez + offset, 'LineWidth', 1.6);

    % label angle near the right side
    text(ax, t_ps(end) + 0.8, offset, sprintf('%g^\\circ', theta_list(i)), ...
        'FontSize', 10, 'VerticalAlignment', 'middle');
end

xlabel(ax, 't (ps)');
ylabel(ax, 'E_z(t) [V/m]');
title(ax, 'Far-field THz waveforms at different observation angles');
set(ax, 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'LineWidth', 1.0, 'TickDir', 'out', 'Box', 'on');

grid(ax, 'on');
xlim(ax, [min(t_ps_all{1}), max(t_ps_all{1}) + 3]);

% optional: hide y tick labels to make it look more like a paper figure
ax.YTickLabel = {};

%% --- 6. Save figure ---
% exportgraphics(fig, 'figure_5angles_Ex.png', 'Resolution', 600);
% fprintf('Saved: figure_5angles_Ex.png\n');