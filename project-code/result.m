theta_deg = [18, 20, 22, 25, 28, 30, 35, 40, 45, 50, ...
             55, 60, 65, 70, 75, 80, 85, 90];

Epp = [0.9596, 0.9704, 0.9845, 1.0000, 0.9908, ...
            0.9734, 0.9293, 0.8673, 0.7927, 0.7207, ...
            0.6426, 0.5576, 0.4692, 0.3783, 0.2852, ...
            0.1915, 0.0957, 0.0001];

% Epp = Epp_raw / max(Epp_raw);
theta_rad = deg2rad(theta_deg);



figure;
pax = polaraxes;
hold(pax,'on');

polarplot(pax, theta_rad, Epp, '-o', 'LineWidth', 2);

thetalim(pax,[0 90]);
rlim(pax,[0 1]);

pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';

% show angle labels
pax.ThetaTick = [0 30 60 90];

% hide default radial labels
pax.RTick = [0 0.5 1];
pax.RTickLabel = {};

% manual y-style labels on left
text(-0.03, 0.98, '1.0', ...
    'Units','normalized', 'HorizontalAlignment','right');
text(-0.03, 0.50, '0.5', ...
    'Units','normalized', 'HorizontalAlignment','right');
text(-0.03, 0.02, '0', ...
    'Units','normalized', 'HorizontalAlignment','right');

% y-axis title
text(-0.15, 0.5, 'Amplitude, E_{pp} (norm.)', ...
    'Units','normalized', ...
    'Rotation', 90, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize', 12);

title('(b) \theta_1 (deg)');