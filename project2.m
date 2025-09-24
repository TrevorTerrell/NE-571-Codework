clear

%% Inputs
% fission_x_sec =         [ 0.003320, 0.07537 ];
% neutron_per_fission =   [ 0.008476 / fission_x_sec(1), 0.18514 / fission_x_sec(2)];
% absorption_x_sec =      [ 0.0004,   0.0197  ];
% diffusion_coefficient =  [ 1.13,     0.16    ];
% scatter_x_sec =         [ 0.0494 ];
materials(1).fiss = [ 0.003320; 0.07537 ];
materials(1).nu   = [ 0.008476; 0.18514 ] ./ materials(1).fiss;
materials(1).abs  = [ 0.012070; 0.12100 ];
materials(1).diff = [ 1.262700; 0.35430 ];
materials(1).sct  = [ 0.01412 ];
materials(1).rm   = [ 0.026190; 0.12100 ];

materials(2).fiss = [ 0;      0      ];
materials(2).nu   = [ 0;      0      ] ./ materials(1).fiss;
materials(2).abs  = [ 0.0004; 0.0197 ];
materials(2).diff = [ 1.13;   0.16   ];
materials(2).sct  = [ 0.0494 ];
materials(2).rm   = [ 0.0498; 0.0197 ];

materials4(1).fiss = [ 0.003378; 0.0004850; 0.006970; 0.07527 ];
materials4(1).nu   = [ 0.009572; 0.001193;  0.01768;  0.18514 ] ./ materials4(1).fiss;
materials4(1).abs  = [ 0.004946; 0.002840;  0.03053;  0.1210  ];
materials4(1).diff = [ 2.1623;   1.0867;    0.6318;   0.3543  ];
materials4(1).sct  = [ 0.083004; 0.0584;    0.06453 ];
materials4(1).rm   = [ 0.08795;  0.06124;   0.09506;  0.1210  ];

materials4(2).fiss = [ 0; 0; 0; 0 ];
materials4(2).nu   = [ 0; 0; 0; 0 ];
materials4(2).abs  = [ 3.1e-4; 2.475e-3; 1.45e-3; 9.114e-3 ];
materials4(2).sct  = [ 3.058e-2; 7.42e-2; 1.018e-1 ];
materials4(2).rm   = [ 3.096e-2; 7.707e-2; 1.033e-1; 9.351e-3];
materials4(2).diff = 1 ./ (3 .* ([1.805e-1; 3.327e-1; 3.611e-1; 1.047] + materials4(2).rm));


mats = materials;

%width_in_x = 10;
%width_in_y = 10;

% layout = [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2;
%            2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2;
%            2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2;
%            2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2;
%            2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2;
%            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2;];
% fprintf('Using custom core')

% Quick Square-in-Square
outer = 51;
inner = 51;
inner_width = 83.64;
layout = 2 * ones(outer);
layout((1 + (outer - inner) / 2):(inner + (outer - inner) / 2), (1 + (outer - inner) / 2):(inner + (outer - inner) / 2)) = ones(inner);

width_in_x = inner_width * outer / inner;
width_in_y = width_in_x;

fprintf('Using %.2f by %.2f cm core surrounded by water\n', inner_width, inner_width)

POWER = 3000e6 / 3.2e-11; %Fissions

%% Constants
groups = size(mats(1).fiss, 1);

[nodes_in_x, nodes_in_y] = size(layout);

node_width_x = width_in_x / nodes_in_x;
node_width_y = width_in_y / nodes_in_y;

%% Building Matrices
disp('Building Matrices')
tic;
SPAN = (nodes_in_x - 2) * (nodes_in_y - 2);

M = zeros(groups * SPAN);
S = zeros(groups * SPAN);

flux = ones(size(M,1),1);
flux_old = zeros(size(M,1),1);
crit = 1;
crit_old = 0;
tol = 1e-5;

A_g = zeros(SPAN);
F_g = zeros(SPAN);
S_g = zeros(SPAN);

for g = 1:groups
    mat_area = (1 + SPAN * (g - 1)):(SPAN * g);
    %A_g = MatCreateA(center_const(:, g + 1), x_axis_const(:, g + 1), x_axis_const(:, g + 1), y_axis_const(:, g + 1), y_axis_const(:, g + 1), nodes_in_x, nodes_in_y);
    A_g = CreateLossMat(layout, mats, node_width_x, node_width_y, g);
    F_g = CreateFissMat(layout, mats, g);
    M(mat_area, mat_area) = A_g;
    S(1:SPAN, mat_area) = F_g;

    if g ~= groups
        S_g = CreateSctrMat(layout, mats, g);
        S(SPAN + mat_area, mat_area) = S_g;
    end
end

clear A_g
clear S_g
clear F_g

%M = sparse(M);
%M_inv = inv(M);
S = sparse(S);
matrix = sparse(inv(M) * S);

K = ones(groups * SPAN, 1);
K(1:(SPAN)) = 1/crit;

fprintf('\tFinsihed in %.2f s\n', toc)
clear M

%% Main Loop
fprintf('Entering Main Loop\n')
tic;
n = 0;
while abs(crit - crit_old) >= tol
    n = n + 1;
    flux_old = flux;
    crit_old = crit;
    flux = (matrix .* K) * flux_old;
    crit = crit * sum(S * flux .* K) / sum(S * flux_old .* K);
    K(1:(SPAN)) = 1/crit;
    flux =  flux ./ sum(flux);
end
fprintf('\tMain Loop converged after %i iterations taking %.2f s\n',n, toc)

%% Flux Cleanup and Plotting
x_vals = linspace(0, width_in_x, nodes_in_x);
nodes_per_group = SPAN;

groupData(groups).flux = zeros(nodes_in_x, nodes_in_y);
groupData(groups).fiss_factor = 0;
unorm_power = 0;

figure(100);

for g = 1:groups
    st = 1 + (nodes_per_group * (g - 1));
    sp = nodes_per_group * g;
    g_flux = flux(st:sp);
    g_spacial_flux = SpatialFlux(g_flux, nodes_in_x, nodes_in_y);

    groupData(g).flux = g_spacial_flux .* node_width_x .* node_width_y;
    groupData(g).fiss_factor = GetFissNormFactor(groupData(g).flux, layout, mats, g);
    unorm_power = unorm_power + groupData(g).fiss_factor;
end

for g = 1:groups
    groupData(g).flux = groupData(g).flux .* POWER ./ unorm_power ./ node_width_y ./ node_width_x;
    center_line = groupData(g).flux(1:end, (nodes_in_y + 1) / 2);
    plot(x_vals, center_line, 'DisplayName', sprintf('Group %i', g - 1), 'LineWidth', 3)
    hold("on");
end

xlabel("Y-Centerline X-Position (cm)")
ylabel("Flux (neutrons/cm^2)")
legend("show");
fprintf('Solution Found for 2D slab with %.3f x %.3f cm, %i x %i nodes, and k = %.5f', width_in_x, width_in_y, nodes_in_x, nodes_in_y, crit)