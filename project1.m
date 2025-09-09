%% Setup
clear
clc

%% Constants
% macroscopic cross sections
neutron_production = 0.1570;    %cm^-1
neutron_transfer = 3.62e-2;     %cm^-1
neutron_absorption = 0.1532;    %cm^-1
power = 1;                      %fractional
energy_per_fission_neutron = 1; %fractional

% spacial data
nodes_x = 256;
nodes_y = 256;
width_x = 512;
width_y = 512;

node_width_x = width_x/(nodes_x - 1);
node_width_y = width_y/(nodes_y - 1);

nodes_x = nodes_x - 2;
nodes_y = nodes_y - 2;

%% Prebuild A matrix
% cell values
diffusion_coefficent = (3 * (neutron_absorption + neutron_transfer)) ^ -1;

center = (2 * diffusion_coefficent) / (node_width_x ^ 2) + (2 * diffusion_coefficent) / (node_width_y ^ 2) + neutron_absorption;
north = - diffusion_coefficent / (node_width_y ^ 2);
south = - diffusion_coefficent / (node_width_y ^ 2);
east = - diffusion_coefficent / (node_width_x ^ 2);
west = - diffusion_coefficent / (node_width_x ^ 2);

% cells to vectors
A_only_x = zeros(nodes_x);
    % only_x is comprised of c, e, and w in a tridiag array
A_only_c = ones(nodes_x, 1);
A_only_e = ones(nodes_x - 1, 1); % -1 [by fencepost]
A_only_w = ones(nodes_x - 1, 1); % ''
    % no "A_only_y" b/c y's contributions are discrete
A_only_n = ones((nodes_y - 1) * nodes_x, 1); % -1 for same as above, * nodes_x because outside
A_only_s = ones((nodes_y - 1) * nodes_x, 1); % ''

A_only_c = center.*A_only_c;
A_only_e = east.*A_only_e;
A_only_w = west.*A_only_w;
A_only_n = north.*A_only_n;
A_only_s = south.*A_only_s;

%% Build A matrix
A = zeros(nodes_x * nodes_y);
A_only_x = A_only_x + diag(A_only_c) + diag(A_only_e, 1) + diag(A_only_w, -1);
for i = 0:(nodes_y - 1)
    ij = i * nodes_x + 1;
    A(ij:(ij + nodes_x - 1), ij:(ij + nodes_x - 1)) = A_only_x;
end
A = A + diag(A_only_n, nodes_x) + diag(A_only_s, - nodes_x);

%% Build B matrix
B = eye(nodes_x * nodes_y) * neutron_production;

%% Final Loop Prep
C = inv(A) * B;
flux = ones(nodes_x * nodes_y, 1);
flux_old = 0;
criticality = 1;
tolerance = 1e-10;

%% Iteration Loop
while abs(min(flux) - min(flux_old)) > tolerance
    flux_old = flux;
    flux = C * flux_old / criticality;
    criticality = criticality * sum(flux) / sum(flux_old);
    flux = flux / sum(flux);
end

%% Flux Cleanup
flux = flux .* (width_x * width_y) ./ sum(flux);
flux = flux ./ max(flux);
spatial_flux = zeros(nodes_y + 2, nodes_x + 2);
index = 1;
for i = 2:(1 + nodes_y)
    for j = 2:(1 + nodes_x)
        spatial_flux(i,j) = flux(index);
        index = index + 1;
    end
end

%% Plotting
data = spatial_flux;
xvals = 0:node_width_x:width_x;
yvals = 0:node_width_y:width_y;

xvals(mod(0:nodes_x, 8) ~= 0) = "";
yvals(mod(0:nodes_y, 8) ~= 0) = "";

figure(1);
h = heatmap(data);
h.Colormap = hot;
h.XDisplayLabels = xvals;
h.YDisplayLabels = yvals;
h.GridVisible = 'off';
h.Title = 'Flux Heatmap (Numerical)';
disp(criticality)
%% Analytical Solution
anal_flux = zeros(nodes_y + 2, nodes_x + 2);
for j = 0:(nodes_y + 1)
    for i = 0:(nodes_x + 1)
        anal_flux(j + 1, i + 1) = cos(pi * i * node_width_x / width_x - pi/2) * cos(pi * j * node_width_y / width_y - pi/2);
    end
end

%anal_flux = anal_flux .* (width_x * width_y) ./ sum(sum(anal_flux));
%anal_flux = anal_flux ./ max(max(anal_flux));

data = anal_flux;
figure(2);
h = heatmap(data);
h.Colormap = hot;
h.XDisplayLabels = xvals;
h.YDisplayLabels = yvals;
h.GridVisible = 'off';
h.Title = 'Flux Heatmap (Analytical)';