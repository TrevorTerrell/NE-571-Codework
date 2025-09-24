function S = CreateSctrMat(layout, materials, g)
% CreateLossMat Creates the F matrix for the A * phi = F * phi / k equation.
%
%   layout: 2D matrix holding the material ID for each node
%   materials: struct with the following components:
%       sct: downscattering cross section
%   g: number of energy groups

    % Shrink layout b/c 0 boundary condition
    layout = layout(2:(end - 1), 2:(end - 1));
    [Nx, Ny] = size(layout);

    Sc = zeros(Nx*Ny, 1);

    for i = 1:Nx
        for j = 1:Ny
            Sc((i - 1) * Ny + j) = materials(layout(i, j)).sct(g);
        end
    end

    S = diag(Sc);
end