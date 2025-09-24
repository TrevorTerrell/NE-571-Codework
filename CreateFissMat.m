function F = CreateFissMat(layout, materials, g)
% CreateLossMat Creates the F matrix for the A * phi = F * phi / k equation.
%
%   layout: 2D matrix holding the material ID for each node
%   materials: struct with the following components:
%       fiss: fission cross section
%       nu: number of neutrons per fission
%   g: number of energy groups

    % Shrink layout b/c 0 boundary condition
    layout = layout(2:(end - 1), 2:(end - 1));
    [Nx, Ny] = size(layout);

    nuSigF = zeros(Nx*Ny, 1);

    for i = 1:Nx
        for j = 1:Ny
            nuSigF((i - 1) * Ny + j) = materials(layout(i, j)).fiss(g) * materials(layout(i, j)).nu(g);
        end
    end

    F = diag(nuSigF);%eye(Nx * Ny) .* nuSigF;
end