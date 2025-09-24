function A = CreateLossMat(layout, materials, dx, dy, g)
% CreateLossMat Creates the A matrix for the A * phi = F * phi / k equation.
%
%   layout: 2D matrix holding the material ID for each node
%   materials: struct with the following components:
%       diff: diffusion coefficient
%       rm: removal cross section
%   dx: width of 1 node in the x direction
%   dx: width of 1 node in the y direction
%   g: number of energy groups

    %Shrink layout b/c 0 boundary condition
    layout = layout(2:(end - 1), 2:(end - 1));
    [Nx, Ny] = size(layout);

    A_only_c = zeros(Nx * Ny, 1);
    A_only_e = zeros(Nx * Ny - 1, 1);
    A_only_w = zeros(Nx * Ny - 1, 1);
    A_only_n = zeros(Nx * (Ny - 1), 1);
    A_only_s = zeros(Nx * (Ny - 1), 1);

    for i = 1:Nx
        for j = 1:Ny
            A_only_c((i - 1) * Ny + j) = 2 * materials(layout(i,j)).diff(g) * (dx ^ 2 + dy ^ 2) / (dx ^ 2 * dy ^ 2) + materials(layout(i,j)).rm(g);
        end
    end

    for i = 1:Nx
        for j = 1:(Ny - 1)
            A_only_e((i - 1) * Ny + j) = - materials(layout(i, j + 1)).diff(g) / (dx ^ 2);
            A_only_w((i - 1) * Ny + j) = - materials(layout(i, j    )).diff(g) / (dx ^ 2);
        end
    end

    for j = 1:Ny
        for i = 1:(Nx - 1)
            A_only_n((i - 1) * Ny + j) = - materials(layout(i,     j)).diff(g) / (dy ^ 2);
            A_only_s((i - 1) * Ny + j) = - materials(layout(i + 1, j)).diff(g) / (dy ^ 2);
        end
    end

    % c = diag(A_only_c) 
    % e = diag(A_only_e, 1) 
    % w = diag(A_only_w, -1) 
    % n = diag(A_only_n, Nx)
    % s = diag(A_only_s, - Nx)

    A = diag(A_only_c) + diag(A_only_e, 1) + diag(A_only_w, -1) + diag(A_only_n, Nx) + diag(A_only_s, - Nx);

end
