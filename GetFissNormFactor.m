function factor = GetFissNormFactor(flux, layout, materials, group)

    [Nx, Ny] = size(layout);
    fiss = zeros(Ny,Ny);
    for i = 1:Nx
        for j = 1:Ny
            fiss(i, j) = materials(layout(i, j)).fiss(group);
        end
    end

    factor = sum(sum(fiss .* flux));
end
