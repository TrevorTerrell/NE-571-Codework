function spatial_flux = SpatialFlux(flux, Nx, Ny)
% SpatialFlux   Converts the linear form of a 2D flux to a square form.
%   Varibles Nx and Ny are the *total* number of nodes in the x and y axes.
    spatial_flux = zeros(Ny, Nx);
    index = 1;
    for i = 2:(Ny - 1)
        for j = 2:(Nx - 1)
            spatial_flux(i,j) = flux(index);
            index = index + 1;
        end
    end
end