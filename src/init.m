

% create the computational grid
PML_size = 20;          % size of the PML in grid points
Nx = 60 - 2*PML_size;  % number of grid points in the x (row) direction
Ny = 60 - 2*PML_size;  % number of grid points in the y (column) direction
dx = 0.1e-3;            % grid point spacing in the x direction [m]
dy = 0.1e-3;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

      
% define the properties of the propagation medium
medium.sound_speed = 1500;           % [m/s]

detectors = double.empty(0,2);
for i = 1 : size(kgrid.x)
    for j = 1 :  size(kgrid.x)
        temp = {kgrid.x(i,j), kgrid.y(i,j)};
        detectors = [detectors;temp];
    end
end

save('variables');


