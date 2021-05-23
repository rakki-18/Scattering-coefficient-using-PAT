
clearvars;
% load the initial pressure distribution from an image and scale the
% magnitude
GruneisenParameter = 1;
p0 = GruneisenParameter * loadImage('./Results/energy_distribution_new.jpg');

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction  [m]
dy = 0.1e-3;        % grid point spacing in the y direction  [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% resize the image to match the size of the computational grid and assign
% to the source input structure
source.p0 = resize(p0, [Nx, Ny]);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% define a centered circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, source.p0 + cart2grid(kgrid, sensor.mask), [-GruneisenParameter, GruneisenParameter]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-GruneisenParameter, GruneisenParameter]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;