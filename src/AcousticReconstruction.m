
clearvars;
% create the computational grid
PML_size = 20;              % size of the PML in grid points
Nx = 128 - 2 * PML_size;    % number of grid points in the x direction
Ny = 256 - 2 * PML_size;    % number of grid points in the y direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% load the initial pressure distribution from an image and scale the
% magnitude
GruneisenParameter = 1;
p0 = GruneisenParameter * loadImage('./Results/energy_distribution_H.jpg');

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(p0, true);
% resize the image to match the size of the computational grid and assign
% to the source input structure
p0 = resize(p0, [Nx, Ny]);

% assign to the source structure
source.p0 = p0;


% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add first order compensation for only recording over a half plane
p0_recon = 2 * p0_recon;


% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + sensor.mask * GruneisenParameter, [-GruneisenParameter, GruneisenParameter]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% plot the reconstructed initial pressure 
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon, [-GruneisenParameter, GruneisenParameter]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% apply a positivity condition
p0_recon(p0_recon < 0) = 0;

% plot the reconstructed initial pressure with positivity condition
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon, [-GruneisenParameter, GruneisenParameter]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% plot a profile for comparison
figure;
profile = 30;
plot(kgrid.y_vec * 1e3, p0(profile, :), 'k-', ...
     kgrid.y_vec * 1e3, p0_recon(profile, :), 'b:');
xlabel('y-position [mm]');
ylabel('Pressure');
legend('Initial Pressure', 'Time Reversal');
axis tight;
set(gca, 'YLim', [0, 5.1]);