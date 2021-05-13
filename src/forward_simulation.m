% Creating the 2D grid for computation
Nx = 128;  % Number of points in the x-direction
Ny = 256;  % Number of points in the y-direction
%Spacing between the points
dx = 50e-6;
dy = 50e-6;
%Make the grid
kgrid = makeGrid(Nx,dx,Ny,dy); 

%Properties of the medium in which the waves are propogating
medium.sound_speed= 1500*ones(Nx,Ny);  % Speed varies spatially
medium.sound_speed(1:50, :) = 1800;
medium.density = 1040;

% Properties of the source
% source is a disc in this simulation
disc_x_pos = 75;
disc_y_pos = 120;
disc_radius = 8;
disc_mag = 3;
source.p0 = disc_mag*makeDisc(Nx,Ny, disc_x_pos, disc_y_pos, disc_radius);

% Specifying the number of sensors, their size, etc.
sensor_radius = 2.5e-3;
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% Simulating the PAT experiment
sensor_data = kspaceFirstOrder2D(kgrid,medium,source,sensor);

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
