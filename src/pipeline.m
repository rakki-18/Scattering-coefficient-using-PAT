meshName = 'meshSpectral';

% Erasing the contents of the log file
fid = fopen(fullfile('', 'conditionMatrix.log'), 'w');

% run the optical propogation on the mesh
optical_propogation(meshName);

% Load variables from the workspace
load('variables','H','Mesh');

if(Mesh.type == 'spec')
    % Only one wavelength is considered as of now
    H = H(:,1);
end
% Convert H vector to a 2D matrix
mesh2matrix(H,Mesh.nodes);

% energy distribution matrix is now saved as "image" in the
% workspace
load('variables.mat','image');
% Simulate the Acoustic Reconstruction on the initial pressure source
AcousticReconstruction(image);

% Load the reconstructed pressure distribution value
load('variables','p0_recon');
% Get the reconstructed H vector from the 2D reconstructed pressure source
inversemesh2matrix(p0_recon,Mesh.nodes);

load('variables','H_recon');
%Plot the reconstructed energy distribution
figure;
plotim(Mesh,H_recon);
title('Reconstructed Energy distribution','FontSize',10);
colorbar('horiz');

% Provide the initial homogeneous guess for mua and mus values
% This will be the known mua and mus value for the background.
initial_value = [0.01,1.00];


% Reconstruct the values of mua and mus from the reconstructed energy
% distribution vector
mua_mus_reconstruction(H_recon,initial_value,Mesh);