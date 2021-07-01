meshName = 'crosstalk2';

% Erasing the contents of the log file and variables workspace
fid = fopen(fullfile('', 'conditionMatrix.log'), 'w');
save('variables');

% run the optical propogation on the mesh
optical_propogation(meshName);

% Load variables from the workspace
load('variables','H','Mesh','Nodes');

% Final reconstructed H value
H_recon = [];
% Iterate over each wavelength and do acoustic reconstruction
for i = 1: size(Mesh.wv)
    mesh2matrix(H((i-1)*Nodes+1:i*Nodes),Mesh.nodes);

    % energy distribution matrix is now saved as "image" in the
    % workspace
    load('variables.mat','image');
    % Simulate the Acoustic Reconstruction on the initial pressure source
    AcousticReconstruction(image);

    % Load the reconstructed pressure distribution value
    load('variables','p0_recon');
    % Get the reconstructed H vector from the 2D reconstructed pressure source
    inversemesh2matrix(p0_recon,Mesh.nodes);

    % Load the reconstructed H for the particular wavelength
    load('variables.mat','H_partial');
%     % Plot the reconstructed energy distribution
    figure;
    plotim(Mesh,H_partial);
    title('Reconstructed Energy distribution','FontSize',10);
    colorbar('horiz');
    H_recon = [H_recon; H_partial];
end




% Provide the initial homogeneous guess for mua and mus values
% This would be calculated using the known concentrations of chromophores
% of the background and the known scattering power and scattering amplitude
% of the background.
% [mua, mus] = calc_initial_value([0.01, 0.01,0.4],[1,1],Mesh);
% save('variables.mat','mua','mus','-append');
% 
% % Reconstruct the values of mua and mus from the reconstructed energy
% % distribution vector
% mua_mus_reconstruction(H_recon,mua,mus,Mesh);