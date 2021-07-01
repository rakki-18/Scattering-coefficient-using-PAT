% Runs the optical propagation on the spectral mesh, "meshName"
% Calculates the fluence and saves the energy distribution as a mat file
function optical_propogation(meshName)

% Load the mesh
meshLoc = "./Meshes/" + string(meshName) + "/" + string(meshName);
meshLoc = char(meshLoc);
Mesh = load_mesh(meshLoc);

Nodes = size(Mesh.nodes,1);

% Calculate the fluence
[fluence_data,mua,mus] = femdata_spectral(Mesh,0,Mesh.wv);
% Adding the mua and mus information to the mesh
Mesh.mua = []; Mesh.mua = mua;
Mesh.mus = []; Mesh.mus = mus;

% Calculating energy distribution
H = fluence_data.phi.*Mesh.mua;


% Iterate through all the wavelengths
% for i = 1: size(Mesh.wv)
for i = 1: 1
    
    % % %
    % PLOTTING THE DISTRIBUTION
    % % %


    figure;
    plotim(Mesh,Mesh.mua((i-1)*Nodes+1:i*Nodes));
    title('\mu_a','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,Mesh.conc(:,1));
    title('Actual HbO','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,Mesh.conc(:,2));
    title('Actual deoxyHbO','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,fluence_data.phi((i-1)*Nodes+1:i*Nodes));
    title('Fluence','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,H((i-1)*Nodes+1:i*Nodes));
    title('Energy distribution','FontSize',20);
    colorbar('horiz');
end

save('variables','H','Mesh','Nodes','-append');
    

end
