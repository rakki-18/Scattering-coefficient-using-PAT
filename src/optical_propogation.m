% Runs the optical propagation on the mesh, "meshName"
% Calculates the fluence and saves the energy distribution as a mat file
function optical_propogation(meshName)

% Load the mesh
meshLoc = "./MeshSample/" + string(meshName) + "/" + string(meshName);
meshLoc = char(meshLoc);
Mesh = load_mesh(meshLoc);


% Calculate the fluence
fluence_data = femdata(meshLoc,0);

%Convert sparse matrix to full matrix 
fluence_data.phi = full(fluence_data.phi);

% Taking sum along the rows as there are more than one source
fluence_data.phi = sum(fluence_data.phi, 2)/size(fluence_data.phi,2);

% Calculating energy distribution
H = fluence_data.phi.*Mesh.mua;

% % %
% PLOTTING THE DISTRIBUTION
% % %


figure;
plotim(Mesh,Mesh.mua);
title('\mu_a','FontSize',20);
colorbar('horiz');


figure;
plotim(Mesh,Mesh.mus);
title('\mu_s','FontSize',20);
colorbar('horiz');

figure;
plotim(Mesh,fluence_data.phi);
title('Fluence','FontSize',20);
colorbar('horiz');

figure;
plotim(Mesh,H);
title('Energy distribution','FontSize',20);
colorbar('horiz');

save('variables','H','Mesh');


end
