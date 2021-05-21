

Mesh = load_mesh('./MeshSample/mesh5.node');

fluence_data = femdata('./MeshSample/mesh5.node',0);

%Convert sparse matrix to full matrix 
fluence_data.phi = full(fluence_data.phi);

H = fluence_data.phi.*Mesh.mua;

% subplot(1,4,1);
figure;
plotim(Mesh,Mesh.mua);
title('\mu_a','FontSize',20);
colorbar('horiz');

% subplot(1,4,2);
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



%% plot image function
function plotim(mesh,val)

       h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val);

shading interp;
view(2);
axis equal; 
axis off;
colormap hot;

end