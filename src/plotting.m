figure;
plotim(Mesh,mua);
title('Fluence','FontSize',20);
colorbar('horiz');

figure;
plotim(Mesh,Mesh.mua);
title('Fluence','FontSize',20);
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