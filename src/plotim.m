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