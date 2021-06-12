function mesh2matrix(H, node)

% Input - 1D vector reresenting the energy distribution at each node
%         mesh.node - the coordinate points of the nodes in the mesh
% Output - a 2D matrix source where source(i,j) represents the energy distribution of the node located at the cell (i,j).
% Note that this is for a rectangular mesh 

offset = min(node(:,1));
length = sqrt(size(H,1));
node_dist = (max(node(:,1)) - offset)/(length - 1);
image = zeros(4*length,4*length);

for i = 1:size(H,1)
    x = (node(i,1) - offset)/node_dist; % finding X-coordinate
    y = (node(i,2) - offset)/node_dist; % finding Y-coordinate
    for j = x*4+1 : (x+1)*4
        for k = y*4+1 : (y+1)*4
            image(j,k) = 100*H(i);     %scaling up the values by a factor of 100.
        end
    end 
    
end
save('variables','image','-append');
end
