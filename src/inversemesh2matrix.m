inverseMesh2matrix(p0_recon,Mesh.nodes);

function inverseMesh2matrix(p0,node)

% Input - p0 : a 2D matrix source where source(i,j) represents the energy distribution of the node located at the cell (i,j) 
%         node : the coordinate points of the nodes in the mesh.
% Output - a 1D vector representing the energy distribution at each node

offset = min(node(:,1));
length = sqrt(size(node,1));
node_dist = (max(node(:,1)) - offset)/(length - 1);
H_recon = zeros(size(node,1),1);

for i = 1:size(node,1)
    x = (node(i,1) - offset)/node_dist;  % finding X-coordinate
    y = (node(i,2) - offset)/node_dist;  % finding Y-coordinate
    average = 0;
    for j = x*4 + 1: (x+1)*4
        for k = y*4 + 1: (y+1)*4
            average = average + p0(j,k);
        end
    end
    average = average/16;
    average = average/100;
    H_recon(i) = average;
end
save('variables','H_recon','-append');
end