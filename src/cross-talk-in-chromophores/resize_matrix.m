function resized = resize_matrix(matrix,nodes,wv)
% Reshapes the mua matrix ( nodes*wv X 1) to (wv X nodes)

resized = zeros(wv,nodes);

for i = 1: wv
    for j = 1 : nodes
        resized(i,j) = matrix((i-1)*nodes+j,1);
    end
end
