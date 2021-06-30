% Reconstructs the value of mua given the value of H , initial mua and mus values
function mua = mua_reconstruction(H,mua,mus,Mesh)

%hyper parameters
max_iterations = 3;
regularisation_parameter = 0.1;


nodes = size(H,1)/size(Mesh.wv,1);

kappa = find_kappa(mua,mus);
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");

% Storing error in each iterations
error_list = [;];
error_list_mua = [];
% iterative update of mua values
for i = 1:max_iterations
    fprintf("%d iterations started\n",i);
    % Calculate fluence at each wavelength for the current mesh
    fluence_new = new_femdata_spectral(mesh_new,0);


    G = find_jacobian(mua,mus,kappa,mesh_new,fluence_new,i);

    delta_t = find_delta_t(G,regularisation_parameter,H, fluence_new.phi.*mua);

    [mua,mus,kappa,mesh_new] = update(delta_t, mua,mus,kappa,Mesh);

    % Find error in calculation
    error_H = sum((fluence_new.phi.*mua - H).*(fluence_new.phi.*mua - H),1);
    error_list = [error_list [error_H;sum(abs(Mesh.mua - mua),1)]];
    error_mua = abs(Mesh.mua - mua);
    error_list_mua = [error_list_mua error_mua];

    save('variables','error_list','error_list_mua','-append');

    %% PLOTTING RESULTS for first wavelength
    figure;
    plotim(Mesh,mua(1:nodes));
    title('mua obtained','FontSize',20);
    colorbar('horiz');
      
    
end





   
end

% Creates a new mesh, "mesh_name" having the mua, mus and kappa values
function [mesh_New] = new_mesh(Mesh,mua,mus,kappa,~)
mesh_New = Mesh;
mesh_New.mua = mua;
mesh_New.mus = mus;
mesh_New.kappa = kappa;
% mesh_loc = "./Meshes/" + mesh_name + "/" + mesh_name;
% mesh_loc = char(mesh_loc);
% save_mesh(mesh_new,mesh_loc);
end

% finding the update of optical parameters from the jacobian matrix,
% calculated energy distribution and the actual energy distribution
function [delta_t] = find_delta_t(G, regularisation_parameter,H, Hcal)

% (GtG + lambda*I)delta_t = Gt*(Em - Ec)
hessian = transpose(G)*G;
delta_t = hessian + regularisation_parameter*eye(size(hessian,1));
clear hessian;
delta_t = delta_t\transpose(G);
clear G;
delta_t = delta_t*(H - Hcal);

end

% helper function to find kappa from mua and mus values
function [kappa] = find_kappa(mua,mus)
kappa = 1./(3*(mus+mua));
end

% updates the mua, mus and kappa values based on delta_t
function [mua,mus,kappa,mesh_new] = update(delta_t, mua, mus, kappa,Mesh)
learning_rate_mua = 1;
mua = mua + learning_rate_mua*delta_t;
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");
end

% finding the jacobian matrix which is the derivative of energy distributin
% wrt mua values
function [G] = find_jacobian(mua,mus,kappa,mesh_new,fluence,iter)

G = zeros(size(mesh_new.mua,1), size(mesh_new.mua,1));
delta = 0.0001;

for i = 1: size(mua,1)
% for i = 1: 10 
    fprintf("%d iteration in Jacobian mua started. %d th iteration\n",i,iter);
    temp_mua = mua;
    temp_mua(i) = mua(i) + delta;
    mesh_jacobian = new_mesh(mesh_new,temp_mua,mus,kappa,"mesh_jacobian");
    fluence_data = new_femdata_spectral(mesh_jacobian,0);
    

    G(:,i) = (fluence_data.phi - fluence.phi)/delta;
    G(:,i) = G(:,i).*mesh_new.mua;
    G(i,i) = G(i,i) + fluence.phi(i);
    clear mesh_jacobian temp_mua fluence_data;
end


% Log the condition number of Jacobian matrix 

fid = fopen(fullfile('', 'conditionMatrix.log'), 'a');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid, 'Condition number of Jacobian matrix at %d iteration: %d\n', iter, cond(G));
fclose(fid);

end
