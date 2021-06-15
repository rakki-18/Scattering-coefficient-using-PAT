% Reconstructs the value of mua and mus given the value of H , initial mua and mus values
function mua_mus_reconstruction(H,mua,mus,Mesh)

%hyper parameters
max_iterations = 4;
regularisation_parameter = 0.01;

nodes = size(H,1)/size(Mesh.wv,1);

kappa = find_kappa(mua,mus);
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");

% Storing error in each iterations
error_list = [;;];
error_list_mua = [];
error_list_mus = [];
% iterative update of mua and mus values
for i = 1:max_iterations
    fprintf("%d iterations started\n",i);
    fluence_new = new_femdata_spectral(mesh_new,0);
    
    G = find_jacobian(mua,mus,kappa,mesh_new,fluence_new,i);
    
    delta_t = find_delta_t(G,regularisation_parameter,H, fluence_new.phi.*mua);
    
    
    % Find error in calculation
    error_H = sum((fluence_new.phi.*mua - H).*(fluence_new.phi.*mua - H),1);
    
    error_list = [error_list [error_H;sum(abs(Mesh.mua - mua),1);sum(abs(Mesh.mus - mus),1)]];
    error_mua = abs(Mesh.mua - mua);
    error_mus = abs(Mesh.mus - mus);
    error_list_mua = [error_list_mua error_mua];
    error_list_mus = [error_list_mus error_mus];

    [mua,mus,kappa,mesh_new] = update(delta_t, mua,mus,kappa,Mesh);
    save('variables','-append');
    %% PLOTTING RESULTS for first wavelength
    figure;
    plotim(Mesh,mua(1:nodes));
    title('mua obtained','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,mus(1:nodes));
    title('mus obtained','FontSize',20);
    colorbar('horiz');
    save('variables','-append');

    
    
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
delta_t = zeros(2*size(H,1),1);
hessian = transpose(G)*G;
delta_t = hessian + regularisation_parameter*eye(size(hessian,1));
delta_t = delta_t\transpose(G);
delta_t = delta_t*(H - Hcal);
end

% helper function to find kappa from mua and mus values
function [kappa] = find_kappa(mua,mus)
kappa = 1./(3*(mus+mua));
end

% updates the mua, mus and kappa values based on delta_t
function [mua,mus,kappa,mesh_new] = update(delta_t, mua, mus, kappa,Mesh)
mua = mua + delta_t(end/2+1:end);
kappa = kappa + delta_t(1:end/2);
mus = (1./(3.*kappa))-mua;
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");
end

% finding the jacobian matrix which is the derivative of energy distributin
% wrt kappa and mua values
function [G] = find_jacobian(mua,mus,kappa,mesh_new,fluence,iter)

G = zeros(size(mesh_new.mua,1), 2*size(mesh_new.mua,1));
delta = 0.0001;

for i = 1: size(kappa,1)
% for i = 1: 10
    fprintf("%d iteration in Jacobian started. %d th iteration\n",i,iter);
    temp_kappa = kappa;
    temp_kappa(i) = kappa(i) + delta;
    temp_mus = (1./(3.*temp_kappa))-mua;
    mesh_jacobian = new_mesh(mesh_new,mua,temp_mus, temp_kappa,"mesh_jacobian");
    
    
    fluence_data = new_femdata_spectral(mesh_jacobian,0);

    
    G(:,i) = (fluence_data.phi - fluence.phi)/delta;
    G(:,i) = G(:,i).*mua;

end

for i = 1: size(mua,1)
% for i = 1: 10 
    fprintf("%d iteration in Jacobian mua started. %d th iteration\n",i,iter);
    temp_mua = mua;
    temp_mua(i) = mua(i) + delta;
    mesh_jacobian = new_mesh(mesh_new,temp_mua,mus,kappa,"mesh_jacobian");
    fluence_data = new_femdata_spectral(mesh_jacobian,0);
    

    G(:,size(kappa,1)+i) = (fluence_data.phi - fluence.phi)/delta;
    G(:,size(kappa,1)+i) = G(:,size(kappa,1)+i).*mesh_new.mua;
    G(i,size(kappa,1)+i) = G(i,size(kappa,1)+i) + fluence.phi(i);
  
end


% Log the condition number of Jacobian matrix 

fid = fopen(fullfile('', 'conditionMatrix.log'), 'a');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid, 'Condition number of Jacobian matrix at %d iteration: %d\n', iter, cond(G));
fclose(fid);

end
