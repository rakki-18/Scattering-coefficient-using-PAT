% Reconstructs the value of mua and mus given the value of H and the initial mua and mus values
function mua_mus_reconstruction(H,initial_value,Mesh)

%hyper parameters
max_iterations = 4;
regularisation_parameter = 0.01;

nodes = size(H,1);

% initial guess of mua and mus
mua = zeros(nodes,1) + initial_value(1);
mus = zeros(nodes,1) + initial_value(2);
kappa = find_kappa(mua,mus);
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");

% iterative update of mua and mus values
error_list = [;;];
error_list_mua = [];
error_list_mus = [];
for i = 1:max_iterations
    fprintf("%d iterations started\n",i);
    fluence_new = femdata('./MeshSample/mesh_new/mesh_new',0);
    
    % Taking average along the rows as  there are more than one source
    fluence_new.phi = sum(fluence_new.phi, 2)/size(fluence_new.phi,2);
    
    G = find_jacobian(mua,mus,kappa,mesh_new,fluence_new,i);
    delta_t = find_delta_t(G,regularisation_parameter,H, fluence_new.phi.*mua);
    
    [mua,mus,kappa,mesh_new] = update(delta_t, mua,mus,kappa,Mesh);
    % Find error in calculation
    error_H = sum((fluence_new.phi.*mua - H).*(fluence_new.phi.*mua - H),1);
    save('variables','-append');
    error_list = [error_list [error_H;sum(abs(Mesh.mua - mua),1);sum(abs(Mesh.mus - mus),1)]];
    error_mua = abs(Mesh.mua - mua);
    error_mus = abs(Mesh.mus - mus);
    error_list_mua = [error_list_mua error_mua];
    error_list_mus = [error_list_mus error_mus];
    %% PLOTTING RESULTS
    figure;
    plotim(Mesh,mua);
    title('mua obtained','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,mus);
    title('mus obtained','FontSize',20);
    colorbar('horiz');
    save('variables','-append');
    
    
    
end





   
end

% Creates a new mesh, "mesh_name" having the mua, mus and kappa values
function [mesh_new] = new_mesh(Mesh,mua,mus,kappa,mesh_name)
mesh_new = Mesh;
mesh_new.mua = mua;
mesh_new.mus = mus;
mesh_new.kappa = kappa;
mesh_loc = "./MeshSample/" + mesh_name + "/" + mesh_name;
mesh_loc = char(mesh_loc);
save_mesh(mesh_new,mesh_loc);
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
kappa = zeros(size(mua));
for i = 1: size(mua,1)
    kappa(i) = 3*(mua(i) + mus(i));
    kappa(i) = 1/kappa(i);
end
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

G = zeros(size(mua,1), 2*size(mua,1));
delta = 0.0001;
for i = 1: size(kappa,1)
% for i = 1: 10
    fprintf("%d iteration in Jacobian started. %d th iteration\n",i,iter);
    temp_mus = mus;
    temp_mus(i) = mus(i) + delta;
    temp_kappa = find_kappa(mua,temp_mus);
    new_mesh(mesh_new,mua,temp_mus, temp_kappa,"mesh_jacobian");
    fluence_data = femdata('./MeshSample/mesh_jacobian/mesh_jacobian',0);
    % Taking average along the rows as  there are more than one source
    fluence_data.phi = sum(fluence_data.phi, 2)/size(fluence_data.phi,2);
    
    delta_kappa = temp_kappa(i) - kappa(i);
    
    for j = 1: size(kappa,1)
        dphi = (fluence_data.phi(j) - fluence.phi(j))/delta_kappa;
        G(j,i) = mua(j)*dphi;
    end
end

for i = 1: size(mua,1)
% for i = 1: 10   
    fprintf("%d iteration in Jacobian mua started. %d th iteration\n",i,iter);
    temp_mua = mua;
    temp_mua(i) = mua(i) + delta;
    new_mesh(mesh_new,temp_mua,mus,kappa,"mesh_jacobian");
    fluence_data = femdata('./MeshSample/mesh_jacobian/mesh_jacobian',0);
    % Taking average along the rows as  there are more than one source
    fluence_data.phi = sum(fluence_data.phi, 2)/size(fluence_data.phi,2);
    for j = 1: size(mua,1)
        if(i~=j)
            dphi = (fluence_data.phi(j) - fluence.phi(j))/delta;
            G(j,size(kappa,1)+i) = mua(j)*dphi;
        else
            dphi = (fluence_data.phi(j) - fluence.phi(j))/delta;
            G(j,size(kappa,1)+i) = fluence.phi(j) + mua(i)*dphi;
        end
    end
end

% Log the condition number of Jacobian matrix 

fid = fopen(fullfile('', 'conditionMatrix.log'), 'a');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid, 'Condition number of Jacobian matrix at %d iteration: %d\n', iter, cond(G));
fclose(fid);

end
