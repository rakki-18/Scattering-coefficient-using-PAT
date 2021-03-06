% Reconstructs the value of mua and mus given the value of H , initial mua and mus values
function mua_mus_reconstruction(H,mua,mus,Mesh)

% Take logarithm of the energy distribution as the range is big
H = log(H);
save('variables','H','-append');
%hyper parameters
max_iterations = 3;
regularisation_parameter = 2;


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
    % Reducing the regularisation parameter by half in each iteration
    regularisation_parameter = regularisation_parameter/2;
    % Find error in calculation
%     error_H = sum((fluence_new.phi.*mua - H).*(fluence_new.phi.*mua - H),1);
%     error_list = [error_list [error_H;sum(abs(Mesh.mua - mua),1);sum(abs(Mesh.mus - mus),1)]];
    error_mua = abs(Mesh.mua - mua);
    error_mus = abs(Mesh.mus - mus);
    error_list_mua = [error_list_mua error_mua];
    error_list_mus = [error_list_mus error_mus];
    G_log = find_jacobian(mua,mus,kappa,mesh_new,fluence_new,i);
%     G_mua = find_jacobian_mua(mua,mus,kappa,mesh_new,fluence_new,i);
    fprintf("Finding G done\n");
    Hcal = log(fluence_new.phi.*mua);
    fprintf("Finding Hcal done\n");
    delta_t = find_delta_t(G_log,regularisation_parameter,H, Hcal);
    fprintf("Finding deltaT done\n");
%     [mua,mus,kappa,mesh_new] = update_mua(delta_t, mua,mus,kappa,Mesh);
    [mua,mus,kappa,mesh_new] = update(delta_t, mua,mus,kappa,Mesh);
    fprintf("Finding updates done\n");
    % Updating Mus and kappa values now.
%     fluence_new = new_femdata_spectral(mesh_new,0);
%     G_kappa = find_jacobian_kappa(mua,mus,kappa,mesh_new,fluence_new,i);
%     delta_t = find_delta_t(G_kappa,regularisation_parameter,H, fluence_new.phi.*mua);
%     [mua,mus,kappa,mesh_new] = update_mus(delta_t, mua,mus,kappa,Mesh);
%     clear delta_t G_kappa G_mus fluence_new;
    save('variables','G_log','Hcal','delta_t','error_list','error_list_mua','error_list_mus','-append');
    %% PLOTTING RESULTS for first wavelength
    figure;
    plotim(Mesh,mua(1:nodes));
    title('mua obtained','FontSize',20);
    colorbar('horiz');

    figure;
    plotim(Mesh,mus(1:nodes));
    title('mus obtained','FontSize',20);
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
% delta_t = zeros(2*size(H,1),1);
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
learning_rate_kappa = 1;
mua = mua + learning_rate_mua*delta_t(end/2+1:end);
kappa = kappa + learning_rate_kappa*delta_t(1:end/2);
% mus = mus + learning_rate_kappa*delta_t(1:end/2);
mus = (1./(3.*kappa))-mua;
% kappa = find_kappa(mua,mus);
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");
end
% updates the mua values based on delta_t
function [mua,mus,kappa,mesh_new] = update_mus(delta_t, mua, mus, kappa,Mesh)
    % learning_rate_mua = 1;
    learning_rate_kappa = 1;
    % mua = mua + learning_rate_mua*delta_t(end/2+1:end);
    kappa = kappa + learning_rate_kappa*delta_t;
    % mus = mus + learning_rate_kappa*delta_t(1:end/2);
    mus = (1./(3.*kappa))-mua;
    % kappa = find_kappa(mua,mus);
    mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");
end
% updates the mus and kappa values based on delta_t
function [mua,mus,kappa,mesh_new] = update_mua(delta_t, mua, mus, kappa,Mesh)
    learning_rate_mua = 1;
    % learning_rate_kappa = 1;
    mua = mua + learning_rate_mua*delta_t;
    % kappa = kappa + learning_rate_kappa*delta_t(1:end/2);
    % mus = mus + learning_rate_kappa*delta_t(1:end/2);
    % mus = (1./(3.*kappa))-mua;
    % kappa = find_kappa(mua,mus);
    mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");
end

% finding the jacobian matrix which is the derivative of energy distributin
% wrt kappa and mua values
function [G_log] = find_jacobian(mua,mus,kappa,mesh_new,fluence,iter)

G_log = zeros(size(mesh_new.mua,1), 2*size(mesh_new.mua,1));
delta = 0.0001;

for i = 1: size(kappa,1)
% for i = 1: 10
    fprintf("%d iteration in Jacobian started. %d th iteration\n",i,iter);
    temp_kappa = kappa;
    temp_kappa(i) = kappa(i) + delta;
    temp_mus = (1./(3.*temp_kappa))-mua;
%     temp_mus = mus;
%     temp_mus(i) = mus(i) + delta;
%     temp_kappa = find_kappa(mua,temp_mus);
    mesh_jacobian = new_mesh(mesh_new,mua,temp_mus, temp_kappa,"mesh_jacobian");
    
    
    fluence_data = new_femdata_spectral(mesh_jacobian,0);

    
    G_log(:,i) = (fluence_data.phi - fluence.phi)/delta;
    phiInverse = 1./fluence.phi;
    G_log(:,i) = G_log(:,i).*phiInverse;
    clear temp_kappa temp_mus fluence_data mesh_jacobian;
end

for i = 1: size(mua,1)
% for i = 1: 10 
    fprintf("%d iteration in Jacobian mua started. %d th iteration\n",i,iter);
    temp_mua = mua;
    temp_mua(i) = mua(i) + delta;
    mesh_jacobian = new_mesh(mesh_new,temp_mua,mus,kappa,"mesh_jacobian");
    fluence_data = new_femdata_spectral(mesh_jacobian,0);
    

    G_log(:,size(kappa,1)+i) = (fluence_data.phi - fluence.phi)/delta;
    phiInverse = 1./fluence.phi;
    G_log(:,size(kappa,1)+i) = G_log(:,size(kappa,1)+i).*phiInverse;
    muaInverse = 1./mua;
    G_log(i,size(kappa,1)+i) = G_log(i,size(kappa,1)+i) + muaInverse(i);
    clear mesh_jacobian temp_mua fluence_data;
end


% Log the condition number of Jacobian matrix 

fid = fopen(fullfile('', 'conditionMatrix.log'), 'a');
if fid == -1
  error('Cannot open log file.');
end
% fprintf(fid, 'Condition number of Jacobian matrix at %d iteration: %d\n', iter, cond(G_log));
fprintf("Condition number done");
fclose(fid);

end
function [G] = find_jacobian_mua(mua,mus,kappa,mesh_new,fluence,iter)
    G = zeros(size(mesh_new.mua,1), size(mesh_new.mua,1));
    delta = 0.0001;
    
    for i = 1: size(mua,1)
%     for i = 1: 10 
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

function [G] = find_jacobian_kappa(mua,mus,kappa,mesh_new,fluence,iter)

G = zeros(size(mesh_new.mua,1), size(mesh_new.mua,1));
delta = 0.0001;

for i = 1: size(kappa,1)
% for i = 1: 10
    fprintf("%d iteration in Jacobian started. %d th iteration\n",i,iter);
    temp_kappa = kappa;
    temp_kappa(i) = kappa(i) + delta;
    temp_mus = (1./(3.*temp_kappa))-mua;
    % temp_mus = mus;
    % temp_mus(i) = mus(i) + delta;
    % temp_kappa = find_kappa(mua,temp_mus);
    mesh_jacobian = new_mesh(mesh_new,mua,temp_mus, temp_kappa,"mesh_jacobian");
    
    
    fluence_data = new_femdata_spectral(mesh_jacobian,0);

    
    G(:,i) = (fluence_data.phi - fluence.phi)/delta;
    G(:,i) = G(:,i).*mua;
    clear temp_kappa temp_mus fluence_data mesh_jacobian;
end

% Log the condition number of Jacobian matrix 

fid = fopen(fullfile('', 'conditionMatrix.log'), 'a');
if fid == -1
  error('Cannot open log file.');
end
fprintf(fid, 'Condition number of Jacobian matrix at %d iteration: %d\n', iter, cond(G));
fclose(fid);
end
