% Reconstructs the value of mua and mus given the value of H and the initial mua and mus values
initial_value = [0.01,1.00];
mua_reconstruction(H,initial_value,Mesh);
function mua_reconstruction(H,initial_value,Mesh)

%hyper parameters
max_iterations = 100;
regularisation_parameter = 40;

nodes = size(H,1);

% initial guess of mua and mus
mua = zeros(nodes,1) + initial_value(1);
mus = zeros(nodes,1) + initial_value(2);
kappa = find_kappa(mua,mus);
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");


for i = 1:max_iterations
    fprintf("%d iterations started",i);
    fluence_new = femdata('./MeshSample/mesh_new',0);
    G = find_jacobian(mua,mus,kappa,mesh_new,fluence_new);
    save('variables');
    delta_t = find_delta_t(G,regularisation_parameter,H, fluence_new.phi.*mua);
    [mua,mus,kappa,mesh_new] = update(delta_t, mua,mus,kappa,Mesh);
    
end

end

function [mesh_new] = new_mesh(Mesh,mua,mus,kappa,mesh_name)
mesh_new = Mesh;
mesh_new.mua = mua;
mesh_new.mus = mus;
mesh_new.kappa = kappa;
mesh_loc = "./MeshSample/" + mesh_name;
mesh_loc = char(mesh_loc);
save_mesh(mesh_new,mesh_loc);
end

function [delta_t] = find_delta_t(G, regularisation_parameter,H, Hcal)
hessian = transpose(G)*G;
delta_t = hessian + regularisation_parameter;
delta_t = inv(delta_t);
delta_t = delta_t*transpose(G)*(H - Hcal);
end

function [kappa] = find_kappa(mua,mus)
kappa = zeros(size(mua));
for i = 1: size(mua,1)
    kappa(i) = 3*(mua(i) + mus(i));
    kappa(i) = 1/kappa(i);
end
end

function [mua,mus,kappa,mesh_new] = update(delta_t, mua, mus, kappa,Mesh)
mua = mua + delta_t(end/2+1:end);
kappa = kappa + delta_t(1:end/2);
mus = (1./(3.*kappa))-mua;
mesh_new = new_mesh(Mesh,mua,mus,kappa,"mesh_new");
end

function [G] = find_jacobian(mua,mus,kappa,mesh_new,fluence)
G = zeros(size(mua,1), 2*size(mua,1));
delta = 0.0001;
% for i = 1: size(kappa,1)
for i = 1: 10
    fprintf("%d iteration in Jacobian started\n",i);
    temp_mua = mua;
    temp_mua(i) = mua(i) + delta;
    temp_mus = mus;
    temp_mus(i) = mus(i) + delta;
    temp_kappa = find_kappa(temp_mua,temp_mus);
    new_mesh(mesh_new,temp_mua,temp_mus, temp_kappa,"mesh_jacobian");
    fluence_data = femdata('./MeshSample/mesh_jacobian',0);
    delta_kappa = temp_kappa(i) - kappa(i);
    
    for j = 1: size(kappa,1)
        dphi = (fluence_data.phi(j) - fluence.phi(j))/delta_kappa;
        G(j,i) = mua(j)*dphi;
    end
end

% for i = 1: size(mua,1)
for i = 1: 10
    fprintf("%d iteration in Jacobian mua started\n",i);
    temp_mua = mua;
    temp_mua(i) = mua(i) + delta;
    new_mesh(mesh_new,temp_mua,mus,kappa,"mesh_jacobian");
    fluence_data = femdata('./MeshSample/mesh_jacobian',0);
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
    
end