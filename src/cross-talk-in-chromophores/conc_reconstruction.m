% Reconstructs the value of concentrations of each chromophore given the value of H and initial concentrations.
function conc = mua_reconstruction(H,conc,Mesh)

    %hyper parameters
    max_iterations = 3;
    regularisation_parameter = 0.1;
    
    
    
    
    mesh_new = new_mesh(Mesh,conc,"mesh_new");
    
    % Storing error in each iterations
    error_list = [;;];
    error_list_conc_HbO = [];
    error_list_conc_deoxy = [];
    % iterative update of mua values
    for i = 1:max_iterations
        fprintf("%d iterations started\n",i);
        % Calculate fluence at each wavelength for the current mesh
        [fluence_new,mua,mus] = femdata_spectral(mesh_new,0);
    
    
        G = find_jacobian(conc,mesh_new,fluence_new.phi.*mua,i);
    
        delta_t = find_delta_t(G,regularisation_parameter,H, fluence_new.phi.*mua);
    
        [conc,mesh_new] = update(delta_t,conc,Mesh);
    
        % Find error in calculation
        error_H = sum((fluence_new.phi.*mua - H).*(fluence_new.phi.*mua - H),1);
        error_list = [error_list [error_H;sum(abs(Mesh.conc(:,1) - conc(:,1)),1);sum(abs(Mesh.conc(:,2) - conc(:,2)),1)]];
        error_list_conc_HbO = [error_list_conc_HbO abs(Mesh.conc(:,1) - conc(:,1))];
        error_list_conc_deoxy = [error_list_conc_deoxy abs(Mesh.conc(:,2)-conc(:,2))];
    
        save('variables','error_list','error_list_conc_HbO', 'error_list_conc_deoxy','G','-append');
    
        
        % Plotting the concentrations of HbO and deoxyHbO
        figure;
        plotim(Mesh,conc(:,1));
        title('reconstructed HbO','FontSize',10);
        colorbar('horiz');

        figure;
        plotim(Mesh,conc(:,2));
        title('reconstructed deoxyHbO','FontSize',10);
        colorbar('horiz');
          
        
    end
    
    
    
    
    
       
    end
    
    % Creates a new mesh, "mesh_name" having the mua, mus and kappa values
    function [mesh_New] = new_mesh(Mesh,conc,~)
    mesh_New = Mesh;
    mesh_New.conc = conc;
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
    function [conc,mesh_new] = update(delta_t, conc,Mesh)
    % parameters
    learning_rate = 1;
    nodes = size(conc,1);

    for i = 1 : size(conc,2) - 1
        conc(:,i) = conc(:,i) + learning_rate*delta_t((i-1)*nodes+1:i*nodes);
    end

    mesh_new = new_mesh(Mesh,conc,"mesh_new");
    end

    % finding the jacobian matrix which is the derivative of energy distributin
    % wrt the conc of chromophores
    function [G] = find_jacobian(conc,mesh_new,H,iter)
    nodes = size(conc,1);
    G = zeros(size(H,1), nodes*(size(conc,2)-1)); % subtracting 1 because we are not estimating the concentration of water as that is already correct
    delta = 0.0001;

    %Iterate through all the chromophores except water
    for i = 1: size(conc,2)-1
        for j = 1: nodes
%         for j = 1: 10
            fprintf("%d iteration in Jacobian mua started. %d th iteration\n",j+(i-1)*nodes,iter);
            temp_conc = conc;
            temp_conc(j,i) = conc(j,i) + delta;
            mesh_jacobian = new_mesh(mesh_new,temp_conc,"mesh_jacobian");
            [fluence_data,mua,mus] = femdata_spectral(mesh_jacobian,0);


            G(:,j+(i-1)*nodes) = (fluence_data.phi.*mua - H)/delta;
            clear mesh_jacobian temp_conc fluence_data;
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
    