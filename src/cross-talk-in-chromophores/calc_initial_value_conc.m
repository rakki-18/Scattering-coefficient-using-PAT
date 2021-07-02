% Calculates a homogeneous initial value of concentrations at each node
% given the initial concentration of HbO, deoxyHbR, Water

function conc = calc_initial_value(conc_abs, Mesh)

    conc = zeros(size(Mesh.conc));
    for i = 1: size(conc_abs)
        conc(:,i) = conc_abs(i);
    end
    
    
end