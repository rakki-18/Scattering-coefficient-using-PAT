% Calculates initial value of mua and mus,
% given the initial concentration of HbO, deoxyHbR, Water, scatter power, scatter amplitude and the wavelength array

function [mua, mus] = calc_initial_value(conc, scatter_parameters, Mesh)

%****************************************************************
% calculate absorption coefficients
E = []; mus = [];

wv_array = Mesh.wv;
sa = zeros(size(Mesh.nodes,1),1) + scatter_parameters(1);
sp = zeros(size(Mesh.nodes,1),1) + scatter_parameters(2);
concentration = zeros(size(Mesh.nodes,1),size(conc,2));
for i = 1: size(conc,2)
    concentration(:,i) = conc(1,i);
end
% for mus, wavelength must be in micrometers.
wv_act = wv_array/1000;

for i = 1:length(wv_array)
    index = find(Mesh.wv == wv_array(i));
    E = [E ;Mesh.excoef(index,:)];
    % Calculate mus
    mus = [mus; sa.*wv_act(i).^(-sp)];
end

% Calculate mua
mua = E*concentration';
mua = mua';
% Stack mua of all the wavelengths below each other
mua = mua(:);
kappa = 1./(3*(mus+mua));

end