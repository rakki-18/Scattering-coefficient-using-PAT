% Calculates initial value of mua and mus,
% given the initial concentration of HbO, deoxyHbR, Water, scatter power, scatter amplitude and the wavelength array

function initial_value = calc_initial_value(conc, scatter_parameters, wv,)

%****************************************************************
% calculate absorption coefficients
E = []; mus = [];

% for mus, wavelength must be in micrometers.
wv_act = wv_array/1000;

for i = 1:length(wv_array)
    index = find(mesh.wv == wv_array(i));
    E = [E;mesh.excoef(index,:)];
    mus = [mus; scatter_parameters(1).*wv_act(i).^(-scatter_parameters(2))];
end

mua = E*conc;
mua = mua';
kappa = 1./(3*(mus+mua));
    
% Iterate through wavelength
for
end