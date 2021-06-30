function conc = spectral_decomposition(mua, Mesh)
% Calculates the concentration of each chromophore at each node
% Given the mua values at each node for every wavelength and the spectral composition matrix
% mua = M*conc

conc = Mesh.excoef\mua;

conc = conc';

disp("about to display the plottings of HbO and deoxyHbO");
% Plotting the concentrations of HbO and deoxyHbO
figure;
plotim(Mesh,conc(:,1));
title('reconstructed HbO','FontSize',10);
colorbar('horiz');

figure;
plotim(Mesh,conc(:,2));
title('reconstructed deoxyHbO','FontSize',10);
colorbar('horiz');