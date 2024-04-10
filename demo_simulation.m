%%%  Demo of CMF algorithms. Reconstruction of the simulation of two axially distributed
%%%  bubbles
%%%  The reconstruction is done with frequency domain DAS and RCB,
%%%  CMF-ElNet and CMF-spTV
clear
warning('off','all')
addpath(genpath("functions/"));
addpath(genpath("data/"));
load("data/RFsimu.mat")

%--  field of view for the reconstruction
x = -10:0.2:10;
z = 60:0.2:80;
[X, Z] = meshgrid(x*1e-3, z*1e-3);

%-- Frequency for the reconstruction
f_recon = 4e6;

%-- CSM parameters
K = 130;
overlap = 0.9;

%% Optimization parameters
lambdaL1 = 0.5;
lambdaL2 = 0.3;
lambdaTVlat = 0.25;
lambdaTVax = 0.03;
lambdaL1TV = 0.3;
rho = 4;
inv_opTV = [];
inv_op_enet = [];

%% Beamforming
disp('Compute CSM and steering')
[CSM, axf] = compute_CSM(RF, param.fs, K, overlap, f_recon);
[A, C] = getAC(RF, f_recon, param, X, Z, K, overlap);

disp('Compute PCI')
map_DAS = do_map('pci', CSM, axf, X, Z, param);

disp('Compute RCB')
map_RCB = do_map('rcb', CSM, axf, X, Z, param, 'eps', 5);

disp('Compute CMF-Elestic net')
[map_enet, history_enet, inv_op_enet] = CMF(A, C, X, 0, 0, lambdaL1, lambdaL2, rho, inv_op_enet);  

disp('Compute CMF-sparse TV')
[map_spTV, history, inv_opTV] = CMF(A, C, X, lambdaTVlat, lambdaTVax, lambdaL1TV, 0, rho, inv_opTV);  
   
%% Display the results
figure;
            
subplot(2,2,1);
imagesc(x, z, 10*log10(rescale(map_DAS)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('DAS')

subplot(2,2,2);
imagesc(x, z, 10*log10(rescale(map_RCB)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('RCB')

subplot(2,2,3);
imagesc(x, z, 10*log10(rescale(map_enet)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('CMF-Elastic net')

subplot(2,2,4);
imagesc(x, z, 10*log10(rescale(map_spTV)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('CMF-spTV')

%-- plot the correct position
for i=1:4
    subplot(2,2,i); 
    hold on; plot(hc1.p(:,1)*1e3, hc1.p(:,3)*1e3, 'r+', 'linewidth', 2); hold off
    legend("True position", 'location', 'south');
end
set(gcf,'color','w');