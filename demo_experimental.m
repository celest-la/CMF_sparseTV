%%%  Demo of CMF algorithms. Reconstruction from experimental signals
%%%  bubbles
%%%  The reconstruction is done with frequency domain DAS and RCB,
%%%  CMF-ElNet and CMF-spTV
clear
warning('off','all')
addpath(genpath("functions/"));
addpath(genpath("data/"));
load("data/RFexpe.mat")

%--  field of view for the reconstruction
x = -2.0:0.1:5.3;
z = 31:0.1:37;
[X, Z] = meshgrid(x*1e-3, z*1e-3);

%-- Frequency for the reconstruction
f_recon = 4.2e6;

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
[CSM,axf] = compute_CSM(RF, param.fs, K, overlap, f_recon);
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
            
nexttile;
imagesc(x, z, 10*log10(rescale(map_DAS)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('DAS')

nexttile;
imagesc(x, z, 10*log10(rescale(map_RCB)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('RCB')

nexttile;
imagesc(x, z, 10*log10(rescale(map_enet)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('CMF-Elastic net')

nexttile;
imagesc(x, z, 10*log10(rescale(map_spTV)));
colormap('parula'); axis image; clim([-20,0]);
subtitle('CMF-spTV')
set(gcf,'color','w');