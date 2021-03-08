clear all

addpath('./csm');
addpath('./Utils');
addpath('./nufft_toolbox_cpu');
addpath(genpath('./gpuNUFFT'));

%% Reconstruction parameters
spiralsToDelete=100;
ninterleavesPerFrame=5;
N = 340;
nChannelsToChoose=8;
numFramesToKeep = 500;
useGPU = 1;
SHRINK_FACTOR = 1.3;
nBasis = 20;
lambdaSmoothness = 0.05;
cRatioI=1:nChannelsToChoose;
NsamplesToKeep=600;
sigma=[4.5];
lam=[0.1];
eig_csm=0.01;
%%
% % ==============================================================
% % Load the data
% % ============================================================== 
load('data/spiral_ungated_cardiac_cine_series11.mat'); 
kdata=kdata(:,:,:,1); % Fourth Slice Data
kdata=kdata/max(abs(kdata(:)));

 %% =========================================
 % -------------Preprocessing Data-------------%
 %===========================================
[nFreqEncoding,nCh,numberSpirals]=size(kdata);
numFrames=floor((numberSpirals-spiralsToDelete)/ninterleavesPerFrame);
kdata=kdata(:,cRatioI(1:nChannelsToChoose),spiralsToDelete+1:numberSpirals);
k = k(:,spiralsToDelete+1:numberSpirals);
kdata=kdata(:,:,1:numFrames*ninterleavesPerFrame);
k=k(:,1:numFrames*ninterleavesPerFrame);

kdata = permute(kdata,[1,3,2]);
kdata = reshape(kdata,[nFreqEncoding,ninterleavesPerFrame,numFrames,nChannelsToChoose]);
ktraj=k;clear k;
ktraj = reshape(ktraj,[nFreqEncoding,ninterleavesPerFrame,numFrames]);

% Keeping only numFramesToKeep

kdata = kdata(:,:,1:numFramesToKeep,cRatioI(1:nChannelsToChoose));
ktraj = ktraj(:,:,1:numFramesToKeep);
%save data kdata ktraj dcf
%% ==============================================================
% Scaling trajectory
% ==============================================================
ktraj_scaled =  SHRINK_FACTOR*ktraj;clear ktraj;
%% ==============================================================
% Compute coil compresession and Compute the coil sensitivity map
% ============================================================== 
ktraj_scaled=reshape(ktraj_scaled,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep]);
kdata=reshape(kdata,[nFreqEncoding*ninterleavesPerFrame*numFramesToKeep,nChannelsToChoose]);
[vkdata] = combine_coils(kdata,0.75);
nChannelsToChoose=size(vkdata,2);
kdata=reshape(vkdata,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
[coilImages] = coil_sens_map_NUFFT(kdata,ktraj_scaled,N,useGPU);
csm=giveEspiritMaps(reshape(coilImages,[size(coilImages,1), size(coilImages,2), nChannelsToChoose]),eig_csm);
ktraj_scaled=reshape(ktraj_scaled,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep]);
kdata=reshape(kdata,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
clear vkdata;
%% ==============================================================
% % Compute the weight matrix
% % ============================================================= 
kdata_com = reshape(kdata(1:NsamplesToKeep,:,:,:),[NsamplesToKeep,ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
ktraj_com = reshape(ktraj_scaled(1:NsamplesToKeep,:,:,:),[NsamplesToKeep,ninterleavesPerFrame,numFramesToKeep]);
N1 = 64;
csm_lowRes = giveEspiritMapsSmall(coilImages,N1,N1);
ktraj_com = reshape(ktraj_com/N1,[size(ktraj_com,1),size(ktraj_com,2),numFramesToKeep]);
kdata_com = reshape(kdata_com, [size(kdata_com,1),size(ktraj_com,2),numFramesToKeep,nChannelsToChoose]);
FT_LR= NUFFT(ktraj_com,1,0,0,[N1,N1]);
tic;lowResRecons = l2Recont(kdata_com,FT_LR,csm_lowRes,0.1,N1);toc
lowResRecons=reshape(lowResRecons,[N1*N1,numFramesToKeep]);
[~,~,L]=estimateLapKernelLR(lowResRecons,sigma,lam);
[ tmp, L] = iterative_est_laplacian(L,FT_LR,kdata_com,csm_lowRes, N1,sigma, lam);
[~,Sbasis,V]=svd(L);
V=V(:,end-nBasis+1:end);
Sbasis=Sbasis(end-nBasis+1:end,end-nBasis+1:end);
%% ==============================================================
% % Final Reconstruction
% % ============================================================= 
ktraj_scaled=reshape(ktraj_scaled,[nFreqEncoding*ninterleavesPerFrame,numFramesToKeep]);
kdata=reshape(kdata,[nFreqEncoding*ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
tic; x = solveUV(ktraj_scaled,kdata,csm, V, N, 60,lambdaSmoothness*Sbasis,useGPU);toc
y = reshape(reshape(x,[N*N,nBasis])*V',[N,N,numFramesToKeep]);

%% ==============================================================
% % Save and Display results
% % ============================================================= 

for i=1:size(y,3);imagesc((rot90(abs(y(:,:,i)),-1))); pause(0.1); colormap gray;end
%for i=1:250;imagesc(((abs(y(:,:,i))))); pause(0.1); colormap gray; end
%clear kdata csm V ;
%save(strcat('res_iter_',num2str(lambdaSmoothness),'_',num2str(sigma(ii)),'_',num2str(sl),'.mat'), 'y','-v7.3');

