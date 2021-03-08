clear all

addpath('./csm');
addpath('./Utils');
addpath('./nufft_toolbox_cpu');
addpath(genpath('./gpuNUFFT'));

%% Reconstruction parameters
spiralsToDelete=60;
ninterleavesPerFrame=6;
N = 340;
nChannelsToChoose=4;
numFramesToKeep = 200;
useGPU = 1;
SHRINK_FACTOR = 1.0;
nBasis = 30;
lambdaSmoothness = 0.05;
cRatioI=1:nChannelsToChoose;
sigma=[4.5];
lam=[0.3];
eig_csm=0.15;
%%
% % ==============================================================
% % Load the data
% % ============================================================== 
load('data/spiral_ungated_cardiac_cine_Series17.mat'); 
kdata=permute(kdata,[1,3,2]); % Fourth Slice Data
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
ktraj_scaled =  SHRINK_FACTOR*ktraj*N;clear ktraj;
%% ==============================================================
% Compute coil compresession and Compute the coil sensitivity map
% ============================================================== 
ktraj_scaled=reshape(ktraj_scaled,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep]);
kdata=reshape(kdata,[nFreqEncoding*ninterleavesPerFrame*numFramesToKeep,nChannelsToChoose]);
[vkdata] = combine_coils(kdata,0.9);
nChannelsToChoose=size(vkdata,2);
kdata=reshape(vkdata,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
[vcoilImages] = coil_sens_map_NUFFT(kdata,ktraj_scaled,N,useGPU);
csm=giveEspiritMaps(reshape(vcoilImages,[size(vcoilImages,1), size(vcoilImages,2), nChannelsToChoose]),eig_csm*max(abs(vcoilImages(:))));
ktraj_scaled=reshape(ktraj_scaled,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep]);
kdata=reshape(kdata,[nFreqEncoding,ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
clear vcoilImages vkdata;
%% ==============================================================
% % Compute the weight matrix
% % ============================================================= 
no_ch=size(csm,3);
Nav=permute((kdata(:,1,:,:)),[1,2,4,3]);
[~,x1,L]=estimateLapKernelLR(reshape(Nav,[nFreqEncoding*no_ch,numFramesToKeep]),sigma,lam);
[~,Sbasis,V]=svd(L);
V=V(:,end-nBasis+1:end);
Sbasis=Sbasis(end-nBasis+1:end,end-nBasis+1:end);
clear L Nav;
%% ==============================================================
% % Final Reconstruction
% % ============================================================= 
ktraj_scaled=reshape(ktraj_scaled,[nFreqEncoding*ninterleavesPerFrame,numFramesToKeep]);
kdata=reshape(kdata,[nFreqEncoding*ninterleavesPerFrame,numFramesToKeep,nChannelsToChoose]);
tic; x = solveUV(ktraj_scaled,kdata,csm, V, N, 50,lambdaSmoothness*Sbasis,useGPU);toc
y = reshape(reshape(x,[N*N,nBasis])*V',[N,N,numFramesToKeep]);

%% ==============================================================
% % Save and Display results
% % ============================================================= 
for i=1:size(y,3);imagesc(fliplr(abs(y(:,:,i)))); pause(0.1); colormap gray;end


