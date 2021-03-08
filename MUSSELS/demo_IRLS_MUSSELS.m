close all
clearvars

Nint=4;
addpath(genpath('utils/'));


%% Load coil sensitivity maps, k-space data
csm= load('data/CSM_diffusion.mat','S1');
csm=(csm.S1);

tmp=load('data/diffusion_raw_data_MS4.mat');
raw_data=(tmp.resampledKSpace);
raw_data=raw_data./max(raw_data(:));

N=size(raw_data);
[N1,N2,Nch]=deal(N(1),N(2),N(3));

for ph=1:61
    %% Separate k-space data from each shot
    k_dwi_ms=unchop(squeeze(raw_data(:,:,:,ph)));
    
    kdata=zeros(N1,256,Nch,Nint);
    for int=1:Nint
        kdata(:,int:Nint:N2,:,int)=k_dwi_ms(:,int:Nint:N2,:);
    end
    mask=squeeze(sum(abs(kdata),3))>0;
    
    %% IRLS MUSSELS Recon
    
    tic;
    %[rec] =mussels_cs(kdata,csm,mask,[6,6],10,8,.005,1,0);
    [rec] =mussels_cs(kdata,csm,mask,[4,4],5,5,.005,1,0); % faster
    toc;
    
    dwi(:,:,ph)=sos(ifft2c(rec));
    figure();imagesc(cat(2,rot90(sos(ifft2c(sum(kdata,4))),1),rot90(dwi(:,:,ph),1)));colormap(gray);title(['without phase compensation  ','  MUSSELS Recon'])
    
    
end

