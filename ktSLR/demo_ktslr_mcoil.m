tic
clear;
%% load the fully sampled cartesian data
 load ('data/FB_ungated_cardiac_perfusion.mat');

 addpath('utils/');
 
for c = 1:size(kspacecoils,4), 
    for t = 1:size(kspacecoils,3)
      kspacecoils(:,:,t,c)=  fftshift(fft2(ifftshift(ifft2(kspacecoils(:,:,t,c)))));
    end
end

kspacecoils = kspacecoils(:,:,1:end,:);
[n1,n2,n3,n4]=size(kspacecoils);
% scaling factor to ensure the below choice of regularization parameters
% approximately holds.
kspacecoils = kspacecoils*(2.5067e-5)./max(abs(kspacecoils(:)));
k=kspacecoils;

imfor_coilmaps = sqrt(n1*n2)*ifft2(kspacecoils);

imfor_coilmaps = squeeze(mean((imfor_coilmaps),3));

smoothness = 15;
 coil_sens = double(ismrm_estimate_csm_walsh_modified(imfor_coilmaps, smoothness,0.0));
%[coil_sens] = gdouble(give_smooth_csm_jsense(imfor_coilmaps, 1,3));
%load coilmaps_192_jsense.mat;
%% load Undersampled 18 golden ratio lines
%load ('raw_datasets/Pt1_rest_moreframes_21.mat');
[n1,n2,n3,n4]=size(kspacecoils);
[mask] = strucrand(n1,n2,n3,16);
mask = fftshift(fftshift(mask,1),2);
mask = repmat(mask,[1 1 1 n4]);
 
%coil_sens_rep = gdouble(coil_sens);

 coil_sens_rep =zeros(n1,n2,n3,n4,'double');
 for nc = 1:n4,
    for nt = 1:n3,
         coil_sens_rep(:,:,nt,nc)=coil_sens(:,:,nc);
     end
 end

%coil_sens_rep = fftshift(fftshift(coil_sens_rep,1),2);

S=double(find(logical(mask)));

b = double(kspacecoils(S)); %% the Fourier data (radial regridded on Cartesian grid)

clear  maskcoils coil_sens gdatcoil ;
%% Define the forward and backward operators A and Atranspose (At)
% A = @(z)A_fhp3D(z,S);
% At=@(z)At_fhp3D(z,S,n1,n2,n3);
smaps=coil_sens_rep;
A = @(z)A_fhp3D_mcoil(z, S,n1,n2,n3,n4,smaps);
At = @(z)At_fhp3D_mcoil(z, S, n1,n2,n3,n4,smaps);


%% Define the forward and backward gradient operators.
% This function has an option of changing the step sizes of the gradients along the x, y and t dimensions
% For ex: step_size = [1,1,0.303] implies Grad_x(U) = [U(x+1,y,t)-U(x,y,t)]/1,
% Grad_y(U) = [U(x,y+1,t)-U(x,y,t)]/1, Grad_t(U) = [U(x,y,t+1)-U(x,y,t)]/0.303
step_size = [1,1,0.303];
[D,Dt] = defDDt(step_size);

%% First guess, direct IFFT
x_init = At(b);

%% Call kt slr using augmented Lagrangian with continuation (refer: S.G.Lingala et al, ISBI 2011)


 mu1 =1e-10; % Regularization parameter for schatten p-norm
 mu2 =4e-9; % Reg. parameter for spatiotemporal TV norm * Note: temporal wt weighted 10 times higher than spatial weight (see DefDDtmod.m)
 opts.mu1 = mu1;
 opts.mu2 = mu2;
 opts.p=0.1; % The value of p in Schatten p-norm; p=0.1: non convex; p = 1: convex
 [~,sq,~]=givefastSVD(reshape(x_init, n1*n2,n3)); % find the singular values of the initial guess
 opts.beta1=10./max(sq(:));% The continuation parameter for low rank norm; initialize it as 1./max(singular value of x_init)
 opts.beta2=10./max(abs(x_init(:))); % The continuation parameter for the TV norm; Initialize it as 1./max((x_init(:)))
 opts.beta1rate = 50; % The continuation parametr increment for low rank norm
 opts.beta2rate = 25; % similar increment for TV norm
 opts.outer_iter =15; % no of outer iterations - INCREASE THIS TO BE CONSERVATIVE
 opts.inner_iter = 50; % no of inner iterations

 
 [Recon,cost,opts] = minSNandTV(A,At,D,Dt, x_init,b,1,opts);
 %%
 close all; 
 
 
 x = sum(ifft2(kspacecoils).*conj(coil_sens_rep),4);
 
 [Error,Recon]= RMSE(Recon,x);

figure(1); colormap(gray);
subplot(3,3,1); imagesc(abs(x(:,:,24))); title('Fully sampled gold standard, a spatial frame'); 
subplot(3,3,2); imagesc(abs(x_init(:,:,24))); title('Direct IFFT, a spatial frame'); 
subplot(3,3,3); imagesc(abs(Recon(:,:,24))); title('k-t SLR, a spatial frame'); 
subplot(3,3,4); imagesc(abs(squeeze(x(130,:,:)))); title('Gold standard, image time profile'); 
subplot(3,3,5); imagesc(abs(squeeze(x_init(130,:,:)))); title('Direct IFFT, image time profile'); 
subplot(3,3,6); imagesc(abs(squeeze(Recon(130,:,:)))); title('k-t SLR, image time profile'); 
subplot(3,3,8); imagesc(abs(fftshift(mask(:,:,5)))); title('Radial sampling: one frame'); 
subplot(3,3,9); plot(abs(cost),'linewidth',2); title('Cost v/s iteration');  

 
