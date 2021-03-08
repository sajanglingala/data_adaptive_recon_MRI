% This program is free software: you can use or modify for research uses.  
%--------------------------------------
% This is the main file. It recovers a 4D free-breathing ECG-gated MRI cardiac CINE dataset in undersampled kspace form.  
% check Readme.pdf file for more details about this dataset and the parameters used in this code.

% Please see : Y. Mohsin, S.G Lingala, E. DiBella, M.Jacob,
% Accelerated dynamic MRI Using Patch Regularization for Implicit motion CompEnsation (PRICE),
% Magnetic Resonance in Medicine, 2016
%
% Last Edit: 11/27/2016
% Contact: Y. Mohsin (yasir-mohsin@uiowa.edu)
%          M. Jacob (mathews-jacob@uiowa.edu)
%
% Copyright 2016, CBIG Lab, College of Engineering, The University of Iowa.
% https://research.engineering.uiowa.edu/cbig/content/price
%--------------------------------------

clear all;close all;clc;
addpath('utils/');
load('data/freebreathing_ungated_Cardiac_cine_Cartesian.mat'); % laoding the data
subplot(321);imagesc(abs(double((kdata(:,:,1,1)))));colormap gray
mask = zeros(size(kdata)); % extracting the mask for the kspace data
mask(kdata~=0)=1;
subplot(322);imagesc(abs(double((mask(:,:,6,1)))));

% fft shift the kspace data 
for c = 1:size(kdata,4)
    for t =1:size(kdata,3)
        knew(:,:,t,c) =  fft2(ifftshift(ifft2(fftshift(kdata(:,:,t,c))))).*fftshift(mask(:,:,t,c));
    end
end
subplot(323);imagesc(abs(double(knew(:,:,1,1))));
clear kdata mask c t;
% flipud and fliplr the coil sensitivities provided to make it
% consistent to the data
for c = 1:size(b1,3)
    b1(:,:,c) = fliplr(flipud(b1(:,:,c))); 
end
subplot(324);imagesc(abs(double(b1(:,:,1))));

S=double(find(knew));
b=double((knew(S)));

n1 = size(knew,1); n2 = size(knew,2); n3 =size(knew,3); n4 = size(knew,4);

coil_sens_rep =double(zeros(n1,n2,n3,n4));% make the 4D version of the CSM
for nc = 1:n4,
    for nt = 1:n3,
    coil_sens_rep(:,:,nt,nc)=b1(:,:,nc);
    end
end
clear knew b1 nt c nc

%%  Define the forward and backward operators (A) and Atranspose (At)
A = @(z)A_fhp3Dmulticoil(z, S,n1,n2,n3,n4,coil_sens_rep);
At = @(z)Ah_fhp3Dmulticoil(z, S, n1,n2,n3,n4,coil_sens_rep);

%%
x_init = (At(b));
scale=255/max(x_init(:)); % scale the data from 0-255
x_init=x_init*scale;
b=b*scale;
acc=n1*n2*n3*n4/size(S,1) % acceleration factor
subplot(325); for t = 1:20, imagesc(abs(x_init(:,:,t))); title('initial guess');pause(0.1); end
subplot(326); imagesc(abs(mean(x_init,3))); title('averaged image: initial guess');
%%
%%%%%   Set up input %%%%%%
clear input;
input.measurements = b;
input.initialguess = x_init;
input.A = A;
input.At = At;
input.lambda = 1e-1;
input.beta = 1e-2;
input.innerIterations = 10;
input.MaxOuterIterations = 15;
input.np = 1;
input.nwin = 1;
input.nt =3;
input.w=0.1;
input.T=50;
input.p=0.6;
input.INNERTHRESHOLD = 10^(-6);
input.OUTERTHRESHOLD = 10^(-5);
input.earray  = [];
input.globalcost = [];
input.sigmaPatch = 1/3;
input.distancefn = 'L1';
input.debug = true;
input.doublesided = 1;
input = initializeVariables(input);

input.phi = @(z,v,ip) L1_sat(z,v,ip);
input.psi = @(z,ip,v) shrinkage_L1_sat(z,ip,v);

% Call PRICE using CG with continuation 
tic;[out,input]=nlm_simple2(input);toc

for i=1:20
pause(0.1);colormap(gray)
imagesc(double(abs(out(:,:,i))));title('recon'); 
end