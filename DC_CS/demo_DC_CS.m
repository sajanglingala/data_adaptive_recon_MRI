% Author- Sajan Lingala
% Example code for: 
% S.G.Lingala, E.DiBella, M.Jacob, "Deformation corrected compressed sensing: a novel framework for
% accelerated dynamic MRI",IEEE Transactions on Medical Imaging, 2014

%%
clear;clc;close all; 

addpath(genpath('utils/demon_registration_version_8f/'));
addpath('utils/');
% load the fully sampled cardiac perfusion data
load ('data/Cardiac_perfusion_with_fake_motion.mat');


x=r1;
% initialize the deformation maps
% Dx, Dy are the forward maps
% Fx, Fy are the backward maps
Dx = zeros(size(x)); Dy= Dx; Fx = Dx; Fy = Dx;
% Assign the data size
n1=size(x,1); n2=size(x,2); n3=size(x,3);

% define a mask that contains regions of heart (used for MSE calculations
% during the iterations)
maskroi = zeros(size(Dx));
maskroi(75:125, 30:70,:)=1;



% Specify the no. of radial lines for under-sampling the k space; The
% Acceleration corresponds to the total no. of phase encodes/ no of radial lines 
% acquired
lines = [24];

% Call the function that performs golden angle radial sampling
% The radial grid is rotated by the golden angle in subsequent frames 
[T3D] = double(goldenratio_samp(n1,n2,n3,lines,0));
mask = fftshift(fftshift(T3D,1),2);
S=find(mask~=0);

%% Define the forward and inverse Fourier sampling operators
A = @(z)A_fhp3D(z,S);
At=@(z)At_fhp3D(z,S,n1,n2,n3);

% Perform undersampling on the fullysampled data to simulate measured k-t space
% data, b
b=A(x);
%% 



%% Define parameters for the DC-CS algorithm
% The algorithm solves the optimization problem: 
% min_{f,\theta) \|A(f)-b\|_{2}^{2} + \lambda \|\psi(T_{\theta(f)})
% \|_{1}


% load the initial guess (eg: Spatial TV reconstruction)
opts.temp=4e-09;
fn = sprintf('spatialTV_lines_lam_%s_%s.mat',num2str(lines),num2str(opts.temp))
load (fn); %initial guess
f_init =(Recon2(:,:,1:n3));

% Define the forward and backward image warping operators
T = @(z)theta_forward(z,Dx,Dy);
Tt = @(z)theta_backward(z,Fx,Fy);

opts.lambda=1e-12;%Regularization parameter;

% Initialize the cost function variables
earray=[];datacons=[];
regn=[];% Cost is stored in this

[D,Dt] = defDDtmodGPU; %Define the Forward and Backward Finite difference operator

f=f_init; % Initialize the reconstruction 
iter=0;

opts.beta =double(5*abs(1/max(double(f(:)))));

[g,earrayg,opts] = g_motionsub_mod_boundarycorrection(f,Dx,Dy,Fx,Fy,opts);

tf=T(f);
e=tf-g;out=0;
norm_tf_g=[]; mse_error=[];

% calculate the cost for initial guess
dc = A(f)-b;
cost = sum(abs(dc(:)).^2); 
tempq = D(T(f));
tempq=tempq{1};
earray = [earray,cost + opts.lambda*sum(abs(tempq(:)))];
datacons = [datacons, cost];
regn = [regn, opts.lambda*sum(abs(tempq(:)))];

% evaluate MSE for initial guess
[recon,mse]= RMSE_cal(double(x).*maskroi, double(f).*maskroi); 
mse_error=[mse_error,mse];


% Demons force strength parameter 
opts.alphademon=4; 







%%  START THE ITERATIVE ALGORITHM 


for out = 1:10 
   
for in = 1:1000

%% gsub
%% Spatio-temporal denoising

 [g,earrayg,opts] = g_motionsub_mod_boundarycorrection(f,Dx,Dy,Fx,Fy,opts);
 


%% Image registration
if in ==1 
    
 
   if out ==1 && in ==1
     [perfout,Dx,Dy,Fx,Fy,~]=demon_function_dynamic_nofig(double(abs((f))),double(abs((g))),((Dx)),((Dy)),10,opts.alphademon);
   end
    if out ==2 && in ==1
     [perfout,Dx,Dy,Fx,Fy,~]=demon_function_dynamic_nofig(double(abs((f))),double(abs((g))),((Dx)),((Dy)),10,opts.alphademon);
   end
    if out ==3 && in ==1
     [perfout,Dx,Dy,Fx,Fy,~]=demon_function_dynamic_nofig(double(abs((f))),double(abs((g))),((Dx)),((Dy)),10,opts.alphademon);
    end
    
end

% Update the warping operator
T = @(z)theta_forward(z,Dx,Dy);
Tt = @(z)theta_backward(z,Fx,Fy);


%% Image reconstruction 
[f,earrayf] = fsub_CG(b,A,At,T,Tt,g,opts, ((f)),1e-6,130);



%% Cost calculations
dc = A(f)-b;

cost = sum(abs(dc(:)).^2); 


tempq = D(T(f));
tempq=tempq{1};
earray = [earray,cost + opts.lambda*sum(abs(tempq(:)))];


datacons = [datacons, cost];
regn = [regn, opts.lambda*sum(abs(tempq(:)))];

[recon,mse]= RMSE_cal(double(x).*maskroi, double(f).*maskroi); 
%erecon = x(:,:,4)-recon;
mse_error=[mse_error,mse];
  
             tf=T(f);
e=tf-g;
norm_tf_g=[norm_tf_g,(norm(e(:))/norm(tf(:)))];


%% Plotting
figure(1); 
brighten(0.5);
subplot(2,3,1); plot(double(earray));  hold on; plot(double(datacons),'r'); hold on; plot(double(regn),'k'); hold off; legend('cost', 'data consistency','\lambda*reg');
subplot(2,3,2); imagesc(abs(f(:,:,5))); colormap(gray);
subplot(2,3,3); plot(mse_error);title('mse error');%imagesc(abs(erecon)); colormap(gray);
subplot(2,3,4); imagesc(abs(squeeze(double(g(100,:,:))))); colormap(gray);title('g');
subplot(2,3,5); imagesc(abs(squeeze(double(tf(100,:,:))))); colormap(gray);title('f')
subplot(2,3,6); plot(norm_tf_g); title('norm tf-g');
pause(0.1);
%%
         if in > 1
             if abs(earray(end) - earray(end-1))/abs(earray(end)) < 1e-3
            opts.gamma = opts.gamma*2;
         break;
             end
           

          end


end

opts.beta = opts.beta*10;
opts.alphademon=opts.alphademon*2;

end
earray=double(earray);
recon=double(recon);
mse_error=double(mse_error);
Dx=double(Dx);
Dy=double(Dy);
Fx=double(Fx);
Fy=double(Fy);
f=double(f);
tf=double(tf);
datacons=double(datacons);
regn=double(regn);


 
 
 

