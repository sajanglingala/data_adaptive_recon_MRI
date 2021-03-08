%  Basic demon registration code. (To easy understand the algorithm)

% Clean

% Compile the mex files

% Read two images
function [perfout,Dx,Dy,Fx,Fy,Hsmooth]=demon_function_dynamic_nofig(f,g,Dx,Dy,sigma_gauss,noise_const);
perfout = zeros(size(f));
Fx = zeros(size(f)); Fy=Fx;
    [g]= scale(f, g);
addpath('/nfs/s-iib52/users/slingala/Research/demon_registration_version_8f/functions');
for nframe = 1:size(f,3)
I1 = f(:,:,nframe); I2=g(:,:,nframe);
%I1 = 256*(I1-min(I1(:)))/range(I1(:));
%I2 = 256*(I2-min(I2(:)))/range(I2(:));
% Set static and moving image
S=I2; M=I1;

% Alpha (noise) constant
alpha=noise_const; %default-2.5

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],sigma_gauss);
%Hsmooth=1;
% The transformation fields
Tx=Dx(:,:,nframe); Ty=Dy(:,:,nframe);%zeros(size(M));
%Tx=zeros(size(I1)); Ty=Tx;
[Sy,Sx] = gradient(S);
cost=[];
for itt=1:130
	    % Difference image between moving and static image
        Idiff=M-S;

        % Default demon force, (Thirion 1998)
       % Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
        %Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);

        % Extended demon force. With forces from the gradients from both
        % moving as static image. (Cachier 1999, He Wang 2005)
       [My,Mx] = gradient(M);
        Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
 
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

      %  Uxs = 3*imgaussian(Ux,sigma_gauss);
      %  Uys = 3*imgaussian(Uy,sigma_gauss);
        
        
        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        
        %[Tx,Ty] = expfield(Tx,Ty);
        M=movepixels(I1,Tx,Ty); 
        
        % cost =[cost,energy(S,M,Tx,Ty,alpha)];
         %if itt>1 && abs(cost(end) - cost(end-1))/abs(cost(end)) < 1e-10
         %   break;
        %end
        %figure(1);subplot(1,3,1); imagesc(abs(M(:,:,1))); subplot(1,3,2); imagesc(abs(S(:,:,1)));pause(0.1);
end
[Fx(:,:,nframe),Fy(:,:,nframe)]=backwards2forwards(Tx,Ty);
Dx(:,:,nframe)=Tx;
Dy(:,:,nframe)=Ty;
perfout(:,:,nframe)=M;
 subplot(1,4,1); imagesc(abs(M));subplot(1,4,2); imagesc(abs(I1)); subplot(1,4,3); imagesc(abs(I2));  title(num2str(nframe));
 subplot(1,4,4);plot(cost);pause(0.1);
end

	