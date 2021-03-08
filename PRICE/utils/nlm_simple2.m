function [out,input]=nlm_simple2(input)
%--------------------------------------
% This function does the PRICE recon using CG algorithm.
% input -
%        input - input is a struct
%
% output -
%        out - reconstructed images using PRICE scheme 
%        input - input is a struct
%
% Please see : Y. Mohsin, S.G Lingala, E. DiBella, M.Jacob,
% Accelerated dynamic MRI Using Patch Regularization for Implicit motion CompEnsation (PRICE),
% Magnetic Resonance in Medicine, 2016.
%
% Last Edit: 11/27/2016
% Contact: Y. Mohsin (yasir-mohsin@uiowa.edu)
%          M. Jacob (mathews-jacob@uiowa.edu)
%
% Copyright 2016, CBIG Lab, College of Engineering, The University of Iowa.
% https://research.engineering.uiowa.edu/cbig/content/price
%--------------------------------------

out = input.initialguess;
input.globalenergy  = [];
%---------------------
DtD = computeDtD(input);
Atb=input.At(input.measurements);
energy_old=10^15;
for i=1:input.MaxOuterIterations, 
    i
    input.beta = 1.5*input.beta;% incrementing the continuation paremeter 
    input.lambeta = input.beta*input.lambda;
     for j=1:input.innerIterations, 
%% right hand side in image domain 
        [Dph,Regenergy,input] =giveH(out,input);    
        right = 2*Atb - input.lambeta*Dph ;
%% left hand side in image domain 
        left = @(z)lefthand(z,input,DtD);
        %  Solve for the recon: Conjugate gradient update
        [out,~,~,~,~] = pcg(left,right(:),1e-8,[],[],[],out(:)); 
        out = reshape(out,input.nx,input.ny,input.nz);
        % cost calculation
        err = input.A(out)-input.measurements;
        dataconsistency = sum(abs(err(:)).^2);
        input.globalenergy = [input.globalenergy,dataconsistency/2+ input.lambda*Regenergy];
     end
       test = abs(energy_old-input.globalenergy(end))/energy_old;
        energy_old = input.globalenergy(end);
        if(test<1e-6 && i>1)
            break;
        end
        input.T=input.T/1.1; % % decreasing the thresholding paremeter 
        if input.T<2.5;
           input.T=2.5;
        end        
        
%         figure(20);
%         subplot(1,2,1);plot(abs(double(input.globalenergy)));colormap(gray);
%         subplot(1,2,2);imagesc(abs(double(out(:,:,1))));colormap(gray);
end
end
