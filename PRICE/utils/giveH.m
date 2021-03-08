function [Dph,Regenergy,input]=giveH(uh1,input)
%--------------------------------------
% This function computes the finite difference operator DpH. 
% Also it computes the energy used for plotting the convergence. 

% This function computes the finite difference operator of the Dirac function.
% input -
%        uh1 - input matrix/vector
%        input - input is a struct
%
% output -
%        Dph - finite difference operator 
%        Regenergy - energy used for the cost fucntion
%        input - input is a struct

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

index = 1;
Regenergy = 0;
T=input.T;
Dph = double(zeros(input.nx,input.ny,input.nz));
for k = input.Wint,
    for i=input.Winx,
        for j=input.Winy,
            if i==0&&j==0
                w=10;
            else
                w=input.w;
            end
            if( index ~= input.midIndex)
              diff = circshift(uh1,[i,j,k]) - uh1;% f(y)-f(x)               
              PXPY = sqrt(ifft2(fft2(abs(diff).^2).*input.GaussImage));%||Px-Py||_heta               
              Regenergy = Regenergy + sum(input.phi(PXPY(:),T,input)/w).^(1/input.p);    
              %% hxy
              wxy = input.psi(PXPY,input,T)*w;
              wxy = ifft2(fft2(abs(wxy)).*input.GaussImage);
              h = wxy.*diff;
              %Dph 
              Dph = Dph + (h - circshift(h,[-i,-j,-k]));
            end
            index = index + 1;
        end
    end
end
clear index diff diff2 diff3 wxy h  