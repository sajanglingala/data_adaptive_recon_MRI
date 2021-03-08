function  DtD = computeDtD(ipStruct)
%--------------------------------------
% %This function computes the finite difference operator of the Dirac function.
% input -
%        ipStruct - input is a struct
%
% output -
%        DTD - finite difference operator.
%
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

DtD = zeros(ipStruct.nx,ipStruct.ny,ipStruct.nz); 
Dirac = zeros(ipStruct.nx,ipStruct.ny,ipStruct.nz); 
Dirac(1,1,1) = 1;
%% Calculating DtD
for k = ipStruct.Wint,
    for i=ipStruct.Winx,
        for j=ipStruct.Winy,
            if i==0&&j==0
                w=10;
                
            else
                w=ipStruct.w;
            end
          DDirac = circshift(Dirac,[i,j,k])-Dirac;
          DDtDirac = DDirac-circshift(DDirac,[-i,-j,-k]);
          DtD = DtD + DDtDirac*w;
         end
    end
end
