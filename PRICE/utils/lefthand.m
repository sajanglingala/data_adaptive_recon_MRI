function Acc = lefthand(out,input,DtD) 
%--------------------------------------
% This fucntion solves the left hand side of the PRICE scheme.
% input -
%        out - input matrix/vector
%        DtD - finite difference operator of the Dirac function
%        input - input is a struct
% output -
%        Acc - output matrix/vector
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

out = reshape(out,input.nx,input.ny,input.nz);
AtAf = input.At(input.A(out));
DTDf = ifftn(fftn(DtD).*fftn(out)); 
Acc = 2*AtAf-input.lambeta*DTDf; 
Acc = Acc(:);
end