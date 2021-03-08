function At = Ah_fhp3Dmulticoil(z, S, n1,n2,n3,n4,coil_sens)
%--------------------------------------
% This function Defines the backward operator Atranspose (At)
%
% input -
%        z - input matrix/vector
%        S - Sampling indices
%        n1,n2,n3,n4 - data size
%        coil_sens - coil sensetivities map 
%
% output -
%        At - backword operator At applied on z: At(z)
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

z=double(z); S = double(S);
p=double(zeros(n1,n2,n3,n4));
p(S)=z;
for nc = 1:n4
    p(:,:,:,nc) = sqrt(n1*n2)*ifft2(p(:,:,:,nc)).*conj(coil_sens(:,:,:,nc));
end
At = sum(p,4);