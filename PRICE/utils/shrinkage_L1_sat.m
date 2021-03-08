function out =  shrinkage_L1_sat(in,input,T)
%--------------------------------------
% This function solves the shrinkage part of the scheme and provides the result obtained by applying
% the l_p threshold to in 
% input -
%        in - input matrix/vector
%        T - threshold parameter
%        input.p - p-value for the l_p penalty (input is a struct)
%
% output -
%        out - thresholded matrix/vector
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

mask = abs(in) < T;
wt = max(abs(in) - (1/input.beta).*((abs(in)).^(input.p-1)),0)./(abs(in)+(in==0));
out = wt.*mask + not(mask);
end



              

             

             
