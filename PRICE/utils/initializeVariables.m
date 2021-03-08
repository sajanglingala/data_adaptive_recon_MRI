function ip=initializeVariables(ip)
%--------------------------------------
% This function initializes the Gaussian kernel weight 
% input -
%        ip - ip is a struct
%
% output -
%        ip - ip is a struct
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

[ip.nx,ip.ny,ip.nz] = size(ip.initialguess);

[x,y] = meshgrid([0:ip.ny/2,-ip.ny/2+1:-1],[0:ip.nx/2,-ip.nx/2+1:-1]);
if(ip.np==0)
%     Gauss = 0*x;
%     Gauss(1,1)=1;
Gauss = ones(ip.nx,ip.ny);
else
    Gauss = exp(-(x.^2+y.^2)/(ip.sigmaPatch*ip.np^2));
%     Gauss = Gauss.*(abs(x)<=1);
%     Gauss = Gauss.*(abs(y)<=1);
    Gauss = Gauss/(sum(sum(Gauss)));
    Gauss = fft2(Gauss);
end
ip.GaussImage = repmat(Gauss,[1,1,ip.nz]);

% Search window
%---------------

ip.Winx = [-ip.nwin:ip.nwin];
ip.Winy = [-ip.nwin:ip.nwin];
if(ip.doublesided)
    ip.Wint = [-ip.nt:ip.nt]; % compare the frame before and the frame after
    ip.sizewts = length(ip.Winx)*length(ip.Winy)*length(ip.Wint);
    ip.midIndex = floor(ip.sizewts/2)+1;
else
    ip.Wint = [0:ip.nt];% compare the frame after only
    ip.sizewts = length(ip.Winx)*length(ip.Winy)*length(ip.Wint);
    ip.midIndex = floor(length(ip.Winx)*length(ip.Winy)/2)+1;
end
% ip.w = zeros(ip.nx,ip.ny,ip.nz,ip.sizewts);
% Specifying the distance function
%----------------------------------
end
