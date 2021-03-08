function data = l2Recont(kdata,FT,csm,lambda,N)
%function data = l2Recont(kdata,FT,csm,lambda,N)
%
% Arguments:
%   kdata   [np nv nf nc]   complex
%   FT      Fourier operator
%   csm     coil sens map
%   lambda      threshold
%   N               reconstucted image size
%
% Optional Arguments:

%
% Outputs:
%   data     [N N N Nt]      complex
%
% Ahmed, Abdul Haseeb <abdul-ahmed@uiowa.edu>

[~,~,nFrames,nCh] = size(kdata); 
Atb = zeros(N,N,nFrames);

for ii=1:nCh
    Atb = Atb + bsxfun(@times,FT'*kdata(:,:,:,ii),conj(csm(:,:,ii)));
end
lambda=max(abs(Atb(:)))*lambda;

ATA = @(x) SenseATA(x,FT,csm,N,nFrames,nCh) + lambda*x;
data = pcg(ATA,Atb(:),1e-5,60);
data = reshape(data,[N,N,nFrames]);

end