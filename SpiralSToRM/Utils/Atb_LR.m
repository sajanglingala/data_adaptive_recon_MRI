function Atb = Atb_LR(FT,kdata,csm,useGPU)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[N,~,nChannelsToChoose] = size(csm);
nf = size(kdata,3);
Atb = zeros(N,N,nf);


if(useGPU)
    for j=1:nChannelsToChoose  
         Atb = Atb + bsxfun(@times,(FT'*kdata(:,:,:,j)),conj(csm(:,:,j)));   
    end
end

