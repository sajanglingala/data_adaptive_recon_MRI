function Atb = Atb_UV(FT,kdata,V,csm,N,useGPU)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,~,nChannelsToChoose] = size(csm);
nBasis = size(V,2);
Atb = zeros(N,N,nBasis);


if(useGPU)
    for i=1:nBasis
      for j=1:nChannelsToChoose
         temp = squeeze(kdata(:,:,j))*diag(V(:,i));
         Atb(:,:,i) = Atb(:,:,i) + (FT'*temp(:)).*conj(csm(:,:,j));  
      end
    end
else
     for i=1:nBasis
      for j=1:nChannelsToChoose
         temp = squeeze(kdata(:,:,j))*diag(V(:,i));
         %temp=reshape(temp,[2496,4,450]);
         Atb(:,:,i) = Atb(:,:,i) + (FT'*temp).*conj(csm(:,:,j));  
      end
    end

end

