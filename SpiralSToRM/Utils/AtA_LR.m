function y = AtA_LR(FT,x,csm,nf,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,~,nChannelsToChoose] = size(csm);
x = reshape(x,[N,N,nf]);
y = zeros(N,N,nf);

for j=1:nChannelsToChoose
        temp = FT'*(FT*(bsxfun(@times,x,csm(:,:,j)))); 
        y = y + bsxfun(@times,temp,conj(csm(:,:,j)));  
end

y = y(:);

end

