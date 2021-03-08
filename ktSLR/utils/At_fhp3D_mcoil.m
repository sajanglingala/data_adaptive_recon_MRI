function At = At_fhp3D_mcoil(z, S, n1,n2,n3,n4,smaps)

z=double(z); S = double(S);

p=zeros(n1,n2,n3,n4);
p(S)=z;

for nc = 1:n4
 
p(:,:,:,nc)= sqrt(n1*n2)*ifft2(squeeze(p(:,:,:,nc))).*conj(squeeze(smaps(:,:,:,nc)));
 
end

At = squeeze(sum(p,4));
