function A = A_fhp3D_mcoil(z, S,n1,n2,n3,n4,smaps)
z=double(z); S = double(S);
%[n1,n2,n3] = size(z);
p = zeros(n1,n2,n3,n4); 
for i = 1: n4, 
    p(:,:,:,i)=1/sqrt(n1*n2)*fft2(z.*smaps(:,:,:,i));
end
A = p(S);
