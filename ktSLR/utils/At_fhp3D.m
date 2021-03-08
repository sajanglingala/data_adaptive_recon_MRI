function At = At_fhp3D(z, S, n1,n2,n3)
z=double(z); S = single(S);
p=zeros(n1,n2,n3);
p(S)=z(S);

At=sqrt(n1*n2)*ifft2(p);
