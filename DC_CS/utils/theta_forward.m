function A = theta_forward(z, Dx,Dy)

z=double(z);
A =  zeros(size(z));
for t = 1:size(z,3); A(:,:,t) = (movepixels((z(:,:,t)),Dx(:,:,t),Dy(:,:,t))); end
A = double(A); 