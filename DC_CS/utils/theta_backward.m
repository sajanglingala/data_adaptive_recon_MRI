function A = theta_backward(z, Fx,Fy)
z=double(z);
A =  zeros(size(z));
for t = 1:size(z,3); A(:,:,t) = (movepixels((z(:,:,t)),Fx(:,:,t),Fy(:,:,t))); end
A=double(A);

