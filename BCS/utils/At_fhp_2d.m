% Backward Fourier transform sampling operator

function x = At_fhp_2d(y,s, OMEGA,n)
s=(s);

nc = n(3);
x = zeros(n(1),n(2),n(4));
OMEGA=double(OMEGA);

for ii=1:nc,
    fx = zeros(n(1),n(2),n(4));
    fx(OMEGA) = y(:,ii);
    x = x + (sqrt(n(1)*n(2))*ifft2(fx).*conj(squeeze(s(:,:,ii,:))));
end

x=reshape(x, n(1)*n(2),n(4));
