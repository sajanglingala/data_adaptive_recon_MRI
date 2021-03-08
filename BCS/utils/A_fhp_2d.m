% Forward Fourier transform sampling operator

function y = A_fhp_2d(x, s, OMEGA,n,nc)

s=(s);

x=reshape(x,n(1),n(2),n(3));
y=zeros(size(OMEGA,1),nc);

for ii=1:nc,
    yc = (1/sqrt(n(1)*n(2)))*fft2(squeeze(s(:,:,ii,:)).*x);
    y(:,ii) = yc(OMEGA);
end


