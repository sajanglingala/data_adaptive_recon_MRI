
function [U]= scale(ideal, recon)
 x = ideal; x_init = recon;

I=(double(x));U=(double(x_init));
 
 for ii=1:size(U,3)
alpha = sum(dot(U(:,:,ii),I(:,:,ii)))/(sum(dot(U(:,:,ii),U(:,:,ii))));

U(:,:,ii)=(alpha)*U(:,:,ii);

end
