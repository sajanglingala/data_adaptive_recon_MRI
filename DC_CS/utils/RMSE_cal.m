
function [U,rmse]= RMSE_cal(ideal, recon)
 x = ideal; x_init = recon;

I=double(abs(x));U=(double(abs(x_init)));
 
 for ii=1:size(U,3)
alpha = sum(dot(U(:,:,ii),I(:,:,ii)))/(sum(dot(U(:,:,ii),U(:,:,ii))));

U(:,:,ii)=(alpha)*U(:,:,ii);

E=sum(sum(abs((I(:,:,ii)-U(:,:,ii))).^2));         

E=E*1/(sum(sum(abs((I(:,:,ii))).^2)));

Error(ii)=(E);
end

rmse=mean((Error))