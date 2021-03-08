function[Error,U]= RMSE(U,I)
U = (U); I = (I);
for ii=1:size(U,3)
alpha = sum(dot(U(:,:,ii),I(:,:,ii)))/(sum(dot(U(:,:,ii),U(:,:,ii))));

U(:,:,ii)=(alpha)*U(:,:,ii);

E=sum(sum(abs((I(:,:,ii)-U(:,:,ii))).^2));         

E=E*1/(sum(sum(abs((I(:,:,ii))).^2)));

Error(ii)=(E);
end


Error=mean((Error));
