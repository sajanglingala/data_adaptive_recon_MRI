function kdata=unchop(x);

S = size(x);
x = reshape(x,S(1),S(2),prod(S(3:end)));

%  x(1:2:end,:,:)=x(1:2:end,:,:)*-1;

x(:,1:2:end,:)=x(:,1:2:end,:)*-1;

kdata = reshape(x,S);
