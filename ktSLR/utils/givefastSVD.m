function [u,s,v] = givefastSVD(X);

[v,s] = eig(X'*X);s = sqrt(abs(s));

s = diag(s);
index = find(s >1e-6*max(s));
v = v(:,index);
s = diag(s(index));
u = X*v;

for i=1:size(s,2),
        u(:,i) = u(:,i)./sqrt(sum(abs(u(:,i)).^2));
end

