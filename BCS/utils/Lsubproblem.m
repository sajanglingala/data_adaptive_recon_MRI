% Shrinkage function 

function L = Lsubproblem(opts, phi, U)
if strcmp(opts.phi, 'TV');
    phifwd = phi(U);
    TV = sqrt(abs(phifwd{1}).^2 + abs(phifwd{2}).^2) ;
    TV(TV==0) = 1;
    TV = max(TV - (1/opts.beta1), 0)./TV;
    L{1} = double(phifwd{1}.*TV);
    L{2} = double(phifwd{2}.*TV);
elseif strcmp(opts.phi, 'l1');
    Z = U;
    L_1 = (abs(Z)-1/opts.beta1); 
    L_1 = L_1.*(L_1>0); 
    L_2 = abs(Z)+(abs(Z)<1e-12);
    L = L_1.*Z./L_2;
    L=(L);
    
elseif strcmp(opts.phi,'lp')
    Z = U ;
    L_1 = (abs(Z)-1/opts.beta1*((abs(Z)+0.00001).^(opts.p-1))); 
    L_1 = L_1.*(L_1>0); 
    L_2 = abs(Z)+(abs(Z)<1e-12);
    L = L_1.*Z./L_2;
    L=double(L);
end
