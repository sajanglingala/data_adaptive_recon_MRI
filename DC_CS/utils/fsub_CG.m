%--------------------------------------------------------------------------
% CG solution to region limited problem
% [X,Potential,ErrorArray,ErrorIndex,efinal] = xupdateHOTV(A,b,baseline,mask,kappa,lambda,mu,Niter,Potential)
% Solves {X*} = arg min_{X} ||Af-b||^2 + mu ||Rf||_{l_1}
%--------------------------------------------------------------------------

function [X,earray1] = fsub_CG(b,A,At,T,Tt,g,C, X, THRESHOLD,Niter)

oldcost = 0;
earray1 = [];
for i=1:Niter,
    
    Ax =A(X);
    resY = (Ax - b);
    eY = sum(abs(resY(:)).^2);
    
    resQ = T(X)-g; 
    eQ = sum (abs(resQ(:)).^2);
    
    
    
    cost1 = eY + (C.lambda*C.beta/2)*eQ;
    
    earray1 = [earray1,cost1];
    
    if(abs(cost1-oldcost)/abs(cost1) < THRESHOLD)
      %  i
        break;
    end
    oldcost = cost1;
    
  %  conjugate gradient direction
   % ------------------------------
    
    % gradient: gn
    gn = At(resY) + (C.lambda*C.beta/2)*Tt(resQ) ;%+ lam2*Dt(resTV);
    gn = 2*gn;%+C.mu1*Lam4 + C.mu2*Dt(LamTV);
    
    % search direction: sn  
    if(i==1)
        sn = gn;                                          
        oldgn = gn;
    else
        gamma = abs(sum(sum(conj(gn(:)).*gn(:)))/sum(sum(conj(oldgn(:)).*oldgn(:))));
        sn = gn + gamma*sn; 
        oldgn = gn;
    end
    
    % line search
    %-------------
    
    Asn = A(sn); 
    Tsn = T(sn);
    
    
    numer = Asn(:)'*resY(:)  +   (C.beta*C.lambda/2)*Tsn(:)'*resQ(:);
  
    denom = Asn(:)'*Asn(:)  +   (C.lambda*C.beta/2)*Tsn(:)'*Tsn(:);
  
  
    if(denom < 1e-18)
       % break;
    end
    alpha = -real(numer)/real(denom);
   
    % updating
    %-------------
    
    X = (X + alpha*sn);
end

    
