%--------------------------------------------------------------------------
% CG solution to region limited problem
% [X,Potential,ErrorArray,ErrorIndex,efinal] = xupdateHOTV(A,b,baseline,mask,kappa,lambda,mu,Niter,Potential)
% Solves {X*} = arg min_{X} ||Af-b||^2 + mu ||Rf||_{l_1}
%--------------------------------------------------------------------------

function [X,earray1] = CG_solver(b,A, At,D,Dt,w,Lam1,Lam2,Lam3,Lam4,Y1,Y2,Y3,C, X, THRESHOLD,Niter)

oldcost = 0;
earray1 = [];
lam1 = 0.5*C.mu1*C.beta1;
%lam1 = 10000;
lam2 = 0.5*C.mu2*C.beta2;

LamTV = cell(3,1);
LamTV{1} = Lam1; 
LamTV{2} = Lam2;
LamTV{3} = Lam3;
Lam4=Lam4+1i*1e-18;
eTV=double(0);eNN=double(0);
for i=1:Niter,
    
    resY = (A(X) - b);
    eY = sum(abs(resY(:)).^2);
    
    resw = X-w; 
    eNN = lam1*sum(abs(resw(:)).^2);
     eNN = eNN + abs(C.mu1*(Lam4(:)'*resw(:)));
    
    Dx = D(X);
    resTV = cell(3,1);
    resTV{1} = Dx{1}-Y1; resTV{2} = Dx{2}-Y2; resTV{3} = Dx{3} - Y3;
% 
    
    %eTV = eTV + C.mu2*abs(conj(LamTV{1}).*resTV{1} + conj(LamTV{2}).*resTV{2} + conj(LamTV{3}).*resTV{3});
    %eTV = sum(eTV(:));
    LamTV1=LamTV{1}+1i*1e-18;LamTV2=LamTV{2}+1i*1e-18;LamTV3=LamTV{3}+1i*1e-18;
    resTV1=resTV{1}; resTV2=resTV{2};resTV3=resTV{3};
    
    eTV = lam2*(abs(resTV1(:)'*resTV1(:)).^2+abs(resTV2(:)'*resTV2(:)).^2+abs(resTV3(:)'*resTV3(:)).^2);
    
    eTV = eTV + C.mu2 * abs(conj(LamTV1(:)'*resTV1(:)) + conj(LamTV2(:)'*resTV2(:)) + conj(LamTV3(:)'*resTV3(:))); 
    
    
    cost1 = eY + eNN + eTV;
    
    earray1 = [earray1,cost1];
    
    if(abs(cost1-oldcost)/abs(cost1) < THRESHOLD)
       % i
        break;
    end
    oldcost = cost1;
    
  %  conjugate gradient direction
   % ------------------------------
    
    % gradient: gn
    
    gn = At(A(X)-b) + lam1*(X-w) + lam2*Dt(resTV);
    gn = 2*gn+C.mu1*Lam4 + C.mu2*Dt(LamTV);
    
    % search direction: sn  
    if(i==1)
        sn = gn;                                          
        oldgn = gn;
    else
        gamma = abs(sum(gn(:)'*gn(:))/sum(oldgn(:)'*oldgn(:)));
        sn = gn + gamma*sn; 
        oldgn = gn;
    end
    
    % line search
    %-------------
    Asn = A(sn);  
    Dsn = D(sn);
    
    numer = Asn(:)'*resY(:) + lam1*sn(:)'*resw(:) + 0.5* C.mu1*sn(:)'*Lam4(:);
    for index = 1:3
         numer = numer + lam2*sum(sum(sum(conj(resTV{index}).*Dsn{index})))   +  0.5*C.mu2*sum(sum(sum(conj(LamTV{index}).*Dsn{index}))); 
    end
    
    denom = Asn(:)'*Asn(:) + lam1*sn(:)'*sn(:); 
    for index = 1:3 
       denom = denom + lam2*sum(sum(sum(conj(Dsn{index}).*Dsn{index}))); 
    end
    if(denom < 1e-18)
        break;
    end
    alpha = -real(numer)/real(denom);
   
    % updating
    %-------------
    
    X = (X + alpha*sn);
end

    
