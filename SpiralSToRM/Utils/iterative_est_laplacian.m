
function [ X, A] = iterative_est_laplacian(A,FT,kdata,csm_lowRes,N1, sigSq, lam)
%function [ X, A] = iterative_est_laplacian(A,FT,kdata,csm_lowRes,N1, sigSq, lam)
%
% Arguments:
%   A       Laplacian matrix [nf nf]   complex
%   FT      Fourier operator
%   kdata   [np nv nf nc]   complex
%   csm_lowRes     coil sens map
%   lam      threshold
%   N1               reconstucted image size
%   sigSq            Sigma to control standard deviation

%
%
% Outputs:
%   X     [N N Nf]      low res images
%   A     [Nf Nf]      improved laplacian matrix
%
% Ahmed, Abdul Haseeb <abdul-ahmed@uiowa.edu>

nf = size(kdata,3);
W=(circshift(eye(nf),[0,1])+circshift(eye(nf),[1,0]));
W=(W+W')/2;
Lnn = (diag(sum(W,1))-W);
lambda2=0.0001;
Xpre=zeros(N1*N1,nf);

for i=1:30
  
    Atb = Atb_LR(FT,kdata,csm_lowRes,true);
    lambda1=lam*max(abs(Atb(:)));
    Reg1 = @(x) reshape(reshape(x,[N1*N1,nf])*(lambda1*A+lambda2*Lnn),[N1*N1*nf,1]);
    AtA = @(x) AtA_LR(FT,x,csm_lowRes,nf,N1)+Reg1(x);
    tic; [x1,~,~,it,res] = pcg(AtA,Atb(:),9e-4,50,[],[],Xpre(:));toc;%dataset virg 9e-4;
    X =(reshape(x1,[N1*N1,nf]));
   
    [~,~,A]=estimateLapKernelLR(X,sigSq,lam);
        
    if (norm(X(:)-Xpre(:))<1e-10) 
        break; 
    end 
    Xpre=reshape(X,[N1*N1,nf]);
end
end

