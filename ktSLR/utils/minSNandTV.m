function [U,cost,opts] = minSNandTV(A,At,D,Dt, x_init,b,T3D,opts);

U=x_init; [m n d] = size(U); 

%% Lam1, Lam2, Lam3 and Lam4 are respectively the Lagrange multipliers. 
% Lam4 is for low rank norm while Lam1, Lam2 and Lam3 are for the x, y and t gradients
% refer SG Lingala et al, ISBI 2011 for details
Lam4 = zeros(m,n,d);Lam1=Lam4;Lam2=Lam4;Lam3=Lam4;

%%
Dfwd=D(U);

o=0;cost=[];
for out = single(1:opts.outer_iter),
    o=o+1
    for in = single(1:opts.inner_iter)
        
     
 % ================================
        %  Begin Alternating Minimization
        % ----------------
        %   w-subprolem (multi dimensional TV shrinkaeg)
        % ----------------
     
    Z1 = Dfwd{1} + Lam1/opts.beta2;
    Z2 = Dfwd{2} + Lam2/opts.beta2;
    Z3 = Dfwd{3} + Lam3/opts.beta2;
    V = abs(Z1).^2 + abs(Z2).^2 + abs(Z3).^2;
    V = sqrt(V);
    V(V==0) = 1;
    V = max(V - 1/opts.beta2, 0)./V;
    Wx = Z1.*V;
    Wy = Z2.*V;
    Wt = Z3.*V;
      % --------------------
        %  Lambda - subproblem (singular value shrnkage)
        % --------------------
         [u,sigma,v] = givefastSVD(reshape((U+Lam4/opts.beta1),m*n,d));
         s=diag(sigma);
      
       
        thres=(1/opts.beta1).*(s.^(opts.p-1));
        
        s=(s-thres);
        s = s.*(s>0);
        Lambda=u*(diag(s))*v';
        Lambda=reshape(Lambda,m,n,d);
     %   imagesc(abs(double(Lambda(:,:,10)))); 
        
        % ----------------
        %  Solve for the recon: Conjugate gradient update
        % ----------------
      

 [U,earray1] = CG_solver(b,A, At,D,Dt,Lambda,Lam1,Lam2,Lam3,Lam4,Wx,Wy,Wt,opts, U, 1e-5,30);
%% plot a frame of the recon while it iterates
   figure(3); 
    subplot(1,2,1);imagesc(abs(double((U(:,:,31))))); title('A frame of the reconstruction'); colormap(gray);
    
 %% cost calculations
 e = A(U) - b;     
 [uq,sigmaq,vq] = givefastSVD(reshape((U),m*n,d));
 V1 = sum( sum(sum( abs(sqrt(Dfwd{1}.^2 + Dfwd{2}.^2 + Dfwd{3}.^2)))));
      
 cost = [cost, sum(abs(e(:)).^2)  +  sum(abs(sigmaq(:)).^(opts.p)./(opts.p))*opts.mu1 + V1*opts.mu2 ];
 figure(3);subplot(1,2,2);plot(double(cost)); hold on; pause(0.1); hold off; title('Cost');
 
 
 
 if in>1
        if abs(cost(end) - cost(end-1))/abs(cost(end)) < 1e-3
            break;
  end
        end

Dfwd=D(U);

%%        Update rules for the Lagrange multipliers
%     Lam4 = Lam4 - 1.618*opts.beta1*(Lambda - U);
%     
%     Lam1 = Lam1 - 1.618*opts.beta2*(Wx - Dfwd{1});
%     
%     Lam2 = Lam2 - 1.618*opts.beta2*(Wy - Dfwd{2});
%     
%     Lam3 = Lam3 - 1.618*opts.beta2*(Wt - Dfwd{3});
    
    
    
    end 
    
    %% increment beta 1 and beta2
    opts.beta1=opts.beta1*opts.beta1rate;

    opts.beta2=opts.beta2*opts.beta2rate;
end
end 