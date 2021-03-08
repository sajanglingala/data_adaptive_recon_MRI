function [g,earray,opts] = g_motionsub_mod_boundarycorrection(f,Dx,Dy,Fx,Fy,opts)

T = @(z)theta_forward(z,Dx,Dy);
Tt = @(z)theta_backward(z,Fx,Fy);

q = T(f);

Q=double(fft(double(q),[],3));



[D,Dt] = defDDtmodGPU;


[n1,n2,n3]=size(q);
Omega = double(conj(psf2otf([1,-1],[n1*n2 n3])));% Omegasquare
Omegasq = double(abs(psf2otf([1,-1],[n1*n2 n3])).^2);
Omega=reshape(Omega,n1,n2,n3);
Omegasq=reshape(Omegasq,n1,n2,n3);


[n1,n2,n3]=size(q);
earray=[];
opts.outer_alpha = 20; 
opts.inner_alpha = 40;

w=zeros(size(q));datac=[]; penalty=[];g=zeros(size(q));


% 
% G = opts.beta*Q+ opts.gamma*Omega.*fft(w,[],3);
%        G=G./(opts.beta + opts.gamma*(Omegasq));
%        g = ifft(G,[],3); temp = D(g);
opts.gamma = abs(1/max(double(q(:)))); %opts.beta=1e-3;



for out_alpha = single(1:opts.outer_alpha)
    
   for in_alpha = single(1:opts.inner_alpha)
       
       %g
       %[g,earray_g] = gsub_CG(D,Dt,w,q,opts, g, 1e-6,30);

       G = opts.beta*Q+ opts.gamma*Omega.*double(fft(double(w),[],3));
       
       G=G./(opts.beta + opts.gamma*(Omegasq));
       
       
       g = double(ifft(double(G),[],3));
       
      % plot(earray_g);
       % w
       
           w = D(g);w=w{1};grad_g=w;
    V = abs(w).^2;% + abs(Z2).^2 + abs(Z3).^2;
    V = sqrt(V);
    V(V==0) = 1;
    thresh_TV = (1/opts.gamma).*(V.^(1-1));
    
    V = max(V - thresh_TV, 0)./V;
    w = w.*V;
       
       
       const=w-grad_g;
       const = sum(abs(const(:)).^2);
       
       
%          w = D(g);w=w{1};grad_g = w;
%          w_1 = (abs(w)-1/opts.alpha); w_1 = w_1.*(w_1>0); 
%          w_2 = abs(w)+(abs(w)<1e-12);
%          w = w_1.*w./w_2;
        
       
         % check cost;
         dc = g-q; dc = sum(abs(dc(:)).^2);
         earray = [earray,( sum(abs(grad_g(:))) + (opts.beta/2)*(dc))];
         datac=[datac,(opts.beta/2)*dc];
         penalty=[penalty,sum(abs(grad_g(:)))];
% %          
%         figure(1);
%        subplot(2,4,1);plot(earray); 
%         subplot(2,4,2); imagesc(abs(squeeze(g(100,:,:))));
%         subplot(2,4,3);plot(datac); 
%         subplot(2,4,4);plot(penalty);
%         pause(0.1);
         
         if in_alpha > 1
             if abs(earray(end) - earray(end-1))/abs(earray(end)) < 1e-3
           % opts.gamma = opts.gamma*2;
         break;
             end
          end
             
            % if in_alpha>10
            % if abs(penalty(end) - penalty(end-1))/abs(penalty(end)) < 1e-10
             % break;
            %end
        % end
       
       
       
       
       
       
       
       
       
   end
  
   opts.gamma=opts.gamma*30;
     %if const < 1e-5 && abs(earray(end) - earray(end-1))/abs(earray(end)) < 1e-7
       %      break;
            end
    
end







