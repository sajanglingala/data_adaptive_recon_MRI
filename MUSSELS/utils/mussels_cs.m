function [x] = mussels_cs(b,Sen,SM,ksize,Niter,iter,Lambda,gamma,r);

%   INPUTS:
%   b: k-space data of size Nx x Npf x Nch X Nsh
%   Sen: coil sensitivity maps of size Nx x Ny x Nch x MB
%   SM : Sampling pattern of the shots Nx x Ny x Nch x MB
%   iter : number of outer iterations
%   Niter  : number of inner iterations of pcg
%   ksize: filter size
%   rank: the threshold of hard thresholding
%   Eg: tic;recon= mussels_cs(b,Sen,SM,[6,6],5,5,5e4,1,0);toc

[N1,N2,Nch,Nsh]=size(b);
k=ksize(1);


%%% define the Hankel operator and its inverse
T2 = @ (z) getHankel(z,ksize);
T2t= @ (z) invHankel(z,N1,N2,Nsh,ksize);
Nsz=ksize(1)*ksize(2)*Nsh;


C = @ (z,tflag) afun_pcg(z,Sen,SM,tflag);
x0=C(b,'transp'); clear b
x=reshape(x0,N1,N2,Nsh);


for int=1:Nsh
    x4(:,:,int)=conj(rot90(squeeze(x(:,:,int)),2));
end

D=cat(2,T2(x),T2(x4));
[U,S,~] = svd(D'*D,'econ');
eps=1*S(r+1,r+1)^0.5;


for in =1: iter
    C = @ (z,tflag) afun_pcg(z,Sen,SM,tflag);
    
    % ================================
    % Nuclear norm implemented as IRLS
    % --------------------
    %  ||Ax-b||^2 + Lambda* ||T(x)||*
    % =||Ax-b||^2 + Lambda* ||T(x)Wn^0.5||^2_F
    % W=(T2'T2)^-0.5
    % --------------------
    
    s=diag(S+eps).^(-0.5);
    W=(U*diag(s));
    eps=min(eps,gamma*s(r+1));
    %  figure(41); plot(diag(S).^-0.5); hold on;plot((s),'r');pause(.01);hold off
    %
    %   figure(42); plot(diag(S),'b'); hold on; plot(diag(S+eps),'r');hold off
    
    % ================================
    %  PCG: CG Sense update
    % 2A'Ax+2Lambda*T2'(T2(x)Wn^0.5Wn'^0.5)=2A'b
    % 2A'Ax+2Lambda*G'(W)G(W)x=2A'b
    % 2A'Ax+2lambda*deconv(conv(f,x),f')=2A'b
    % 2A'Ax+2lambda*IFFT(FFT(IFFT(FFT(f')*FFT(x'))),FFT(f'))=2A'b
    % ================================
    
    W=W*W';
    
    WWt = @ (z) getWWt(z,W,T2,T2t,Nsz);
    
    H1=@(x)C(x,'notransp')+(Lambda*reshape(WWt(reshape((x),N1,N2,Nsh)),[],1));
    
    
    [x,~,~,~,~] = pcg(H1,2.*x0(:),1e-10,Niter,[],[],x(:));
    
    x=reshape((x),N1,N2,Nsh);
    %     figure(53);imshow(rot90(sos(ifft2c(x)),0),[],'InitialMagnification','fit');title(['after CG, iter : ',num2str(in)]);pause(.01);
    
    for int=1:Nsh
        x4(:,:,int)=conj(rot90(squeeze(x(:,:,int)),2));
    end
    
    D=cat(2,T2(x),T2(x4));
    [U,S,~] = svd(D'*D,'econ');
    
    
end
end

%------------------- pcg ----------------%
function y = afun_pcg(data,Sen,SM,transp_flag)
CSen=conj(Sen);
N=size(Sen);

if strcmp(transp_flag,'transp')
    
    for j=1:size(data,4) % Nint
        for i=1:N(3) % Nchannel
            idata_tmp1(:,:,i)=sqrt((N(1)*N(2))).*CSen(:,:,i).*fftshift(ifft2(fftshift(data(:,:,i,j))));%.*mask(:,:,j)
        end
        y(:,:,j)=fftshift(fft2(fftshift(sum(idata_tmp1(:,:,:),3))))./sqrt((N(1)*N(2)));
    end
    
elseif strcmp(transp_flag,'notransp')
    
    data=reshape(data,N(1),N(2),[]);
    
    for j=1:size(data,3) % Nint
        gdata = sqrt((N(1)*N(2))).*fftshift(ifft2(fftshift(data(:,:,j))));
        for i= 1:size(Sen,3)% Nchannel
            gdata1=fftshift(fft2(fftshift(gdata.*Sen(:,:,i)))).*SM(:,:,j)./sqrt((N(1)*N(2)));%
            gdata2(:,:,i)=sqrt((N(1)*N(2))).*fftshift(ifft2(fftshift(gdata1))).*CSen(:,:,i);
        end
        gdata4(:,:,j)=fftshift(fft2(fftshift(sum(gdata2,3))))./sqrt((N(1)*N(2)));
    end
    y=2.*gdata4(:);
    
end
end



function fwd_Im = getWWt(x,W,T2,T2t,Nsz);

for int=1:size(x,3)
    x4(:,:,int)=conj(rot90(squeeze(x(:,:,int)),2));
end
D=cat(2,T2(x),T2(x4));
D=D*W;
x1=T2t(D(:,1:Nsz));
x2=T2t(D(:,Nsz+1:end));
for int=1:size(x,3)
    fwd_Im(:,:,int)=x1(:,:,int)+conj(rot90(squeeze(x2(:,:,int)),2));
end

end

%-----------------------------------%
