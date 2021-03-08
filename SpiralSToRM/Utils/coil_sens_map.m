function [csm,coilimages] = coil_sens_map(kdat,k,w,N)
numFramesforcsm=100;
kdata=kdat(:,:,1:numFramesforcsm,:);
[n,nint,nf,nCh]=size(kdata);
coilimages=zeros(N,N,nCh);

w = w(:,:,1:numFramesforcsm);
k=k(:,:,1:numFramesforcsm);
kdata=reshape(kdata,[n*nint*nf,nCh]);
[A1,B1] = creatA(1.5*transpose(k(:)),N);
At_DFT = @(x) INFT_new(x,A1,B1,sqrt(w(:)));

for i=1:nCh
    coilimages(:,:,i)=double(At_DFT(kdata(:,i).*sqrt(w(:))));
end

[A2, B2] = creatA_zeropadded(transpose(k(:)),N,2*N);
Q = ifftshift(fft2(fftshift(INFT_new(ones(n*nint*nf,1),A2,B2,(w(:)))/N^2)));
M=@(u)AhA(reshape(u,[size(coilimages,1), size(coilimages,2), nCh]),Q,N,nCh);
coilImages = pcg(M,coilimages(:),1e-6,100);

csm=giveEspiritMaps(reshape(coilimages,[size(coilimages,1), size(coilimages,2), nCh]));
csm=fftshift(fftshift(csm,1),2);
end

 function res= AhA(X,Q,n,nc)

res = zeros(n,n,nc);
tmp = zeros(2*n,2*n);

    for k=1:nc
        tmp(n/2+1:end-n/2,n/2+1:end-n/2)=X(:,:,k);
        tmp1=ifftshift(fft2(fftshift(tmp))).*Q;
        tmp1 = ifftshift(ifft2(fftshift(tmp1)));       
        res(:,:,k)=tmp1(n/2+1:end-n/2,n/2+1:end-n/2); 
        tmp(:)=0;
    end
    

res = res(:);
 end