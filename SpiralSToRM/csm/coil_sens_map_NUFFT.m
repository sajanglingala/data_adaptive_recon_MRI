function [coilimages] = coil_sens_map_NUFFT(kdata,trajectory,N,useGPU)

[nReadouts,nInterleaves,nFrames,nCh] = size(kdata); 
coilimages=zeros(N,N,nCh);

if(useGPU)
    osf = 2; wg = 3; sw = 8;
    %w=ones(nReadouts*nInterleaves*nFrames,1);
    %w=repmat(dcf,[1 nInterleaves*nFrames]);
    ktraj_gpu = [real(trajectory(:)),imag(trajectory(:))]';
    FT = gpuNUFFT(ktraj_gpu/N,ones(nReadouts*nInterleaves*nFrames,1),osf,wg,sw,[N,N],[],true);
else
    
    %w=repmat(dcf,[1 nInterleaves]);
    FT= NUFFT(trajectory(:)/N,ones(nReadouts*nInterleaves*nFrames,1),0,0,[N,N]);
    %Atb = Atb_UV(FT,kdata,V,csm,false);
    %AtA = @(x) AtA_UV(FT,x,V,csm,nFreqEncoding*ninterleavesPerFrame);
end

ATA = @(x) reshape(FT'*(FT*reshape(x,[N,N,1])),[N*N,1]);

kdata = reshape(kdata,[nReadouts*nInterleaves*nFrames,nCh]);
for i=1:nCh
    temp = FT'*kdata(:,i);
    coilimages(:,:,i)= reshape(pcg(ATA,temp(:),1e-6,70),[N,N]);
end
%csm=giveEspiritMaps(reshape(coilimages,[size(coilimages,1), size(coilimages,2), nCh]));

end