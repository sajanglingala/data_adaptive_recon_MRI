function ress = mtimes(a,bb)
% performs the normal nufft
nf=size(bb,3);
for n=1:nf
b = bb(:,:,n);
if a.adjoint
	b = b(:).*a.w(:);
    if size(a.st)==1
	res = nufft_adj(b, a.st{1})/sqrt(prod(a.imSize));
    else
        res = nufft_adj(b, a.st{n})/sqrt(prod(a.imSize));
    end
	res = reshape(res, a.imSize(1), a.imSize(2));
% 	res = res.*conj(a.phase);
% 	if a.mode==1
% 		res = real(res);
% 	end

else
	b = reshape(b,a.imSize(1),a.imSize(2));
% 	if a.mode==1
% 		b = real(b);
% 	end
% 	b = b.*a.phase;
if size(a.st)==1
	res = nufft(b, a.st{1})/sqrt(prod(a.imSize));
else
    res = nufft(b, a.st{n})/sqrt(prod(a.imSize));
end
	res = reshape(res,a.dataSize(1),a.dataSize(2));
end
ress(:,:,n) = res;
end

