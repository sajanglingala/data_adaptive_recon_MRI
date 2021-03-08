function out = SenseATA(x,FT,csm,N,nFrames,nCh)

  x = reshape(x,[N,N,nFrames]);
  out = zeros(N,N,nFrames);
  for ii=1:nCh,
      xnew = FT'*(FT*bsxfun(@times,x,csm(:,:,ii)));
      out = out + bsxfun(@times,xnew,conj(csm(:,:,ii)));
  end
  out = out(:);
    
end