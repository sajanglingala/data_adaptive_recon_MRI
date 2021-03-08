function f = fft2(x, mrows, ncols)
%FFT2 Two-dimensional discrete Fourier Transform.
%   FFT2(X) returns the two-dimensional Fourier transform of matrix X.
%   If X is a vector, the result will have the same orientation.
%
%   FFT2(X,MROWS,NCOLS) pads matrix X with zeros to size MROWS-by-NCOLS
%   before transforming.
%
%   Class support for input X: 
%      float: double, single
%
%   See also FFT, FFTN, FFTSHIFT, FFTW, IFFT, IFFT2, IFFTN.

%   Copyright 1984-2010 The MathWorks, Inc. 

if ismatrix(x)
    if nargin==1
        f = fftn(x);
    else
        f = fftn(x,[mrows ncols]);
    end
else
    if nargin==1
        f = fft(fft(x,[],2),[],1);
    else
        f = fft(fft(x,ncols,2),mrows,1);
    end
end   
