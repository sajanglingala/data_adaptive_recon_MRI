function [ img ] = invHankel( H,sx,sy,nc,ksize )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%     H1=transpose(H);
    tsx=size(H,1);
    tsy=size(H,2)/nc;
    H2 = reshape(H,tsx,tsy,nc);
% tsx=size(H,2);
% tsy=size(H,1)/nc;
% H = reshape(permute(H,[2,1]),tsx,tsy,nc);

img = row2im_tpose(H2,[sx,sy,nc],ksize);
end

