clear ; 
load /Users/slingala/Downloads/phantoms/vtPhantom/x_im.mat; 
x_im=abs(x_im); 
x_reg=x_im;

for i = 2:140
[ x_reg(:,:,i) ] = demon_reg( x_reg(:,:,i-1),x_reg(:,:,i) );
end