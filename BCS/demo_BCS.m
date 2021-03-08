% Main file to execute the Blind compressed sensing (BCS) scheme
% BCS optimization problem ||A(UV)-b||^2 + lambda*||U||_l1 s.t. ||V|^2 < r
% V - dictionary and U are the coefficients, r - number of dictionary atoms

% Author : Sampada Bhave
% Refer to "Accelerated whole brain multi-parameter mapping using
% blind compresssed sensing", MRM, 2016 for furthur details  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

load data/brain_T2_T1rho.mat

addpath('utils/');

data_orig=squeeze(rssq(img1,3));  % fully sampled coil combined data

[nx,ny,nc,nd] = size(img1);       % dimensions of the image
csm=repmat(csm,[1 1 1 nd]);       % coil sensitivities
m=nx*ny;
n=nd;

%%% sampling pattern 
OMEGA = find(samp);

%%% A and At operators are the forward and backward Fourier sampling operators respectively
A = @(x)A_fhp_2d(x, csm, OMEGA,[nx,ny,nd],nc);
At = @(y)At_fhp_2d(y, csm, OMEGA,[nx,ny,nc,nd]);

%%% undersampled raw data
y=zeros(size(OMEGA,1),nc);
for i=1:nc
    junk = fft2(img1(:,:,i,:));
    y(:,i) = junk(OMEGA);
end
b=y;

%%% first guess
x_init = At(y);
x_init1=reshape((x_init),nx,ny,nd);
RMSE((x_init1),(data_orig))

% r denotes the number of temporal basis functions in the dictionary
r =40;

%%% The BCS algorithm 
% parameters are all specified in opts.

opts.phi = 'l1' ;   % choices of  'TV' - TV on U
                    %             'l1' - l1 norm on U
                    %             'lp' - lp norm on U 
opts.con='Fro';     % choice of contraint 'Fro' - Frobenius norm on V
                    %                     'Col' - Column norm on V

                
lambda1=0.05;    % regularization parameter
opts.outer =200; % The iterations of the outer loop
opts.inner =70;  % The iterations of the inner loop
opts.lambda1=lambda1;
opts.r=r;
opts.n1=nx;
opts.n2=ny;
opts.n3 = nd;
opts.beta1 =0.01;
opts.beta2 =10;
opts.p=1; 
earrayout=[];rmse_un=[];
earrayin=[];

if strcmp(opts.phi, 'TV')
    [phi,phit] = def_phi_phit(opts);
    RU = @(z)sp_TV_eval(z);
    
elseif strcmp(opts.phi, 'l1')
    [phi,phit] = def_phi_phit_ID;
    RU = @(z)(sum(abs(z(:))));
    
elseif strcmp(opts.phi, 'lp')
    [phi,phit] = def_phi_phit_ID;
    RU = @(z)(sum(abs(z(:)).^opts.p))^(1/opts.p);
end

U=zeros(m,r,'double');

% random dictionary initialization
V = (rand(r,n));
Lam2=zeros(size(V));

% Initial U sub problem - CG
if strcmp(opts.phi, 'TV'), L{1}=zeros(size(U),'double'); L{2}=zeros(size(U),'double'); else  L=0; end;
gradU_left = @(z)Eval_gradUleft(z,A,At,phi, phit, V,opts);
gradU_right = At(b)*V' +(opts.lambda1*opts.beta1/2)*phit(L);%gradU_right = gradU_right(:);
x0= zeros(size(U)); x0=x0(:);
[U,~,~,~,resvec] = pcg(gradU_left,gradU_right(:),1e-5,30,[],[],x0(:));
U=reshape(U,opts.n1*opts.n2, opts.r);

opts.beta1=1/max(abs(U(:)));earray_out=[];

reconerr=[];
for out =1:opts.outer
    out
    tic
    opts.beta2=1/max(abs(V(:)));
    for in = 1:opts.inner
        
        
        % L sub problem - Shrinkage
        L = Lsubproblem(opts, phi, U);
        
        
        % U sub problem - CG
        gradU_left = @(z)Eval_gradUleft(z,A,At,phi, phit, V,opts);
        gradU_right = (At(b)*V' +(opts.lambda1*opts.beta1/2)*phit(L));gradU_right = gradU_right(:);
        x0=(U(:));
        [U,~,~,~,earray_u] = pcg(gradU_left,gradU_right,1e-4,30,[],[],x0(:));
        U=reshape(U,opts.n1*opts.n2, opts.r);
        
        % update Q
        KK=V+1/opts.beta2*Lam2;
        if strcmp(opts.con,'Col')
            C=1;
            Q=KK*diag(min(sum(KK.^2).^(-.5),1));
        elseif strcmp(opts.con,'Fro')
            C=r;
            Kfn = sum(sum(KK.^2));
            Q = min(1,sqrt(C/Kfn))*KK;
        end
        
        % V sub problem - CG
        gradV_left=@(z)Eval_Vleft(z,A,At,U,opts);
        gradV_right= 2*(U'*At(b))-Lam2+opts.beta2*Q;
        if out ==1 && in ==1
            vx0=zeros(size(V)); vx0=vx0(:);
        else
            vx0 = V(:);
        end
        [V,~,~,~,earray_v] = pcg(gradV_left,gradV_right(:),1e-4,30,[],[],vx0(:));
        V = reshape(V, opts.r, opts.n3);
        
        % COST Calculations
        dc = A(U*V)-b; % data consistency
        reg_U = RU(phi(U)); % l1 norm on the spatial weights
        
        cost = sum(abs(dc(:)).^2) +opts.lambda1*reg_U;
        earrayin = [earrayin,cost];
        con2=abs(sum(sum((V-Q).^2)));
        
        Lam2 = Lam2 -opts.beta2*(Q - V);
        
        if in>1
            if (abs(earrayin(end)-earrayin(end-1))/abs(earrayin(end))) < 1e-2
                opts.beta2=opts.beta2*10;
            end
            
            if con2 < 1e-4
                in
                break;
            end
            
        end
    end
    dc = A(U*V)-b; % data consistency
    reg_U = RU(phi(U)); % l1 norm on the spatial weights
    
    cost = sum(abs(dc(:)).^2) +opts.lambda1*reg_U;
    earrayout = [earrayout,cost];
    
    % plot the results during while the algorithm is converging
    X_est = U*V;
    X_est = reshape(X_est,nx,ny,nd);
    figure(10);
    subplot(2,3,1); imagesc(abs((U))); title('Example spatial weight')
    subplot(2,3,2); imagesc(abs((V))); title('Dictionary V')
    subplot(2,3,4); plot(abs((V(3,:)))); title('Example temporal basis: 3');
    subplot(2,3,3); imagesc(abs((X_est(:,:,2)))); title('A spatial frame: reconstruction');%reshape(L(:,1),n1,n2,1))); title('L');
    subplot(2,3,6); plot((earrayout)); title('Decreasing cost');
    subplot(2,3,5); imagesc(abs(squeeze((X_est(100,:,:)))));
    pause(0.01);
    
    
    if out>1
        if (abs(earrayout(end)-earrayout(end-1))/abs(earrayout(end))) < 1e-3;
            opts.beta1=opts.beta1*40;
        end
        if (abs(earrayout(end)-earrayout(end-1))/abs(earrayout(end))) < 1e-4;  % stop rule for the outer loop
            out
            break;
        end
    end
    err=RMSE(gather(X_est),data_orig);
    reconerr=[reconerr,err];
    TIME(out)=toc;
end
X_corr=X_est;
[MSE,X_corr] = RMSE(X_corr, data_orig);
MSE

%%% parameter map estimation 
% S0Maps, t1rhoMaps and t2Maps are maps estimated from fullysampled data

% TE and TSL parameters for parameter map estimation
TEs=10:10:120;TEs=TEs';
TSLs=10:10:120;TSLs=TSLs';

[t1rho_recon,t2_recon,S0_recon]=t1rho_t2MapsCalc_ver2(X_corr.*repmat(mapmask,[1,1,nd]),TSLs,TEs);
[MSE_t1rho, t1rho_recon]=RMSE(t1rho_recon.*mask,t1rhoMaps.*mask);
[MSE_t2, t2_recon]=RMSE(t2_recon.*mask,t2Maps.*mask);
[MSE_S0, S0_recon]=RMSE(S0_recon.*mask,S0Maps.*mask);


figure,colormap(jet);
subplot(2,2,1),imagesc(imrotate(t1rhoMaps.*mask,270), [0 100]); axis image; axis off;  title('T1_rho map original');
subplot(2,2,2),imagesc(imrotate(t1rho_recon,270), [0 100]); axis image; axis off;  title('T1_rho map estimated');
subplot(2,2,3),imagesc(imrotate(t2Maps.*mask,270), [0 120]); axis image; axis off;  title('T2 map original');
subplot(2,2,4),imagesc(imrotate(t2_recon,270), [0 120]); axis image; axis off;  title('T2 map estimated');

