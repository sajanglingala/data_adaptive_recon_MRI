function [phi,phit] = def_phi_phit(opts)

  
        phi = @(U) ForwardD((U));warning off all;
        phit = @(DU) Dive((DU));
        
         function DU = ForwardD(U)
            U=double(U);
             U = reshape(U, opts.n1, opts.n2,opts.r);
            DU = cell(2,1);
            % Forward finite difference operator
            xshift = zeros(size(U),'double');
            xshift(:,1:end-1,:) = U(:,2:end,:);
            xshift(:,end,:) = U(:,1,:);
            DU{1} = xshift - U; clear xshift;
            
            
            yshift = zeros(size(U),'double'); 
            yshift(1:end-1,:,:) = U(2:end,:,:);
            yshift(end,:,:) = U(1,:,:);
            DU{2} = yshift - U; clear yshift;
           
         end
        
         function DtXY = Dive(DU)
            X = double(DU{1});
            Y = double(DU{2});
           
            warning off all;
            % Transpose of the forward finite difference operator
            xshift = zeros(size(X),'double');
            xshift(:,1:end-1,:) = X(:,2:end,:);
            xshift(:,end,:) = X(:,1,:);
            tempx = -(xshift - X); clear xshift;
            diffX = tempx(:,1:end-1,:);
            
            DtXY = zeros(size(X),'double');
            DtXY(:,1,:) = X(:,end,:) - X(:, 1,:);
            DtXY(:,2:end,:) = diffX;
            
            
            yshift = zeros(size(Y),'double'); 
            yshift(1:end-1,:,:) = Y(2:end,:,:);
            yshift(end,:,:) = Y(1,:,:);
            tempy = -(yshift - Y); clear yshift;
            diffY = tempy(1:end-1,:,:);
            
            uy = zeros(size(Y),'double');
            uy(1,:,:) = Y(end,:,:) - Y(1,:,:); 
            uy(2:end,:,:) = diffY;
            DtXY = DtXY + uy; 
            
            DtXY = reshape(DtXY, opts.n1*opts.n2, opts.r);
            
         end
end
        