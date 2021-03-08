function [D,Dt] = defDDt(step_size)
        % defines finite difference operator D
        % and its transpose operator
        
        D = @(U) ForwardD((U));warning off all;
        Dt = @(DU) Dive((DU));
        
         function DU = ForwardD(U)
            DU = cell(3,1);
            % Forward finite difference operator
            xshift = zeros(size(U));
            xshift(:,1:end-1,:) = U(:,2:end,:);
            xshift(:,end,:) = U(:,1,:);
            DU{1} = xshift - U; clear xshift;
            DU{1} = DU{1}./step_size(1); 
            
            yshift = zeros(size(U)); 
            yshift(1:end-1,:,:) = U(2:end,:,:);
            yshift(end,:,:) = U(1,:,:);
            DU{2} = yshift - U; clear yshift;
            DU{2} = DU{2}./step_size(2); 
            
            
            tshift = zeros(size(U));
            tshift(:,:,1:end-1) = U(:,:,2:end);
            tshift(:,:,end) = U(:,:,1);
            DU{3} = tshift - U; clear tshift;
            DU{3} = DU{3}./step_size(3); 
            
         end
        
         function DtXY = Dive(DU)
            X = double(DU{1});
            Y = double(DU{2});
            t = double(DU{3});
            
            warning off all;
            % Transpose of the forward finite difference operator
            xshift = zeros(size(X));
            xshift(:,1:end-1,:) = X(:,2:end,:);
            xshift(:,end,:) = X(:,1,:);
            tempx = -(xshift - X); clear xshift;
            diffX = tempx(:,1:end-1,:);
            
            DtXY = zeros(size(X));
            DtXY(:,1,:) = X(:,end,:) - X(:, 1,:);
            DtXY(:,2:end,:) = diffX;
            DtXY = DtXY./step_size(1); 
            
            yshift = zeros(size(Y)); 
            yshift(1:end-1,:,:) = Y(2:end,:,:);
            yshift(end,:,:) = Y(1,:,:);
            tempy = -(yshift - Y); clear yshift;
            diffY = tempy(1:end-1,:,:);
            
            uy = zeros(size(Y));
            uy(1,:,:) = Y(end,:,:) - Y(1,:,:); 
            uy(2:end,:,:) = diffY;
            DtXY = DtXY + (uy)./step_size(2); 
            
            tshift = zeros(size(t));
            tshift(:,:,1:end-1) = t(:,:,2:end);
            tshift(:,:,end) = t(:,:,1);
            tempt = -(tshift - t); clear tshift;
            difft = tempt(:,:,1:end-1);
            
            
            
            temp = zeros(size(t));
            temp(:,:,2:end) = difft;
             temp(:,:,1) = t(:,:,end) - t(:,:,1);
             DtXY = DtXY + (temp)./step_size(3);
        end
        
    end