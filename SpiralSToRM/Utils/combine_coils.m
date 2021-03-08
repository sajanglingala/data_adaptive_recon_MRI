function [vKSpaceData] = combine_coils(kSpaceData,threshold)

    %temp = permute(kSpaceData,[1,3,2]);
    %temp = reshape(temp,param.n*param.nf*param.lines_per_SG_block,param.NcoilSelect);
    %[N,~,NcoilSelect]=size(coilimages);
    temp=kSpaceData;
    Rs = real(temp'*temp);


    [v,s] = eig(Rs);
    s = diag(s);[s,i] = sort(s,'descend');
    v = v(:,i);
    s=s./sum(s);s = cumsum(s);
    nvchannels = min(find(s>threshold));

    vKSpaceData = temp*v(:,1:nvchannels);
    
    %vCoilImages = reshape(coilimages,N^2,NcoilSelect)*v(:,1:nvchannels);
    %vCoilImages = reshape(vCoilImages,N,N,nvchannels);
%     vKSpaceData = reshape(vKSpaceData,[param.n,param.lines_per_SG_block*param.nf,nvchannels]);
%     vKSpaceData = permute(vKSpaceData,[1,3,2]);