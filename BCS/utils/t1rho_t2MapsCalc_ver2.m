% funtion to calculate parameter maps
% mono exponential fitting model
% M(p)= S0 * exp(-tsl(p)/t1rho) * exp(-te(p)/t2)

% Author: Sampada Bhave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t1rhoMaps,t2Maps,S0Maps] = t1rho_t2MapsCalc_ver2(vol,TSLs,TEs)

[Nx, Ny, Ns, Nz] = size(vol);
numTSLs = size(TSLs,1);
numTEs=size(TEs,1);
sets = floor(Ns/(numTSLs+numTEs));
disp('T1rho map generation...');
disp(['--> imaging matrix = ' num2str(Nx) ' x ' num2str(Ny) ' x ' num2str(Nz)]);
disp(['--> TSL values = ' num2str(TSLs') ' ms']);
disp(['--> TSL sets = ' num2str(sets)]);
disp(['--> TE values = ' num2str(TEs') ' ms']);

tic;
toc0 = toc;

% calculate T1rho maps
t1rhoMaps = zeros(Nx,Ny,sets,Nz);  % one map for each set of TSL values
t2Maps = zeros(Nx,Ny,sets,Nz);  % one map for each set of TE values
S0Maps = zeros(Nx,Ny,sets,Nz); 
NxNy = Nx*Ny;
NxNyTSLs = NxNy*numTSLs;   %actually not needed... multiplied by 0 ahead..

stepInc=(0:99)'*NxNy;
te=[TEs,zeros(length(TEs),1)]; te=te(:);
tsl=[zeros(length(TSLs),1),TSLs];tsl=tsl(:);
A=[ones(length(te),1),tsl,te];

for nz = 1:Nz
    disp(['Analyzing slice ' num2str(nz) ' / ' num2str(Nz) '... (t = ' num2str(toc-toc0) ' sec)']);
    % work with smaller volumes for increased performance
    volBuffer = vol(:,:,:,nz);
    volMaskBuffer = squeeze(vol(:,:,1,nz));
    t1rhoMapsBuffer = zeros(Nx,Ny,sets);
    t2MapsBuffer = zeros(Nx,Ny,sets);
    S0Buffer=zeros(Nx,Ny,sets);
    nonZeroIndices = find(volMaskBuffer(:)>0);
    for ns = 0:sets-1  % index shifted by 1 for faster processing
        steps = stepInc + NxNyTSLs * ns;
        NxNyns = NxNy * ns;
        for ni = 1:size(nonZeroIndices,1)
            Svalues = volBuffer(nonZeroIndices(ni) + steps(1:numTSLs+numTEs));
            P=A\log(Svalues);
            t1rhoMapsBuffer(nonZeroIndices(ni) + NxNyns) = -1/P(2);
            t2MapsBuffer(nonZeroIndices(ni) + NxNyns)=-1/P(3);
            S0Buffer(nonZeroIndices(ni) + NxNyns)=exp(P(1));
        end
    end
S0Maps(:,:,:,nz)=S0Buffer;
t1rhoMaps(:,:,:,nz) = t1rhoMapsBuffer;
t2Maps(:,:,:,nz) = t2MapsBuffer;
end

% clean up T1rho maps
t1rhoMaps(~isfinite(t1rhoMaps)) = 0;  % set any NANs or INFs to zero
t1rhoMaps(t1rhoMaps<0) = 0;  % set any negative T1rho values to zero
t1rhoMaps(t1rhoMaps>(median(t1rhoMaps(t1rhoMaps(:)>0))*10)) = 0;  % set high noise values to zero

% clean up T2 maps
t2Maps(~isfinite(t2Maps)) = 0;  % set any NANs or INFs to zero
t2Maps(t2Maps<0) = 0;  % set any negative T1rho values to zero
t2Maps(t2Maps>(median(t2Maps(t2Maps(:)>0))*10)) = 0;  % set high noise values to zero

% clean up S0 maps
S0Maps(~isfinite(S0Maps)) = 0;  % set any NANs or INFs to zero
S0Maps(S0Maps<0) = 0;  % set any negative T1rho values to zero
S0Maps(S0Maps>(median(S0Maps(S0Maps(:)>0))*10)) = 0;  % set high noise values to zero
