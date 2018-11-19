%% setup

Nx = 101;
Ny = 101;
Nz = 150;

numSim = 20;

simpleSource        = zeros(Nx*Ny*Nz,numSim);
simpleSourceError   = zeros(Nx*Ny*Nz,numSim);

%% read in all simulations

for i = 1:numSim
    
    fname = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\',sprintf('kfMethod_getSig_simpleSource_seed%d.dos',i));
    
    fid = fopen(fname,'r');
    
    Header      = fread(fid,1,'int32');
    no_regions  = Header(1);
    
    temp       = fread(fid, no_regions, 'float32');
    tempError  = fread(fid, no_regions, 'float32');
    
    temp = reshape(temp,[Nx Ny Nz]);
    temp = permute(temp,[2 1 3]);
    
    tempError = reshape(tempError,[Nx Ny Nz]);
    tempError = permute(tempError,[2 1 3]);
    
    fclose(fid);
    
    simpleSource(:,i)       = temp(:);
    simpleSourceError(:,i)  = tempError(:);
    
end

%% calculate means, stds

simpleSource_mean   = mean(simpleSource,2);
simpleSource_std    = std(simpleSource,0,2);

simpleSourceError_mean  = mean(simpleSourceError,2);
simpleSourceError_std   = std(simpleSourceError,0,2);

%% find > 50%

maxDose = max(simpleSource_mean);

deleteInd = simpleSource_mean < 0.5*maxDose;

simpleSource_mean(deleteInd)        = [];
simpleSource_std(deleteInd)         = [];
simpleSourceError_mean(deleteInd)   = [];
simpleSourceError_std(deleteInd)    = [];

ratioStd        = simpleSourceError_mean./simpleSource_std;
ratioStd_mean   = mean(ratioStd);

%% make plots

figure
hold on
plot(simpleSource_mean,simpleSourceError_mean,'k.')
plot(simpleSource_mean,simpleSource_std,'r.')
xlabel('dose')
ylabel('\sigma')
legend({'vmc++' sprintf('estimated from %d runs',numSim)},'location','best')
title('standard deviation in D > 50%')

figure
hold on
plot(simpleSource_mean,ratioStd,'k.')
xlabel('dose')
ylabel('\sigma_{vmc} / \sigma_{est}')
title(sprintf('ratio of standard deviations in D > 50%%, mean = %f',ratioStd_mean))
