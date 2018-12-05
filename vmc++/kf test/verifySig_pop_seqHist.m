%% setup

Nx = 101;
Ny = 101;
Nz = 150;

numSim = 50;

dose_nMillion            = zeros(Nx*Ny*Nz,numSim);

dose_nthMillion         = zeros(Nx*Ny*Nz,numSim);
doseError_1stMillion    = zeros(Nx*Ny*Nz,1);

%% read in all simulations

for i = 1:numSim
    
    fname = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\',sprintf('kfMethod_getSig_%dmillion.dos',i));
    
    fid = fopen(fname,'r');
    
    Header      = fread(fid,1,'int32');
    no_regions  = Header(1);
    
    tempDose        = fread(fid, no_regions, 'float32');
    tempDoseError   = fread(fid, no_regions, 'float32');
    
    tempDose = reshape(tempDose,[Nx Ny Nz]);
    tempDose = permute(tempDose,[2 1 3]);
    
    tempDoseError = reshape(tempDoseError,[Nx Ny Nz]);
    tempDoseError = permute(tempDoseError,[2 1 3]);
    
    fclose(fid);
    
    if i == 1
        doseError_1stMillion    = tempDoseError(:);
    end
    
    dose_nMillion(:,i) = tempDose(:);
    
end

%% calculate the dose for the n millionth histories

for i = 1:numSim
    
    if i == 1
        
        % dose in the first 1 million is just the dose for 1 million
        % particles
        dose_nthMillion(:,i) = dose_nMillion(:,i);
        
    else
        
        % first determine the total dose delivered by all n million and
        % (n-1) million histories
        totalDose_nMillion      = dose_nMillion(:,i).*i.*10^6;
        totalDose_nm1Million    = dose_nMillion(:,i-1).*(i-1).*10^6;
        
        % now subtract the two to get the total dose delivered by the nth
        % million histories
        totalDose_nthMillion = totalDose_nMillion-totalDose_nm1Million;
        
        % now divide by a million histories to get the average dose of the
        % nth million
        dose_nthMillion(:,i) = totalDose_nthMillion./(10^6);
        
    end
    
end

%% calculate means, stds

dose_nthMillion_mean   = mean(dose_nthMillion,2);
dose_nthMillion_std    = std(dose_nthMillion,0,2);

%% find > 50%

maxDose = max(dose_nthMillion_mean);

deleteInd = dose_nthMillion_mean < 0.5*maxDose;

dose_nthMillion_mean(deleteInd)     = [];
dose_nthMillion_std(deleteInd)      = [];
doseError_1stMillion(deleteInd)     = [];

ratioStd        = doseError_1stMillion./dose_nthMillion_std;
ratioStd_mean   = mean(ratioStd);

%% make plots

figure
hold on
plot(dose_nthMillion_mean,doseError_1stMillion,'k.','MarkerSize',1)
plot(dose_nthMillion_mean,dose_nthMillion_std,'r.','MarkerSize',1)
xlabel('dose')
ylabel('\sigma')
legend({'vmc++' sprintf('estimated from %d runs',numSim)},'location','best')
title('standard deviation in D > 50%')

figure
hold on
plot(dose_nthMillion_mean,ratioStd,'k.')
xlabel('dose')
ylabel('\sigma_{vmc} / \sigma_{est}')
title(sprintf('ratio of standard deviations in D > 50%%, mean = %f',ratioStd_mean))
