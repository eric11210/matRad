%% calculate mean and variance (analytical and MC)

% setup options
options.numOfScenarios  = resultGUI.apertureInfo.numPhases;
options.bioOpt          = 'none';

resultGUI.apertureInfo = matRad_daoVec2ApertureInfo(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector);

% analytical
dMean       = matRad_backProjection(resultGUI.apertureInfo.bixelWeights,dij,options);
% prep cuts
dMax        = max(dMean);
voxelMask   = dMean == dMax;
% continue analytical
dVar        = matRad_doseVariance(resultGUI.apertureInfo,dij,voxelMask);

% Monte Carlo
nHistoriesPerTrial  = 1;
nTrials             = 10000;

dMean_MC    = zeros(nTrials,1);
dVar_MC     = zeros(nTrials,1);

for trial = 1:nTrials
    [dMean_trial,dVar_trial]  = matRad_doseVarianceMC(resultGUI.apertureInfo,dij,nHistoriesPerTrial);
    
    dMean_MC(trial) = dMean_trial(voxelMask);
    dVar_MC(trial) = dVar_trial(voxelMask);
end

%% analysis

% calculate delta, the difference of MC to analytical divided by the SD of
% the mean (which is equal to SD divided by sqrt(nHistories)
delta       = (dMean(voxelMask)-dMean_MC)./sqrt(dVar(voxelMask)./nHistoriesPerTrial);

% calculate chi square
chiSquare       = delta'*delta;

% calculate p values
p       = chi2cdf(chiSquare,numel(delta),'upper');

%% graphing deltas

figure
hold on
histogram(delta,'Normalization','pdf')
plot(-4:0.01:4,normpdf(-4:0.01:4,0,1))
xlabel('\Delta')
ylabel('probability density')
title(sprintf('Analytical variance \\chi^2 = %.0f, p = %.2f',chiSquare,p))
savefig('Analytical variance oneVoxel')

% save
save('analysis_oneVoxel','dMean*','dVar*','delta*','chiSquare*','p*','voxelMask')