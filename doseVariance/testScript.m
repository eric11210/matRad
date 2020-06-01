%% calculate mean and variance (analytical and MC)

% lungPatient0

load('postSequencing.mat')

% setup options
options.numOfScenarios  = resultGUI.apertureInfo.numPhases;
options.bioOpt          = 'none';

resultGUI.apertureInfo = matRad_daoVec2ApertureInfo(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector);

% analytical
dMean       = matRad_backProjection(resultGUI.apertureInfo.bixelWeights,dij,options);
% prep cuts
dMax        = max(dMean);
voxelMask   = dMean > 0.5*dMax;
% continue analytical
dVar        = matRad_doseVariance(resultGUI.apertureInfo,dij,voxelMask);
dVar_summed = sum(dVar(dij.targetVox));
[dVarSum,~] = matRad_doseVarianceSum(resultGUI.apertureInfo,dij);

% Monte Carlo
nHistories          = 50000;
[dMean_MC,dVar_MC]  = matRad_doseVarianceMC(resultGUI.apertureInfo,dij,nHistories);

%% analysis

% calculate delta, the difference of MC to analytical divided by the SD of
% the mean (which is equal to SD divided by sqrt(nHistories)
delta       = (dMean(voxelMask)-dMean_MC(voxelMask))./sqrt(dVar(voxelMask)./nHistories);
delta_varMC = (dMean(voxelMask)-dMean_MC(voxelMask))./sqrt(dVar_MC(voxelMask)./nHistories);

% calculate chi square
chiSquare       = delta'*delta;
chiSquare_varMC = delta_varMC'*delta_varMC;

% calculate p values
p       = chi2cdf(chiSquare,numel(delta),'upper');
p_varMC = chi2cdf(chiSquare_varMC,numel(delta_varMC),'upper');

% save
save('analysis','dMean*','dVar*','delta*','chiSquare*','p*','voxelMask')

%% graphing deltas

figure
hold on
histogram(delta,'Normalization','pdf')
plot(-4.575:0.01:4.575,normpdf(-4.575:0.01:4.575,0,1))
xlabel('\Delta')
xlim([-4.575 4.575])
ylabel('probability density')
title(sprintf('Analytical variance \\chi^2 = %.0f, p = %.2f',chiSquare,p))
savefig('Analytical variance')

figure
hold on
histogram(delta_varMC,'Normalization','pdf')
plot(-4:0.01:4,normpdf(-4:0.01:4,0,1))
xlabel('\Delta')
ylabel('probability density')
title(sprintf('Monte Carlo variance, \\chi^2 = %.0f p = %.2f',chiSquare_varMC,p_varMC))
savefig('Monte Carlo variance')

%% graphing difference in analytical vs MC

figure
hold on
plot(dVar(dij.structVox))
plot(dVar_MC(dij.structVox))
xlabel('voxel index')
ylabel('\sigma^2')
legend({'analytical' 'Monte Carlo'},'Location','Best')
title('Analytical and Monte Carlo variance')
savefig('Both variances')

figure
hold on
plot(dMean(dij.structVox))
plot(dMean_MC(dij.structVox))
xlabel('voxel index')
ylabel('\mu')
legend({'analytical' 'Monte Carlo'},'Location','Best')
title('Analytical and Monte Carlo mean')
savefig('Both means')

figure
plot((dVar_MC(dij.structVox)-dVar(dij.structVox))./dVar(dij.structVox))
xlabel('voxel index')
ylabel('(\sigma^2_{MC}-\sigma^2_{analytical})/\sigma^2_{analytical}')
title('Comparing analytical to Monte Carlo variance')
savefig('Variance comparison')

figure
plot((dMean_MC(dij.structVox)-dMean(dij.structVox))./dMean(dij.structVox))
xlabel('voxel index')
ylabel('(\mu_{MC}-\mu_{analytical})/\mu_{analytical}')
title('Comparing analytical to Monte Carlo mean')
savefig('Mean comparison')