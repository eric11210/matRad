%% calculate mean and variance (analytical and MC)

% setup options
options.numOfScenarios  = resultGUI.apertureInfo.numPhases;
options.bioOpt          = 'none';

% analytical
dMean       = matRad_backProjection(resultGUI.apertureInfo.bixelWeights,dij,options);
dVar        = matRad_doseVariance(resultGUI.apertureInfo,dij);
dVar_summed = sum(dVar(dij.targetVox));
[dVarSum,~] = matRad_doseVarianceSum(resultGUI.apertureInfo,dij);

% Monte Carlo
nHistories          = 10000;
[dMean_MC,dVar_MC]  = matRad_doseVarianceMC(resultGUI.apertureInfo,dij,nHistories);

%% analysis

% prep 50% cuts
dMax        = max(dMean(dij.structVox));
deleteInd   = dMean(dij.structVox) < 0.01*dMax;

% calculate delta, the difference of MC to analytical divided by the SD of
% the mean (which is equal to SD divided by sqrt(nHistories)
delta       = (dMean(dij.structVox)-dMean_MC(dij.structVox))./sqrt(dVar./nHistories);
delta_varMC = (dMean(dij.structVox)-dMean_MC(dij.structVox))./sqrt(dVar_MC(dij.structVox)./nHistories);

% perform 50% cuts
delta(deleteInd)        = [];
delta_varMC(deleteInd)  = [];

% calculate chi square
chiSquare       = delta'*delta;
chiSquare_varMC = delta_varMC'*delta_varMC;

% calculate p values
p       = chi2cdf(chiSquare,numel(delta),'upper');
p_varMC = chi2cdf(chiSquare_varMC,numel(delta),'upper');

%% graphing deltas

figure
hold on
histogram(delta,'Normalization','pdf')
plot(-4:0.01:4,normpdf(-4:0.01:4,0,1))
xlabel('\Delta')
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

% save
save('analysis','dMean*','dVar*','delta*','chiSquare*','p*')

%% graphing difference in analytical vs MC

figure
hold on
plot(dVar)
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
plot((dVar_MC(dij.structVox)-dVar)./dVar)
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