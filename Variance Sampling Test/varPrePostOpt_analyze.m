%% analysis and plotting

% load variances, PDVHs
load('varPrePostOpt_oneFrac.mat')

% calculate standard deviations
dStd_preOpt         = sqrt(dVar_preOpt);
dStd_postProbOpt    = sqrt(dVar_postProbOpt);
dStd_postConvOpt    = sqrt(dVar_postConvOpt);

% calculate relative standard deviations
dRelStd_preOpt  = dStd_preOpt./d_preOpt;
dRelStd_postProbOpt = dStd_postProbOpt./d_postProbOpt;
dRelStd_postConvOpt = dStd_postConvOpt./d_postConvOpt;

% plot variances
figure
hold on
plot(dVar_preOpt(dij.targetVox))
plot(dVar_postProbOpt(dij.targetVox))
plot(dVar_postConvOpt(dij.targetVox))
xlabel('target voxel index')
ylabel('variance')
legend({'preOpt' 'postProbOpt' 'postConvOpt'},'location','best')

% plot standard deviations
figure
hold on
plot(dStd_preOpt(dij.targetVox))
plot(dStd_postProbOpt(dij.targetVox))
plot(dStd_postConvOpt(dij.targetVox))
xlabel('target voxel index')
ylabel('standard deviation')
legend({'preOpt' 'postProbOpt' 'postConvOpt'},'location','best')

% plot relative standard deviations
figure
hold on
plot(dRelStd_preOpt(dij.targetVox))
plot(dRelStd_postProbOpt(dij.targetVox))
plot(dRelStd_postConvOpt(dij.targetVox))
xlabel('target voxel index')
ylabel('relative standard deviation')
legend({'preOpt' 'postProbOpt' 'postConvOpt'},'location','best')

% plot ratio of variances
figure
hold on
plot(dVar_postProbOpt(dij.targetVox)./dVar_preOpt(dij.targetVox))
xlabel('target voxel index')
ylabel('ratio of variance')

% plot ratio of standard deviations
figure
hold on
plot(dStd_postProbOpt(dij.targetVox)./dStd_preOpt(dij.targetVox))
xlabel('target voxel index')
ylabel('ratio of standard deviations')

% plot ratio of relative standard deviations
figure
hold on
plot(dRelStd_postProbOpt(dij.targetVox)./dRelStd_preOpt(dij.targetVox))
xlabel('target voxel index')
ylabel('ratio of relative standard deviations')