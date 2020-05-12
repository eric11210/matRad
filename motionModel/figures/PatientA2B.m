%% get model for patient A

% file
options_A.data.origin             = 'file';
options_A.data.fileInfo.f         = 19;
options_A.data.fileInfo.m         = 1;
options_A.data.fileInfo.processed = true;
options_A.data.fileInfo.t_tot     = 600; % seconds

% processing options

% for binning
options_A.processing.percExtTarg  = 1;
options_A.processing.nPosPhases   = 5;
options_A.processing.nVelPhases   = 1;    % is this necessary?
options_A.processing.velBinning   = true; % is this necessary? just set nVelSubPhases to 1

% for resampling
options_A.processing.fResample = 3; % Hz

% for FSM
options_A.processing.doFSM                = true;
options_A.processing.FSM.timeLd           = 0.5; % s
options_A.processing.FSM.velocityRangeInh = [-inf -0.5]; % mm/s
options_A.processing.FSM.velocityRangeExh = [0.5 inf]; % mm/s 0.3
options_A.processing.FSM.velocityRangeEOE = [0 inf]; % mm/s
%options.processing.FSM.cTheta           = 4.5; % mm/s
%options.processing.FSM.cLambda          = 4; % mm
options_A.processing.FSM.cSLength         = 0.133; % s

% training ratio
options_A.processing.trainRatio = 0.5;

% FFT options

options_A.FFT.doFFT       = true;
options_A.FFT.doWindowing = true;

% histogram options

options_A.hist.doHist     = false;
options_A.hist.timePoints = 0:0.04:90; %s

% convergence time options

options_A.convTime.doConvTime = true;
options_A.convTime.l2_targ    = 0.01;

% get model

% best fit parameters
options_A.processing.nSubPerPosPhase  = 1;
options_A.processing.nSubPerVelPhase  = 2;
options_A.processing.FSM.nTimeFracs   = 6;

% extract model
model_A = matRad_generateMotionModel(options_A);

% absolute path to dumping folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

%% get data for patient B

% file
options_B.data.origin             = 'file';
options_B.data.fileInfo.f         = 21;
options_B.data.fileInfo.m         = 3;
options_B.data.fileInfo.processed = true;
options_B.data.fileInfo.t_tot     = 600; % seconds

% processing options

% for binning
options_B.processing.percExtTarg  = 1;
options_B.processing.nPosPhases   = 5;
options_B.processing.nVelPhases   = 1;    % is this necessary?
options_B.processing.velBinning   = true; % is this necessary? just set nVelSubPhases to 1

% for resampling
options_B.processing.fResample = 6; % Hz

% for FSM
options_B.processing.doFSM                = true;
options_B.processing.FSM.timeLd           = 0.5; % s
options_B.processing.FSM.velocityRangeInh = [-inf -0.5]; % mm/s
options_B.processing.FSM.velocityRangeExh = [0.5 inf]; % mm/s 0.3
options_B.processing.FSM.velocityRangeEOE = [0 inf]; % mm/s
%options.processing.FSM.cTheta           = 4.5; % mm/s
%options.processing.FSM.cLambda          = 4; % mm
options_B.processing.FSM.cSLength         = 0.133; % s

% training ratio
options_B.processing.trainRatio = 0.5;

% FFT options

options_B.FFT.doFFT       = true;
options_B.FFT.doWindowing = true;

% histogram options

options_B.hist.doHist     = true;
options_B.hist.timePoints = 0:0.04:90; %s

% convergence time options

options_B.convTime.doConvTime = true;
options_B.convTime.l2_targ    = 0.01;

% get model

% best fit parameters
options_B.processing.nSubPerPosPhase  = 10;
options_B.processing.nSubPerVelPhase  = 10;
options_B.processing.FSM.nTimeFracs   = 10;

% get data

data_B = matRad_readMotionData(options_B.data.fileInfo,options_B.processing.percExtTarg);

[data_B_train,data_B_test,parametersNoGood] = matRad_processMotionData(data_B,options_B.processing,options_B.data.fileInfo);

%% do histogram stuff

model_A = matRad_calcPosHist(model_A,data_B_test,options_A.hist);

%% ChiSquaredSigmaDistribution

%display it
figure
hold on
histogram(sum(model_A.chiSquaresMC,2),'Normalization','pdf')
yl = ylim;
plot(repmat(sum(model_A.chiSquares),1,2),yl,'k')
xlabel('\chi^2_{help}')
ylabel('probability density')
fname = 'ChiSquaredSigmaDistribution_A2B_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

fprintf('The chi squared value is %f and the p-value is %f.',sum(model_A.chiSquares),model_A.pSum);

%% chiSquaredDeltaT, pDeltaT

%display chiSquaredDeltaT
figure
plot(options_B.hist.timePoints,model_A.chiSquares)
xlabel('\Delta t / s')
ylabel('\chi^2_{help}')
fname = 'ChiSquaredDeltaT_A2B_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%display pDeltaT
figure
plot(options_B.hist.timePoints,model_A.p)
xlabel('\Delta t / s')
ylabel('p_{help}')
fname = 'PDeltaT_A2B_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%% PhaseHistogramsDeltaT

% find dynamic range
minVal = min([model_A.hist_obs(:); model_A.hist_pred(:)]);
maxVal = max([model_A.hist_obs(:); model_A.hist_pred(:)]);

maxVal = 3000;

% display observed
figure
imagesc([options_B.hist.timePoints(1) options_B.hist.timePoints(end)],[1 options_B.processing.nPosPhases],model_A.hist_obs)
caxis([minVal maxVal])
xlabel('\Delta t / s')
ylabel('position bin $l$')
yticks(1:options_B.processing.nPosPhases)
colorbar
fname = 'ObservedPhaseHistogramsDeltaT_A2B_raw';
set(gca,'XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

% display expected
figure
imagesc([options_B.hist.timePoints(1) options_B.hist.timePoints(end)],[1 options_B.processing.nPosPhases],model_A.hist_pred)
caxis([minVal maxVal])
xlabel('\Delta t / s')
ylabel('position bin $l$')
yticks(1:options_B.processing.nPosPhases)
colorbar
fname = 'ExpectedPhaseHistogramsDeltaT_A2B_raw';
set(gca,'XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%% SampleAndObservedSequences PREP

% stop in calcPosHist, right after generating a sample trajectory
model_A = matRad_generateMotionModel(options);

%% SampleAndObservedSequences

% absolute path to dumping folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

% sample sequence
t_simulated = (0:(numStepsModel_simulate-1)).*model_A.deltaT_sample;
l_simulated = matRad_runMarkovChain_P(model_A,numStepsModel_simulate,false);
x_simulated = model_A.indices.posSubPhase2Pos(model_A.indices.subPhase2PosSubPhase(l_simulated));

% display observed
figure
plot(data.t_sample-data.t_sample(1),data.x_sample,'k')
xlim([0 20])
xlabel('time / s')
ylabel('principal component of motion / mm')
fname = 'ObservedSequence_A2B_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

% display sampled
figure
plot(t_simulated,x_simulated,'k')
xlim([0 20])
xlabel('time / s')
ylabel('principal component of motion / mm')
fname = 'SampleSequence_A2B_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

