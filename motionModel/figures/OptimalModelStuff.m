%% model options

% file
options.data.origin             = 'file';
options.data.fileInfo.f         = 19;
options.data.fileInfo.m         = 1;
options.data.fileInfo.processed = true;
options.data.fileInfo.t_tot     = 600; % seconds

% processing options

% for binning
options.processing.percExtTarg  = 1;
options.processing.nPosPhases   = 5;
options.processing.nVelPhases   = 1;    % is this necessary?
options.processing.velBinning   = true; % is this necessary? just set nVelSubPhases to 1

% for resampling
options.processing.fResample = 3; % Hz

% for FSM
options.processing.doFSM                = true;
options.processing.FSM.timeLd           = 0.5; % s
options.processing.FSM.velocityRangeInh = [-inf -0.5]; % mm/s
options.processing.FSM.velocityRangeExh = [0.5 inf]; % mm/s 0.3
options.processing.FSM.velocityRangeEOE = [0 inf]; % mm/s
%options.processing.FSM.cTheta           = 4.5; % mm/s
%options.processing.FSM.cLambda          = 4; % mm
options.processing.FSM.cSLength         = 0.133; % s

% training ratio
options.processing.trainRatio = 0.5;

% FFT options

options.FFT.doFFT       = true;
options.FFT.doWindowing = true;

% histogram options

options.hist.doHist     = true;
options.hist.timePoints = 0:0.04:80; %s

% convergence time options

options.convTime.doConvTime = true;
options.convTime.l2_targ    = 0.01;

% get model

% best fit parameters
options.processing.nSubPerPosPhase  = 1;
options.processing.nSubPerVelPhase  = 2;
options.processing.FSM.nTimeFracs   = 6;

% extract model
model = matRad_generateMotionModel(options);

% absolute path to dumping folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

% find index of deltaT at 30 s
ind = find(options.hist.timePoints == 30);

%% ObservedandExpectedPositionHistograms

% display it
figure
bar([model.hist_obs(1:options.processing.nPosPhases,ind) model.hist_pred(1:options.processing.nPosPhases,ind)])
xlabel('position bin')
ylabel('frequency')
legend({'observed' 'expected'},'Location','Best')
fname = 'ObservedandExpectedPositionHistograms_raw';
grid on
set(gca,'YMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%% ChiSquaredDistribution

%display it
figure
hold on
histogram(model.chiSquaresMC(:,ind),'Normalization','pdf')
yl = ylim;
plot(repmat(model.chiSquares(ind),1,2),yl,'k')
xlabel('\chi^2_{help}')
ylabel('probability density')
fname = 'ChiSquaredDistribution_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

fprintf('The chi squared value is %f and the p-value is %f.',model.chiSquares(ind),model.p(ind));

%% ChiSquaredSigmaDistribution

%display it
figure
hold on
histogram(sum(model.chiSquaresMC,2),'Normalization','pdf')
yl = ylim;
plot(repmat(sum(model.chiSquares),1,2),yl,'k')
xlabel('\chi^2_{help}')
ylabel('probability density')
fname = 'ChiSquaredSigmaDistribution_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

fprintf('The chi squared value is %f and the p-value is %f.',sum(model.chiSquares),model.pSum);

%% chiSquaredDeltaT, pDeltaT

%display chiSquaredDeltaT
figure
plot(options.hist.timePoints,model.chiSquares)
xlabel('\Delta t / s')
ylabel('\chi^2_{help}')
fname = 'ChiSquaredDeltaT_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%display pDeltaT
figure
plot(options.hist.timePoints,model.p)
xlabel('\Delta t / s')
ylabel('p_{help}')
fname = 'PDeltaT_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%% PhaseHistogramsDeltaT

% find dynamic range
minVal = min([model.hist_obs(:); model.hist_pred(:)]);
maxVal = max([model.hist_obs(:); model.hist_pred(:)]);

maxVal = 260;

% display observed
figure
imagesc([options.hist.timePoints(1) options.hist.timePoints(end)],[1 options.processing.nPosPhases],model.hist_obs)
caxis([minVal maxVal])
xlabel('\Delta t / s')
ylabel('position bin $l$')
yticks(1:options.processing.nPosPhases)
colorbar
fname = 'ObservedPhaseHistogramsDeltaT_raw';
set(gca,'XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

% display expected
figure
imagesc([options.hist.timePoints(1) options.hist.timePoints(end)],[1 options.processing.nPosPhases],model.hist_pred)
caxis([minVal maxVal])
xlabel('\Delta t / s')
ylabel('position bin $l$')
yticks(1:options.processing.nPosPhases)
colorbar
fname = 'ExpectedPhaseHistogramsDeltaT_raw';
set(gca,'XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

%% SampleAndObservedSequences PREP

% stop in calcPosHist, right after generating a sample trajectory
model = matRad_generateMotionModel(options);

%% SampleAndObservedSequences

% absolute path to dumping folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

% sample sequence
t_simulated = (0:(numStepsModel_simulate-1)).*model.deltaT_sample;
l_simulated = matRad_runMarkovChain_P(model,numStepsModel_simulate,false);
x_simulated = model.indices.posSubPhase2Pos(model.indices.subPhase2PosSubPhase(l_simulated));

% display observed
figure
plot(data.t_sample-data.t_sample(1),data.x_sample,'k')
xlim([0 20])
xlabel('time / s')
ylabel('principal component of motion / mm')
fname = 'ObservedSequence_raw';
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
fname = 'SampleSequence_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);


