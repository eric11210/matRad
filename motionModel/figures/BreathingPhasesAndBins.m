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
options.processing.nSubPerPosPhase  = 3;
options.processing.nSubPerVelPhase  = 7;
options.processing.FSM.nTimeFracs   = 10;


% get raw data
data = matRad_readMotionData(options.data.fileInfo,options.processing.percExtTarg);


%% start the processing, pause at the end

[data_train,data_test] = matRad_processMotionData(data,options.processing,options.data.fileInfo);

%% perform this code

% absolute path to figure folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

% modify FS_sample to cut off first/last cycles
FS_sample = FS_sample(indFirstCycle_dec_train:indLastCycle_dec_train);

% get indices where breathing phase changes
newFS           = false(size(FS_sample));
newFS(2:end)    = logical(diff(FS_sample) ~= 0);

% get indices for bin dividers
binDividers     = 1:nSubPerPosPhase:(nPosBins+1);
subBinDividers  = setdiff(1:(nPosBins+1),binDividers);

subBinDividers(subBinDividers > max(binDividers)) = [];

% display it
figure
hold on
plot(data_train.t_sample-data_train.t_sample(1),data_train.x_sample,'k')
plot(data_train.t_sample(newFS)-data_train.t_sample(1),data_train.x_sample(newFS),'kp','MarkerSize',10,'MarkerFaceColor','yellow')
plot([0 20],repmat(xBounds(binDividers),2,1),'b')
plot([0 20],repmat(xBounds(subBinDividers),2,1),'r-.')
xlim([0 20])
xlabel('time / s')
ylabel('principal component of motion / mm')
fname = 'BreathingPhasesAndBins_raw';
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);
