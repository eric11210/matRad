%% data origin

% file
options.data.origin             = 'file';
options.data.fileInfo.f         = 19;
options.data.fileInfo.m         = 1;
options.data.fileInfo.processed = true;
options.data.fileInfo.t_tot     = 600; % seconds

%% processing options

% for binning
options.processing.percExtTarg  = 1;
options.processing.nPosPhases   = 5;
options.processing.nVelPhases   = 1;    % is this necessary?
options.processing.velBinning   = true; % is this necessary? just set nVelSubPhases to 1

% for resampling
options.processing.fResample = 3; % Hz, 6 deg/s * 1 deg

% for FSM
options.processing.doFSM                = true;
options.processing.FSM.timeLd           = 0.5; % s
options.processing.FSM.velocityRangeInh = [-inf -0.5]; % mm/s
options.processing.FSM.velocityRangeExh = [0.5 inf]; % mm/s
options.processing.FSM.velocityRangeEOE = [0 10]; % mm/s
%options.processing.FSM.cTheta           = 4.5; % mm/s
%options.processing.FSM.cLambda          = 4; % mm
options.processing.FSM.cSLength         = 0.133; % s

% training ratio
options.processing.trainRatio = 0.5;

%% FFT options

options.FFT.doFFT       = true;
options.FFT.doWindowing = true;

%% histogram options

options.hist.doHist     = true;
options.hist.timePoints = 0:0.04:80; %s

%% convergence time options

options.convTime.doConvTime = true;
options.convTime.l2_targ    = 0.01;

%% get model

% best fit parameters
options.processing.nSubPerPosPhase  = 1;
options.processing.nSubPerVelPhase  = 7;
options.processing.FSM.nTimeFracs   = 10;

% extract model
model = matRad_generateMotionModel(options);



