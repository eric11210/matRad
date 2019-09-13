%% data origin

% file
options.data.origin             = 'file';
options.data.fileInfo.p         = 9;
options.data.fileInfo.f         = 19;
options.data.fileInfo.m         = 1;
options.data.fileInfo.processed = true;

%% processing options

% for binning
options.processing.percExtTarg  = 1;
options.processing.nPosPhases   = 5;
options.processing.nVelPhases   = 1;    % is this necessary?
options.processing.velBinning   = true; % is this necessary? just set nVelSubPhases to 1

% for resampling
options.processing.fResample = 6; % Hz, 6 deg/s * 1 deg

% for FSM
options.processing.doFSM                = true;
options.processing.FSM.timeLd           = 0.5; % s
options.processing.FSM.velocityRangeInh = [-inf -0.5]; % mm/s
options.processing.FSM.velocityRangeExh = [0.5 inf]; % mm/s
options.processing.FSM.velocityRangeEOE = [0 10]; % mm/s
%options.processing.FSM.cTheta           = 4.5; % mm/s
%options.processing.FSM.cLambda          = 4; % mm
options.processing.FSM.cSLength         = 0.133; % s

%% deprecated?

options.convTime.percRMSD_targ  = 1;
options.stdANDfft.nHistories    = 10;
options.stdANDfft.doWindowing   = true;

%% histogram options

options.hist.timePoints = [1 4 6 22.5 45 90 180 270 360]/6; %s

%% get model

% best fit parameters
options.processing.nSubPerPosPhase  = 3;
options.processing.nSubPerVelPhase  = 10;
options.processing.FSM.nTimeFracs   = 10;

% extract model
model = matRad_generateMotionModel(options);



