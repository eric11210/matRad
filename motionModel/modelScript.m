%% data origin

% file
options.data.origin             = 'file';
options.data.fileInfo.p         = 9;
options.data.fileInfo.f         = 19;
options.data.fileInfo.m         = 1;
options.data.fileInfo.processed = true;

% function
%options.data.origin             = 'function';
options.data.function           = 'sawtooth';
options.data.samplingFreq       = 26;
options.data.totTime            = 300;
options.data.period             = 5;
options.data.amplitude          = 1.5;
%options.data.function           = 'square';

%% processing options

% for binning
options.processing.percExtTarg  = 1;
options.processing.nPosPhases   = 10;
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

%% parameter options

nSubPerPosPhaseVec  = 1:10;
nSubPerVelPhaseVec  = 1:10;
nTimeFracsVec       = 2:10; %%

totNumPosSubPhases  = numel(nSubPerPosPhaseVec);
totNumVelSubPhases  = numel(nSubPerVelPhaseVec);
totNumTimeFracs     = numel(nTimeFracsVec);

nRuns = totNumPosSubPhases.*totNumVelSubPhases.*totNumTimeFracs;

% deprecated?
Pij_meanOfStdVec = zeros(size(nSubPerPosPhaseVec));
Pij_relMeanOfStdVec = zeros(size(nSubPerPosPhaseVec));
convergeTVec = zeros(size(nSubPerPosPhaseVec));

%% set up outputs

pSumArr         = zeros(nRuns,1);
chiSquareSumArr = zeros(nRuns,1);

pArr            = zeros(nRuns,numel(options.hist.timePoints));
chiSquaresArr   = zeros(nRuns,numel(options.hist.timePoints));

nSubPerPosPhaseArr  = zeros(nRuns,1);
nSubPerVelPhaseArr  = zeros(nRuns,1);
nTimeFracsArr       = zeros(nRuns,1);

%% loop over parameters

i = 101; %%

for nTimeFracs = nTimeFracsVec
    
    fprintf('Number of time fractions: %d.\n',nTimeFracs);
    
    for nSubPerVelPhase = nSubPerVelPhaseVec
        
        fprintf('Number of vel subphases per vel phase: %d.\n',nSubPerVelPhase);
        
        for nSubPerPosPhase = nSubPerPosPhaseVec
            
            fprintf('Number of pos subphases per pos phase: %d.\n',nSubPerPosPhase);
            
            fprintf('Run %d of %d.\n',i,nRuns);
            
            % variable parameters
            options.processing.nSubPerPosPhase  = nSubPerPosPhase;
            options.processing.nSubPerVelPhase  = nSubPerVelPhase;
            options.processing.FSM.nTimeFracs   = nTimeFracs;
            
            % extract model
            model = matRad_generateMotionModel(options);
            
            % get outputs
            pSumArr(i)          = model.pSum;
            chiSquareSumArr(i)  = sum(model.chiSquares);
            
            pArr(i,:)           = model.p;
            chiSquaresArr(i,:)  = model.chiSquares;
            
            nSubPerPosPhaseArr(i)   = nSubPerPosPhase;
            nSubPerVelPhaseArr(i)   = nSubPerVelPhase;
            nTimeFracsArr(i)        = nTimeFracs;
            
            % put data in vectors
            %{
                Pij_meanOfStdVec(i) = mean(Pij_std(Pij_std ~= 0));
                Pij_relMeanOfStdVec(i) = mean(Pij_relStd(Pij_relStd ~= 0));
                convergeTVec(i) = convergeT;
            %}
            i = i+1;
            
            % save after every iteration
            save('stats','*Arr','i')
        end
        
    end
    
end

