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
options.processing.nPosPhases   = 5;
options.processing.nVelPhases   = 1;    % is this necessary?
options.processing.velBinning   = true; % is this necessary? just set nVelSubPhases to 1

% for FSM
options.processing.doFSM                = true;
options.processing.FSM.timeLd           = 0.5; % s
options.processing.FSM.velocityRangeInh = [-inf -0.5]; % mm/s
options.processing.FSM.velocityRangeExh = [0.5 inf]; % mm/s
options.processing.FSM.velocityRangeEOE = [0 10]; % mm/s
%options.processing.FSM.cTheta           = 4.5; % mm/s
%options.processing.FSM.cLambda          = 4; % mm
options.processing.FSM.cSLength         = 0.133; % s

%% FFT options

options.FFT.doFFT       = true;
options.FFT.doWindowing = true;

%% histogram options

options.hist.doHist     = true;
options.hist.timePoints = 0:0.04:80; %s

%% convergence time options

options.convTime.doConvTime = true;
options.convTime.l2_targ    = 0.01;

%% parameter options

nSubPerPosPhaseVec  = 1:10;
nSubPerVelPhaseVec  = 1:10;
nTimeFracsVec       = 1:10;
fResampleVec        = [1 3 6 9 12 24];

totNumPosSubPhases  = numel(nSubPerPosPhaseVec);
totNumVelSubPhases  = numel(nSubPerVelPhaseVec);
totNumTimeFracs     = numel(nTimeFracsVec);
totNumFResample     = numel(fResampleVec);

nRuns = totNumPosSubPhases.*totNumVelSubPhases.*totNumTimeFracs.*totNumFResample;

%% set up outputs

pSumArr         = zeros(nRuns,1);
chiSquareSumArr = zeros(nRuns,1);

pArr            = zeros(nRuns,numel(options.hist.timePoints));
chiSquaresArr   = zeros(nRuns,numel(options.hist.timePoints));

nSubPerPosPhaseArr  = zeros(nRuns,1);
nSubPerVelPhaseArr  = zeros(nRuns,1);
nTimeFracsArr       = zeros(nRuns,1);
fResampleArr        = zeros(nRuns,1);
convergeTArr        = zeros(nRuns,1);

%% loop over parameters

i = 1;

for fResample = fResampleVec
    
    fprintf('Resampling frequency: %d Hz.\n',fResample);
    
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
                options.processing.fResample        = fResample;
                
                % extract model
                model = matRad_generateMotionModel(options);
                
                % get outputs
                if options.hist.doHist
                    pSumArr(i)          = model.pSum;
                    chiSquareSumArr(i)  = sum(model.chiSquares);
                    pArr(i,:)           = model.p;
                    chiSquaresArr(i,:)  = model.chiSquares;
                end
                if options.convTime.doConvTime
                    convergeTArr(i)     = model.convergeT;
                end
                
                nSubPerPosPhaseArr(i)   = nSubPerPosPhase;
                nSubPerVelPhaseArr(i)   = nSubPerVelPhase;
                nTimeFracsArr(i)        = nTimeFracs;
                fResampleArr(i)         = fResample;
                
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
    
end

