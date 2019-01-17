options.data.origin             = 'file';
options.data.fileInfo.p         = 6;
options.data.fileInfo.f         = 19;
options.data.fileInfo.m         = 1;
options.data.fileInfo.processed = true;

options.data.origin             = 'file';
options.data.function           = 'sawtooth';
options.data.samplingFreq       = 26;
options.data.totTime            = 600;
options.data.period             = 5;
options.data.amplitude          = 1.5;
%options.data.function           = 'square';

options.processing.percExtTarg  = 0;
options.processing.nPosPhases   = 10;
options.processing.nVelPhases   = 10;
options.processing.velBinning   = false;
options.convTime.percRMSD_targ  = 1;
options.std.nHistories          = 10;

nSubPerPhaseVec = [1:10 20:10:100 200:100:500];
totNumSubPhases = numel(nSubPerPhaseVec);
Pij_meanOfStdVec = zeros(size(nSubPerPhaseVec));
Pij_relMeanOfStdVec = zeros(size(nSubPerPhaseVec));
convergeTVec = zeros(size(nSubPerPhaseVec));

i = 1;
convergeT = 0;
for nSubPerPhase = nSubPerPhaseVec
    
    fprintf('Number of subphases per phase: %d.\n',nSubPerPhase);
    fprintf('Run %d of %d.\n',i,totNumSubPhases);
    
    options.processing.nSubPerPosPhase = 100;
    options.processing.nSubPerVelPhase = 1;
    
    model = matRad_generateMotionModel(options);
    
    % put data in vectors
    Pij_meanOfStdVec(i) = mean(Pij_std(Pij_std ~= 0));
    Pij_relMeanOfStdVec(i) = mean(Pij_relStd(Pij_relStd ~= 0));
    convergeTVec(i) = convergeT;
    i = i+1;
end


save('Convergence times and errors.mat')
