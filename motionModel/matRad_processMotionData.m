function [data_train,data_test] = matRad_processMotionData(data,options)

%% get variables from inputs
% data
t_cut   = data.t_cut;
x_cut   = data.x_cut;

% options
percExtTarg         = options.percExtTarg;
nPosPhases          = options.nPosPhases;
nVelPhases          = options.nVelPhases;
nSubPerPosPhase     = options.nSubPerPosPhase;
nSubPerVelPhase     = options.nSubPerVelPhase;

%% process data, prepare for binning

% resample data to have consistent sampling interval
deltaT_sample = mean(diff(t_cut));
t_sample = (min(t_cut):deltaT_sample:max(t_cut))';
x_sample = interp1(t_cut,x_cut,t_sample);

% remove mean so data is centred around 0
x_sample = x_sample-mean(x_sample);

% decimate data to desired frequency
decFactor           = round(1./(deltaT_sample.*options.fResample));
%x_sample        = lowpass(x_sample,0.8./decFactor);
deltaT_sample_dec   = decFactor.*deltaT_sample;
x_sample_dec        = decimate(x_sample,decFactor,8);
t_sample_dec        = flipud((max(t_sample):(-deltaT_sample_dec):min(t_sample))');

% determine velocity
v_sample_dec = gradient(x_sample_dec,deltaT_sample_dec);

% determine number of bins for the motion data (not incl. inhale/exhale
% information)
xBoundsMax = max(x_sample_dec);
xBoundsMin = min(x_sample_dec);
vBoundsMax = max(v_sample_dec);
vBoundsMin = min(v_sample_dec);
vBoundsM = max(abs([vBoundsMax vBoundsMin]));

% find min and max thresholds for the extreme phases of motion

% position thresholds, assume asymmetric
nPosSubPhasesMaxExt = 0;
nPosSubPhasesMinExt = 0;
nPosSubPhasesMaxExt_final = 0;
nPosSubPhasesMinExt_final = 0;

percAbove = 0;
percBelow = 0;

while percAbove < percExtTarg || percBelow < percExtTarg
    
    nPosBins = nPosPhases.*nSubPerPosPhase/2+nPosSubPhasesMaxExt+nPosSubPhasesMinExt;
    deltaL = (xBoundsMax-xBoundsMin)./nPosBins;
    lMaxThres = xBoundsMax-nPosSubPhasesMaxExt.*deltaL;
    lMinThres = xBoundsMin+nPosSubPhasesMinExt.*deltaL;
    
    percAbove = 100.*nnz(x_sample_dec > lMaxThres)./numel(x_sample_dec);
    percBelow = 100*nnz(x_sample_dec < lMinThres)./numel(x_sample_dec);
    
    if percAbove < percExtTarg
        nPosSubPhasesMaxExt_final = nPosSubPhasesMaxExt;
        nPosSubPhasesMaxExt = nPosSubPhasesMaxExt+1;
    end
    
    if percBelow < percExtTarg
        nPosSubPhasesMinExt_final = nPosSubPhasesMinExt;
        nPosSubPhasesMinExt = nPosSubPhasesMinExt+1;
    end
end

nPosSubPhasesMaxExt = nPosSubPhasesMaxExt_final;
nPosSubPhasesMinExt = nPosSubPhasesMinExt_final;

nPosBins = nPosPhases.*nSubPerPosPhase/2+nPosSubPhasesMaxExt+nPosSubPhasesMinExt;
deltaL = (xBoundsMax-xBoundsMin)./nPosBins;
lMaxThres = xBoundsMax-nPosSubPhasesMaxExt.*deltaL;
lMinThres = xBoundsMin+nPosSubPhasesMinExt.*deltaL;
percAbove = 100.*nnz(x_sample_dec > lMaxThres)./numel(x_sample_dec);
percBelow = 100*nnz(x_sample_dec < lMinThres)./numel(x_sample_dec);

if options.doFSM
    nPosSubPhases = nPosBins;
else
    nPosSubPhases = 2.*nPosBins;
end

xBounds = linspace(xBoundsMax,xBoundsMin,nPosBins+1);

%% bin data

% do position binning
l_sample_dec = zeros(size(x_sample_dec));
for posBin = 1:nPosBins
    l_sample_dec(xBounds(posBin+1) <= x_sample_dec & x_sample_dec <= xBounds(posBin)) = posBin;
end

if options.doFSM
    % FSM takes care of distinction between exhaling, inhaling, EOE, and
    % IRR
    
    % determine threshold velocity and amplitude values
    % this is 1/4 the average velocity, using a typical breathing period of
    % 4 s
    options.FSM.cTheta           = 2*2*std(x_sample)./(4*4); % mm/s
    % this is 1/3 the average amplitude
    options.FSM.cLambda          = 2*std(x_sample)/3; % mm
    
    % do the FSM on the non-decimated trace
    FS_sample = matRad_doFSM(x_sample,deltaT_sample,options.FSM);
    
    % now decimate it
    endInd = numel(x_sample);
    startInd = decFactor-(decFactor*ceil(endInd./decFactor)-endInd);
    FS_sample_dec = FS_sample(startInd:decFactor:endInd);
    
    % cut off the signal until we begin and end with a complete cycle
    % (1>2>3>1>...)
    % also getting starting and stopping indices for the training and
    % testing data
    [FS_sample_dec,t_sample_dec,x_sample_dec,v_sample_dec,l_sample_dec,indFirstCycle_train,indLastCycle_train,indFirstCycle_test,indLastCycle_test] = matRad_splitTrainTest(FS_sample_dec,t_sample_dec,x_sample_dec,v_sample_dec,l_sample_dec,0.5);
    
    % do time split
    FS_sample_dec = matRad_FSMfracTime(FS_sample_dec,deltaT_sample_dec,options.FSM);
    
    % each non-IRR is split into time fractions; the IRR state is not
    nStates = 3.*options.FSM.nTimeFracs+1;
    
    % include finite states in l_sample
    l_sample_dec = (FS_sample_dec-1).*nPosSubPhases + l_sample_dec;
    
else
    % only do this part if not doing FSM
    
    % number of states is 1
    nStates = 1;
    
    % include inhale/exhale information:
    % breathing in is the beginning of the cycle, so it should have
    % lower phase number
    % breathing out is the end of the cycle, so it should have higher phase
    % number
    
    % define breathing in/out by peaks: points after a min and before a max are
    % exhale, points after a max and before a min are inhale.
    [~,ind_maxPeaks] = findpeaks(x_sample_dec,'MinPeakProminence',0.5);
    [~,ind_minPeaks] = findpeaks(-x_sample_dec,'MinPeakProminence',0.5);
    
    numMinPeaks = numel(ind_minPeaks);
    numMaxPeaks = numel(ind_maxPeaks);
    
    if abs(numMinPeaks-numMaxPeaks) > 1
        error('Error in peak-finding algorithm: difference in min and max peaks greater than 1.');
    end
    
    if ind_maxPeaks(1) < ind_minPeaks(1)
        % first peak is a max
        l_sample_dec(1:ind_maxPeaks(1)) = nPosSubPhases+1-l_sample_dec(1:ind_maxPeaks(1));
        
        for i = 1:numMinPeaks
            if i+1 > numMaxPeaks
                l_sample_dec(ind_minPeaks(i):end) = nPosSubPhases+1-l_sample_dec(ind_minPeaks(i):end);
                break
            end
            if ind_minPeaks(i) > ind_maxPeaks(i+1)
                error('Error in peak-finding algorithm: min and max peaks not arranged correctly.');
            end
            l_sample_dec(ind_minPeaks(i):ind_maxPeaks(i+1)) = nPosSubPhases+1-l_sample_dec(ind_minPeaks(i):ind_maxPeaks(i+1));
        end
        
    else
        % first peak is a min
        
        for i = 1:numMinPeaks
            if i > numMaxPeaks
                l_sample_dec(ind_minPeaks(i):end) = nPosSubPhases+1-l_sample_dec(ind_minPeaks(i):end);
                break
            end
            if ind_minPeaks(i) > ind_maxPeaks(i)
                error('Error in peak-finding algorithm: min and max peaks not arranged correctly.');
            end
            l_sample_dec(ind_minPeaks(i):ind_maxPeaks(i)) = nPosSubPhases+1-l_sample_dec(ind_minPeaks(i):ind_maxPeaks(i));
        end
        
    end
    
end


if options.velBinning
    % velocity thresholds, assume symmetric around 0
    nVelSubPhasesExt        = 0;
    nVelSubPhasesExt_final  = 0;
    
    percAboveAndBelow = 0;
    
    %{
    while percAboveAndBelow < 2*percExtTarg
        
        nVelBins = nVelPhases.*nSubPerVelPhase+2*nVelSubPhasesExt;
        if mod(nVelBins,2) ~= 0
            % we want an even number of bins, so that they split nicely into
            % positive and negative
            if percAboveAndBelow < 2*percExtTarg
                nVelSubPhasesExt = nVelSubPhasesExt+1;
            end
            continue
        end
        deltaL = 2*vBoundsM./nVelBins;
        lMaxThres = vBoundsM-nVelSubPhasesExt.*deltaL;
        lMinThres = -vBoundsM+nVelSubPhasesExt.*deltaL;
        
        percAboveAndBelow = 100.*nnz(v_sample > lMaxThres)./numel(v_sample)+100*nnz(v_sample < lMinThres)./numel(v_sample);
        
        if percAboveAndBelow < 2*percExtTarg
            nVelSubPhasesExt_final = nVelSubPhasesExt;
            nVelSubPhasesExt = nVelSubPhasesExt+1;
        end 
    end
    %}
    
    nVelSubPhasesExt = nVelSubPhasesExt_final;
    
    nVelBins = nVelPhases.*nSubPerVelPhase+2*nVelSubPhasesExt;
    deltaL = 2*vBoundsM./nVelBins;
    lMaxThres = vBoundsM-nVelSubPhasesExt.*deltaL;
    lMinThres = -vBoundsM+nVelSubPhasesExt.*deltaL;
    percAboveAndBelow = 100.*nnz(v_sample_dec > lMaxThres)./numel(v_sample_dec)+100*nnz(v_sample_dec < lMinThres)./numel(v_sample_dec);
    
    nVelSubPhases = nVelBins;
    vBounds = linspace(-vBoundsM,vBoundsM,nVelBins+1);
    
    % now do velocity binning
    for velBin = 1:nVelBins
        l_sample_dec(vBounds(velBin) <= v_sample_dec & v_sample_dec <= vBounds(velBin+1)) = l_sample_dec(vBounds(velBin) <= v_sample_dec & v_sample_dec <= vBounds(velBin+1))+nStates.*nPosSubPhases.*(velBin-1);
    end
    
    
else
    
    nVelPhases          = 0;
    nVelSubPhases       = 1;
    nVelSubPhasesExt    = 0;
end

nSubPhases = nPosSubPhases*nStates*nVelSubPhases;

%% conversion indices

% convert from sub phase to position-state mix
subPhase2PosStateSubPhase = mod((1:nSubPhases)',nStates.*nPosSubPhases);
subPhaseChangeInd = subPhase2PosStateSubPhase == 0;
subPhase2PosStateSubPhase(subPhaseChangeInd) = nStates.*nPosSubPhases;

% convert from position-state subphase to position only
posStateSubPhase2PosSubPhase = mod((1:(nPosSubPhases*nStates))',nPosSubPhases);
posSubPhaseChangeInd = posStateSubPhase2PosSubPhase == 0;
posStateSubPhase2PosSubPhase(posSubPhaseChangeInd) = nPosSubPhases;

% convert from position-state subphase to state only
posStateSubPhase2State = floor(((1:(nPosSubPhases*nStates))-1)'./nPosSubPhases)+1;

% convert from sub phase to position only
subPhase2PosSubPhase = posStateSubPhase2PosSubPhase(subPhase2PosStateSubPhase);

% convert from sub phase to state only
subPhase2State = posStateSubPhase2State(subPhase2PosStateSubPhase);

% add the exhale and EOE states if doing FSM
if options.doFSM
    
    % first convert from state to FS (1,2,3,4)
    state2FS = floor(((1:(nStates))-1)'./options.FSM.nTimeFracs)+1;
    
    % then convert from sub phase to FS
    subphase2FS = state2FS(subPhase2State);
    
    % if state is exhale or EOE, modify subPhase2PosSubPhase
    subPhase2PosSubPhase(subphase2FS == 2 | subphase2FS == 3) = 2.*nPosSubPhases+1-subPhase2PosSubPhase(subphase2FS == 2 | subphase2FS == 3);
    
end

% convert from sub phase to velocity only
subPhase2VelSubPhase = floor((1:nSubPhases)'./(nStates.*nPosSubPhases))+1;
subPhase2VelSubPhase(subPhaseChangeInd) = subPhase2VelSubPhase(subPhaseChangeInd)-1;

% convert from pos sub phase to pos phase
posSubPhase2PosPhase = [repmat(nPosPhases+1,nPosSubPhasesMaxExt,1); repelem((1:nPosPhases/2)',nSubPerPosPhase,1); repmat(nPosPhases+2,nPosSubPhasesMinExt*2,1); repelem(((nPosPhases/2+1):nPosPhases)',nSubPerPosPhase,1); repmat(nPosPhases+1,nPosSubPhasesMaxExt,1)];

% convert from vel sub phase to vel phase
velSubPhase2VelPhase = [repmat(nVelPhases+1,nVelSubPhasesExt,1); repelem((1:nVelPhases)',nSubPerVelPhase,1); repmat(nVelPhases+2,nVelSubPhasesExt,1)];
if isempty(velSubPhase2VelPhase)
    velSubPhase2VelPhase = 1;
end

% convert from sub phase to pos phase
subPhase2PosPhase = posSubPhase2PosPhase(subPhase2PosSubPhase);
nSubPhasePerPosPhase = accumarray(subPhase2PosPhase,1);

% convert from sub phase to vel phase
subPhase2VelPhase = velSubPhase2VelPhase(subPhase2VelSubPhase);
nSubPhasePerVelPhase = accumarray(subPhase2VelPhase,1);

% convert from sub phase to phase
subPhase2Phase = subPhase2PosPhase+(nPosPhases+2).*(subPhase2VelPhase-1);
nSubPhasePerPhase = accumarray(subPhase2Phase,1);

% convert from position subphase to position
% remember to duplicate this to have both inhale and exhale
posSubPhase2Pos = (diff(xBounds)./2+xBounds(1:nPosBins))';
posSubPhase2Pos = [posSubPhase2Pos; flipud(posSubPhase2Pos)];

% convert from position phase to position
posPhase2Pos = zeros(nPosPhases+2,1);
for posPhase = 1:(nPosPhases+2)
    posPhase2Pos(posPhase) = mean(posSubPhase2Pos(posSubPhase2PosPhase == posPhase));
end

if options.velBinning
    % convert from velocity subphase to velocity
    velSubPhase2Vel = diff(vBounds)./2+vBounds(1:nVelBins);
    velSubPhase2Vel = velSubPhase2Vel';
    
    indices.velSubPhase2Vel         = velSubPhase2Vel;
end

%% put variables in struct
data.l_sample                   = l_sample_dec;
data.t_sample                   = t_sample_dec;
data.x_sample                   = x_sample_dec;
data.deltaT_sample              = deltaT_sample_dec;

indices.subPhase2PosSubPhase    = subPhase2PosSubPhase;
indices.subPhase2VelSubPhase    = subPhase2VelSubPhase;
indices.subPhase2State          = subPhase2State;

indices.posSubPhase2PosPhase    = posSubPhase2PosPhase;
indices.subPhase2PosPhase       = subPhase2PosPhase;
indices.nSubPhasePerPosPhase    = nSubPhasePerPosPhase;

indices.velSubPhase2VelPhase    = velSubPhase2VelPhase;
indices.subPhase2VelPhase       = subPhase2VelPhase;
indices.nSubPhasePerVelPhase    = nSubPhasePerVelPhase;

indices.subPhase2Phase          = subPhase2Phase;
indices.nSubPhasePerPhase       = nSubPhasePerPhase;

indices.nSubPhases              = nSubPhases;
indices.nPosSubPhases           = nPosSubPhases;
indices.nVelSubPhases           = nVelSubPhases;
%indices.nPhases                = nPhases;

indices.nPosPhases              = max(subPhase2PosPhase);

indices.posSubPhase2Pos         = posSubPhase2Pos;
indices.posPhase2Pos            = posPhase2Pos;

data.indices    = indices;

%% split data into training and testing
data_train = data;
data_test = data;

% indFirstCycle_train indLastCycle_train indFirstCycle_test indLastCycle_test
data_train.l_sample = data.l_sample(indFirstCycle_train:indLastCycle_train);
data_train.t_sample = data.t_sample(indFirstCycle_train:indLastCycle_train);
data_train.x_sample = data.x_sample(indFirstCycle_train:indLastCycle_train);

data_test.l_sample = data.l_sample(indFirstCycle_test:indLastCycle_test);
data_test.t_sample = data.t_sample(indFirstCycle_test:indLastCycle_test);
data_test.x_sample = data.x_sample(indFirstCycle_test:indLastCycle_test);

end
