function data = matRad_processMotionData(data,options)

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

% determine velocity
v_sample = gradient(x_sample,deltaT_sample);

% determine number of bins for the motion data (not incl. inhale/exhale
% information)
xBoundsMax = max(x_sample);
xBoundsMin = min(x_sample);
vBoundsMax = max(v_sample);
vBoundsMin = min(v_sample);
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
    
    percAbove = 100.*nnz(x_sample > lMaxThres)./numel(x_sample);
    percBelow = 100*nnz(x_sample < lMinThres)./numel(x_sample);
    
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
percAbove = 100.*nnz(x_sample > lMaxThres)./numel(x_sample);
percBelow = 100*nnz(x_sample < lMinThres)./numel(x_sample);

nPosSubPhases = 2*nPosBins;
xBounds = linspace(xBoundsMax,xBoundsMin,nPosBins+1);

%% bin data

% do position binning
l_sample = zeros(size(x_sample));
for posBin = 1:nPosBins
    l_sample(xBounds(posBin+1) <= x_sample & x_sample <= xBounds(posBin)) = posBin;
end

% include inhale/exhale information:
% breathing in is the beginning of the cycle, so it should have
% lower phase number
% breathing out is the end of the cycle, so it should have higher phase
% number

% define breathing in/out by peaks: points after a min and before a max are
% exhale, points after a max and before a min are inhale.
[~,ind_maxPeaks] = findpeaks(x_sample,'MinPeakProminence',0.5);
[~,ind_minPeaks] = findpeaks(-x_sample,'MinPeakProminence',0.5);

numMinPeaks = numel(ind_minPeaks);
numMaxPeaks = numel(ind_maxPeaks);

if abs(numMinPeaks-numMaxPeaks) > 1
    error('Error in peak-finding algorithm: difference in min and max peaks greater than 1.');
end

if ind_maxPeaks(1) < ind_minPeaks(1)
    % first peak is a max
    l_sample(1:ind_maxPeaks(1)) = nPosSubPhases+1-l_sample(1:ind_maxPeaks(1));
    
    for i = 1:numMinPeaks
        if i+1 > numMaxPeaks
            l_sample(ind_minPeaks(i):end) = nPosSubPhases+1-l_sample(ind_minPeaks(i):end);
            break
        end
        if ind_minPeaks(i) > ind_maxPeaks(i+1)
            error('Error in peak-finding algorithm: min and max peaks not arranged correctly.');
        end
        l_sample(ind_minPeaks(i):ind_maxPeaks(i+1)) = nPosSubPhases+1-l_sample(ind_minPeaks(i):ind_maxPeaks(i+1));
    end
    
else
    % first peak is a min
    
    for i = 1:numMinPeaks
        if i > numMaxPeaks
            l_sample(ind_minPeaks(i):end) = nPosSubPhases+1-l_sample(ind_minPeaks(i):end);
            break
        end
        if ind_minPeaks(i) > ind_maxPeaks(i)
            error('Error in peak-finding algorithm: min and max peaks not arranged correctly.');
        end
        l_sample(ind_minPeaks(i):ind_maxPeaks(i)) = nPosSubPhases+1-l_sample(ind_minPeaks(i):ind_maxPeaks(i));
    end
    
end


if options.velBinning
    % velocity thresholds, assume symmetric around 0
    nVelSubPhasesExt = 0;
    
    percAboveAndBelow = 0;
    
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
    
    nVelSubPhasesExt = nVelSubPhasesExt_final;
    
    nVelBins = nVelPhases.*nSubPerVelPhase+2*nVelSubPhasesExt;
    deltaL = 2*vBoundsM./nVelBins;
    lMaxThres = vBoundsM-nVelSubPhasesExt.*deltaL;
    lMinThres = -vBoundsM+nVelSubPhasesExt.*deltaL;
    percAboveAndBelow = 100.*nnz(v_sample > lMaxThres)./numel(v_sample)+100*nnz(v_sample < lMinThres)./numel(v_sample);
    
    nVelSubPhases = nVelBins;
    vBounds = linspace(-vBoundsM,vBoundsM,nVelBins+1);
    
    % now do velocity binning
    for velBin = 1:nVelBins
        l_sample(vBounds(velBin) <= v_sample & v_sample <= vBounds(velBin+1)) = l_sample(vBounds(velBin) <= v_sample & v_sample <= vBounds(velBin+1))+nPosSubPhases.*(velBin-1);
    end
    
    
else
    
    nVelPhases          = 0;
    nVelSubPhases       = 1;
    nVelSubPhasesExt    = 0;
end

nSubPhases = nPosSubPhases*nVelSubPhases;

%% conversion indices

% convert from sub phase to position only
subPhase2PosSubPhase = mod((1:nSubPhases)',nPosSubPhases);
changeInd = subPhase2PosSubPhase == 0;
subPhase2PosSubPhase(changeInd) = nPosSubPhases;

% convert from sub phase to velocity only
subPhase2VelSubPhase = floor((1:nSubPhases)'./nPosSubPhases)+1;
subPhase2VelSubPhase(changeInd) = subPhase2VelSubPhase(changeInd)-1;

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

%convert from sub phase to phase
subPhase2Phase = subPhase2PosPhase+(nPosPhases+2).*(subPhase2VelPhase-1);
nSubPhasePerPhase = accumarray(subPhase2Phase,1);

% convert from position subphase to position
% remember to duplicate this to have both inhale and exhale
posSubPhase2Pos = diff(xBounds)./2+xBounds(1:nPosBins);
posSubPhase2Pos = [posSubPhase2Pos fliplr(posSubPhase2Pos)]';

if options.velBinning
    % convert from velocity subphase to velocity
    velSubPhase2Vel = diff(vBounds)./2+vBounds(1:nVelBins);
    velSubPhase2Vel = velSubPhase2Vel';
    
    indices.velSubPhase2Vel         = velSubPhase2Vel;
end

%% put variables in struct
data.l_sample                   = l_sample;
data.t_sample                   = t_sample;
data.x_sample                   = x_sample;
data.deltaT_sample              = deltaT_sample;

indices.subPhase2PosSubPhase    = subPhase2PosSubPhase;
indices.subPhase2VelSubPhase    = subPhase2VelSubPhase;

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

indices.posSubPhase2Pos         = posSubPhase2Pos;

data.indices    = indices;

end

