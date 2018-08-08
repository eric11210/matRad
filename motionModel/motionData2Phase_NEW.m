function [l_sample,t_sample,deltaT_sample,subPhase2Phase,nSubPhases] = motionData2Phase_NEW(x_cut,t_cut,nPhases,nSubPerPhase)

% resample data to have consistent sampling interval
deltaT_sample = mean(diff(t_cut));
t_sample = (min(t_cut):deltaT_sample:max(t_cut))';
x_sample = interp1(t_cut,x_cut,t_sample);

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

percTarg = 1;
percAbove = 0;
percBelow = 0;

while percAbove < percTarg || percBelow < percTarg
    
    if percAbove < percTarg
        nPosSubPhasesMaxExt = nPosSubPhasesMaxExt+1;
    end
    
    if percBelow < percTarg
        nPosSubPhasesMinExt = nPosSubPhasesMinExt+1;
    end
    
    nPosBins = nPhases.*nSubPerPhase+nPosSubPhasesMaxExt+nPosSubPhasesMinExt;
    deltaL = (xBoundsMax-xBoundsMin)./nPosBins;
    lMaxThres = xBoundsMax-nPosSubPhasesMaxExt.*deltaL;
    lMinThres = xBoundsMin+nPosSubPhasesMinExt.*deltaL;
    
    percAbove = 100.*nnz(x_sample > lMaxThres)./numel(x_sample);
    percBelow = 100*nnz(x_sample < lMinThres)./numel(x_sample);
    
end

% velocity thresholds, assume symmetric around 0
nVelSubPhasesExt = 0;

percTarg = 1;
percAboveAndBelow = 0;

while percAboveAndBelow < 2*percTarg
    
    nVelSubPhasesExt = nVelSubPhasesExt+1;
    
    nVelBins = nPhases.*nSubPerPhase+nVelSubPhasesExt;
    if mod(nVelBins,2) ~= 0
        % we want an even number of bins, so that they split nicely into
        % positive and negative
        continue
    end
    deltaL = 2*vBoundsM./nVelBins;
    lMaxThres = vBoundsM-nVelSubPhasesExt.*deltaL;
    lMinThres = -vBoundsM+nVelSubPhasesExt.*deltaL;
    
    percAboveAndBelow = 100.*nnz(v_sample > lMaxThres)./numel(v_sample)+100*nnz(v_sample < lMinThres)./numel(v_sample);
end

nFrames = nPosBins*2;
subPhase2Phase = [repmat(nPhases+1,nPosSubPhasesMaxExt,1); repelem((1:nPhases/2)',nSubPerPhase,1); repmat(nPhases+2,nPosSubPhasesMinExt*2,1); repelem(((nPhases/2+1):nPhases)',nSubPerPhase,1); repmat(nPhases+1,nPosSubPhasesMaxExt,1)];
nSubPhases = numel(subPhase2Phase);

% determine bounds on phases
% phase 1 is max tumour position, min velocity (beginning of inhale)
xBounds = linspace(xBoundsMax,xBoundsMin,nPosBins+1);
vBounds = linspace(-vBoundsM,vBoundsM,nVelBins+1);

% do position binning
l_sample = zeros(size(x_sample));
for posBin = 1:nPosBins
    l_sample(xBounds(posBin+1) <= x_sample & x_sample <= xBounds(posBin)) = posBin;
end

% now do velocity binning
for velBin = 1:nVelBins
    l_sample(vBounds(velBin) <= v_sample & v_sample <= vBounds(velBin+1)) = l_sample(vBounds(velBin) <= v_sample & v_sample <= vBounds(velBin+1))+nPosBins.*(velBin-1);
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
    l_sample(1:ind_maxPeaks(1)) = nFrames+1-l_sample(1:ind_maxPeaks(1));
    
    for i = 1:numMinPeaks
        if i+1 > numMaxPeaks
            l_sample(ind_minPeaks(i):end) = nPhases_sample+1-l_sample(ind_minPeaks(i):end);
            break
        end
        if ind_minPeaks(i) > ind_maxPeaks(i+1)
            error('Error in peak-finding algorithm: min and max peaks not arranged correctly.');
        end
        l_sample(ind_minPeaks(i):ind_maxPeaks(i+1)) = nFrames+1-l_sample(ind_minPeaks(i):ind_maxPeaks(i+1));
    end
    
else
    % first peak is a min
    
    for i = 1:numMinPeaks
        if i > numMaxPeaks
            l_sample(ind_minPeaks(i):end) = nFrames+1-l_sample(ind_minPeaks(i):end);
            break
        end
        if ind_minPeaks(i) > ind_maxPeaks(i)
            error('Error in peak-finding algorithm: min and max peaks not arranged correctly.');
        end
        l_sample(ind_minPeaks(i):ind_maxPeaks(i)) = nFrames+1-l_sample(ind_minPeaks(i):ind_maxPeaks(i));
    end
    
end

end