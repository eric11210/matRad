function [l_sample,t_sample,deltaT_sample,subPhase2Phase,nSubPhases] = motionData2Phase(x_cut,t_cut,nPhases,nSubPerPhase)

% resample data to have consistent sampling interval
deltaT_sample = mean(diff(t_cut));
t_sample = (min(t_cut):deltaT_sample:max(t_cut))';
x_sample = interp1(t_cut,x_cut,t_sample);

% determine number of bins for the motion data (not incl. inhale/exhale
% information)
lBoundsMax = max(x_sample);
lBoundsMin = min(x_sample);

% find min and max thresholds for the extreme phases of motion
nSubPhasesMaxExt = 0;
nSubPhasesMinExt = 0;

percTarg = 1;
percAbove = 0;
percBelow = 0;

while percAbove < percTarg || percBelow < percTarg
    
    if percAbove < percTarg
        nSubPhasesMaxExt = nSubPhasesMaxExt+1;
    end
    
    if percBelow < percTarg
        nSubPhasesMinExt = nSubPhasesMinExt+1;
    end
    
    nBins = nPhases.*nSubPerPhase./2+nSubPhasesMaxExt+nSubPhasesMinExt;
    deltaL = (lBoundsMax-lBoundsMin)./nBins;
    lMaxThres = lBoundsMax-nSubPhasesMaxExt.*deltaL;
    lMinThres = lBoundsMin+nSubPhasesMinExt.*deltaL;
    
    percAbove = 100.*nnz(x_sample > lMaxThres)./numel(x_sample);
    percBelow = 100*nnz(x_sample < lMinThres)./numel(x_sample);
    
end

nFrames = nBins*2;
subPhase2Phase = [repmat(nPhases+1,nSubPhasesMaxExt,1); repelem((1:nPhases/2)',nSubPerPhase,1); repmat(nPhases+2,nSubPhasesMinExt*2,1); repelem(((nPhases/2+1):nPhases)',nSubPerPhase,1); repmat(nPhases+1,nSubPhasesMaxExt,1)];
nSubPhases = numel(subPhase2Phase);

% preliminary binning of data
lBounds = linspace(lBoundsMax,lBoundsMin,nBins+1);

l_sample = zeros(size(x_sample));
for bin = 1:nBins
    l_sample(lBounds(bin+1) <= x_sample & x_sample <= lBounds(bin)) = bin;
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