function ct = matRad_binFrames2Phases(ct,numPhases)

% get data from ct
tumourPos = ct.tumourMotion.coordsVox;
t = ct.tumourMotion.t;

% extract principal component
[~,components] = pca(tumourPos);
x_XCAT = components(:,1);

% get boundaries on bins
lBoundsMax = max(x_XCAT);
lBoundsMin = min(x_XCAT);
nBins = numPhases/2;
lBounds = linspace(lBoundsMax,lBoundsMin,nBins+1);

% preliminary binning of data
l_XCAT = zeros(size(x_XCAT));
for bin = 1:nBins
    l_XCAT(lBounds(bin+1) <= x_XCAT & x_XCAT <= lBounds(bin)) = bin;
end

% include inhale/exhale information:
% breathing in is the beginning of the cycle, so it should have
% lower phase number
% breathing out is the end of the cycle, so it should have higher phase
% number

% define breathing in/out by peaks: points after a min and before a max are
% exhale, points after a max and before a min are inhale.
%[~,ind_maxPeaks] = findpeaks(x_XCAT,'MinPeakProminence',0.5);
[~,ind_minPeaks] = findpeaks(-x_XCAT,'MinPeakProminence',0.5);

l_XCAT(ind_minPeaks:end) = numPhases+1-l_XCAT(ind_minPeaks:end);

ct.tumourMotion.frames2Phases = l_XCAT(1:(end-1));
ct.tumourMotion.nFramesPerPhase = accumarray(ct.tumourMotion.frames2Phases,1);
ct.tumourMotion.nFramesPerPhase = repelem(ct.tumourMotion.nFramesPerPhase,ct.tumourMotion.nFramesPerPhase);
ct.tumourMotion.numFrames = numel(l_XCAT)-1;
ct.tumourMotion.numPhases = numPhases;


end