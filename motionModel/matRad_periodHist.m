function histT = matRad_periodHist(l,deltaT,lPeriodThresh,binEdgesMin,binWidths)

% initialize up the histogram
histT = zeros(numel(binEdgesMin),1);

% we define a new cycle to be when the patient first begins inhaling
% i.e. when the phase resets to a value less than the threshold

% note that the first part of the signal is probably an incomplete cycle,
% so just ignore it. i.e. look for the first complete cycle

% first find which is the first index above the threshold
startInd = find(l > lPeriodThresh,1,'first');

% then find the next index which is below the threshold. this is the
% beginning of the first complete cycle. add 1 to start at next index
startInd = find(l(startInd:end) <= lPeriodThresh,1,'first')+1;

% set the half-cycle flag to false. this indicates that a half-cycle has
% been completed, i.e. that the patient has started expiration
halfCycle = false;

% set the period to 0
period = 0;

% loop through the entire phase array
for currentInd = startInd:numel(l)
    
    % add deltaT to the period
    period = period+deltaT;
    
    if ~halfCycle
        % half-cycle has been completed as soon as phase exceeds threshold
        halfCycle = l(currentInd) > lPeriodThresh;
    end
    
    if halfCycle && l(currentInd) <= lPeriodThresh
        % a new cycle has begun as soon both a half-cycle has been completed,
        % and when the phase dips below the threshold
        
        % find the correct bin
        bin = binEdgesMin <= period & period < binEdgesMin + binWidths;
        
        % increment bin by 1
        histT(bin) = histT(bin)+1;
        
        % then reset the period to 0
        period = 0;
        
        % reset halfCycle to false
        halfCycle = false;
    end
end


end

