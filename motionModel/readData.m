p = 6;
f = 19;
m = 1;

data = dlmread(sprintf('p%d_f%d_m%d_proc2_pca.csv',p,f,m));

t = data(:,1);
z = data(:,2);

t_cut = t(57135:end);
z_cut = z(57135:end);

deltaTSample = mean(diff(t_cut));
t_sample = (min(t_cut):deltaTSample:max(t_cut))';
z_sample = interp1(t_cut,z_cut,t_sample);

z_filter = movmean(z_sample,5);
z_filterGrad = gradient(z_filter,mean(diff(t_sample)));

nPhases = 10;
nBins = nPhases/2;

pBounds_global = linspace(max(z_sample),min(z_sample),nBins+1);

z_filterMax = findpeaks(z_sample,t_sample,'MinPeakDistance',2);
z_filterMin = -findpeaks(-z_sample,t_sample,'MinPeakDistance',2);
pBounds_local = linspace(mean(z_filterMax),mean(z_filterMin),nBins+1);

delta = (mean(z_filterMax)-mean(z_filterMin))./(nBins-1);
pBounds_localMod = linspace(mean(z_filterMax)+delta/2,mean(z_filterMin)-delta/2,nBins+1);

%{
p = zeros(size(z_filter));
for bin = 1:nBins
    p(pBounds_localMod(bin+1) <= z_filter & z_filter <= pBounds_localMod(bin)) = bin;
end
p(z_filterGrad > 0) = 11-p(z_filterGrad > 0);
%}

percAbove_global = 100*nnz(z_sample > max(pBounds_global))./numel(z_sample);
percBelow_global = 100*nnz(z_sample < min(pBounds_global))./numel(z_sample);
figure
hold on
plot(t_sample,z_sample)
plot(t_sample,repmat(pBounds_global,numel(t_sample),1),'r')
title(sprintf('Global: %.2f above, %.2f below',percAbove_global,percBelow_global))
xlabel('time / s')
ylabel('Primary component of target / mm')

percAbove_local = 100*nnz(z_sample > max(pBounds_local))./numel(z_sample);
percBelow_local = 100*nnz(z_sample < min(pBounds_local))./numel(z_sample);
figure
hold on
plot(t_sample,z_sample)
plot(t_sample,repmat(pBounds_local,numel(t_sample),1),'r')
title(sprintf('Local: %.2f above, %.2f below',percAbove_local,percBelow_local))
xlabel('time / s')
ylabel('Primary component of target / mm')

percAbove_localMod = 100*nnz(z_sample > max(pBounds_localMod))./numel(z_sample);
percBelow_percTarg = 100*nnz(z_sample < min(pBounds_localMod))./numel(z_sample);
figure
hold on
plot(t_sample,z_sample)
plot(t_sample,repmat(pBounds_localMod,numel(t_sample),1),'r')
title(sprintf('Local (modified): %.2f above, %.2f below',percAbove_localMod,percBelow_percTarg))
xlabel('time / s')
ylabel('Primary component of target / mm')



pBoundsMax_percTarg = max(pBounds_global);
pBoundsMin_percTarg = min(pBounds_global);
stepSize = (pBoundsMax_percTarg-pBoundsMin_percTarg)./10^4;

percAbove_percTarg = 100.*nnz(z_sample > pBoundsMax_percTarg)./numel(z_sample);
percTarg = 2;
while percAbove_percTarg < percTarg
    pBoundsMax_percTarg = pBoundsMax_percTarg-stepSize;
    
    percAbove_percTarg = 100.*nnz(z_sample > pBoundsMax_percTarg)./numel(z_sample);
end

percBelow_percTarg = 100*nnz(z_sample < pBoundsMin_percTarg)./numel(z_sample);
while percBelow_percTarg < percTarg
    pBoundsMin_percTarg = pBoundsMin_percTarg+stepSize;
    
    percBelow_percTarg = 100*nnz(z_sample < pBoundsMin_percTarg)./numel(z_sample);
end

pBounds_percTarg = linspace(pBoundsMax_percTarg,pBoundsMin_percTarg,nBins+1);


figure
hold on
plot(t_sample,z_sample)
plot(t_sample,repmat(pBounds_percTarg,numel(t_sample),1),'r')
title(sprintf('Targeted percentage: %.2f above, %.2f below',percAbove_percTarg,percBelow_percTarg))
xlabel('time / s')
ylabel('Primary component of target / mm')




pBounds = pBounds_percTarg;

p = zeros(size(z_filter));
for bin = 1:nBins
    p(pBounds(bin+1) <= z_filter & z_filter <= pBounds(bin)) = bin;
end

p(z_filterGrad > 0) = nPhases+1-p(z_filterGrad > 0);

p_fixed = p;

notFinished = 1;
i = 1;
while notFinished
    if p_fixed(i) == nPhases
        % patient is at the end of exhale, frequent oscillation between
        % phase nPhases and 1 due to "bumpiness" of motion trace
        
        % set all phases to nPhases until patient starts inhaling for real
        % this is when the phase following 1 is 2
        
        index_2after1 = find(p_fixed(i:end) == 2,1,'first')+i-1;
        index_lastnPhases = find(p_fixed(i:index_2after1) == nPhases,1,'last')+i-1;
        
        p_fixed(i:index_lastnPhases) = nPhases;
        
        i = index_lastnPhases+1;
    else
        
        i = i+1;
    end
    
    if i == numel(p_fixed)+1
        notFinished = 0;
    end
end

p_fixed(z_filter > pBounds(1)) = nPhases+1;
p_fixed(z_filter < pBounds(nBins+1)) = nPhases+2;



figure
hold on
plot(t_sample(1:1000),z_filter(1:1000))
plot(t_sample(1:1000),repmat(pBounds_percTarg,numel(t_sample(1:1000)),1),'r')
title(sprintf('Targeted percentage: %.2f above, %.2f below',percAbove_percTarg,percBelow_percTarg))
xlabel('time / s')
ylabel('Primary component of target / mm')

figure
hold on
plot(t_sample(1:1000),p_fixed(1:1000))
title(sprintf('Targeted percentage: %.2f above, %.2f below',percAbove_percTarg,percBelow_percTarg))
xlabel('time / s')
ylabel('Phase')

timeInPhase = zeros(nPhases+2,1);
hitsInPhase = zeros(nPhases+2,1);
numTransP2P = zeros(nPhases+2,nPhases+2);
numTransP2P_inclii = zeros(nPhases+2,nPhases+2);


for i = 1:numel(p_fixed)
    
    if i == 1       
        currPhase = p_fixed(i);
        nextPhase = p_fixed(i+1);
    elseif i == numel(p_fixed)
        prevPhase = p_fixed(i-1);
        currPhase = p_fixed(i);
    else
        prevPhase = p_fixed(i-1);
        currPhase = p_fixed(i);
        nextPhase = p_fixed(i+1);
    end
    
    % increment time in the current phase
    timeInPhase(currPhase) = timeInPhase(currPhase)+deltaTSample;
    
    if i ~= 1 && prevPhase ~= currPhase
        % increment the number of times spent in the current phase
        hitsInPhase(currPhase) = hitsInPhase(currPhase)+1;
    elseif i == 1
        % increment the number of times spent in the current phase
        hitsInPhase(currPhase) = hitsInPhase(currPhase)+1;
    end
    
    
    if i ~= numel(p_fixed) && nextPhase ~= currPhase
        % if we are transitioning from one phase to another, increment the
        % number of such transitions
        numTransP2P(currPhase,nextPhase) = numTransP2P(currPhase,nextPhase)+1;
    end
    
    if i ~= numel(p_fixed)
        numTransP2P_inclii(currPhase,nextPhase) = numTransP2P_inclii(currPhase,nextPhase)+1;
    end
end

q = zeros(nPhases+2,nPhases+2);
P_deltaTSample = zeros(nPhases+2,nPhases+2);
P_i = zeros(nPhases+2,1);

for i = 1:nPhases+2
    for j = 1:nPhases+2
        
        if i == j
            
            q(i,j) = -hitsInPhase(i)./timeInPhase(i);
            
        else
            
            q(i,j) = (hitsInPhase(i)./timeInPhase(i)).*numTransP2P(i,j)./sum(numTransP2P(i,:),2);
            
        end
        
        P_deltaTSample(i,j) = numTransP2P_inclii(i,j)./sum(numTransP2P_inclii(i,:),2);
    end
    
    P_i(i) = nnz(p_fixed == i)./numel(p_fixed);
end

qTest = (P_deltaTSample-eye(size(P_deltaTSample)))./deltaTSample;


deltaTModel = 0.5;
totalTModel = 2;
startingPhase = 1;
numPoints = round(totalTModel/deltaTModel)+1;

P_deltaTModel = expm(deltaTModel*q);
possibleTransitions = P_deltaTModel >= 0.05;
P_deltaTModel(~possibleTransitions) = 0;
P_deltaTModel = P_deltaTModel./repmat(sum(P_deltaTModel,2),1,size(P_deltaTModel,2));
P_mod = eye(size(q))+deltaTModel*q;

numPossibleTransitions = sum(possibleTransitions,2);
maxPossibleTransitions = max(numPossibleTransitions);

maxNumPaths = maxPossibleTransitions^(numPoints-1);

tracks.sequence = zeros(maxNumPaths,numPoints);
tracks.probability = zeros(maxNumPaths,1);

numPaths = 1;
numPaths_new = numPaths;
currPhase = startingPhase;
tracks.sequence(numPaths,1) = currPhase;
tracks.probability(numPaths) = 1;

for point = 1:(numPoints-1)
    
    currTotalNumPossibleTransitions = 0;
    
    for path = 1:numPaths
        
        currPhase = tracks.sequence(path,point);
        currPossibleTransitions = find(possibleTransitions(currPhase,:));
        currNumPossibleTransitions = numel(currPossibleTransitions);
        
        % do the first transition outside of the loop first
        currPath = path;
        nextPhase = currPossibleTransitions(1);
        
        tracks.sequence(currPath,1:point) = tracks.sequence(path,1:point);
        tracks.sequence(currPath,point+1) = nextPhase;
        
        tempProbability = tracks.probability(path)*P_deltaTModel(currPhase,nextPhase);
        
        for transition = 2:currNumPossibleTransitions
            
            currPath = numPaths_new+transition-1;
            
            nextPhase = currPossibleTransitions(transition);
            
            tracks.sequence(currPath,1:point) = tracks.sequence(path,1:point);
            tracks.sequence(currPath,point+1) = nextPhase;
            
            tracks.probability(currPath) = tracks.probability(path)*P_deltaTModel(currPhase,nextPhase);
            
        end
        
        tracks.probability(path) = tempProbability;
        
        numPaths_new = numPaths_new+currNumPossibleTransitions-1;
        
    end
    
    numPaths = numPaths_new;
    
    matRad_progress(point,numPoints-1);
end


tracks.probability((numPaths+1):end) = [];
tracks.sequence((numPaths+1):end,:) = [];

[tracks.probability,ind] = sort(tracks.probability,'descend');
tracks.sequence = tracks.sequence(ind,:);


cumProb = cumsum(tracks.probability);



P_pt = zeros(nPhases+2,numPoints);
P_pt_fromPModel = zeros(nPhases+2,numPoints);
initPVec = zeros(1,nPhases+2);
initPVec(startingPhase) = 1;
for t= 1:numPoints
    for p = 1:nPhases+2
        
        P_pt(p,t) = sum(tracks.probability(tracks.sequence(:,t) == p));
        
        P_pt_fromPModel(:,t) = (initPVec*P_deltaTModel^(t-1))';
        
    end
end
