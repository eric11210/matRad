close all
clear variables
load data

deltaT = mean(diff(t));
deltaTSamp = deltaT./10;
tSamp = min(t):deltaTSamp:max(t);
tumourPosSamp = interp1(t,tumourPos,tSamp);

nPhases = 10;

[~,components] = pca(tumourPosSamp);
x_XCAT = components(:,1);

lBoundsMax = max(x_XCAT);
lBoundsMin = min(x_XCAT);
nBins = nPhases/2;
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
[~,ind_maxPeaks] = findpeaks(x_XCAT,'MinPeakProminence',0.5);
[~,ind_minPeaks] = findpeaks(-x_XCAT,'MinPeakProminence',0.5);

l_XCAT(ind_minPeaks:end) = nPhases+1-l_XCAT(ind_minPeaks:end);

% right
xl_XCAT = zeros(nPhases,1);
t_xl_XCAT = zeros(nPhases,1);
t_lBounds = zeros(nPhases+1,1);

% wrong
tl = zeros(nPhases,1);


for i = 1:nPhases
    
    xl_XCAT(i) = mean(x_XCAT(l_XCAT == i));
    t_xl_XCAT(i) = interp1(x_XCAT(l_XCAT == i),tSamp(l_XCAT == i),xl_XCAT(i));
    
    tl(i) = mean(tSamp(l_XCAT == i));
end

t_lBounds(1) = 0;
for i = 2:numel(lBounds)
    
    prev_l = i-1;
    next_l = i;
    t_lBounds(i) = interp1(x_XCAT(l_XCAT == prev_l | l_XCAT == next_l),tSamp(l_XCAT == prev_l | l_XCAT == next_l),lBounds(i));
    
    prev_l = nPhases+1-i;
    next_l = nPhases+2-i;
    t_lBounds(nPhases+2-i) = interp1(x_XCAT(l_XCAT == prev_l | l_XCAT == next_l),tSamp(l_XCAT == prev_l | l_XCAT == next_l),lBounds(i));
end
t_lBounds(end) = max(tSamp);

x_XCAT_tl = interp1(tSamp,x_XCAT,tl);

figure
plot(tSamp,x_XCAT)
hold
plot(t_xl_XCAT,xl_XCAT,'k.')
plot(tl,x_XCAT_tl,'r.')
plot(tSamp,repmat(lBounds,numel(tSamp),1),'r--')
plot(repmat(t_lBounds',numel(x_XCAT),1),x_XCAT,'r--')
xlabel('time / s')
ylabel('principal component')
legend({'XCAT motion' 'average position in bin' 'average time in bin' 'bin edges'},'location','best')
title('position binning')


% now try time binning
nBoundsMax = max(tSamp);
nBoundsMin = min(tSamp);
nBounds = linspace(nBoundsMin,nBoundsMax,nPhases+1);

n_XCAT = zeros(size(x_XCAT));
for phase = 1:nPhases
    n_XCAT(nBounds(phase) <= tSamp & tSamp <= nBounds(phase+1)) = phase;
end

% right
xn_XCAT = zeros(nPhases,1);
t_xn_XCAT = zeros(nPhases,1);

% wrong
tn = zeros(nPhases,1);
for i = 1:nPhases
    
    xn_XCAT(i) = mean(x_XCAT(n_XCAT == i));
    t_xn_XCAT(i) = interp1(x_XCAT(n_XCAT == i),tSamp(n_XCAT == i),xn_XCAT(i));
    
    tn(i) = mean(tSamp(n_XCAT == i));
end

x_XCAT_tn = interp1(tSamp,x_XCAT,tn);

figure
plot(tSamp,x_XCAT)
hold
plot(t_xn_XCAT,xn_XCAT,'k.')
plot(tn,x_XCAT_tn,'r.')
plot(repmat(nBounds,numel(x_XCAT),1),x_XCAT,'r--')
xlabel('time / s')
ylabel('principal component')
legend({'XCAT motion' 'average position in bin' 'average time in bin' 'bin edges'},'location','best')
title('time binning')