%% load in CT and setup

load lungPatient0_5mm_rep

% threshold factor
threshFac = 0.25;

% this is the reference plan, the most accurate way of calculating dose
fname = sprintf('10 phases CH, repeat.mat');
load(fname)
refDose = recalc.resultGUI.physicalDose;
refDoseError = recalc.resultGUI.physicalDoseError;

% adjust overlap priorities
cst_Over = matRad_setOverlapPriorities(cst);

V_TargAndNorm = false(size(refDose));
% adjust objectives _and_ constraints internally for fractionation
for i = 1:size(cst_Over,1)
    for j = 1:size(cst_Over{i,6},1)
        cst_Over{i,6}(j).dose = cst_Over{i,6}(j).dose/pln.numOfFractions;
    end
    if ~isempty(cst_Over{i,6}) && ~strcmp(cst_Over{i,2},'BODY') && ~strcmp(cst_Over{i,2},'External')
        [x, y, z] = ind2sub(size(refDose),cst_Over{i,4}{1});
        for k = 1:numel(x)
            V_TargAndNorm(x(k),y(k),z(k)) = true;
        end
    end
end

%percentage of the volume with at least a x% error relative to the
%reference dose
numPhases_vec = 1:10;

percVErr1 = zeros(size(numPhases_vec));
percVErr3 = zeros(size(numPhases_vec));
percVErr5 = zeros(size(numPhases_vec));
percVErr10 = zeros(size(numPhases_vec));
fluence = zeros(size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,2), numel(numPhases_vec));
weight = zeros(size(numPhases_vec));
alpha1 = zeros(size(numPhases_vec));
delta1 = zeros(size(numPhases_vec));
alpha2 = zeros(size(numPhases_vec));
delta2 = zeros(size(numPhases_vec));

numPVHPoints = 1e3;

PVH = zeros(numel(numPhases_vec),numPVHPoints);

PVHPoints = zeros(numel(numPhases_vec),numPVHPoints);

sigma = zeros(size(numPhases_vec));


%% determine the correlation factor for each gantry angle resolution

% determine which voxels to delete
deleteInd = refDose < threshFac*max(refDose(:));
%deleteInd = ~V_TargAndNorm;

% setup optimization
global D;
x0 = 1;

i = 1;
for numPhases = numPhases_vec
    
    % load two calcs of same plan
    fname = sprintf('%d phases CH',numPhases);
    load(fname)
    dose1 = recalc.resultGUI.physicalDose;
    doseError1 = recalc.resultGUI.physicalDoseError;
    
    fname = sprintf('%d phases CH, repeat.mat',numPhases);
    load(fname)
    dose2 = recalc.resultGUI.physicalDose;
    doseError2 = recalc.resultGUI.physicalDoseError;
    
    % calculate deltas
    doseDiff = dose1-dose2;
    doseDiffError = sqrt(doseError1.^2+doseError2.^2);
    delta = doseDiff./doseDiffError;
    delta(doseDiffError == 0) = 0;
    delta(deleteInd) = [];
    delta = delta(:);
    
    D = delta;
    
    % perform fit to K-F
    x = fminsearch(@mll_KF_noDelta,x0);
    
    % extract parameters
    sigma(i) = x(1);
    alpha1 = 0;
    delta1 = 0;
    alpha2 = 0;
    delta2 = 0;
    
    % graph results
    figure
    histogram(delta,'Normalization','pdf');
    hold on
    
    delta_limits    = xlim;
    delta_fit       = linspace(delta_limits(1),delta_limits(2),1000);
    pdf_fit         = pdf_KF(delta_fit,sigma(i),alpha1,delta1,alpha2,delta2);
    
    plot(delta_fit,pdf_fit,'r-')
    
    xlabel('\Delta')
    ylabel('pdf')
    title(sprintf('%d phases, RNG seed change only, \\sigma = %f',numPhases,sigma(i)))
    savefig(sprintf('Sigma for %d phases CH.fig',numPhases))
    
    i = i+1;
end

%% determine K-F parameters when comparing doses to reference dose

% setup optimization
x0  = [0    0       0   0]';
lb  = [0    -inf    0   -inf]';
ub  = [1    inf     1   inf]';
A   = [1    0       1   0];
b   = 1;

% scale the reference dose error by sigma(end) to account for correlations
refDoseError = sigma(end).*refDoseError;

i = 1;
for numPhases = numPhases_vec
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('%d phases CH.mat',numPhases);
    load(fname);
    dose = recalc.resultGUI.physicalDose;
    % scale the dose error by sigma(i) to account for correlations
    doseError = sigma(i).*recalc.resultGUI.physicalDoseError;
    
    doseDiff = dose-refDose;
    doseDiffError = sqrt(doseError.^2+refDoseError.^2);
    delta = doseDiff./doseDiffError;
    delta(doseDiffError == 0) = 0;
    delta(deleteInd) = [];
    delta = delta(:);
    
    D = delta;
    
    x = fmincon(@mll_KF_noSigma,x0,A,b,[],[],lb,ub);
    alpha1(i) = x(1);
    delta1(i) = x(2);
    alpha2(i) = x(3);
    delta2(i) = x(4);
    
    % graph pdf
    figure
    subplot(1,3,1)
    histogram(delta,'Normalization','pdf');
    hold on
    
    delta_limits    = xlim;
    delta_fit       = linspace(delta_limits(1),delta_limits(2),1000);
    pdf_fit         = pdf_KF(delta_fit,1,alpha1(i),delta1(i),alpha2(i),delta2(i));
    
    plot(delta_fit,pdf_fit,'r-')
    
    xlabel('\Delta')
    ylabel('pdf')
    title(sprintf('%d phases: K-F pdf',numPhases))
    
    % graph systematic differences
    relDiff = doseDiffError./refDose;
    relDiff(deleteInd) = [];
    
    subplot(1,3,2)
    histogram(abs(delta1(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_1 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('%d phases, \\alpha_1 = %f, \\delta_1 = %f',numPhases,alpha1(i),delta1(i)))
    
    subplot(1,3,3)
    histogram(abs(delta2(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_2 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('%d phases, \\alpha_2 = %f, \\delta_2 = %f',numPhases,alpha2(i),delta2(i)))
    
    savefig(sprintf('K-F for %d phases CH.fig',numPhases))
    
    
    % percent difference stuff
    percDiff = 100*abs(dose-refDose)./refDose;
    percDiff(refDose == 0) = 0;
    percDiff(deleteInd) = [];
    percDiff = percDiff(:);
    numVox = numel(percDiff);
    
    PVHPoints(i,:) = linspace(0,max(percDiff)*1.05,numPVHPoints);
    for j = 1:numPVHPoints
        PVH(i,j) = sum(percDiff > PVHPoints(i,j));
    end
    PVH(i,:) = 100.*PVH(i,:)./numVox;
    
    percVErr1(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    
    for i_flu = 1:numel(recalc.apertureInfo.beam)
        for phase_flu = 1:recalc.apertureInfo.numPhases
            fluence(:,:,i) = fluence(:,:,i)+recalc.apertureInfo.beam(i_flu).shape{phase_flu}(1).shapeMap;
            weight(i) = weight(i)+recalc.apertureInfo.beam(i_flu).shape{phase_flu}(1).weight;
        end
    end
    
    i = i+1;
end


save('Phase Res Results CH')
%{
h = figure;
hold on

for i = 1:numel(numPhases_vec)
    figure(h)
    semilogy(PVHPoints(i,:),PVH(i,:))

end

path = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Papers\Dynamic Fluence Calculation\Figures\';

figure(h)
xlabel('Relative dose difference (\%)')
ylabel('Volume (\%)')
ylim([1 100])
fname = 'Lung';
%title(fname)
legend({'$\Delta\theta_{\mathrm{dose}} = \SI{0.5}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{1}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{2}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{4}{\degree}$'},'Location','Best')
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [path fname '.tex'];
%matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);
%}
