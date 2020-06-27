%% load in CT and setup

%load lungPatient0_5mm5p_rep

% threshold factor
threshFac = 0.5;

% set of fluence gantry angle spacings
maxFluGantryAngleSpacingS = pln.propOpt.VMAToptions.fluGantryAngleSpacing./[1 round(pln.propOpt.prop4D.motionModel.deltaT_sample/0.04)];

% this is the reference plan, the most accurate way of calculating dose
fname = sprintf('convFrac max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacingS(end));
load(fname)
refDose = recalc.resultGUI.dMean_MC;
refDoseError = sqrt(recalc.resultGUI.dVar_MC);

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

percVErr1 = zeros(size(maxFluGantryAngleSpacingS));
percVErr3 = zeros(size(maxFluGantryAngleSpacingS));
percVErr5 = zeros(size(maxFluGantryAngleSpacingS));
percVErr10 = zeros(size(maxFluGantryAngleSpacingS));
fluence = zeros(size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,2), numel(maxFluGantryAngleSpacingS));
weight = zeros(size(maxFluGantryAngleSpacingS));
alpha1 = zeros(size(maxFluGantryAngleSpacingS));
delta1 = zeros(size(maxFluGantryAngleSpacingS));
alpha2 = zeros(size(maxFluGantryAngleSpacingS));
delta2 = zeros(size(maxFluGantryAngleSpacingS));

numPVHPoints = 1e3;

PVH = zeros(numel(maxFluGantryAngleSpacingS),numPVHPoints);

PVHPoints = zeros(numel(maxFluGantryAngleSpacingS),numPVHPoints);

sigma = zeros(size(maxFluGantryAngleSpacingS));


%% determine the correlation factor for each gantry angle resolution

% determine which voxels to delete
deleteInd = refDose < threshFac*max(refDose(:));
%deleteInd = ~V_TargAndNorm;

% setup optimization
global D;
x0 = 1;

i = 1;
for maxFluGantryAngleSpacing = maxFluGantryAngleSpacingS
    
    % load two calcs of same plan
    fname = sprintf('convFrac max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacing);
    load(fname)
    if isfield(recalc.resultGUI,'dMean_MC')
        dose1 = recalc.resultGUI.dMean_MC;
        doseError1 = sqrt(recalc.resultGUI.dVar_MC);
    else
        dose1 = recalc.resultGUI.physicalDose(:);
        doseError1 = zeros(size(dose1));
    end
    
    fname = sprintf('convFrac repeat max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacing);
    load(fname)
    if isfield(recalc.resultGUI,'dMean_MC')
        dose2 = recalc.resultGUI.dMean_MC;
        doseError2 = sqrt(recalc.resultGUI.dVar_MC);
    else
        dose2 = recalc.resultGUI.physicalDose(:);
        doseError2 = zeros(size(dose2));
    end
    
    
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
    title(sprintf('\\Delta\\theta_{flu} = %f, RNG seed change only, \\sigma = %f',maxFluGantryAngleSpacing,sigma(i)))
    savefig(sprintf('Sigma for flu gantry angle spacing = %f.fig',maxFluGantryAngleSpacing))
    
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
for maxFluGantryAngleSpacing = maxFluGantryAngleSpacingS(1:end-1)
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('convFrac max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacing);
    load(fname);
    if isfield(recalc.resultGUI,'dMean_MC')
    dose = recalc.resultGUI.dMean_MC;
    % scale the dose error by sigma(i) to account for correlations
    doseError = sigma(i).*sqrt(recalc.resultGUI.dVar_MC);
    else
        dose = recalc.resultGUI.physicalDose(:);
        doseError = zeros(size(dose1));
    end
    
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
    title(sprintf('\\Delta\\theta_{flu} = %f: K-F pdf',maxFluGantryAngleSpacing))
    
    % graph systematic differences
    relDiff = doseDiffError./refDose;
    relDiff(deleteInd) = [];
    
    subplot(1,3,2)
    histogram(abs(delta1(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_1 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('\\Delta\\theta_{flu} = %f, \\alpha_1 = %f, \\delta_1 = %f',maxFluGantryAngleSpacing,alpha1(i),delta1(i)))
    
    subplot(1,3,3)
    histogram(abs(delta2(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_2 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('\\Delta\\theta_{flu} = %f, \\alpha_2 = %f, \\delta_2 = %f',maxFluGantryAngleSpacing,alpha2(i),delta2(i)))
    
    savefig(sprintf('K-F for flu gantry angle spacing = %f.fig',maxFluGantryAngleSpacing))
    
    
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


save('Fluence Res Results')
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
