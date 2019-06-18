%% load in CT and setup

load lungPatient0_5mm_rep

% threshold factor
threshFac = 0.25;

%this is the reference plan, the most accurate way of calculating dose
fname = sprintf('0.5 degrees, dyn + interp.mat');
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
angularResS = [0.5 1 2 4 8];

percVErr1_NY = zeros(size(angularResS));
percVErr3_NY = zeros(size(angularResS));
percVErr5_NY = zeros(size(angularResS));
percVErr10_NY = zeros(size(angularResS));
fluence_NY = zeros(size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,2), numel(angularResS));
weight_NY = zeros(size(angularResS));
alpha1_NY = zeros(size(angularResS));
delta1_NY = zeros(size(angularResS));
alpha2_NY = zeros(size(angularResS));
delta2_NY = zeros(size(angularResS));

percVErr1_YY = zeros(size(angularResS));
percVErr3_YY = zeros(size(angularResS));
percVErr5_YY = zeros(size(angularResS));
percVErr10_YY = zeros(size(angularResS));
fluence_YY = zeros(size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,2), numel(angularResS));
weight_YY = zeros(size(angularResS));
alpha1_YY = zeros(size(angularResS));
delta1_YY = zeros(size(angularResS));
alpha2_YY = zeros(size(angularResS));
delta2_YY = zeros(size(angularResS));

numPVHPoints = 1e3;

PVH_YY = zeros(numel(angularResS),numPVHPoints);
PVH_NY = zeros(numel(angularResS),numPVHPoints);

PVHPoints_YY = zeros(numel(angularResS),numPVHPoints);
PVHPoints_NY = zeros(numel(angularResS),numPVHPoints);

sigma = zeros(size(angularResS));


%% determine the correlation factor for each gantry angle resolution

% determine which voxels to delete
deleteInd = refDose < threshFac*max(refDose(:));
%deleteInd = ~V_TargAndNorm;

% setup optimization
global D;
x0 = 1;

i = 1;
for angularRes = angularResS
    
    % load two calcs of same plan
    fname = sprintf('%.1f degrees, dyn + interp.mat',angularRes);
    load(fname)
    dose1 = recalc.resultGUI.physicalDose;
    doseError1 = recalc.resultGUI.physicalDoseError;
    
    fname = sprintf('%.1f degrees, dyn + interp, repeat.mat',angularRes);
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
    title(sprintf('%.1f^\\circ resolution, RNG seed change only, \\sigma = %f',angularRes,sigma(i)))
    savefig(sprintf('Sigma for %.1f degrees.fig',angularRes))
    
    i = i+1;
end

%% determine K-F parameters when comparing doses to reference dose

% setup optimization
x0  = [0    0       0   0]';
lb  = [0    -inf    0   -inf]';
ub  = [1    inf     1   inf]';
A   = [1    0       1   0];
b   = 1;

% scale the reference dose error by sigma(1) to account for correlations
refDoseError = sigma(1).*refDoseError;

i = 1;
for angularRes = angularResS
    %for each angular resolution, proceed from the best approximation to
    %the worst
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('%.1f degrees, dyn + interp.mat',angularRes);
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
    alpha1_YY(i) = x(1);
    delta1_YY(i) = x(2);
    alpha2_YY(i) = x(3);
    delta2_YY(i) = x(4);
    
    % graph pdf
    figure
    subplot(1,3,1)
    histogram(delta,'Normalization','pdf');
    hold on
    
    delta_limits    = xlim;
    delta_fit       = linspace(delta_limits(1),delta_limits(2),1000);
    pdf_fit         = pdf_KF(delta_fit,1,alpha1_YY(i),delta1_YY(i),alpha2_YY(i),delta2_YY(i));
    
    plot(delta_fit,pdf_fit,'r-')
    
    xlabel('\Delta')
    ylabel('pdf')
    title(sprintf('%.1f^\\circ resolution, CA: K-F pdf',angularRes))
    
    % graph systematic differences
    relDiff = doseDiffError./refDose;
    relDiff(deleteInd) = [];
    
    subplot(1,3,2)
    histogram(abs(delta1_YY(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_1 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('%.1f^\\circ resolution, CA, \\alpha_1 = %f, \\delta_1 = %f',angularRes,alpha1_YY(i),delta1_YY(i)))
    
    subplot(1,3,3)
    histogram(abs(delta2_YY(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_2 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('%.1f^\\circ resolution, CA, \\alpha_2 = %f, \\delta_2 = %f',angularRes,alpha2_YY(i),delta2_YY(i)))
    
    savefig(sprintf('K-F for CA, %.1f degrees.fig',angularRes))
    
    
    % percent difference stuff
    percDiff = 100*abs(dose-refDose)./refDose;
    percDiff(refDose == 0) = 0;
    percDiff(deleteInd) = [];
    percDiff = percDiff(:);
    numVox = numel(percDiff);
    
    PVHPoints_YY(i,:) = linspace(0,max(percDiff)*1.05,numPVHPoints);
    for j = 1:numPVHPoints
        PVH_YY(i,j) = sum(percDiff > PVHPoints_YY(i,j));
    end
    PVH_YY(i,:) = 100.*PVH_YY(i,:)./numVox;
    
    percVErr1_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_YY(:,:,i) = fluence_YY(:,:,i)+recalc.apertureInfo.beam(j).shape{1}(1).shapeMap;
        weight_YY(i) = weight_YY(i)+recalc.apertureInfo.beam(j).shape{1}(1).weight;
    end
    
    %next, do interpolation but no dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + interp.mat',angularRes);
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
    alpha1_NY(i) = x(1);
    delta1_NY(i) = x(2);
    alpha2_NY(i) = x(3);
    delta2_NY(i) = x(4);
    
    % graph results
    figure
    subplot(1,3,1)
    histogram(delta,'Normalization','pdf');
    hold on
    
    delta_limits    = xlim;
    delta_fit       = linspace(delta_limits(1),delta_limits(2),1000);
    pdf_fit         = pdf_KF(delta_fit,1,alpha1_NY(i),delta1_NY(i),alpha2_NY(i),delta2_NY(i));
    
    plot(delta_fit,pdf_fit,'r-')
    
    xlabel('\Delta')
    ylabel('pdf')
    title(sprintf('%.1f^\\circ resolution, DA: K-F pdf',angularRes))
    
    % graph systematic differences
    relDiff = doseDiffError./refDose;
    relDiff(deleteInd) = [];
    
    subplot(1,3,2)
    histogram(abs(delta1_NY(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_1 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('%.1f^\\circ resolution, DA, \\alpha_1 = %f, \\delta_1 = %f',angularRes,alpha1_NY(i),delta1_NY(i)))
    
    subplot(1,3,3)
    histogram(abs(delta2_NY(i).*relDiff),'Normalization','pdf')
    xlabel('\delta_2 \cdot RCSU')
    ylabel('pdf')
    title(sprintf('%.1f^\\circ resolution, DA, \\alpha_2 = %f, \\delta_2 = %f',angularRes,alpha2_NY(i),delta2_NY(i)))
    
    savefig(sprintf('K-F for DA, %.1f degrees.fig',angularRes))
    
    
    percDiff = 100*abs(dose-refDose)./refDose;
    percDiff(refDose == 0) = 0;
    percDiff(deleteInd) = [];
    percDiff = percDiff(:);
    numVox = numel(percDiff);
    
    PVHPoints_NY(i,:) = linspace(0,max(percDiff)*1.05,numPVHPoints);
    for j = 1:numPVHPoints
        PVH_NY(i,j) = sum(percDiff > PVHPoints_NY(i,j));
    end
    PVH_NY(i,:) = 100.*PVH_NY(i,:)./numVox;
    
    percVErr1_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_NY(:,:,i) = fluence_NY(:,:,i)+recalc.apertureInfo.beam(j).shape{1}(1).shapeMap;
        weight_NY(i) = weight_NY(i)+recalc.apertureInfo.beam(j).shape{1}(1).weight;
    end
    
    i = i+1;
end


save('Res Results', '*_NY', '*_YY','sigma')
%{
figure
hold
plot(angularResS,percVErr1_YY)
plot(angularResS,percVErr3_YY)
plot(angularResS,percVErr5_YY)
plot(angularResS,percVErr10_YY)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Dynamic, interpolated';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_YY_oldDij)
plot(angularResS,percVErr3_YY_oldDij)
plot(angularResS,percVErr5_YY_oldDij)
plot(angularResS,percVErr10_YY_oldDij)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Dynamic, interpolated, oldDij';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_YY_oldDij_itself)
plot(angularResS,percVErr3_YY_oldDij_itself)
plot(angularResS,percVErr5_YY_oldDij_itself)
plot(angularResS,percVErr10_YY_oldDij_itself)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Dynamic, interpolated, oldDij cf itself';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_NY)
plot(angularResS,percVErr3_NY)
plot(angularResS,percVErr5_NY)
plot(angularResS,percVErr10_NY)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Not dynamic, interpolated';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_NN)
plot(angularResS,percVErr3_NN)
plot(angularResS,percVErr5_NN)
plot(angularResS,percVErr10_NN)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Not dynamic, not interpolated';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_NN_itself)
plot(angularResS,percVErr3_NN_itself)
plot(angularResS,percVErr5_NN_itself)
plot(angularResS,percVErr10_NN_itself)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Not dynamic, not interpolated cf itself';
title(fname)
grid
savefig(fname)
%}
h_YY = figure;
hold on
h_NY = figure;
hold on

for i = 1:numel(angularResS)
    figure(h_YY)
    semilogy(PVHPoints_YY(i,:),PVH_YY(i,:))
    
    figure(h_NY)
    semilogy(PVHPoints_NY(i,:),PVH_NY(i,:))
end

path = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Papers\Dynamic Fluence Calculation\Figures\';

figure(h_YY)
xlabel('Relative dose difference (\%)')
ylabel('Volume (\%)')
ylim([1 100])
fname = 'Lung_Dynamic_interpolated';
%title(fname)
legend({'$\Delta\theta_{\mathrm{dose}} = \SI{0.5}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{1}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{2}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{4}{\degree}$'},'Location','Best')
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [path fname '.tex'];
%matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

figure(h_NY)
xlabel('Relative dose difference (\%)')
ylabel('Volume (\%)')
ylim([1 100])
fname = 'Lung_Notdynamic_interpolated';
%title(fname)
legend({'$\Delta\theta_{\mathrm{dose}} = \SI{0.5}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{1}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{2}{\degree}$','$\Delta\theta_{\mathrm{dose}} = \SI{4}{\degree}$'},'Location','Best')
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [path fname '.tex'];
%matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

