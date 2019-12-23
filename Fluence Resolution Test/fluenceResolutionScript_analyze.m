%% load in CT and setup

load lungPatient2_3mm5p_rep

% threshold factor
threshFac = 0.5;

% set of fluence gantry angle spacings
maxFluGantryAngleSpacingS = [4 2 1 0.5 0.125 0.0625];

% this is the reference plan, the most accurate way of calculating dose
fname = sprintf('SBRT max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacingS(end));
load(fname)
refDose = recalc.resultGUI.physicalDose;
%refDoseError = recalc.resultGUI.physicalDoseError;

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
obj = zeros(size(maxFluGantryAngleSpacingS));

numPVHPoints = 1e3;

PVH = zeros(numel(maxFluGantryAngleSpacingS),numPVHPoints);

PVHPoints = zeros(numel(maxFluGantryAngleSpacingS),numPVHPoints);

%% determine PVHs

% determine which voxels to delete
deleteInd = refDose < threshFac*max(refDose(:));
%deleteInd = ~V_TargAndNorm;

i = 1;
for maxFluGantryAngleSpacing = maxFluGantryAngleSpacingS(1:(end-1))
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('SBRT max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacing);
    load(fname);
    dose = recalc.resultGUI.physicalDose;
    
    % percent difference stuff
    percDiff = 100*abs(dose-refDose)./max(refDose(:));
    percDiff(refDose == 0) = 0;
    percDiff(deleteInd) = [];
    percDiff = percDiff(:);
    numVox = numel(percDiff);
    
    PVHPoints(i,:) = linspace(0,max(percDiff)*1.05,numPVHPoints);
    for j = 1:numPVHPoints
        PVH(i,j) = sum(percDiff > PVHPoints(i,j));
    end
    PVH(i,:) = 100.*PVH(i,:)./numVox;
    
    percVErr1(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & ~deleteInd)./nnz(~deleteInd);
    percVErr3(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & ~deleteInd)./nnz(~deleteInd);
    percVErr5(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & ~deleteInd)./nnz(~deleteInd);
    percVErr10(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & ~deleteInd)./nnz(~deleteInd);
    
    for i_flu = 1:numel(recalc.apertureInfo.beam)
        for phase_flu = 1:recalc.apertureInfo.numPhases
            fluence(:,:,i) = fluence(:,:,i)+recalc.apertureInfo.beam(i_flu).shape{phase_flu}(1).shapeMap;
            weight(i) = weight(i)+recalc.apertureInfo.beam(i_flu).shape{phase_flu}(1).weight;
        end
    end
    
    obj(i) = recalc.f;
    
    i = i+1;
end


save('SBRT Conventional Fluence Res Results','fluence','weight','percVErr*','PVH*')

%%

h = figure;

for i = 1:(numel(maxFluGantryAngleSpacingS)-1)
    figure(h)
    semilogy(PVHPoints(i,:),PVH(i,:))
    hold on
end

path = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Papers\Dynamic Fluence Calculation\Figures\';

figure(h)
xlabel('Relative dose difference (\%)')
ylabel('Volume (\%)')
ylim([0.1 100])
fname = 'SBRT Fluence Res Results';
%title(fname)
legend({'$\Delta\theta_{\mathrm{flu}} = \SI{4}{\degree}$','$\Delta\theta_{\mathrm{flu}} = \SI{2}{\degree}$','$\Delta\theta_{\mathrm{flu}} = \SI{1}{\degree}$','$\Delta\theta_{\mathrm{flu}} = \SI{0.5}{\degree}$','$\Delta\theta_{\mathrm{flu}} = \SI{0.25}{\degree}$','$\Delta\theta_{\mathrm{flu}} = \SI{0.125}{\degree}$'},'Location','Best')
grid on
set(gca,'YMinorTick','on','XMinorTick','on')
savefig(fname)
fullpath = [path fname '.tex'];
%matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=4,minor x tick num=4,width=\linewidth,height=\linewidth','showInfo',false);

