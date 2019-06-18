% we want to find a sigma for each gantry angle... put this in the analyze
% script?

load lungPatient0_5mm_rep

% load two calcs of same plan
fname = sprintf('8.0 degrees, dyn + interp.mat');
load(fname)
dose1 = recalc.resultGUI.physicalDose;
doseError1 = recalc.resultGUI.physicalDoseError;

fname = sprintf('8.0 degrees, dyn + interp, repeat.mat');
load(fname)
dose2 = recalc.resultGUI.physicalDose;
doseError2 = recalc.resultGUI.physicalDoseError;

% adjust overlap priorities
cst_Over = matRad_setOverlapPriorities(cst);

V_TargAndNorm = false(size(dose1));
% adjust objectives _and_ constraints internally for fractionation
for i = 1:size(cst_Over,1)
    for j = 1:size(cst_Over{i,6},1)
        cst_Over{i,6}(j).dose = cst_Over{i,6}(j).dose/pln.numOfFractions;
    end
    if ~isempty(cst_Over{i,6}) && ~strcmp(cst_Over{i,2},'BODY') && ~strcmp(cst_Over{i,2},'External')
        [x, y, z] = ind2sub(size(dose1),cst_Over{i,4}{1});
        for k = 1:numel(x)
            V_TargAndNorm(x(k),y(k),z(k)) = true;
        end
    end
end

% setup optimization
global D;
x0 = 1;

% delete all voxels with dose less than 50% of dMax
dThresh = 0.25*max([dose1(:); dose2(:)]);
deleteInd = dose1 <= dThresh | dose2 <= dThresh;

% calculate deltas
doseDiff = dose1-dose2;
doseDiffError = sqrt(doseError1.^2+doseError2.^2);
delta = doseDiff./doseDiffError;
delta(doseDiffError == 0) = 0;
delta(deleteInd) = [];
delta = delta(:);

D = delta;

% perform fit to K-F
%x = fmincon(@mll_KF_noDelta,x0,[],[],[],[],lb,ub);
x = fminsearch(@mll_KF_noDelta,x0);

% extract parameters
sigma = x(1);
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
pdf_fit         = pdf_KF(delta_fit,sigma,alpha1,delta1,alpha2,delta2);

plot(delta_fit,pdf_fit,'r-')

xlabel('\Delta')
ylabel('pdf')
title(sprintf('8^\\circ resolution, RNG seed change only, \\sigma = %f',sigma))


