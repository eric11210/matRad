load lungPatient0_5mm_rep

% load two calcs of same plan
fname = sprintf('8.0 degrees, dyn + interp.mat');
load(fname)
refDose1 = recalc.resultGUI.physicalDose;
refDoseError1 = recalc.resultGUI.physicalDoseError;

fname = sprintf('8.0 degrees, dyn + interp, repeat.mat');
load(fname)
refDose2 = recalc.resultGUI.physicalDose;
refDoseError2 = recalc.resultGUI.physicalDoseError;

% adjust overlap priorities
cst_Over = matRad_setOverlapPriorities(cst);

V_TargAndNorm = false(size(refDose1));
% adjust objectives _and_ constraints internally for fractionation
for i = 1:size(cst_Over,1)
    for j = 1:size(cst_Over{i,6},1)
        cst_Over{i,6}(j).dose = cst_Over{i,6}(j).dose/pln.numOfFractions;
    end
    if ~isempty(cst_Over{i,6}) && ~strcmp(cst_Over{i,2},'BODY') && ~strcmp(cst_Over{i,2},'External')
        [x, y, z] = ind2sub(size(refDose1),cst_Over{i,4}{1});
        for k = 1:numel(x)
            V_TargAndNorm(x(k),y(k),z(k)) = true;
        end
    end
end

% setup optimization
global D;
x0 = [1     0   0       0   0]';
lb = [0     0   -inf    0   -inf]';
ub = [inf   1   inf     1   inf]';
A = [0      1   0       1   0];
b = 1;

% calculate deltas
doseDiff = dose-refDose;
doseDiffError = sqrt(doseError.^2+refDoseError.^2);
delta = doseDiff./doseDiffError;
delta(doseDiffError == 0) = 0;
delta(deleteInd) = [];
delta = delta(:);

D = delta;

% perform fit to K-F
x = fmincon(@minusLogLikelihood,x0,A,b,[],[],lb,ub);

% extract parameters
sigma   = x(1);
alpha1  = x(2);
delta1  = x(3);
alpha2  = x(4);
delta2  = x(5);

