%% setup

load lungPatient3_3mm5p_rep

% meta information for treatment plan

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth = 5;

% dose calculation settings
pln.propDoseCalc.marginOptions.addMargin    = true; % margin around targets for determining which bixels to calculate dose
pln.propDoseCalc.marginOptions.margin.x     = 10; % margin size in mm
pln.propDoseCalc.marginOptions.margin.y     = 10; % margin size in mm
pln.propDoseCalc.marginOptions.margin.z     = 10; % margin size in mm
pln.propDoseCalc.memorySaverPhoton          = false;
pln.propDoseCalc.vmc                        = true;
pln.propDoseCalc.vmcOptions.keepError       = false;
pln.propDoseCalc.vmcOptions.source          = 'phsp';
pln.propDoseCalc.vmcOptions.phspBaseName    = '5cmx5cm_SSD50cm';
pln.propDoseCalc.vmcOptions.SCD             = 500;
pln.propDoseCalc.vmcOptions.dumpDose        = 1;
pln.propDoseCalc.vmcOptions.version         = 'Carleton';
pln.propDoseCalc.vmcOptions.nCasePerBixel   = 500;
pln.propDoseCalc.vmcOptions.numOfParMCSim   = 16;
pln.propDoseCalc.sampleTargetProb           = 1;

% optimization settings
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runVMAT         = true;
pln.propOpt.runDAO          = true;
pln.propOpt.runSequencing   = true;
pln.propOpt.preconditioner  = true;
pln.propOpt.numLevels       = 7;

pln.propOpt.VMAToptions.machineConstraintFile   = [pln.radiationMode '_' pln.machine];
pln.propOpt.VMAToptions.continuousAperture      = true;
pln.propOpt.VMAToptions.fixedGantrySpeed        = true;
pln.propOpt.VMAToptions.deliveryTime            = 70;

pln.propOpt.VMAToptions.startingAngle               = -180;
pln.propOpt.VMAToptions.finishingAngle              = 180;
pln.propOpt.VMAToptions.maxGantryAngleSpacing       = 4;      % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxFluGantryAngleSpacing    = 4;    % Max gantry angle spacing for fluence calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing    = 8;      % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing    = 32;      % Max gantry angle spacing for FMO

pln.propOpt.run4D                   = true;
pln.propOpt.varOpt                  = false;
pln.propOpt.prop4D.singlePhaseFMO   = true;
% multi-phase FMO hasn't been implemented fully (would have to do changes in FMO and leaf
% sequencing - probably better only for fluence, not DAO).

pln = matRad_VMATGantryAngles(pln,cst,ct);

numOfFractions = pln.numOfFractions;

% absolute path to figure folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

%% CO

load('CO')

V66_5       = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(1,:),66);
V66_50      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(3,:),66);
V66_95      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(5,:),66);
V66_mean    = interp1(pln.numOfFractions.*dvh(25).doseGrid,dvh(25).volumePoints(:),66);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D98_5           = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D98_50          = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D98_95          = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D98_mean        = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),98);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D2_5            = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D2_50           = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D2_95           = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D2_mean         = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),2);

fprintf('\n\nCO\n\n')

fprintf('V66: %.1f (%.1f, %.1f, %.1f)\n',V66_mean,V66_5,V66_50,V66_95);
fprintf('D98: %.1f (%.1f, %.1f, %.1f)\n',D98_mean,D98_5,D98_50,D98_95);
fprintf('D2 : %.1f (%.1f, %.1f, %.1f)\n',D2_mean,D2_5,D2_50,D2_95);

%% PO

load('PO')

V66_5       = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(1,:),66);
V66_50      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(3,:),66);
V66_95      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(5,:),66);
V66_mean    = interp1(pln.numOfFractions.*dvh(25).doseGrid,dvh(25).volumePoints(:),66);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D98_5           = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D98_50          = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D98_95          = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D98_mean        = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),98);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D2_5            = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D2_50           = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D2_95           = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D2_mean         = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),2);

fprintf('\n\nPO\n\n')

fprintf('V66: %.1f (%.1f, %.1f, %.1f)\n',V66_mean,V66_5,V66_50,V66_95);
fprintf('D98: %.1f (%.1f, %.1f, %.1f)\n',D98_mean,D98_5,D98_50,D98_95);
fprintf('D2 : %.1f (%.1f, %.1f, %.1f)\n',D2_mean,D2_5,D2_50,D2_95);

%% STO

load('STO')

V66_5       = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(1,:),66);
V66_50      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(3,:),66);
V66_95      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(5,:),66);
V66_mean    = interp1(pln.numOfFractions.*dvh(25).doseGrid,dvh(25).volumePoints(:),66);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D98_5           = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D98_50          = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D98_95          = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D98_mean        = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),98);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D2_5            = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D2_50           = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D2_95           = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D2_mean         = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),2);

fprintf('\n\nSTO\n\n')

fprintf('V66: %.1f (%.1f, %.1f, %.1f)\n',V66_mean,V66_5,V66_50,V66_95);
fprintf('D98: %.1f (%.1f, %.1f, %.1f)\n',D98_mean,D98_5,D98_50,D98_95);
fprintf('D2 : %.1f (%.1f, %.1f, %.1f)\n',D2_mean,D2_5,D2_50,D2_95);

%% DAD

load('DAD')

V66_5       = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(1,:),66);
V66_50      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(3,:),66);
V66_95      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(5,:),66);
V66_mean    = interp1(pln.numOfFractions.*dvh(25).doseGrid,dvh(25).volumePoints(:),66);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D98_5           = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D98_50          = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D98_95          = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D98_mean        = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),98);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D2_5            = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D2_50           = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D2_95           = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D2_mean         = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),2);

fprintf('\n\nDAD\n\n')

fprintf('V66: %.1f (%.1f, %.1f, %.1f)\n',V66_mean,V66_5,V66_50,V66_95);
fprintf('D98: %.1f (%.1f, %.1f, %.1f)\n',D98_mean,D98_5,D98_50,D98_95);
fprintf('D2 : %.1f (%.1f, %.1f, %.1f)\n',D2_mean,D2_5,D2_50,D2_95);

%% 3DITV

load('3DITV')

V66_5       = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(1,:),66);
V66_50      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(3,:),66);
V66_95      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(5,:),66);
V66_mean    = interp1(pln.numOfFractions.*dvh(25).doseGrid,dvh(25).volumePoints(:),66);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D98_5           = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D98_50          = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D98_95          = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D98_mean        = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),98);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D2_5            = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D2_50           = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D2_95           = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D2_mean         = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),2);

fprintf('\n\n3DITV\n\n')

fprintf('V66: %.1f (%.1f, %.1f, %.1f)\n',V66_mean,V66_5,V66_50,V66_95);
fprintf('D98: %.1f (%.1f, %.1f, %.1f)\n',D98_mean,D98_5,D98_50,D98_95);
fprintf('D2 : %.1f (%.1f, %.1f, %.1f)\n',D2_mean,D2_5,D2_50,D2_95);

%% 3DCTV

load('3DCTV')

V66_5       = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(1,:),66);
V66_50      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(3,:),66);
V66_95      = interp1(pdvh_MC(25).doseGrid,pdvh_MC(25).volumePoints(5,:),66);
V66_mean    = interp1(pln.numOfFractions.*dvh(25).doseGrid,dvh(25).volumePoints(:),66);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D98_5           = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D98_50          = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D98_95          = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),98);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D98_mean        = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),98);

[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(1,:),'stable');
D2_5            = interp1(pdvh_MC(25).volumePoints(1,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(3,:),'stable');
D2_50           = interp1(pdvh_MC(25).volumePoints(3,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(pdvh_MC(25).volumePoints(5,:),'stable');
D2_95           = interp1(pdvh_MC(25).volumePoints(5,uniqueInd),pdvh_MC(25).doseGrid(uniqueInd),2);
[~,uniqueInd]   = unique(dvh(25).volumePoints(:),'stable');
D2_mean         = interp1(dvh(25).volumePoints(uniqueInd),pln.numOfFractions.*dvh(25).doseGrid(uniqueInd),2);

fprintf('\n\n3DCTV\n\n')

fprintf('V66: %.1f (%.1f, %.1f, %.1f)\n',V66_mean,V66_5,V66_50,V66_95);
fprintf('D98: %.1f (%.1f, %.1f, %.1f)\n',D98_mean,D98_5,D98_50,D98_95);
fprintf('D2 : %.1f (%.1f, %.1f, %.1f)\n',D2_mean,D2_5,D2_50,D2_95);
