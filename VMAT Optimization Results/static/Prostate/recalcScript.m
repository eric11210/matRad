load PROSTATE.mat

% meta information for treatment plan

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

pln.propDoseCalc.memorySaverPhoton = false;

% beam geometry settings
pln.propStf.bixelWidth = 5;


% optimization settings
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runVMAT = true;
pln.propOpt.runDAO = true;
pln.propOpt.runSequencing = true;
pln.propOpt.preconditioner = true;
pln.propOpt.numLevels = 7;

pln.propOpt.VMAToptions.machineConstraintFile = [pln.radiationMode '_' pln.machine];
pln.propOpt.VMAToptions.continuousAperture = false;

pln.propOpt.VMAToptions.maxGantryAngleSpacing = 4;      % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 4;      % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 28;      % Max gantry angle spacing for FMO

pln = matRad_VMATGantryAngles(pln,cst,ct);

% load results
load('Results_NEWCONTOURS.mat');


angularRes = 4;

recalc.pln = pln;
recalc.pln.propOpt.VMAToptions.maxGantryAngleSpacing = angularRes;

recalc.continuousAperture = true;
recalc.interpNew = true;
recalc.dijNew = true;

recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
save('Results_NEWCONTOURS.mat','resultGUI','recalc');