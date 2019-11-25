load lungPatient2_3mm5p_rep

% meta information for treatment plan

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

% dose calculation settings
pln.propDoseCalc.marginOptions.addMargin    = true; % margin around targets for determining which bixels to calculate dose
pln.propDoseCalc.marginOptions.margin.x     = 10; % margin size in mm
pln.propDoseCalc.marginOptions.margin.y     = 10; % margin size in mm
pln.propDoseCalc.marginOptions.margin.z     = 10; % margin size in mm
pln.propDoseCalc.memorySaverPhoton          = false;
pln.propDoseCalc.vmc                        = true;
pln.propDoseCalc.vmcOptions.keepError       = true;
pln.propDoseCalc.vmcOptions.source          = 'phsp';
pln.propDoseCalc.vmcOptions.phspBaseName    = '5cmx5cm_SSD50cm';
pln.propDoseCalc.vmcOptions.SCD             = 500;
pln.propDoseCalc.vmcOptions.dumpDose        = 1;
pln.propDoseCalc.vmcOptions.version         = 'Carleton';
pln.propDoseCalc.vmcOptions.nCasePerBixel   = 500;
pln.propDoseCalc.vmcOptions.numOfParMCSim   = 16;

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
pln.propOpt.VMAToptions.continuousAperture = true;

pln.propOpt.VMAToptions.startingAngle               = -180;
pln.propOpt.VMAToptions.finishingAngle              = 180;
pln.propOpt.VMAToptions.maxFluGantryAngleSpacing    = 4;    % Max gantry angle spacing for fluence calculation
pln.propOpt.VMAToptions.maxGantryAngleSpacing       = 4;    % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing    = 8;    % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing    = 32;   % Max gantry angle spacing for FMO

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;
pln.propOpt.prop4D.singlePhaseFMO = true;
% multi-phase FMO hasn't been implemented fully (would have to do changes in FMO and leaf
% sequencing - probably better only for fluence, not DAO).

pln = matRad_VMATGantryAngles(pln,cst,ct);

% load results
load('CO.mat');

oldDir = pwd;

%%

maxFluGantryAngleSpacingS = [4 2 1 0.5 0.25];

for maxFluGantryAngleSpacing = maxFluGantryAngleSpacingS
    % calculate dose at original fluence gantry spacing, work our way down
    % to smaller spacings to get a more accurate dose calculation
    
    recalc.pln = pln;
    recalc.pln.propOpt.VMAToptions.maxFluGantryAngleSpacing = maxFluGantryAngleSpacing;
    
    fname = sprintf('Max fluence gantry angle spacing = %.3f.mat',maxFluGantryAngleSpacing);
    fprintf('%s\n',fname);
    recalc.continuousAperture = true;
    recalc.interpNew = true;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo,false,dij,'flu');
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    
end