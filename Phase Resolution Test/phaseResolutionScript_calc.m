% load results
load('Results.mat');

% prep for loop
numPhases_vec   = 1;
oldDir          = pwd;
ctName_base     = 'lungPatient0_5mm';

% these parameters don't change throughout the loop
recalc.continuousAperture = true;
recalc.interpNew = true;
recalc.dijNew = true;

for numPhases = numPhases_vec
    %% prep
    % load ct, cst, pln
    ctName = sprintf('%s%dp_rep.mat',ctName_base,numPhases);
    load(ctName,'ct','cst','pln');
    
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
    pln.propDoseCalc.vmcOptions.nCasePerBixel   = 2000;
    pln.propDoseCalc.vmcOptions.numOfParMCSim   = 128;
    
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
    
    pln.propOpt.VMAToptions.startingAngle = -180;
    pln.propOpt.VMAToptions.finishingAngle = 180;
    pln.propOpt.VMAToptions.maxGantryAngleSpacing = 8;      % Max gantry angle spacing for dose calculation
    pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 8;      % Max gantry angle spacing for DAO
    pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 40;      % Max gantry angle spacing for FMO
    
    pln.propOpt.run4D = false;
    pln.propOpt.prop4D.singlePhaseFMO = false;
    % multi-phase FMO hasn't been implemented fully (would have to do changes in FMO and leaf
    % sequencing - probably better only for fluence, not DAO).
    
    pln = matRad_VMATGantryAngles(pln,cst,ct);
    
    recalc.pln                                              = pln;
    recalc.pln.propOpt.run4D                                = true;
    recalc.pln.propOpt.VMAToptions.maxGantryAngleSpacing    = 4;
    recalc.numPhases    = numPhases;
    
    %% recalc
    
    % determine name
    fname = sprintf('%d phases.mat',numPhases);
    fprintf('%s\n',fname);
    
    % do recalc, save
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    
    % do it again for correlation
    
    % determine name
    fname = sprintf('%d phases, repeat.mat',numPhases);
    fprintf('%s\n',fname);
    
    % do recalc, save
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    
end