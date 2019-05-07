load lungPatient0_3DVMAT.mat

% meta information for treatment plan

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

% dose calculation settings
pln.propDoseCalc.memorySaverPhoton          = false;
pln.propDoseCalc.vmc                        = true;
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

% load results
load('Results.mat');

angularResS = [0.5 1 2 4 8];

oldDir = pwd;

for angularRes = angularResS
    %for each angular resolution, proceed from the best approximation to
    %the worst
    
    recalc.pln = pln;
    recalc.pln.propOpt.VMAToptions.maxGantryAngleSpacing = angularRes;
    
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('%.1f degrees, dyn + interp.mat',angularRes);
    fprintf('%s\n',fname);
    recalc.continuousAperture = true;
    recalc.interpNew = true;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    
    %{
    %next, do dynamic fluence and interpolation, but using old dij matrices
    fname = sprintf('%.1f degrees, dyn + interp oldDij.mat',angularRes);
    fprintf('%s\n',fname);
    recalc.continuousAperture = true;
    recalc.interpNew = true;
    recalc.dijNew = false;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    %}
    
    %{
    %NOT SURE IT MAKES SENSE TO EVER DO THIS
    %next, do dynamic fluence but no interpolation
    fname = sprintf('%.1f degrees, dyn + Ninterp.mat',angularRes);
    recalc.continuousAperture = true;
    recalc.interpNew = false;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    save(fname,'resultGUI','recalc');
    %}
    
    
    %next, do interpolation but no dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + interp.mat',angularRes);
    fprintf('%s\n',fname);
    recalc.continuousAperture = false;
    recalc.interpNew = true;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    
    %{
    %finally, do neither interpolation nor dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + Ninterp.mat',angularRes);
    fprintf('%s\n',fname);
    recalc.continuousAperture = false;
    recalc.interpNew = false;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    cd(oldDir);
    save(fname,'resultGUI','recalc');
    %}
    
end