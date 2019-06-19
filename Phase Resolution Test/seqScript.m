% prep

%clearvars -except *dir
%close all

% meta information for treatment plan
load lungPatient0_5mm_rep

cst{9,6}.penalty    = 50;
cst{25,6}.penalty   = 2000;

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
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 32;      % Max gantry angle spacing for FMO

pln.propOpt.run4D = false;
pln.propOpt.varOpt = false;
pln.propOpt.prop4D.singlePhaseFMO = true;
% multi-phase FMO hasn't been implemented fully (would have to do changes in FMO and leaf
% sequencing - probably better only for fluence, not DAO).

pln = matRad_VMATGantryAngles(pln,cst,ct);

% stf
stf = matRad_generateStf(ct,cst,pln);

% calc Dij
t0 = tic;
dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
tDij = toc(t0);

% inverse planning for imrt
t0 = tic;
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
tFMO = toc(t0);

% DAO
fname = 'Results';
t0 = tic;
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);
tSeq = toc(t0);

t0 = tic;
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
tDAO = toc(t0);
savefig(fname)

% recalculation
%{
angularRes = 4;

recalc.pln = pln;
recalc.pln.propOpt.VMAToptions.maxGantryAngleSpacing = angularRes;

recalc.continuousAperture = true;
recalc.interpNew = true;
recalc.dijNew = true;

recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo,0,dij);
%}
save(fname,'resultGUI','-v7.3')

tOpt = tFMO+tSeq+tDAO;

save('timings','t*')



