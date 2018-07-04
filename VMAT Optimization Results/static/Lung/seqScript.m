% prep

%clearvars -except *dir
%close all

% meta information for treatment plan
%{
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

pln.propOpt.run4D = false;
pln.propOpt.prop4D.singlePhaseFMO = false;

pln = matRad_VMATGantryAngles(pln,cst,ct);

% stf
stf = matRad_generateStf(ct,cst,pln);

% calc Dij
t0 = tic;
dij = matRad_calcPhotonDose_CELL(ct,stf,pln,cst);
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
%}
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
%save(fname,'resultGUI','recalc','-v7.3')

tOpt = tFMO+tSeq+tDAO;

save('timings','t*')



