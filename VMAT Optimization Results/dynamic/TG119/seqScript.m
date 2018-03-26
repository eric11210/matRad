%% prep

clearvars -except *dir
close all

% load patient data, i.e. ct, voi, cst

load TG119_NEW.mat

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
pln.propOpt.VMAToptions.continuousAperture = true;

pln.propOpt.VMAToptions.maxGantryAngleSpacing = 4;      % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 4;      % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 28;      % Max gantry angle spacing for FMO

pln = matRad_VMATGantryAngles(pln,cst,ct);

% stf
stf = matRad_generateStf(ct,cst,pln);

% calc Dij
dij = matRad_calcPhotonDose(ct,stf,pln,cst);


% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);

% DAO
fname = 'Results_NEWCONTOURS';
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

t0_nDij_nJ = tic;
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
t_nDij_nJ = toc(t0_nDij_nJ);
savefig(fname)
save(fname,'resultGUI','-v7.3')

save('timings','t_*')



