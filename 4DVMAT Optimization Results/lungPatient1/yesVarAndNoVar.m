%% setup

load lungPatient1_3mm5p_rep

currentDir = pwd;

% meta information for treatment plan

pln.radiationMode   = 'photons';   % either photons / protons / carbon
pln.machine         = 'Generic';

%pln.numOfFractions  = 30;

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
pln.propOpt.VMAToptions.continuousAperture = true;

pln.propOpt.VMAToptions.startingAngle = -180;
pln.propOpt.VMAToptions.finishingAngle = 180;
pln.propOpt.VMAToptions.maxGantryAngleSpacing = 4;      % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = 8;      % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = 32;      % Max gantry angle spacing for FMO

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;
pln.propOpt.prop4D.singlePhaseFMO = true;
% multi-phase FMO hasn't been implemented fully (would have to do changes in FMO and leaf
% sequencing - probably better only for fluence, not DAO).

pln = matRad_VMATGantryAngles(pln,cst,ct);


stf = matRad_generateStf(ct,cst,pln);


dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);

%% no variance term


pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('CO_FMO')

resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);


resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
cd(currentDir);
savefig('CO_DAO')

cd(currentDir);
save('CO','resultGUI');

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);

save('CO','resultGUI','*dvh*');

clear resultGUI *dvh*


%% yes variance term

pln.propOpt.run4D = true;
pln.propOpt.varOpt = true;

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('PO_FMO')


resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);


resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
cd(currentDir);
savefig('PO_DAO')

[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,50);

cd(currentDir);
save('PO','resultGUI','*dvh*');

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);

save('PO','resultGUI','*dvh*');

clear resultGUI *dvh*

%{
%% 3D optimization on CTV

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('3DCTV_FMO')

pln.propOpt.run4D = false;


resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);


resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
cd(currentDir);
savefig('3DCTV_DAO')

cd(currentDir);
save('3DCTV','resultGUI');

clear resultGUI

%% 3D optimization on ITV

pln.propOpt.run4D = false;
pln.propOpt.varOpt = false;

cst{26,6}       = cst{25,6};
cst{25,6}       = [];
pln.RxStruct    = 26;

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('3DITV_FMO')


resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);


resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
cd(currentDir);
savefig('3DITV_DAO')

cd(currentDir);
save('3DITV','resultGUI');

clear resultGUI
%}

