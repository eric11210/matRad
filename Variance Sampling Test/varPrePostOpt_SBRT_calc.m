%% setup

load lungPatient2_SBRT_3mm5p_rep

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
pln.propDoseCalc.vmcOptions.keepError       = false;
pln.propDoseCalc.vmcOptions.source          = 'phsp';
pln.propDoseCalc.vmcOptions.phspBaseName    = '5cmx5cm_SSD50cm';
pln.propDoseCalc.vmcOptions.SCD             = 500;
pln.propDoseCalc.vmcOptions.dumpDose        = 1;
pln.propDoseCalc.vmcOptions.version         = 'Carleton';
pln.propDoseCalc.vmcOptions.nCasePerBixel   = 500;
pln.propDoseCalc.vmcOptions.numOfParMCSim   = 16;
pln.propDoseCalc.sampleTargetProb           = 1;

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
pln.propOpt.VMAToptions.fixedGantrySpeed = true;
pln.propOpt.VMAToptions.deliveryTime = 300;

pln.propOpt.VMAToptions.startingAngle               = -180;
pln.propOpt.VMAToptions.finishingAngle              = 180;
pln.propOpt.VMAToptions.maxFluGantryAngleSpacing    = 4;    % Max gantry angle spacing for fluence calculation
pln.propOpt.VMAToptions.maxGantryAngleSpacing       = 4;    % Max gantry angle spacing for dose calculation
pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing    = 8;    % Max gantry angle spacing for DAO
pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing    = 32;   % Max gantry angle spacing for FMO

pln.propOpt.run4D = true;
pln.propOpt.varOpt = true;
pln.propOpt.prop4D.singlePhaseFMO = true;
% multi-phase FMO hasn't been implemented fully (would have to do changes in FMO and leaf
% sequencing - probably better only for fluence, not DAO).

pln = matRad_VMATGantryAngles(pln,cst,ct);
stf = matRad_generateStf(ct,cst,pln);

%dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
dij.numOfFractions = pln.numOfFractions;

%% do FMO and sequencing

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);

resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

%% calculate variance, PDVHs before optimization

% setup options
options.numOfScenarios  = resultGUI.apertureInfo.numPhases;
options.bioOpt          = 'none';

% save pre optimization
save('varPrePostOpt_SBRT.mat','resultGUI')

% calculate dose and variance of each voxel before optimization
d_preOpt    = matRad_backProjection(resultGUI.apertureInfo.bixelWeights,dij,options);
dVar_preOpt = matRad_doseVariance(resultGUI.apertureInfo,dij);

% save pre optimization
save('varPrePostOpt_SBRT.mat','resultGUI','*Opt')

% do dvhs
[pdvh_MC_preOpt,dvh_mean_MC_preOpt,dvh_std_MC_preOpt] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh_preOpt,~] = matRad_indicatorWrapper(cst,pln,resultGUI);
close all

% save pre optimization
save('varPrePostOpt_SBRT.mat','resultGUI','*Opt')


%% calculate variance, PDVHs after conventional optimization

% redo leaf sequencing to get aperture library before optimization
pln.propOpt.varOpt = false;
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% redo optimization
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% calculate dose and variance of each voxel after optimization
d_postConvOpt       = matRad_backProjection(resultGUI.apertureInfo.bixelWeights,dij,options);
dVar_postConvOpt    = matRad_doseVariance(resultGUI.apertureInfo,dij);

% save post conventional optimization
save('varPrePostOpt_SBRT.mat','resultGUI','*Opt')

% do dvhs
[pdvh_MC_postConvOpt,dvh_mean_MC_postConvOpt,dvh_std_MC_postConvOpt] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh_postConvOpt,~] = matRad_indicatorWrapper(cst,pln,resultGUI);
close all

% save post conventional optimization
save('varPrePostOpt_SBRT.mat','resultGUI','*Opt')

%% calculate variance, PDVHs after probabilistic optimization

% redo leaf sequencing to get aperture library before optimization
pln.propOpt.varOpt = true;
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% redo optimization
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% calculate dose and variance of each voxel after optimization
d_postProbOpt       = matRad_backProjection(resultGUI.apertureInfo.bixelWeights,dij,options);
dVar_postProbOpt    = matRad_doseVariance(resultGUI.apertureInfo,dij);

% save post conventional optimization
save('varPrePostOpt_SBRT.mat','resultGUI','*Opt')

% do dvhs
[pdvh_MC_postProbOpt,dvh_mean_MC_postProbOpt,dvh_std_MC_postProbOpt] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh_postProbOpt,~] = matRad_indicatorWrapper(cst,pln,resultGUI);
close all

% save post conventional optimization
save('varPrePostOpt_SBRT.mat','resultGUI','*Opt')

clear '*Opt'