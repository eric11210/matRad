%% setup

load lungPatient0_5mm5p_rep

currentDir = pwd;

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

stf = matRad_generateStf(ct,cst,pln);

%% calculate (load) dij

%dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
dij = matRad_loadDij('lungPatient0_5mm5p_rep');
dij.numOfFractions = pln.numOfFractions;

%% conventional optimization

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('CO_FMO')
close all

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('CO_DAO')
save('CO','resultGUI');
close all

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('CO','resultGUI','*dvh*');

% now do dvhs for single fraction
resultGUI.physicalDose  = resultGUI.physicalDose./pln.numOfFractions;
numOfFractions          = pln.numOfFractions;
pln.numOfFractions      = 1;

[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('CO_oneFrac','resultGUI','*dvh*');

clear resultGUI *dvh*

% reset number of fractions
pln.numOfFractions = numOfFractions;


%% 3D optimization on CTV

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('3DCTV_FMO')
close all

% turn off 4d
pln.propOpt.run4D = false;

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('3DCTV_DAO')
save('3DCTV','resultGUI');
close all

% convert sequence to library
resultGUI.apertureInfo = matRad_apertures2Library(resultGUI.apertureInfo,pln,stf,dij.numPhases);

% change initProb
resultGUI.apertureInfo.motionModel.initProb = resultGUI.apertureInfo.motionModel.Pi_deltaTSample';

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('3DCTV','resultGUI','*dvh*');

% now do dvhs for single fraction
resultGUI.physicalDose  = resultGUI.physicalDose./pln.numOfFractions;
numOfFractions          = pln.numOfFractions;
pln.numOfFractions      = 1;

[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('3DCTV_oneFrac','resultGUI','*dvh*');

clear resultGUI *dvh*

% reset number of fractions
pln.numOfFractions = numOfFractions;

%% DAD

% load up 3D-CTV
load('3DCTV','resultGUI');

% do DAD
resultGUI.apertureInfo.numPhases    = ct.tumourMotion.numPhases;
resultGUI.apertureInfo              = matRad_doDAD(resultGUI.apertureInfo,stf);
resultGUI.apertureInfo.run4D        = true;

% prepare motion model
resultGUI.apertureInfo.motionModel = matRad_prepModelForOpt(pln,stf,resultGUI.apertureInfo);

% update aperture vector
[resultGUI.apertureInfo.apertureVector, resultGUI.apertureInfo.mappingMx, resultGUI.apertureInfo.limMx] = matRad_daoApertureInfo2Vec(resultGUI.apertureInfo);

% update apertureInfo
for i = 1:numel(resultGUI.apertureInfo.beam)
    resultGUI.apertureInfo.beam(i).numUniqueVar = (resultGUI.apertureInfo.beam(i).numUniqueVar-resultGUI.apertureInfo.totalNumOfShapes).*resultGUI.apertureInfo.numPhases;
end
resultGUI.apertureInfo = matRad_daoVec2ApertureInfo(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector);

% save results
cd(currentDir);
save('DAD','resultGUI');

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('DAD','resultGUI','*dvh*');

% now do dvhs for single fraction
resultGUI.physicalDose  = resultGUI.physicalDose./pln.numOfFractions;
numOfFractions          = pln.numOfFractions;
pln.numOfFractions      = 1;

[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('DAD_oneFrac','resultGUI','*dvh*');

clear resultGUI *dvh*

% reset number of fractions
pln.numOfFractions = numOfFractions;

%% STO

% load up CO, for the apertureInfo
load('CO','resultGUI');

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% construct effective dij from most probable trajectory
% NOTE that right now this trajectory doesn't go through all of the phases
[dij_STO,trajectory] = matRad_dijSTO(dij,pln,stf,resultGUI.apertureInfo);

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('STO_FMO')
close all

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% turn off 4D
pln.propOpt.run4D = false;

% pick out particular trajectory from library, convert to single sequence
resultGUI.apertureInfo = matRad_library2ST(resultGUI.apertureInfo,pln,stf,trajectory);

% do DAO
resultGUI = matRad_directApertureOptimization(dij_STO,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('STO_DAO')
save('STO','resultGUI');
close all

% convert sequence to library
resultGUI.apertureInfo = matRad_apertures2Library(resultGUI.apertureInfo,pln,stf,dij.numPhases);

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('STO','resultGUI','*dvh*');

% now do dvhs for single fraction
resultGUI.physicalDose  = resultGUI.physicalDose./pln.numOfFractions;
numOfFractions          = pln.numOfFractions;
pln.numOfFractions      = 1;

[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('STO_oneFrac','resultGUI','*dvh*');

clear resultGUI *dvh*

% reset number of fractions
pln.numOfFractions = numOfFractions;

%% 3D optimization on ITV

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% construct dij for ITV by doing a time-average of the dij across each
% phase
dij_ITV                 = dij;
dij_ITV.numPhases       = 1;
dij_ITV.numFrames       = 1;
dij_ITV.physicalDose    = cell(1);
dij_ITV.physicalDose{1} = spalloc(dij_ITV.numOfVoxels,dij_ITV.totalNumOfBixels,nnz(dij.physicalDose{1}));
phaseProb               = accumarray(pln.propOpt.prop4D.motionModel.indices.subPhase2PosPhase,pln.propOpt.prop4D.motionModel.Pi_deltaTSample);

for phase = 1:dij.numPhases
    dij_ITV.physicalDose{1} = dij_ITV.physicalDose{1} + phaseProb(phase).*dij.physicalDose{phase};
end

% turn off 4d
pln.propOpt.run4D = false;

% do FMO
resultGUI = matRad_fluenceOptimization(dij_ITV,cst,pln,stf);
cd(currentDir);
savefig('3DITV_FMO')
close all

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij_ITV,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij_ITV,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('3DITV_DAO')
save('3DITV','resultGUI');
close all

% convert sequence to library
resultGUI.apertureInfo = matRad_apertures2Library(resultGUI.apertureInfo,pln,stf,dij.numPhases);

% change initProb
resultGUI.apertureInfo.motionModel.initProb = resultGUI.apertureInfo.motionModel.Pi_deltaTSample';

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('3DITV','resultGUI','*dvh*');

% now do dvhs for single fraction
resultGUI.physicalDose  = resultGUI.physicalDose./pln.numOfFractions;
numOfFractions          = pln.numOfFractions;
pln.numOfFractions      = 1;

[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('3DITV_oneFrac','resultGUI','*dvh*');

clear resultGUI *dvh*

% reset number of fractions
pln.numOfFractions = numOfFractions;

%% probabilistic optimization

pln.propOpt.run4D = true;
pln.propOpt.varOpt = true;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('PO_FMO')
close all

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('PO_DAO')
save('PO','resultGUI');
close all

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('PO','resultGUI','*dvh*');

clear resultGUI *dvh*

pln.propOpt.run4D = true;
pln.propOpt.varOpt = true;

% now redo optimization, dvhs for single fraction
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        cst{i,6}.dose = cst{i,6}.dose./pln.numOfFractions;
    end
end

pln.DRx             = pln.DRx./pln.numOfFractions;
numOfFractions      = pln.numOfFractions;
pln.numOfFractions  = 1;
dij.numOfFractions  = 1;

% redo FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('PO_oneFrac_FMO')
close all

% redo leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% redo DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('PO_oneFrac_DAO')
save('PO_oneFrac','resultGUI');
close all

% redo dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('PO_oneFrac','resultGUI','*dvh*');

clear resultGUI *dvh*

% reset number of fractions
pln.numOfFractions = numOfFractions;
dij.numOfFractions = numOfFractions;
