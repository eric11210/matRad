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

%% calculate dij

dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);

%% conventional optimization

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('CO_FMO')

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('CO_DAO')
save('CO','resultGUI');

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('CO','resultGUI','*dvh*');

clear resultGUI *dvh*

%% probabilistic optimization

pln.propOpt.run4D = true;
pln.propOpt.varOpt = true;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('PO_FMO')

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('PO_DAO')
save('PO','resultGUI');

% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('PO','resultGUI','*dvh*');

clear resultGUI *dvh*

%% 3D optimization on CTV

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('3DCTV_FMO')

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

%{
% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('3DCTV','resultGUI','*dvh*');
%}

clear resultGUI *dvh*

%% 3D optimization on ITV

pln.propOpt.run4D = true;
pln.propOpt.varOpt = false;

% change obj function goals
cst{26,6}       = cst{25,6};
cst{25,6}       = [];
pln.RxStruct    = 26;

% do FMO
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
cd(currentDir);
savefig('3DITV_FMO')

% turn off 4d
pln.propOpt.run4D = false;

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('3DITV_DAO')
save('3DITV','resultGUI');

%{
% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('3DCTV','resultGUI','*dvh*');
%}

clear resultGUI *dvh*

% change obj function goals
cst{25,6}       = cst{26,6};
cst{26,6}       = [];
pln.RxStruct    = 25;


%% DAD

% load up 3D-CTV
load('3DCTV','resultGUI');

% do DAD
resultGUI.apertureInfo.numPhases    = ct.tumourMotion.numPhases;
resultGUI.apertureInfo              = matRad_doDAD(resultGUI.apertureInfo,stf);

% update aperture vector
[resultGUI.apertureInfo.apertureVector, resultGUI.apertureInfo.mappingMx, resultGUI.apertureInfo.limMx] = matRad_daoApertureInfo2Vec(resultGUI.apertureInfo);

% save results
cd(currentDir);
save('DAD','resultGUI');

%{
% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('DAD','resultGUI','*dvh*');
%}


%% STO

pln.propOpt.run4D = false;
pln.propOpt.varOpt = false;

% select trajectory
singleTrajectory.phase = zeros(pln.propStf.numOfBeams,1);

% construct effective dij
dij_fieldNames = fieldnames(dij);

for i = 1:numel(dij_fieldNames)
    
    if ~strcmp(dij_fieldNames{i},'physicalDose')
        
        dijEff.(dij_fieldNames{i}) = dij.(dij_fieldNames{i});
    end
end

dijEff.physicalDose{1}  = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
dij.numOfScenarios      = 1;

offset = 0;
for i = 1:pln.propStf.numOfBeams
    
    dijEff.physicalDose{1}(:,offset+(1:stf(i).totalNumOfBixels)) = dij.physicalDose{singleTrajectory.phase(i)}(:,offset+(1:stf(i).totalNumOfBixels));
    
    offset = offset+stf(i).totalNumOfBixels;
end

% do FMO
resultGUI = matRad_fluenceOptimization(dijEff,cst,pln,stf);
cd(currentDir);
savefig('STO_FMO')

% turn off 4D
pln.propOpt.run4D = false;

% do leaf sequencing
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dijEff,pln,0);

% do DAO
resultGUI = matRad_directApertureOptimization(dijEff,cst,resultGUI.apertureInfo,resultGUI,pln,stf);

% save results
cd(currentDir);
savefig('STO_DAO')
save('STO','resultGUI');

%{
% do dvhs
[pdvh_MC,dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(resultGUI.apertureInfo,dij,cst,pln,100);
[dvh,~] = matRad_indicatorWrapper(cst,pln,resultGUI);

% save results
save('STO','resultGUI','*dvh*');
%}

clear resultGUI *dvh*