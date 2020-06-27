%% setup

load lungPatient1_3mm5p_rep

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

numOfFractions = pln.numOfFractions;

% absolute path to figure folder
figurePath = 'C:\Users\eric\Carleton University\OneDrive - Carleton University\Carleton\PhD Project\Thesis\figures\generator\';

%% CO

load('CO')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_CO_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% CO_oneFrac

load('CO_oneFrac')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_CO_oneFrac_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% 3DCTV

load('3DCTV')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_3DCTV_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% 3DCTV_oneFrac

load('3DCTV_oneFrac')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_3DCTV_oneFrac_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% DAD

load('DAD')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_DAD_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% DAD_oneFrac

load('DAD_oneFrac')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_DAD_oneFrac_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% STO

load('STO')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_STO_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% STO_oneFrac

load('STO_oneFrac')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_STO_oneFrac_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% 3DITV

load('3DITV')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_3DITV_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% 3DITV_oneFrac

load('3DITV_oneFrac')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_3DITV_oneFrac_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% PO

load('PO')
pln.numOfFractions = numOfFractions;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_PO_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

%% PO_oneFrac

load('PO_oneFrac')
pln.numOfFractions = 1;

matRad_showDVHBands(dvh,pdvh_MC,[1 3 5],cst,pln);

% save figure
fname = 'lungPatient3_PO_oneFrac_DVH_raw';
fullpath = [figurePath fname '.tex'];
matlab2tikz('filename',fullpath,'interpretTickLabelsAsTex',true,'parseStrings',false,'noSize',true,'extraAxisOptions','minor y tick num=3,minor x tick num=3,width=\linewidth,height=\linewidth','showInfo',false);

