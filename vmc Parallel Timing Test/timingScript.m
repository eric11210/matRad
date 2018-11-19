%% setup 

load TG119.mat

% meta information for treatment plan
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.machine         = 'Generic';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = 0:72:359; % [?]
pln.propStf.couchAngles     = [0 0 0 0 0]; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% dose calculation settings
pln.propDoseCalc.memorySaverPhoton          = false;
pln.propDoseCalc.vmc                        = true;
pln.propDoseCalc.vmcOptions.source          = 'phsp';
pln.propDoseCalc.vmcOptions.phspBaseName    = '5cmx5cm_SSD50cm';
pln.propDoseCalc.vmcOptions.SCD             = 500;
pln.propDoseCalc.vmcOptions.dumpDose        = 1;
pln.propDoseCalc.vmcOptions.version         = 'Carleton';
pln.propDoseCalc.vmcOptions.nCasePerBixel   = 5000;
pln.propDoseCalc.vmcOptions.numOfParMCSim   = 8;

% optimization settings
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                            % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.runVMAT         = false;  % 1/true: run VMAT, 0/false: don't
pln.propOpt.run4D           = false;

% generate steering file
stf = matRad_generateStf(ct,cst,pln);

% vector of parallel simulations
numOfParMCSim_vec   = 2.^(8:-1:0)';
calcTime_vec        = zeros(size(numOfParMCSim_vec));

curDir = pwd;

for i = 1:numel(calcTime_vec)
    
    % modify num parallel simulations
    pln.propDoseCalc.vmcOptions.numOfParMCSim   = numOfParMCSim_vec(i);
    
    % calc dose with timer
    tstart = tic;
    
    dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
    
    calcTime_vec(i) = toc(tstart);
    
    % save results
    cd(curDir);
    save('results','*_vec');
end