function dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,nCasePerBixel,numOfParallelMCSimulations,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,nCasePerBixel,numOfParallelMCSimulations)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nCasePerBixel:              number of photons simulated per bixel
%   numOfParallelMCSimulations: number of simultaneously performed simulations (optional) 
%   calcDoseDirect:             boolian switch to bypass dose influence matrix
%                               computation and directly calculate dose; only makes
%                               sense in combination with matRad_calcDoseDirect.m%
% output
%   dij:                        matRad dij struct
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% set output level. 0 = no vmc specific output. 1 = print to matlab cmd.
% 2 = open in terminal(s)
verbose = 0;

if ~isdeployed % only if _not_ running as standalone    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'vmc++'))
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if calcDoseDirect 
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
end

% set environment variables for vmc++
cd(fileparts(mfilename('fullpath')))

if exist(['vmc++' filesep 'bin'],'dir') ~= 7
    error(['Could not locate vmc++ environment. ' ...
        'Please provide the files in the correct folder structure at matRadroot' filesep 'vmc++.']);
else
    VMCPath     = fullfile(pwd , 'vmc++');
    switch pln.propDoseCalc.vmcOptions.version
        case 'Carleton'
            runsPath    = fullfile(VMCPath, 'run');
        case 'dkfz'
            runsPath    = fullfile(VMCPath, 'runs');
    end
    phantomPath = fullfile(runsPath, 'phantoms');

    setenv('vmc_home',VMCPath);
    setenv('vmc_dir',runsPath);
    setenv('xvmc_dir',VMCPath);
    
    if isunix
        system(['chmod a+x ' VMCPath filesep 'bin' filesep 'vmc_Linux.exe']);
    end
    
end

% set consistent random seed (enables reproducibility)
rng(0);

% set number of photons simulated per bixel and number of parallel MC simulations if not specified by user
if nargin < 5
    nCasePerBixel              = 5000;
    numOfParallelMCSimulations = 4;
    
    warning(['Number of photons simulated per bixel (nCasePerBixel) and number of parallel MC simulations (numOfParallelMCSimulations) not specified by user. ',...
        'Use default settings with nCasePerBixel = ',num2str(nCasePerBixel),...
        ' and numOfParallelMCSimulations = ',num2str(numOfParallelMCSimulations),...
        ' in vmc++ calculations.'])
    
elseif nargin < 6
    numOfParallelMCSimulations = 4;
    
    warning(['Number of parallel MC simulations (numOfParallelMCSimulations) not specified by user. ',...
        'Use default settings with numOfParallelMCSimulations = ',num2str(numOfParallelMCSimulations),...
        ' in vmc++ calculations.'])
end

if isunix && numOfParallelMCSimulations > 1
    numOfParallelMCSimulations = 1;
    warning(['Running Unix environment: Number of parallel MC simulations (numOfParallelMCSimulations) set to default settings with numOfParallelMCSimulations = ',num2str(numOfParallelMCSimulations),...
        ' in vmc++ calculations.'])
end

% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        absCalibrationFactorVmc = 1;
    case 'dkfz'
        absCalibrationFactorVmc = 99.818252282632300;
end

% set general vmc++ parameters
VmcOptions.version = pln.propDoseCalc.vmcOptions.version;
% 1 source
VmcOptions.source.myName       = 'some_source';                        % name of source
VmcOptions.source.monitorUnits = 1;
if strcmp(pln.propDoseCalc.vmcOptions.source,'beamlet')
    VmcOptions.source.spectrum     = fullfile(runsPath,'spectra','var_6MV.spectrum');    % energy spectrum source (only used if no mono-Energy given)
    VmcOptions.source.charge       = 0;                                 % charge (-1,0,1
    VmcOptions.source.type         = 'beamlet';
elseif strcmp(pln.propDoseCalc.vmcOptions.source,'phsp')
    VmcOptions.source.particleType  = 2;
    VmcOptions.source.type          = 'phsp';
end
% 2 transport parameter
VmcOptions.McParameter.automatic_parameter  = 'yes';                       % if yes, automatic transport parameters are used
VmcOptions.McParameter.spin                 = 0;                           % 0: spin effects ignored; 1: simplistic; 2: full treatment
% 3 MC control
VmcOptions.McControl.ncase  = nCasePerBixel;                               % number of histories
VmcOptions.McControl.nbatch = 10;                                          % number of batches
% 4 variance reduction
VmcOptions.varianceReduction.repeatHistory      = 0.041;
VmcOptions.varianceReduction.splitPhotons       = 1;   
VmcOptions.varianceReduction.photonSplitFactor = -80;  
% 5 quasi random numbers
VmcOptions.quasi.base      = 2;                                                 
VmcOptions.quasi.dimension = 60;                                             
VmcOptions.quasi.skip      = 1;
% 6 geometry
switch pln.propDoseCalc.vmcOptions.version
    case 'Carleton'
        VmcOptions.geometry.XyzGeometry.methodOfInput = 'MMC-PHANTOM';             % input method ('CT-PHANTOM', 'individual', 'groups')
    case 'dkfz'
        VmcOptions.geometry.XyzGeometry.methodOfInput = 'CT-PHANTOM';             % input method ('CT-PHANTOM', 'individual', 'groups')
end
VmcOptions.geometry.dimensions          = dij.dimensions;
VmcOptions.geometry.XyzGeometry.Ct      = 'CT';                      % name of geometry
% 7 scoring manager
VmcOptions.scoringOptions.startInGeometry               = 'CT';            % geometry in which partciles start their transport
VmcOptions.scoringOptions.doseOptions.scoreInGeometries = 'CT';            % geometry in which dose is recorded
VmcOptions.scoringOptions.doseOptions.scoreDoseToWater  = 'yes';           % if yes output is dose to water
VmcOptions.scoringOptions.outputOptions.name            = 'CT';            % geometry for which dose output is created (geometry has to be scored)
VmcOptions.scoringOptions.outputOptions.dumpDose        = pln.propDoseCalc.vmcOptions.dumpDose;               % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

writeCounter                  = 0;
readCounter                   = 0;
maxNumOfParallelMcSimulations = 0;

% initialize waitbar
figureWait = waitbar(0,'VMC++ photon dose influence matrix calculation..');

fprintf('matRad: VMC++ photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
       
   % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
        dij.bixelNum(i)   = i;
    end
    
    if strcmp(pln.propDoseCalc.vmcOptions.source,'phsp')
        % set angle-specific vmc++ parameters
        
        % export CT cube as binary file for vmc++
        % for the phsp source, rotations are made around the origin
        % need to make the origin the isocenter in the CT
        matRad_exportCtVmc(ct, fullfile(phantomPath, sprintf('matRad_CT_beam%d.ct',i)), stf(i).isoCenter);
        
        % phsp starts off pointed in the +z direction, with source at -z
        
        % phsp source gets translated, then rotated (-z, +y, -x) around
        % origin
        
        % a) correct for the source to collimator distance and change units mm -> cm
        translation = [0 0 pln.propDoseCalc.vmcOptions.SCD + stf(i).sourcePoint_bev(2)]/10;
        
        % b) determine vmc++ rotation angles from gantry and couch
        % angles
        angles = matRad_matRad2vmcSourceAngles(stf(i).gantryAngle,stf(i).couchAngle);
        
        % c) set vmc++ parameters
        VmcOptions.source.translation = translation;
        VmcOptions.source.angles = angles;
        
    else
        % export CT cube as binary file for vmc++
        
        % use dummy isocenter
        matRad_exportCtVmc(ct, fullfile(phantomPath, sprintf('matRad_CT_beam%d.ct',i)), [0 0 0]);
    end
    
    % use beam-specific CT name
    VmcOptions.geometry.XyzGeometry.CtFile  = strrep(fullfile(runsPath,'phantoms',sprintf('matRad_CT_beam%d.ct',i)),'\','/'); % path of density matrix (only needed if input method is 'CT-PHANTOM')
    
    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        writeCounter = writeCounter + 1;

        % create different seeds for every bixel
        VmcOptions.McControl.rngSeeds = [randi(30000),randi(30000)];

        % remember beam and bixel number
        if ~calcDoseDirect
           dij.beamNum(writeCounter)  = i;
           dij.rayNum(writeCounter)   = j;
           dij.bixelNum(writeCounter) = j;
        end
        beamNum(writeCounter)  = i;
        rayNum(writeCounter)   = j;
        bixelNum(writeCounter) = j;
        
        % set ray specific vmc++ parameters
        switch pln.propDoseCalc.vmcOptions.source
            case 'beamlet'
                % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
                rayCorner1 = (stf(i).ray(j).rayCorners_SCD(1,:) + stf(i).isoCenter)/10;
                rayCorner2 = (stf(i).ray(j).rayCorners_SCD(2,:) + stf(i).isoCenter)/10;
                rayCorner3 = (stf(i).ray(j).rayCorners_SCD(3,:) + stf(i).isoCenter)/10; %vmc needs only three corners (counter-clockwise)
                beamSource = (stf(i).sourcePoint + stf(i).isoCenter)/10;
                
                % b) swap x and y (CT-standard = [y,x,z])
                rayCorner1 = rayCorner1([2,1,3]);
                rayCorner2 = rayCorner2([2,1,3]);
                rayCorner3 = rayCorner3([2,1,3]);
                beamSource = beamSource([2,1,3]);
                
                % c) set vmc++ parameters
                VmcOptions.source.monoEnergy                    = stf(i).ray(j).energy;                 % photon energy
                %VmcOptions.source.monoEnergy                   = []                  ;                  % use photon spectrum
                VmcOptions.source.beamletEdges                  = [rayCorner1,rayCorner2,rayCorner3];    % counter-clockwise beamlet edges
                VmcOptions.source.virtualPointSourcePosition    = beamSource;                            % virtual beam source position
                
            case 'phsp'
                % use ray-specific file name for the phsp source (bixelized
                % phsp)
                VmcOptions.source.file_name     = strrep(fullfile(runsPath,'phsp',stf(i).ray(j).phspFileName),'\','/');
        end
        
        
        % create inputfile with vmc++ parameters
        outfile = ['MCpencilbeam_temp_',num2str(mod(writeCounter-1,numOfParallelMCSimulations)+1)];
        matRad_createVmcInput(VmcOptions,fullfile(runsPath, [outfile,'.vmc']));
        
        % parallelization: only run this block for every numOfParallelMCSimulations!!!
        if mod(writeCounter,numOfParallelMCSimulations) == 0 || writeCounter == dij.totalNumOfBixels
            
            % create batch file (enables parallel processes)
            if writeCounter == dij.totalNumOfBixels && mod(writeCounter,numOfParallelMCSimulations) ~= 0
                currNumOfParallelMcSimulations = mod(writeCounter,numOfParallelMCSimulations);
            else
                currNumOfParallelMcSimulations = numOfParallelMCSimulations;
            end
            matRad_createVmcBatchFile(currNumOfParallelMcSimulations,fullfile(VMCPath,'run_parallel_simulations.bat'),verbose);
            
            % save max number of executed parallel simulations
            if currNumOfParallelMcSimulations > maxNumOfParallelMcSimulations 
                maxNumOfParallelMcSimulations = currNumOfParallelMcSimulations;
            end
            
            % perform vmc++ simulation
            current = pwd;
            cd(VMCPath);
            if verbose > 0 % only show output if verbose level > 0
                dos('run_parallel_simulations.bat');
                fprintf(['Completed ' num2str(writeCounter) ' of ' num2str(dij.totalNumOfBixels) ' beamlets...\n']);
            else
                [dummyOut1,dummyOut2] = dos('run_parallel_simulations.bat'); % supress output by assigning dummy output arguments
            end
            cd(current);
            
            for k = 1:currNumOfParallelMcSimulations
                readCounter = readCounter + 1;
                
                % Display progress
                if verbose == 0
                   % matRad_progress(readCounter,dij.totalNumOfBixels);
                end
                
                % update waitbar
                waitbar(writeCounter/dij.totalNumOfBixels);
                
                % import calculated dose
                idx = regexp(outfile,'_');
                switch pln.propDoseCalc.vmcOptions.version
                    case 'Carleton'
                        filename = sprintf('%s%d.dos',outfile(1:idx(2)),k);
                    case 'dkfz'
                        filename = sprintf('%s%d_%s.dos',outfile(1:idx(2)),k,VmcOptions.scoringOptions.outputOptions.name);
                end
                [bixelDose,~] = matRad_readDoseVmc(fullfile(runsPath,filename),VmcOptions);

                % apply relative dose cutoff
                doseCutoff                        = relDoseCutoff*max(bixelDose);
                bixelDose(bixelDose < doseCutoff) = 0;

                % apply absolute calibration factor
                bixelDose = bixelDose*absCalibrationFactorVmc;

                % Save dose for every bixel in cell array
                doseTmpContainer{mod(readCounter-1,numOfBixelsContainer)+1,1} = sparse(V,1,bixelDose(V),dij.numOfVoxels,1);
                
                % save computation time and memory by sequentially filling the 
                % sparse matrix dose.dij from the cell array
                if mod(readCounter,numOfBixelsContainer) == 0 || readCounter == dij.totalNumOfBixels
                    if calcDoseDirect
                        if isfield(stf(beamNum(readCounter)).ray(rayNum(readCounter)),'weight')
                            % score physical dose
                            dij.physicalDose{1}(:,i) = dij.physicalDose{1}(:,i) + stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight * doseTmpContainer{1,1};
                        else
                            error(['No weight available for beam ' num2str(beamNum(readCounter)) ', ray ' num2str(rayNum(readCounter))]);
                        end
                    else
                        % fill entire dose influence matrix
                        dij.physicalDose{1}(:,(ceil(readCounter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:readCounter) = ...
                            [doseTmpContainer{1:mod(readCounter-1,numOfBixelsContainer)+1,1}];
                    end
                end
            end
            
        end
        
    end
end

% delete temporary files
delete(fullfile(VMCPath, 'run_parallel_simulations.bat')); % batch file
delete(fullfile(phantomPath, 'matRad_CT_beam*.ct'));             % phantom file
for j = 1:maxNumOfParallelMcSimulations
    delete(fullfile(runsPath, ['MCpencilbeam_temp_',num2str(mod(j-1,numOfParallelMCSimulations)+1),'.vmc'])); % vmc inputfile
    switch pln.propDoseCalc.vmcOptions.version
        case 'Carleton'
            filename = sprintf('%s%d.dos','MCpencilbeam_temp_',j);
        case 'dkfz'
            filename = sprintf('%s%d_%s.dos','MCpencilbeam_temp_',j,VmcOptions.scoringOptions.outputOptions.name);
    end
    delete(fullfile(runsPath,filename));    % vmc outputfile
end

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end
