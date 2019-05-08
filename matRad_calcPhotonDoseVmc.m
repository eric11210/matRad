function dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
%
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
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

dij.radiationMode = pln.radiationMode;

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
dij.weightToMU         = 100;
dij.scaleFactor        = 1;
dij.memorySaverPhoton  = pln.propDoseCalc.memorySaverPhoton;
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.numOfScenarios     = 1;
if pln.propOpt.run4D
    dij.numPhases          = ct.tumourMotion.numPhases;
    dij.numFrames          = ct.tumourMotion.numFrames;
else
    dij.numPhases          = 1;
    dij.numFrames          = 1;
end

% check if full dose influence data is required
if calcDoseDirect
    numOfColumnsDij         = length(stf);
    numOfBixelsContainer    = 1;
else
    numOfColumnsDij         = dij.totalNumOfBixels;
    numOfBixelsContainer    = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

% Allocate space for dij.physicalDose sparse matrix
for k = 1:dij.numPhases
    dij.physicalDose{k}         = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
    if pln.propDoseCalc.vmcOptions.keepError
        dij.physicalDoseError{k}    = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
    end
end

% Allocate memory for dose_temp cell array
doseTmpContainer        = cell(numOfBixelsContainer,dij.numOfScenarios);
if pln.propDoseCalc.vmcOptions.keepError
    doseTmpContainerError   = cell(numOfBixelsContainer,dij.numOfScenarios);
end

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));
dij.structVox = V;

% find voxels in target
% find all target voxels from cst cell array
VTarget = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        VTarget = [VTarget;vertcat(cst{i,4}{:})];
    end
end
VTarget = unique(VTarget);
dij.targetVox = VTarget;

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
    vectorsPath = fullfile(runsPath, 'vectors');
    
    setenv('vmc_home',VMCPath);
    setenv('vmc_dir',runsPath);
    setenv('xvmc_dir',VMCPath);
    
    if isunix
        system(['chmod a+x ' VMCPath filesep 'bin' filesep 'vmc_Linux.exe']);
    end
    
end

% set consistent random seed (enables reproducibility)
rng(0);

% get default vmc options
VmcOptions = matRad_vmcOptions(pln,ct);

% export CT cube as binary file for vmc++
matRad_exportCtVmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

% save CT name
VmcOptions.geometry.Geometry.CtFile  = strrep(fullfile(phantomPath,'matRad_CT.ct'),'\','/'); % path of density matrix (only needed if input method is 'CT-PHANTOM')

maxNumOfParMCSim    = 0;

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for photons (vmc++)...');
% show busy state
set(figureWait,'pointer','watch');

fprintf('matRad: VMC++ photon dose calculation...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for frame = 1:dij.numFrames
    % loop over all frames first
    % this is because the .vectors file takes the longest to write, so we
    % want to have as few writes as possible
    
    fprintf('Frame %d of %d ...\n',frame,dij.numFrames);
    
    % reset these after every completed frame
    writeCounter        = 0;
    readCounter         = 0;
    
    if pln.propOpt.run4D
        % export vectors cubes as text file for vmc++
        matRad_exportVectorsVmc(ct,frame,fullfile(vectorsPath, 'matRad_MVF.vectors'));
        
        % save vectors name
        VmcOptions.geometry.Geometry.VectorsFile  = strrep(fullfile(vectorsPath,'matRad_MVF.vectors'),'\','/'); % path of density matrix (only needed if input method is 'CT-PHANTOM')
    end
    
    for i = 1:dij.numOfBeams % loop over all beams
        
        fprintf('Beam %d of %d ...\n',i,dij.numOfBeams);
        
        if strcmp(pln.propDoseCalc.vmcOptions.source,'phsp')
            % set angle-specific vmc++ parameters
            
            % phsp starts off pointed in the +z direction, with source at -z
            % phsp source gets translated, then rotated (-z, +y, -x) around
            % 0, then pushed to isocenter
            
            % correct for the source to collimator distance and change units mm -> cm
            translation = stf(i).isoCenter/10+[0 0 pln.propDoseCalc.vmcOptions.SCD + stf(i).sourcePoint_bev(2)]/10;
            
            % enter in isocentre
            isocenter = stf(i).isoCenter/10;
            
            % determine vmc++ rotation angles from gantry and couch
            % angles
            angles = matRad_matRad2vmcSourceAngles(stf(i).gantryAngle,stf(i).couchAngle);
            
            % set vmc++ parameters
            VmcOptions.source.translation   = translation;
            VmcOptions.source.isocenter     = isocenter;
            VmcOptions.source.angles        = angles;
        end
        
        for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
            
            writeCounter = writeCounter + 1;
            
            % create different seeds for every bixel
            VmcOptions.McControl.rngSeeds = [randi(30000),randi(30000)];
            
            % remember beam and bixel number
            if calcDoseDirect
                dij.beamNum(i)    = i;
                dij.rayNum(i)     = i;
                dij.bixelNum(i)   = i;
            else
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
                    VmcOptions.source.file_name     = strrep(stf(i).ray(j).phspFileName,'\','/');
            end
            
            
            %% create input file with vmc++ parameters
            outfile = ['MCpencilbeam_temp_',num2str(mod(writeCounter-1,VmcOptions.run.numOfParMCSim)+1)];
            matRad_createVmcInput(VmcOptions,fullfile(runsPath, [outfile,'.vmc']));
            
            % parallelization: only run this block for every numOfParallelMCSimulations!!!
            if mod(writeCounter,VmcOptions.run.numOfParMCSim) == 0 || writeCounter == dij.totalNumOfBixels
                
                % create batch file (enables parallel processes)
                if writeCounter == dij.totalNumOfBixels && mod(writeCounter,VmcOptions.run.numOfParMCSim) ~= 0
                    currNumOfParMCSim = mod(writeCounter,VmcOptions.run.numOfParMCSim);
                else
                    currNumOfParMCSim = VmcOptions.run.numOfParMCSim;
                end
                matRad_createVmcBatchFile(currNumOfParMCSim,fullfile(VMCPath,'run_parallel_simulations.bat'),verbose);
                
                % save max number of executed parallel simulations
                if currNumOfParMCSim > maxNumOfParMCSim
                    maxNumOfParMCSim = currNumOfParMCSim;
                end
                
                %% perform vmc++ simulation
                current = pwd;
                cd(VMCPath);
                if verbose > 0 % only show output if verbose level > 0
                    dos('run_parallel_simulations.bat');
                    fprintf('Completed %d of %d beamlets...\n',writeCounter,dij.totalNumOfBixels);
                else
                    %[dummyOut1,dummyOut2] = dos('run_parallel_simulations.bat'); % supress output by assigning dummy output arguments
                    dos('run_parallel_simulations.bat &','-echo'); % supress output by assigning dummy output arguments
                    
                    while ~exist('EmptyFile.txt','file')
                        pause(1);
                    end
                    delete('EmptyFile.txt');
                end
                cd(current);
                
                for k = 1:currNumOfParMCSim
                    readCounter     = readCounter+1;
                    
                    % update waitbar
                    waitbar((writeCounter+(frame-1).*dij.totalNumOfBixels)/(dij.numFrames.*dij.totalNumOfBixels));
                    
                    %% import calculated dose
                    idx = regexp(outfile,'_');
                    switch pln.propDoseCalc.vmcOptions.version
                        case 'Carleton'
                            filename = sprintf('%s%d.dos',outfile(1:idx(2)),k);
                        case 'dkfz'
                            filename = sprintf('%s%d_%s.dos',outfile(1:idx(2)),k,VmcOptions.scoringOptions.outputOptions.name);
                    end
                    [bixelDose,bixelDoseError] = matRad_readDoseVmc(fullfile(runsPath,filename),VmcOptions);
                    
                    
                    %%% Don't do any sampling, since the correct error is
                    % difficult to figure out. We also don't really need it on
                    % the Graham cluster.
                    
                    if ~calcDoseDirect
                        % if not calculating dose directly, sample dose
                        
                        % determine cutoff
                        doseCutoff          = VmcOptions.run.relDoseCutoff*max(bixelDose);
                        
                        % determine which voxels to sample
                        indSample = bixelDose < doseCutoff & bixelDose ~= 0;
                        r = rand(nnz(indSample),1);
                        
                        % sample them
                        thresRand = bixelDose(indSample)./doseCutoff;
                        indKeepSampled = r < thresRand;
                        indKeep = indSample;
                        indKeep(indKeep) = indKeepSampled;
                        
                        bixelDose(indKeep) = doseCutoff;
                        bixelDose(indSample & ~indKeep) = 0;
                        
                    end
                    
                    % apply absolute calibration factor
                    bixelDose       = bixelDose*VmcOptions.run.absCalibrationFactorVmc;
                    if pln.propDoseCalc.vmcOptions.keepError
                        bixelDoseError  = sqrt((VmcOptions.run.absCalibrationFactorVmc.*bixelDoseError).^2+(bixelDose.*VmcOptions.run.absCalibrationFactorVmc_err).^2);
                    end
                    
                    % determine the phase and normalization factor
                    if pln.propOpt.run4D
                        phase = ct.tumourMotion.frames2Phases(frame);
                        normFactor = ct.tumourMotion.nFramesPerPhase(phase);
                    else
                        phase = 1;
                        normFactor = 1;
                    end
                    
                    % Save dose for every bixel in cell array
                    % THIS IS CORRECT, FIX NON-VMC
                    % JUST ALWAYS OVER-WRITE THE TEMP ARRAYS
                    doseTmpContainer{mod(readCounter-1,numOfBixelsContainer)+1,1}   = sparse(V,1,bixelDose(V),dij.numOfVoxels,1)./normFactor;
                    if pln.propDoseCalc.vmcOptions.keepError
                        doseTmpContainerError{mod(readCounter-1,numOfBixelsContainer)+1,1}  = sparse(V,1,bixelDoseError(V),dij.numOfVoxels,1)./normFactor;
                    end
                    
                    % save computation time and memory by sequentially filling the
                    % sparse matrix dose.dij from the cell array
                    if mod(readCounter,numOfBixelsContainer) == 0 || readCounter == dij.totalNumOfBixels
                        if calcDoseDirect
                            if isfield(stf(beamNum(readCounter)).ray(rayNum(readCounter)),'weight')
                                % score physical dose
                                dij.physicalDose{phase}(:,i)        = dij.physicalDose{phase}(:,i) + stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight{1} * doseTmpContainer{1,1};
                                dij.physicalDoseError{phase}(:,i)   = sqrt(dij.physicalDoseError{phase}(:,i).^2 + (stf(beamNum(readCounter)).ray(rayNum(readCounter)).weight{1} * doseTmpContainerError{1,1}).^2);
                            else
                                error(['No weight available for beam ' num2str(beamNum(readCounter)) ', ray ' num2str(rayNum(readCounter))]);
                            end
                        else
                            % fill entire dose influence matrix
                            dij.physicalDose{phase}(:,(ceil(readCounter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:readCounter) = ...
                                dij.physicalDose{phase}(:,(ceil(readCounter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:readCounter) + ...
                                [doseTmpContainer{1:mod(readCounter-1,numOfBixelsContainer)+1,1}];
                            
                            if pln.propDoseCalc.vmcOptions.keepError
                                dij.physicalDoseError{phase}(:,(ceil(readCounter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:readCounter) = ...
                                    sqrt( dij.physicalDoseError{phase}(:,(ceil(readCounter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:readCounter).^2 + ...
                                    [doseTmpContainerError{1:mod(readCounter-1,numOfBixelsContainer)+1,1}].^2 );
                            end
                        end
                    end
                end 
            end 
        end
    end
    
    fprintf('Done all beams!\n');
end

fprintf('Done all phases!\n');

%% delete temporary files
delete(fullfile(VMCPath, 'run_parallel_simulations.bat'));  % batch file
delete(fullfile(phantomPath, 'matRad_CT.ct'));              % phantom file
delete(fullfile(vectorsPath, 'matRad_MVF.vectors'));        % vectors file
for j = 1:maxNumOfParMCSim
    delete(fullfile(runsPath, ['MCpencilbeam_temp_',num2str(mod(j-1,VmcOptions.run.numOfParMCSim)+1),'.vmc'])); % vmc inputfile
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
