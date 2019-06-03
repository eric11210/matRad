function apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to generate a shape info struct based on the result of
% multileaf collimator sequencing
%
% call
%   apertureInfo = matRad_sequencing2ApertureInfo(Sequencing,stf)
%
% input
%   Sequencing: matRad sequencing result struct
%   stf:        matRad steering information struct
%
% output
%   apertureInfo: matRad aperture weight and shape info struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MLC parameters:
bixelWidth = stf(1).bixelWidth; % [mm]
numOfMLCLeafPairs = 80;
%     define central leaf pair (here we want the 0mm position to be in the
%     center of a leaf pair (e.g. leaf 41 stretches from -2.5mm to 2.5mm
%     for a bixel/leafWidth of 5mm and 81 leaf pairs)
centralLeafPair = ceil(numOfMLCLeafPairs/2);

% initializing variables
bixelIndOffset = 0; % used for creation of bixel index maps
totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
totalNumOfShapes = sum([sequencing.beam.numOfShapes]);
weightOffset = 1;
vectorOffset = totalNumOfShapes + 1; % used for bookkeeping in the vector for optimization
bixOffset = 1; %used for gradient calculations

if pln.propOpt.runVMAT
    totalNumOfOptBixels = 0;
    totalNumOfLeafPairs = 0;
    
    apertureInfo.propVMAT.jacobT = zeros(sum([sequencing.beam.numOfShapes]),numel(sequencing.beam));
end

apertureInfo.jacobiScale = ones(totalNumOfShapes.*sequencing.numPhases,1);

% loop over all beams
for i=1:size(stf,2)
    
    %% 1. read stf and create maps (Ray & Bixelindex)
    
    % get x- and z-coordinates of bixels
    rayPos_bev = reshape([stf(i).ray.rayPos_bev],3,[]);
    X = rayPos_bev(1,:)';
    Z = rayPos_bev(3,:)';
    
    % create ray-map
    maxX = max(X); minX = min(X);
    maxZ = max(Z); minZ = min(Z);
    
    dimX = (maxX-minX)/stf(i).bixelWidth + 1;
    dimZ = (maxZ-minZ)/stf(i).bixelWidth + 1;
    
    rayMap = zeros(dimZ,dimX);
    
    % get indices for x and z positions
    xPos = (X-minX)/stf(i).bixelWidth + 1;
    zPos = (Z-minZ)/stf(i).bixelWidth + 1;
    
    % get indices in the ray-map
    indInRay = zPos + (xPos-1)*dimZ;
    
    % fill ray-map
    rayMap(indInRay) = 1;
    
    % create map of bixel indices
    bixelIndMap = NaN * ones(dimZ,dimX);
    bixelIndMap(indInRay) = (1:stf(i).numOfRays) + bixelIndOffset;
    bixelIndOffset = bixelIndOffset + stf(i).numOfRays;
    
    % store physical position of first entry in bixelIndMap
    posOfCornerBixel = [minX 0 minZ];
    
    % get leaf limits from the leaf map
    lim_l = NaN * ones(dimZ,1);
    lim_r = NaN * ones(dimZ,1);
    % looping oder leaf pairs
    for l = 1:dimZ
        lim_lInd = find(rayMap(l,:),1,'first');
        lim_rInd = find(rayMap(l,:),1,'last');
        % the physical position [mm] can be calculated from the indices
        lim_l(l) = (lim_lInd-1)*bixelWidth + minX - 1/2*bixelWidth;
        lim_r(l) = (lim_rInd-1)*bixelWidth + minX + 1/2*bixelWidth;
    end
    
    % get leaf positions for all shapes
    % leaf positions can be extracted from the shapes created in Sequencing
    for m = 1:sequencing.beam(i).numOfShapes
        
        if isfield(sequencing.beam(i),'shapes')
            
            % loading shape from Sequencing result
            shapeMap = sequencing.beam(i).shapes(:,:,m);
            % get left and right leaf indices from shapemap
            % initializing limits
            leftLeafPos = NaN * ones(dimZ,1);
            rightLeafPos = NaN * ones(dimZ,1);
            % looping over leaf pairs
            for l = 1:dimZ
                leftLeafPosInd  = find(shapeMap(l,:),1,'first');
                rightLeafPosInd = find(shapeMap(l,:),1,'last');
                
                if isempty(leftLeafPosInd) && isempty(rightLeafPosInd) % if no bixel is open, use limits from Ray positions
                    leftLeafPos(l) = (lim_l(l)+lim_r(l))/2;
                    rightLeafPos(l) = leftLeafPos(l);
                else
                    % the physical position [mm] can be calculated from the indices
                    leftLeafPos(l) = (leftLeafPosInd-1)*bixelWidth...
                        + minX - 1/2*bixelWidth;
                    rightLeafPos(l) = (rightLeafPosInd-1)*bixelWidth...
                        + minX + 1/2*bixelWidth;
                    
                    %Can happen in some cases in SW trajectory sampling
                    if leftLeafPos(l) < lim_l(l)
                        leftLeafPos(l) = lim_l(l);
                    end
                    if rightLeafPos(l) > lim_r(l)
                        rightLeafPos(l) = lim_r(l);
                    end
                    
                end
            end
            
            % save data for each shape of this beam
            apertureInfo.beam(i).shape{1}(m).leftLeafPos = leftLeafPos;
            apertureInfo.beam(i).shape{1}(m).rightLeafPos = rightLeafPos;
            apertureInfo.beam(i).shape{1}(m).weight = sequencing.beam(i).shapesWeight(m);
            apertureInfo.beam(i).shape{1}(m).shapeMap = shapeMap;
            
        elseif isfield(sequencing.beam(i).shape(m),'leftLeafPos')
            % leaf positions already determined
            
            % save data for each shape of this beam
            apertureInfo.beam(i).shape{1}(m).leftLeafPos = sequencing.beam(i).shape(m).leftLeafPos;
            apertureInfo.beam(i).shape{1}(m).rightLeafPos = sequencing.beam(i).shape(m).rightLeafPos;
            apertureInfo.beam(i).shape{1}(m).weight = sequencing.beam(i).shapesWeight(m);
            
            if sequencing.dynamic
                apertureInfo.beam(i).shape{1}(m).leftLeafPos_I = sequencing.beam(i).shape(m).leftLeafPos_I;
                apertureInfo.beam(i).shape{1}(m).rightLeafPos_I = sequencing.beam(i).shape(m).rightLeafPos_I;
                apertureInfo.beam(i).shape{1}(m).leftLeafPos_F = sequencing.beam(i).shape(m).leftLeafPos_F;
                apertureInfo.beam(i).shape{1}(m).rightLeafPos_F = sequencing.beam(i).shape(m).rightLeafPos_F;
            end
            
        end
        
        if pln.propOpt.runVMAT
            apertureInfo.beam(i).shape{1}(m).MURate = sequencing.beam(i).MURate;
        end
        
        apertureInfo.beam(i).shape{1}(m).jacobiScale = 1;
        
        if pln.propOpt.VMAToptions.continuousAperture
            apertureInfo.beam(i).shape{1}(m).vectorOffset = [vectorOffset vectorOffset+dimZ];
            
            % update index for bookkeeping
            vectorOffset = vectorOffset + dimZ*nnz(stf(i).propVMAT.doseAngleDAO);
        else
            apertureInfo.beam(i).shape{1}(m).vectorOffset = vectorOffset;
            
            % update index for bookkeeping
            vectorOffset = vectorOffset + dimZ;
        end
        apertureInfo.beam(i).shape{1}(m).weightOffset = weightOffset;
        weightOffset = weightOffset+1;
    end
    
    % z-coordinates of active leaf pairs
    % get z-coordinates from bixel positions
    leafPairPos = unique(Z);
    
    % find upmost and downmost leaf pair
    topLeafPairPos = maxZ;
    bottomLeafPairPos = minZ;
    
    topLeafPair = centralLeafPair - topLeafPairPos/bixelWidth;
    bottomLeafPair = centralLeafPair - bottomLeafPairPos/bixelWidth;
    
    % create bool map of active leaf pairs
    isActiveLeafPair = zeros(numOfMLCLeafPairs,1);
    isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;
    
    % create MLC window
    % getting the dimensions of the MLC in order to be able to plot the
    % shapes using physical coordinates
    MLCWindow = [minX-bixelWidth/2 maxX+bixelWidth/2 ...
        minZ-bixelWidth/2 maxZ+bixelWidth/2];
    
    % save data for each beam
    apertureInfo.beam(i).numOfShapes = sequencing.beam(i).numOfShapes;
    apertureInfo.beam(i).numOfActiveLeafPairs = dimZ;
    apertureInfo.beam(i).leafPairPos = leafPairPos;
    apertureInfo.beam(i).isActiveLeafPair = isActiveLeafPair;
    apertureInfo.beam(i).centralLeafPair = centralLeafPair;
    apertureInfo.beam(i).lim_l = lim_l;
    apertureInfo.beam(i).lim_r = lim_r;
    apertureInfo.beam(i).bixelIndMap = bixelIndMap;
    apertureInfo.beam(i).posOfCornerBixel = posOfCornerBixel;
    apertureInfo.beam(i).MLCWindow = MLCWindow;
    apertureInfo.beam(i).gantryAngle = stf(i).gantryAngle;
    apertureInfo.beam(i).numBixels = nnz(~isnan(bixelIndMap));
    
    if pln.propOpt.runVMAT
        
        apertureInfo.beam(i).bixOffset = bixOffset;
        bixOffset = bixOffset+apertureInfo.beam(i).numOfActiveLeafPairs;
        
        apertureInfo.propVMAT.beam(i).DAOBeam = stf(i).propVMAT.DAOBeam;
        apertureInfo.propVMAT.beam(i).FMOBeam = stf(i).propVMAT.FMOBeam;
        
        apertureInfo.propVMAT.beam(i).leafDir = sequencing.beam(i).leafDir;
        
        apertureInfo.propVMAT.beam(i).doseAngleBorders = stf(i).propVMAT.doseAngleBorders;
        apertureInfo.propVMAT.beam(i).doseAngleBorderCentreDiff = stf(i).propVMAT.doseAngleBorderCentreDiff;
        apertureInfo.propVMAT.beam(i).doseAngleBordersDiff = stf(i).propVMAT.doseAngleBordersDiff;
        
        if apertureInfo.propVMAT.beam(i).DAOBeam
            
            totalNumOfOptBixels = totalNumOfOptBixels+stf(i).totalNumOfBixels;
            totalNumOfLeafPairs = totalNumOfLeafPairs+apertureInfo.beam(i).numOfShapes*apertureInfo.beam(i).numOfActiveLeafPairs;
            
            apertureInfo.beam(i).gantryRot = sequencing.beam(i).gantryRot;
            
            apertureInfo.propVMAT.beam(i).DAOAngleBorders = stf(i).propVMAT.DAOAngleBorders;
            apertureInfo.propVMAT.beam(i).DAOAngleBorderCentreDiff = stf(i).propVMAT.DAOAngleBorderCentreDiff;
            apertureInfo.propVMAT.beam(i).DAOAngleBordersDiff = stf(i).propVMAT.DAOAngleBordersDiff;
            
            apertureInfo.propVMAT.beam(i).timeFacCurr = stf(i).propVMAT.timeFacCurr;
            apertureInfo.propVMAT.beam(i).timeFac = stf(i).propVMAT.timeFac;
            
            apertureInfo.propVMAT.beam(i).lastDAOIndex = stf(i).propVMAT.lastDAOIndex;
            apertureInfo.propVMAT.beam(i).nextDAOIndex = stf(i).propVMAT.nextDAOIndex;
            apertureInfo.propVMAT.beam(i).DAOIndex = stf(i).propVMAT.DAOIndex;
            
            if apertureInfo.propVMAT.beam(i).FMOBeam
                apertureInfo.propVMAT.beam(i).FMOAngleBorders = stf(i).propVMAT.FMOAngleBorders;
                apertureInfo.propVMAT.beam(i).FMOAngleBorderCentreDiff = stf(i).propVMAT.FMOAngleBorderCentreDiff;
                apertureInfo.propVMAT.beam(i).FMOAngleBordersDiff = stf(i).propVMAT.FMOAngleBordersDiff;
            end
            
            if pln.propOpt.VMAToptions.continuousAperture
                apertureInfo.propVMAT.beam(i).timeFacInd = stf(i).propVMAT.timeFacInd;
                apertureInfo.propVMAT.beam(i).doseAngleDAO = stf(i).propVMAT.doseAngleDAO;
            end
            
            apertureInfo.propVMAT.jacobT(stf(i).propVMAT.DAOIndex,i) = stf(i).propVMAT.timeFacCurr;
            
        else
            apertureInfo.propVMAT.beam(i).fracFromLastDAO = stf(i).propVMAT.fracFromLastDAO;
            apertureInfo.propVMAT.beam(i).timeFracFromLastDAO = stf(i).propVMAT.timeFracFromLastDAO;
            apertureInfo.propVMAT.beam(i).timeFracFromNextDAO = stf(i).propVMAT.timeFracFromNextDAO;
            apertureInfo.propVMAT.beam(i).lastDAOIndex = stf(i).propVMAT.lastDAOIndex;
            apertureInfo.propVMAT.beam(i).nextDAOIndex = stf(i).propVMAT.nextDAOIndex;
            
            if pln.propOpt.VMAToptions.continuousAperture
                apertureInfo.propVMAT.beam(i).fracFromLastDAO_I = stf(i).propVMAT.fracFromLastDAO_I;
                apertureInfo.propVMAT.beam(i).fracFromLastDAO_F = stf(i).propVMAT.fracFromLastDAO_F;
                apertureInfo.propVMAT.beam(i).fracFromNextDAO_I = stf(i).propVMAT.fracFromNextDAO_I;
                apertureInfo.propVMAT.beam(i).fracFromNextDAO_F = stf(i).propVMAT.fracFromNextDAO_F;
            end
            
            apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOIndex,i) = stf(stf(i).propVMAT.lastDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.timeFracFromLastDAO.*stf(i).propVMAT.doseAngleBordersDiff./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBordersDiff;
            apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOIndex,i) = stf(stf(i).propVMAT.nextDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.timeFracFromNextDAO.*stf(i).propVMAT.doseAngleBordersDiff./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBordersDiff;
        end
    end
end

% save global data
apertureInfo.runVMAT            = pln.propOpt.runVMAT;
apertureInfo.preconditioner     = pln.propOpt.preconditioner;
apertureInfo.run4D              = pln.propOpt.run4D;
apertureInfo.varOpt             = pln.propOpt.varOpt;
apertureInfo.numPhases          = sequencing.numPhases;
apertureInfo.bixelWidth         = bixelWidth;
apertureInfo.numOfMLCLeafPairs  = numOfMLCLeafPairs;
apertureInfo.totalNumOfBixels   = totalNumOfBixels;
apertureInfo.totalNumOfShapes   = sum([apertureInfo.beam.numOfShapes]);

% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later
optBixelFactor = 7;
% For optimized beams: 7 = (1 from weights) + (3 from left leaf positions (I, M, and F)) + (3 from
% right leaf positions (I, M, and F))

if isfield(sequencing,'weightToMU')
    apertureInfo.weightToMU = sequencing.weightToMU;
end
if pln.propOpt.runVMAT
    
    tempStruct_beam = apertureInfo.propVMAT.beam;
    tempStruct_jacobT = apertureInfo.propVMAT.jacobT;
    apertureInfo.propVMAT = pln.propOpt.VMAToptions;
    apertureInfo.propVMAT.beam = tempStruct_beam;
    apertureInfo.propVMAT.jacobT = tempStruct_jacobT;
    
    apertureInfo.totalNumOfOptBixels = totalNumOfOptBixels;
    apertureInfo.doseTotalNumOfLeafPairs = sum([apertureInfo.beam(:).numOfActiveLeafPairs]);
    
    % individual size of jacobian vector
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after), then add 2 (influenced by the time of the
    % times before and after)
    intBixelFactor = 2*optBixelFactor+2;
    % for the time (probability) gradients
    optBixelFactor = optBixelFactor+apertureInfo.totalNumOfShapes;
    intBixelFactor = intBixelFactor+apertureInfo.totalNumOfShapes;
    
    if apertureInfo.propVMAT.continuousAperture
        apertureInfo.totalNumOfLeafPairs = sum(reshape([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).doseAngleDAO],2,[]),1)*[apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).numOfActiveLeafPairs]';
        
        shapeInd = 1;
        for i = 1:numel(apertureInfo.beam)
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                %apertureInfo.propVMAT.beam(i).initialLeftLeafInd = apertureInfo.beam(i).shape{1}(1).vectorOffset:(apertureInfo.beam(i).shape{1}(1).vectorOffset+dimZ-1);
                %apertureInfo.propVMAT.beam(i).finalLeftLeafInd = (apertureInfo.beam(i).shape{1}(1).vectorOffset+dimZ):(apertureInfo.beam(i).shape{1}(1).vectorOffset+2*dimZ-1);
                
                %apertureInfo.propVMAT.beam(i).initialRightLeafInd = apertureInfo.propVMAT.beam(i).initialLeftLeafInd+apertureInfo.totalNumOfLeafPairs;
                %apertureInfo.propVMAT.beam(i).finalRightLeafInd = apertureInfo.propVMAT.beam(i).finalLeftLeafInd+apertureInfo.totalNumOfLeafPairs;
                
                apertureInfo.propVMAT.beam(i).timeInd = repmat(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd,1,apertureInfo.beam(1).numOfActiveLeafPairs);
                
                if apertureInfo.propVMAT.beam(i).timeFac(1) ~= 0
                    % some of the time for this optimized angle "leaks"
                    % into a previous dose angle arc (where we need to
                    % constrain leaf speed)
                    
                    % the initial leaves here act as final leaves for the
                    % previous angle
                    %finalLeftLeafInd = [apertureInfo.propVMAT.beam(i).initialLeftLeafInd apertureInfo.propVMAT.beam(i).finalLeftLeafInd];
                    %finalRightLeafInd = [apertureInfo.propVMAT.beam(i).initialRightLeafInd apertureInfo.propVMAT.beam(i).finalRightLeafInd];
                    
                    % part of the optimized time is used for the leaf speed
                    % calc
                    apertureInfo.propVMAT.beam(i).timeInd = [repmat(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd,1,apertureInfo.beam(1).numOfActiveLeafPairs) apertureInfo.propVMAT.beam(i).timeInd];
                else
                    %finalLeftLeafInd = apertureInfo.propVMAT.beam(i).finalLeftLeafInd;
                    %finalRightLeafInd = apertureInfo.propVMAT.beam(i).finalRightLeafInd;
                end
                
                if apertureInfo.propVMAT.beam(i).timeFac(3) ~= 0
                    % some of the time for this optimized angle "leaks"
                    % into a next dose angle arc (where we need to
                    % constrain leaf speed)
                    
                    %initialLeftLeafInd = [apertureInfo.propVMAT.beam(i).initialLeftLeafInd apertureInfo.propVMAT.beam(i).finalLeftLeafInd];
                    %initialRightLeafInd = [apertureInfo.propVMAT.beam(i).initialRightLeafInd apertureInfo.propVMAT.beam(i).finalRightLeafInd];
                    
                    % part of the optimized time is used for the leaf speed
                    % calc
                    apertureInfo.propVMAT.beam(i).timeInd = [apertureInfo.propVMAT.beam(i).timeInd repmat(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd,1,apertureInfo.beam(1).numOfActiveLeafPairs)];
                else
                    %initialLeftLeafInd = apertureInfo.propVMAT.beam(i).initialLeftLeafInd;
                    %initialRightLeafInd = apertureInfo.propVMAT.beam(i).initialRightLeafInd;
                end
                
                %apertureInfo.propVMAT.beam(i).initialLeftLeafInd = initialLeftLeafInd;
                %apertureInfo.propVMAT.beam(i).initialRightLeafInd = initialRightLeafInd;
                %apertureInfo.propVMAT.beam(i).finalLeftLeafInd = finalLeftLeafInd;
                %apertureInfo.propVMAT.beam(i).finalRightLeafInd = finalRightLeafInd;
                
                shapeInd = shapeInd+1;
            end
        end
        
    else
        apertureInfo.totalNumOfLeafPairs = totalNumOfLeafPairs;
    end
    
    % fix instances of leaf touching
    apertureInfo = matRad_leafTouching(apertureInfo);
    
    if pln.propOpt.run4D
        apertureInfo = matRad_doDAD(apertureInfo,stf);
        
        % store transition probabilities in apertureInfo
        if isfield(pln.propOpt.prop4D,'motionModel')
            
            % use model if it exists
            apertureInfo.motionModel        = pln.propOpt.prop4D.motionModel;
            apertureInfo.motionModel.type   = 'Markov';
            
            % strip model of the "out of bounds" phases
            apertureInfo.motionModel = matRad_stripMarkovOOB(apertureInfo.motionModel);
            
            % determine initial position phase
            posPhaseProb = accumarray(apertureInfo.motionModel.indices.subPhase2PosPhase,apertureInfo.motionModel.Pi_deltaTSample);
            initPosPhase = find(posPhaseProb == max(posPhaseProb));
            
            % determine initial probabililty
            apertureInfo.motionModel.initProb = zeros(1,apertureInfo.motionModel.indices.nSubPhases);
            % trigger on first phase and EOE
            initSubPhases   = apertureInfo.motionModel.indices.subPhase2PosPhase == initPosPhase & apertureInfo.motionModel.indices.subPhase2FS == 3;
            % let all subphases corresponding to the initial position phase
            % and FS (EOE) have the same probability
            apertureInfo.motionModel.initProb(initSubPhases) = 1;
            % let all other subphases have the same (much smaller, but
            % nonzero) probability
            apertureInfo.motionModel.initProb(~initSubPhases) = 1e-8;
            % normalize
            apertureInfo.motionModel.initProb = apertureInfo.motionModel.initProb./sum(apertureInfo.motionModel.initProb);
            
        else
            
            apertureInfo.motionModel.type   = 'Markov';
            apertureInfo.motionModel.qij    = rand(apertureInfo.numPhases,apertureInfo.numPhases)/2;
            for i = 1:apertureInfo.numPhases
                apertureInfo.motionModel.qij(i,i) = 0;
                apertureInfo.motionModel.qij(i,i) = -sum(apertureInfo.motionModel.qij(i,:),2);
            end
            apertureInfo.motionModel.initProb  = rand(1,apertureInfo.numPhases);
            apertureInfo.motionModel.initProb  = apertureInfo.motionModel.initProb./sum(apertureInfo.motionModel.initProb);
            
        end
        
    else
        apertureInfo.motionModel.type   = 'Markov';
        apertureInfo.motionModel.qij       = 0;
        apertureInfo.motionModel.initProb  = 1;
        
        apertureInfo.motionModel.indices.subPhase2PosPhase_gridI    = 1;
        apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ    = 1;
        apertureInfo.motionModel.indices.nSubPhasePerPosPhase       = 1;
        apertureInfo.motionModel.indices.nSubPhases                 = 1;
    end
    
    fileName = apertureInfo.propVMAT.machineConstraintFile;
    try
        load(fileName,'machine');
    catch
        error(['Could not find the following machine file: ' fileName ]);
    end
    
    interpGetsTrans = false;
    % prep for gradient offset
    gradOffset      = 0;
    for i = 1:numel(apertureInfo.beam)
        
        maxTime = apertureInfo.propVMAT.beam(i).doseAngleBordersDiff./machine.constraints.gantryRotationSpeed(1);
        [Pij_transT,~,~,~] = matRad_transAndTProb(maxTime,0,apertureInfo.motionModel);
        PIJ_transT = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT(:))./apertureInfo.motionModel.indices.nSubPhasePerPosPhase;
        transMask = PIJ_transT > 0.01;
        
        apertureInfo.propVMAT.beam(i).transMask             = repmat(1:apertureInfo.numPhases,[apertureInfo.numPhases 1]);
        apertureInfo.propVMAT.beam(i).transMask(~transMask) = 0;
        
        if apertureInfo.propVMAT.continuousAperture
            if apertureInfo.propVMAT.beam(i).DAOBeam || interpGetsTrans
                
                apertureInfo.propVMAT.beam(i).leafConstMask             = repmat(1:apertureInfo.numPhases,[apertureInfo.numPhases 1]);
                apertureInfo.propVMAT.beam(i).leafConstMask(~transMask) = 0;
            end
            if apertureInfo.propVMAT.beam(i).DAOBeam
                interpGetsTrans = apertureInfo.propVMAT.beam(i).timeFac(3) ~= 0;
            else
                interpGetsTrans = false;
            end
            
            % count number of transitions
            apertureInfo.propVMAT.numLeafSpeedConstraint      = nnz([apertureInfo.propVMAT.beam.leafConstMask]);
            apertureInfo.propVMAT.numLeafSpeedConstraintDAO   = nnz([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).leafConstMask]);
            
        end
        
        % individual size of jacobian vector
        if apertureInfo.propVMAT.beam(i).DAOBeam
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*optBixelFactor.*apertureInfo.numPhases;
            apertureInfo.beam(i).numUniqueVar   = apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases.*(1+4.*apertureInfo.beam(i).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes;%apertureInfo.beam(i).numOfShapes.*(1+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            
            % conversion from local opt variable to global opt variable
            apertureInfo.beam(i).local2GlobalVar = zeros(apertureInfo.beam(i).numUniqueVar,1);
            
            for phase = 1:apertureInfo.numPhases
                % weight variables
                apertureInfo.beam(i).local2GlobalVar(phase) = apertureInfo.propVMAT.beam(i).DAOIndex+(phase-1)*apertureInfo.totalNumOfShapes;
                % leaf position variables
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(1) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(2) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(1) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(2) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                % time variables
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.beam(i).numUniqueVar-apertureInfo.totalNumOfShapes+(1:apertureInfo.totalNumOfShapes)) = ...
                    (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);
            end
            
            % keep weight and time variables when calculating gradients for
            % d2
            apertureInfo.beam(i).d2KeepVar = true(apertureInfo.beam(i).numUniqueVar,1);
            %apertureInfo.beam(i).d2KeepVar((apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases+1):end) = false;
            %apertureInfo.beam(i).d2KeepVar((apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases+1):(end-apertureInfo.totalNumOfShapes)) = false;
            
        else
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*intBixelFactor.*apertureInfo.numPhases;
            apertureInfo.beam(i).numUniqueVar   = apertureInfo.numPhases.*(2+4.*apertureInfo.beam(i).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes;%(2+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            
            for phase = 1:apertureInfo.numPhases
                % weight variables
                apertureInfo.beam(i).local2GlobalVar(2.*(phase-1)+(1:2)) = [apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex]+(phase-1)*apertureInfo.totalNumOfShapes;
                % leaf position variables
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(2) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(1) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(2) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(1) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                % time variables
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.beam(i).numUniqueVar-apertureInfo.totalNumOfShapes+(1:apertureInfo.totalNumOfShapes)) = ...
                    (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);
            end
            
            % keep weight and time variables when calculating gradients for
            % d2
            apertureInfo.beam(i).d2KeepVar = true(apertureInfo.beam(i).numUniqueVar,1);
            %apertureInfo.beam(i).d2KeepVar((2.*apertureInfo.numPhases+1):end) = false;
            %apertureInfo.beam(i).d2KeepVar((2.*apertureInfo.numPhases+1):(end-apertureInfo.totalNumOfShapes)) = false;
        end
        
        apertureInfo.beam(i).d2KeepVar = find(apertureInfo.beam(i).d2KeepVar);
        
        apertureInfo.beam(i).numKeepVar = numel(apertureInfo.beam(i).d2KeepVar);
        
        apertureInfo.beam(i).gradOffset = gradOffset;
        gradOffset = gradOffset+apertureInfo.beam(i).numKeepVar.*apertureInfo.numPhases.^2;
    end
    
else
    apertureInfo.totalNumOfLeafPairs = sum([apertureInfo.beam.numOfShapes]*[apertureInfo.beam.numOfActiveLeafPairs]');
    apertureInfo.doseTotalNumOfLeafPairs = apertureInfo.totalNumOfLeafPairs;
    
    % individual size of jacobian vector
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after)
    intBixelFactor = 2*optBixelFactor;
    
    for i = 1:numel(apertureInfo.beam)
        if apertureInfo.propVMAT.beam(i).DAOBeam
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*optBixelFactor.*apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases;
        else
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*intBixelFactor.*apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases;
        end
    end
end

% global size of jacobian vector
apertureInfo.bixelJApVec_sz = sum([apertureInfo.beam.bixelJApVec_sz]);
% apertureInfo.bixelJApVec_sz = (apertureInfo.totalNumOfOptBixels*optBixelFactor+(apertureInfo.totalNumOfBixels-apertureInfo.totalNumOfOptBixels)*intBixelFactor)*apertureInfo.numPhases;

% create vectors for optimization
[apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] = matRad_daoApertureInfo2Vec(apertureInfo);

end






