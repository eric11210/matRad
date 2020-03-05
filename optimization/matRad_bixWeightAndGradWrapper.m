function [apertureInfo,results] = matRad_bixWeightAndGradWrapper(apertureInfo,i,results)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad wrapper function to calculate the bixel weights from the aperture vector,
% and also the Jacobian matrix relating these two.
%
% call
%   [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMapW,counters] = matRad_bixWeightAndGradWrapper(apertureInfo,sumGradSq,shapeMapW)
%
% input
%
% output
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

round2 = @(a,b) round(a*10^b)/10^b;

%% precalculate certain quantities
if apertureInfo.runVMAT
    numOfShapes = 1;
else
    numOfShapes = updatedInfo.beam(i).numOfShapes;
end

lastFixelIndMap         = apertureInfo.beam(i).lastFixelIndMap-max(apertureInfo.beam(i).lastFixelIndMap(:))+apertureInfo.beam(i).lastNumBixels;
nextFixelIndMap         = apertureInfo.beam(i).nextFixelIndMap-max(apertureInfo.beam(i).nextFixelIndMap(:))+apertureInfo.beam(i).nextNumBixels;% +apertureInfo.beam(i).lastNumBixels
bixelWidth              = apertureInfo.bixelWidth;
lim_l                   = apertureInfo.beam(i).lim_l;
lim_r                   = apertureInfo.beam(i).lim_r;
edges_l                 = apertureInfo.beam(i).posOfCornerBixel(1)...
    + ((1:size(lastFixelIndMap,2))-1-1/2)*bixelWidth;
edges_r                 = apertureInfo.beam(i).posOfCornerBixel(1)...
    + ((1:size(lastFixelIndMap,2))-1+1/2)*bixelWidth;
centres                 = (edges_l+edges_r)/2;
widths                  = edges_r-edges_l;
numRow                  = apertureInfo.beam(i).numOfActiveLeafPairs;
numCol                  = size(lastFixelIndMap,2);
lastNumBix              = apertureInfo.beam(i).lastNumBixels;
nextNumBix              = apertureInfo.beam(i).nextNumBixels;
effNumBix_lastDose_arcI = apertureInfo.beam(i).effNumBixels_lastDose_arcI;
effNumBix_lastDose_arcF = apertureInfo.beam(i).effNumBixels_lastDose_arcF;
effNumBix_nextDose_arcI = apertureInfo.beam(i).effNumBixels_nextDose_arcI;
effNumBix_nextDose_arcF = apertureInfo.beam(i).effNumBixels_nextDose_arcF;
totalNumOfShapes        = apertureInfo.totalNumOfShapes;
full2UniqueLocalVar     = apertureInfo.beam(i).full2UniqueLocalVar;

% indices for various variables
tIx_Vec     = (1:totalNumOfShapes)+apertureInfo.beam(i).numFullVar-totalNumOfShapes;%(totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:totalNumOfShapes);

if ~apertureInfo.propVMAT.beam(i).DAOBeam
    time_last   = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time;
    time_next   = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time;
    time        = apertureInfo.beam(i).time;
    
    fluAngleBordersDiff         = apertureInfo.propVMAT.beam(i).fluAngleBordersDiff;
    fluAngleBordersDiff_last    = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).fluAngleBordersDiff;
    fluAngleBordersDiff_next    = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).fluAngleBordersDiff;
    timeFacCurr_last            = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr;
    timeFacCurr_next            = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr;
    fracFromLastDAO_MU          = apertureInfo.propVMAT.beam(i).fracFromLastDAO_MU;
    fracFromNextDAO_MU          = apertureInfo.propVMAT.beam(i).fracFromNextDAO_MU;
    fracFromLastDAO_gantryRot   = apertureInfo.propVMAT.beam(i).fracFromLastDAO_gantryRot;
    fracFromNextDAO_gantryRot   = apertureInfo.propVMAT.beam(i).fracFromNextDAO_gantryRot;
    
    DAOindex_last   = 1;%apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
    DAOindex_next   = 2;%apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
    tIx_last        = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex+apertureInfo.beam(i).numFullVar-totalNumOfShapes;%(totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
    tIx_next        = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex+apertureInfo.beam(i).numFullVar-totalNumOfShapes;%(totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
    
    if ~apertureInfo.propVMAT.continuousAperture
        fracFromLastDAOI_leafI = apertureInfo.propVMAT.beam(i).fracFromLastDAO;
        fracFromLastDAOF_leafI = apertureInfo.propVMAT.beam(i).fracFromLastDAO;
        fracFromNextDAOI_leafF = (1-apertureInfo.propVMAT.beam(i).fracFromLastDAO);
        fracFromNextDAOF_leafF = (1-apertureInfo.propVMAT.beam(i).fracFromLastDAO);
    end
end

if apertureInfo.propVMAT.continuousAperture
    fracFromLastDAOI_leafI = apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI;
    fracFromLastDAOF_leafI = apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI;
    fracFromNextDAOI_leafF = apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF;
    fracFromNextDAOF_leafF = apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF;
end

fracFromFluI_leafI_arcI = 1.*ones(numRow,1);
fracFromFluI_leafF_arcI = apertureInfo.propVMAT.beam(i).fluAngleBorderCentreDiff(2)./apertureInfo.propVMAT.beam(i).fluAngleBordersDiff.*ones(numRow,1);
fracFromFluF_leafI_arcI = 0.*ones(numRow,1);
fracFromFluF_leafF_arcI = apertureInfo.propVMAT.beam(i).fluAngleBorderCentreDiff(1)./apertureInfo.propVMAT.beam(i).fluAngleBordersDiff.*ones(numRow,1);
fracFromFluI_leafI_arcF = apertureInfo.propVMAT.beam(i).fluAngleBorderCentreDiff(2)./apertureInfo.propVMAT.beam(i).fluAngleBordersDiff.*ones(numRow,1);
fracFromFluI_leafF_arcF = 0.*ones(numRow,1);
fracFromFluF_leafI_arcF = apertureInfo.propVMAT.beam(i).fluAngleBorderCentreDiff(1)./apertureInfo.propVMAT.beam(i).fluAngleBordersDiff.*ones(numRow,1);
fracFromFluF_leafF_arcF = 1.*ones(numRow,1);

%% prep static

static.lim_r                    = lim_r;
static.edges_l                  = edges_l;
static.edges_r                  = edges_r;
static.centres                  = centres;
static.widths                   = widths;
static.lastBixelIndMap          = cast(lastFixelIndMap,'like',results.arcI.lastDose.bixelJApVec_j{1});
static.nextBixelIndMap          = cast(nextFixelIndMap,'like',results.arcI.nextDose.bixelJApVec_j{1});
static.numRow                   = numRow;
static.numCol                   = numCol;
static.lastNumBix               = lastNumBix;
static.nextNumBix               = nextNumBix;
static.DAOBeam                  = apertureInfo.propVMAT.beam(i).DAOBeam;
static.fixedGantrySpeed         = apertureInfo.propVMAT.fixedGantrySpeed;
static.totalNumOfBixels         = apertureInfo.totalNumOfBixels;
static.lastBixIndVec            = 1:lastNumBix;
static.nextBixIndVec            = 1:nextNumBix;
static.fracToLastDose_arcI      = apertureInfo.propVMAT.beam(i).fracToLastDose_arcI;
static.fracToNextDose_arcI      = apertureInfo.propVMAT.beam(i).fracToNextDose_arcI;
static.fracToLastDose_arcF      = apertureInfo.propVMAT.beam(i).fracToLastDose_arcF;
static.fracToNextDose_arcF      = apertureInfo.propVMAT.beam(i).fracToNextDose_arcF;
static.weightFactor_I           = apertureInfo.propVMAT.beam(i).fluAngleBorderCentreDiff(1)./apertureInfo.propVMAT.beam(i).fluAngleBordersDiff;
static.weightFactor_F           = apertureInfo.propVMAT.beam(i).fluAngleBorderCentreDiff(2)./apertureInfo.propVMAT.beam(i).fluAngleBordersDiff;
static.fracFromLastDAOI_leafI   = fracFromLastDAOI_leafI;
static.fracFromLastDAOF_leafI   = fracFromLastDAOF_leafI;
static.fracFromNextDAOI_leafF   = fracFromNextDAOI_leafF;
static.fracFromNextDAOF_leafF   = fracFromNextDAOF_leafF;

if apertureInfo.propVMAT.beam(i).DAOBeam
    
    numVarMult              = 9;
else
    
    static.fracFromLastDAO_MU           = fracFromLastDAO_MU;
    static.fracFromNextDAO_MU           = fracFromNextDAO_MU;
    static.fracFromLastDAO_gantryRot    = fracFromLastDAO_gantryRot;
    static.fracFromNextDAO_gantryRot    = fracFromNextDAO_gantryRot;
    
    static.time         = time;
    static.time_next    = time_next;
    static.time_last    = time_last;
    
    static.tIx_last = tIx_last;
    static.tIx_next = tIx_next;
    
    static.fluAngleBordersDiff      = fluAngleBordersDiff;
    static.fluAngleBordersDiff_last = fluAngleBordersDiff_last;
    static.fluAngleBordersDiff_next = fluAngleBordersDiff_next;
    
    static.timeFacCurr_last = timeFacCurr_last;
    static.timeFacCurr_next = timeFacCurr_next;
    
    if apertureInfo.propVMAT.fixedGantrySpeed
        numVarMult = 10;
    else
        numVarMult = 12;
    end
end

% determine probabilities
[Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(apertureInfo.beam(i).time,sum([apertureInfo.beam(1:(i-1)).time]),apertureInfo.motionModel);

% determine probability and gradient matrices
pMat_sub    = Pi_T'.*Pij_transT;
pMat        = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],pMat_sub(:));

pGradMat1_sub       = Pi_T_dot'.*Pij_transT;
pGradMat1           = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],pGradMat1_sub(:));
pGradMat1_part      = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) repelem((1:apertureInfo.motionModel.indices.nSubPhases)',apertureInfo.motionModel.indices.nSubPhases)],pGradMat1_sub(:));
pGradMat2_sub       = Pi_T'.*Pij_transT_dot;
pGradMat2           = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],pGradMat2_sub(:));
pGradMat2_part      = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) repelem((1:apertureInfo.motionModel.indices.nSubPhases)',apertureInfo.motionModel.indices.nSubPhases)],pGradMat2_sub(:));
PijGrad_transT      = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT_dot(:));
PijGrad_transT_part = accumarray([repmat((1:apertureInfo.motionModel.indices.nSubPhases)',[apertureInfo.motionModel.indices.nSubPhases 1]) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT_dot(:));

pGradMat                = reshape(sum(apertureInfo.propVMAT.jacobT(:,1:(i-1)),2),[1 1 apertureInfo.totalNumOfShapes]).*pGradMat1 + reshape(apertureInfo.propVMAT.jacobT(:,i),[1 1 apertureInfo.totalNumOfShapes]).*pGradMat2;
pGradMat_part           = reshape(sum(apertureInfo.propVMAT.jacobT(:,1:(i-1)),2),[1 1 apertureInfo.totalNumOfShapes]).*pGradMat1_part + reshape(apertureInfo.propVMAT.jacobT(:,i),[1 1 apertureInfo.totalNumOfShapes]).*pGradMat2_part;
PijGradMat_transT       = reshape(apertureInfo.propVMAT.jacobT(:,i),[1 1 apertureInfo.totalNumOfShapes]).*PijGrad_transT;
PijGradMat_transT_part  = reshape(apertureInfo.propVMAT.jacobT(:,i),[1 1 apertureInfo.totalNumOfShapes]).*PijGrad_transT_part;

% put probability and gradient matrices into struct
apertureInfo.probI_IJ{i}        = pMat;
apertureInfo.probIGrad_IJ{i}    = pGradMat;

apertureInfo.probF_KL{i}        = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT(:));
apertureInfo.probFGrad_KL{i}    = PijGradMat_transT;

apertureInfo.probI_Ij{i}        = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) repelem((1:apertureInfo.motionModel.indices.nSubPhases)',apertureInfo.motionModel.indices.nSubPhases)],pMat_sub(:));
apertureInfo.probIGrad_Ij{i}    = pGradMat_part;

apertureInfo.probF_kL{i}        = accumarray([repmat((1:apertureInfo.motionModel.indices.nSubPhases)',[apertureInfo.motionModel.indices.nSubPhases 1]) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT(:));
apertureInfo.probFGrad_kL{i}    = PijGradMat_transT_part;

% loop over all (initial) phases
for phase_I = 1:apertureInfo.numPhases
    
    % loop over all (final) phases
    for phase_F = 1:apertureInfo.numPhases
        
        % determine cell index
        cellInd = (phase_I-1).*apertureInfo.numPhases+phase_F;
        
        % skip to next phase if probability to transition from phase_I
        % to phase_F and its derivative are both 0
        if pMat(phase_I,phase_F)  == 0 && pGradMat1(phase_I,phase_F) == 0 && pGradMat2(phase_I,phase_F) == 0
            
            % cut out zeros, i.e. all elements after the offset
            results.arcI.lastDose.bixelJApVec_i{cellInd}((results.arcI.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcI.lastDose.bixelJApVec_j{cellInd}((results.arcI.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcI.lastDose.bixelJApVec_vec{cellInd}((results.arcI.lastDose.bixelJApVec_offset{cellInd}+1):end) = [];
            results.arcI.nextDose.bixelJApVec_i{cellInd}((results.arcI.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcI.nextDose.bixelJApVec_j{cellInd}((results.arcI.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcI.nextDose.bixelJApVec_vec{cellInd}((results.arcI.nextDose.bixelJApVec_offset{cellInd}+1):end) = [];
            results.arcF.lastDose.bixelJApVec_i{cellInd}((results.arcF.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcF.lastDose.bixelJApVec_j{cellInd}((results.arcF.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcF.lastDose.bixelJApVec_vec{cellInd}((results.arcF.lastDose.bixelJApVec_offset{cellInd}+1):end) = [];
            results.arcF.nextDose.bixelJApVec_i{cellInd}((results.arcF.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcF.nextDose.bixelJApVec_j{cellInd}((results.arcF.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
            results.arcF.nextDose.bixelJApVec_vec{cellInd}((results.arcF.nextDose.bixelJApVec_offset{cellInd}+1):end) = [];
            
            continue
        end
        
        
        
        
        % loop over all shapes
        for j = 1:numOfShapes
            %% determine variable quantities
            
            weight_I        = apertureInfo.beam(i).shape{phase_I}(j).weight;
            weight_F        = apertureInfo.beam(i).shape{phase_F}(j).weight;
            leftLeafPos_I   = apertureInfo.beam(i).shape{phase_I}(j).leftLeafPos_I;
            leftLeafPos_F   = apertureInfo.beam(i).shape{phase_F}(j).leftLeafPos_F;
            rightLeafPos_I  = apertureInfo.beam(i).shape{phase_I}(j).rightLeafPos_I;
            rightLeafPos_F  = apertureInfo.beam(i).shape{phase_F}(j).rightLeafPos_F;
            
            %% sort out order, set up calculation of bixel weight and gradients
            
            %leftLeafPosI = round2(leftLeafPosI,10);
            %leftLeafPosF = round2(leftLeafPosF,10);
            %rightLeafPosI = round2(rightLeafPosI,10);
            %rightLeafPosF = round2(rightLeafPosF,10);
            
            leftLeafPos_I(leftLeafPos_I <= lim_l) = lim_l(leftLeafPos_I <= lim_l);
            leftLeafPos_F(leftLeafPos_F <= lim_l) = lim_l(leftLeafPos_F <= lim_l);
            rightLeafPos_I(rightLeafPos_I <= lim_l) = lim_l(rightLeafPos_I <= lim_l);
            rightLeafPos_F(rightLeafPos_F <= lim_l) = lim_l(rightLeafPos_F <= lim_l);
            leftLeafPos_I(leftLeafPos_I >= lim_r) = lim_r(leftLeafPos_I >= lim_r);
            leftLeafPos_F(leftLeafPos_F >= lim_r) = lim_r(leftLeafPos_F >= lim_r);
            rightLeafPos_I(rightLeafPos_I >= lim_r) = lim_r(rightLeafPos_I >= lim_r);
            rightLeafPos_F(rightLeafPos_F >= lim_r) = lim_r(rightLeafPos_F >= lim_r);
            
            % determine middle leaf positions
            leftLeafPos_M   = fracFromFluI_leafF_arcI.*leftLeafPos_I+fracFromFluF_leafF_arcI.*leftLeafPos_F;
            rightLeafPos_M  = fracFromFluI_leafF_arcI.*rightLeafPos_I+fracFromFluF_leafF_arcI.*rightLeafPos_F;
            
            % determine initial and final leaf positions for initial and
            % final arcs
            % set the initial leaf positions to the minimum leaf positions
            % always, instead of the leaf positions at the actual beginning
            % of the arc
            % this simplifies the calculation
            % remember which one is actually I and F in left/right MinInd
            [leftLeafPosI_arcI,leftMinInd]      = min([leftLeafPos_I,leftLeafPos_M],[],2);
            leftLeafPosF_arcI                   = max([leftLeafPos_I,leftLeafPos_M],[],2);
            leftLeafPosI_arcF                   = min([leftLeafPos_M,leftLeafPos_F],[],2);
            leftLeafPosF_arcF                   = max([leftLeafPos_M,leftLeafPos_F],[],2);
            [rightLeafPosI_arcI,rightMinInd]    = min([rightLeafPos_I,rightLeafPos_M],[],2);
            rightLeafPosF_arcI                  = max([rightLeafPos_I,rightLeafPos_M],[],2);
            rightLeafPosI_arcF                  = min([rightLeafPos_M,rightLeafPos_F],[],2);
            rightLeafPosF_arcF                  = max([rightLeafPos_M,rightLeafPos_F],[],2);
            
            % find bixel indices where leaves are located
            xPosIndLeftLeafI_arcI   = min(floor((leftLeafPosI_arcI-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndLeftLeafF_arcI   = max(ceil((leftLeafPosF_arcI-edges_r(1))./bixelWidth)+1,1);
            xPosIndLeftLeafI_arcF   = min(floor((leftLeafPosI_arcF-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndLeftLeafF_arcF   = max(ceil((leftLeafPosF_arcF-edges_r(1))./bixelWidth)+1,1);
            xPosIndRightLeafI_arcI  = min(floor((rightLeafPosI_arcI-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndRightLeafF_arcI  = max(ceil((rightLeafPosF_arcI-edges_r(1))./bixelWidth)+1,1);
            xPosIndRightLeafI_arcF  = min(floor((rightLeafPosI_arcF-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndRightLeafF_arcF  = max(ceil((rightLeafPosF_arcF-edges_r(1))./bixelWidth)+1,1);
            %
            xPosLinearIndLeftLeafI_arcI     = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafI_arcI);
            xPosLinearIndLeftLeafF_arcI     = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafF_arcI);
            xPosLinearIndLeftLeafI_arcF     = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafI_arcF);
            xPosLinearIndLeftLeafF_arcF     = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafF_arcF);
            xPosLinearIndRightLeafI_arcI    = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafI_arcI);
            xPosLinearIndRightLeafF_arcI    = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafF_arcI);
            xPosLinearIndRightLeafI_arcF    = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafI_arcF);
            xPosLinearIndRightLeafF_arcF    = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafF_arcF);
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                jacobiScale_I = apertureInfo.beam(i).shape{phase_I}(1).jacobiScale;
                jacobiScale_F = apertureInfo.beam(i).shape{phase_F}(1).jacobiScale;
                
                weightOffset = apertureInfo.numPhases;
            else
                
                weight_last_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).weight;
                weight_last_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).weight;
                weight_next_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).weight;
                weight_next_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).weight;
                
                jacobiScale_last_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(1).jacobiScale;
                jacobiScale_last_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(1).jacobiScale;
                jacobiScale_next_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(1).jacobiScale;
                jacobiScale_next_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(1).jacobiScale;
                
                weightOffset = 2.*apertureInfo.numPhases;
            end
            
            %if apertureInfo.propVMAT.continuousAperture
            variable.vectorIx_L_lastDAOI = weightOffset+4.*numRow.*(phase_I-1)+0.*numRow+(1:numRow);%apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(1).vectorOffset(1) + ((1:numRow)-1);
            variable.vectorIx_L_lastDAOF = weightOffset+4.*numRow.*(phase_I-1)+1.*numRow+(1:numRow);%apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).vectorOffset(2) + ((1:numRow)-1);
            variable.vectorIx_L_nextDAOI = weightOffset+4.*numRow.*(phase_F-1)+2.*numRow+(1:numRow);%apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(1).vectorOffset(1) + ((1:numRow)-1);
            variable.vectorIx_L_nextDAOF = weightOffset+4.*numRow.*(phase_F-1)+3.*numRow+(1:numRow);%apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).vectorOffset(2) + ((1:numRow)-1);
            
            variable.vectorIx_R_lastDAOI = variable.vectorIx_L_lastDAOI+4.*numRow.*apertureInfo.numPhases;%vectorIx_L_lastDAOI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
            variable.vectorIx_R_lastDAOF = variable.vectorIx_L_lastDAOF+4.*numRow.*apertureInfo.numPhases;%vectorIx_L_lastDAOF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
            variable.vectorIx_R_nextDAOI = variable.vectorIx_L_nextDAOI+4.*numRow.*apertureInfo.numPhases;%vectorIx_L_nextDAOI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
            variable.vectorIx_R_nextDAOF = variable.vectorIx_L_nextDAOF+4.*numRow.*apertureInfo.numPhases;%vectorIx_L_nextDAOF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
            %else
            %    vectorIx_LF_last  = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).vectorOffset + ((1:numRow)-1);
            %    vectorIx_LI_next  = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).vectorOffset + ((1:numRow)-1);
            %    vectorIx_RF_last  = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
            %    vectorIx_RI_next  = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
            %end
            
            % switch the arcI and arcF leaf factors where the initial and
            % final leaf positions have been switched
            % we need to remember how the sub-arc leaf positions have been
            % influenced by the leaf positions at the actual beginning and
            % end of the subarc
            variable.fracFromFluI_leftLeafI_arcI = fracFromFluI_leafI_arcI;
            variable.fracFromFluI_leftLeafF_arcI = fracFromFluI_leafF_arcI;
            variable.fracFromFluF_leftLeafI_arcI = fracFromFluF_leafI_arcI;
            variable.fracFromFluF_leftLeafF_arcI = fracFromFluF_leafF_arcI;
            variable.fracFromFluI_leftLeafI_arcF = fracFromFluI_leafI_arcF;
            variable.fracFromFluI_leftLeafF_arcF = fracFromFluI_leafF_arcF;
            variable.fracFromFluF_leftLeafI_arcF = fracFromFluF_leafI_arcF;
            variable.fracFromFluF_leftLeafF_arcF = fracFromFluF_leafF_arcF;
            
            variable.fracFromFluI_rightLeafI_arcI = fracFromFluI_leafI_arcI;
            variable.fracFromFluI_rightLeafF_arcI = fracFromFluI_leafF_arcI;
            variable.fracFromFluF_rightLeafI_arcI = fracFromFluF_leafI_arcI;
            variable.fracFromFluF_rightLeafF_arcI = fracFromFluF_leafF_arcI;
            variable.fracFromFluI_rightLeafI_arcF = fracFromFluI_leafI_arcF;
            variable.fracFromFluI_rightLeafF_arcF = fracFromFluI_leafF_arcF;
            variable.fracFromFluF_rightLeafI_arcF = fracFromFluF_leafI_arcF;
            variable.fracFromFluF_rightLeafF_arcF = fracFromFluF_leafF_arcF;
            
            temp_fracFromFluI_leftLeafI_arcI    = variable.fracFromFluI_leftLeafI_arcI;
            temp_fracFromFluF_leftLeafI_arcI    = variable.fracFromFluF_leftLeafI_arcI;
            temp_fracFromFluI_leftLeafI_arcF    = variable.fracFromFluI_leftLeafI_arcF;
            temp_fracFromFluF_leftLeafI_arcF    = variable.fracFromFluF_leftLeafI_arcF;
            temp_fracFromFluI_rightLeafI_arcI   = variable.fracFromFluI_rightLeafI_arcI;
            temp_fracFromFluF_rightLeafI_arcI   = variable.fracFromFluF_rightLeafI_arcI;
            temp_fracFromFluI_rightLeafI_arcF   = variable.fracFromFluI_rightLeafI_arcF;
            temp_fracFromFluF_rightLeafI_arcF   = variable.fracFromFluF_rightLeafI_arcF;
            
            variable.fracFromFluI_leftLeafI_arcI(leftMinInd == 2) = variable.fracFromFluI_leftLeafF_arcI(leftMinInd == 2);
            variable.fracFromFluI_leftLeafF_arcI(leftMinInd == 2) = temp_fracFromFluI_leftLeafI_arcI(leftMinInd == 2);
            variable.fracFromFluF_leftLeafI_arcI(leftMinInd == 2) = variable.fracFromFluF_leftLeafF_arcI(leftMinInd == 2);
            variable.fracFromFluF_leftLeafF_arcI(leftMinInd == 2) = temp_fracFromFluF_leftLeafI_arcI(leftMinInd == 2);
            variable.fracFromFluI_leftLeafI_arcF(leftMinInd == 2) = variable.fracFromFluI_leftLeafF_arcF(leftMinInd == 2);
            variable.fracFromFluI_leftLeafF_arcF(leftMinInd == 2) = temp_fracFromFluI_leftLeafI_arcF(leftMinInd == 2);
            variable.fracFromFluF_leftLeafI_arcF(leftMinInd == 2) = variable.fracFromFluF_leftLeafF_arcF(leftMinInd == 2);
            variable.fracFromFluF_leftLeafF_arcF(leftMinInd == 2) = temp_fracFromFluF_leftLeafI_arcF(leftMinInd == 2);
            
            variable.fracFromFluI_rightLeafI_arcI(rightMinInd == 2) = variable.fracFromFluI_rightLeafF_arcI(rightMinInd == 2);
            variable.fracFromFluI_rightLeafF_arcI(rightMinInd == 2) = temp_fracFromFluI_rightLeafI_arcI(rightMinInd == 2);
            variable.fracFromFluF_rightLeafI_arcI(rightMinInd == 2) = variable.fracFromFluF_rightLeafF_arcI(rightMinInd == 2);
            variable.fracFromFluF_rightLeafF_arcI(rightMinInd == 2) = temp_fracFromFluF_rightLeafI_arcI(rightMinInd == 2);
            variable.fracFromFluI_rightLeafI_arcF(rightMinInd == 2) = variable.fracFromFluI_rightLeafF_arcF(rightMinInd == 2);
            variable.fracFromFluI_rightLeafF_arcF(rightMinInd == 2) = temp_fracFromFluI_rightLeafI_arcF(rightMinInd == 2);
            variable.fracFromFluF_rightLeafI_arcF(rightMinInd == 2) = variable.fracFromFluF_rightLeafF_arcF(rightMinInd == 2);
            variable.fracFromFluF_rightLeafF_arcF(rightMinInd == 2) = temp_fracFromFluF_rightLeafI_arcF(rightMinInd == 2);
            
            %% determine probabilities and derivatives
            
            %probability         = Pi_T(phase_I).*Pij_transT(phase_I,phase_F);
            variable.probability = pMat(phase_I,phase_F);
            
            if apertureInfo.propVMAT.fixedGantrySpeed
                % for fixed gantry speed, set totalNumOfShapes (the number
                % of time variables) to 0
                variable.totalNumOfShapes       = 0;
                variable.lastNumShapbixIndVec   = 0;
                variable.nextNumShapbixIndVec   = 0;
            else
                % we only need derivatives for variable gantry speed
                
                %probability_dTVec   = sum(apertureInfo.propVMAT.jacobT(:,1:(i-1)),2).*Pi_T_dot(phase_I).*Pij_transT(phase_I,phase_F) + apertureInfo.propVMAT.jacobT(:,i).*Pi_T(phase_I).*Pij_transT_dot(phase_I,phase_F);
                probability_dTVec   = squeeze(pGradMat(phase_I,phase_F,:));
                
                % delete any derivatives with value less than eps
                delInd = abs(probability_dTVec) < eps;
                probability_dTVec(delInd)   = [];
                
                variable.tIx_Vec                = cast(tIx_Vec,'like',results.arcI.lastDose.bixelJApVec_i{1});
                variable.tIx_Vec(delInd)        = [];
                variable.totalNumOfShapes       = numel(probability_dTVec);
                variable.lastNumShapbixIndVec   = 1:(variable.totalNumOfShapes*lastNumBix);
                variable.nextNumShapbixIndVec   = 1:(variable.totalNumOfShapes*nextNumBix);
                
                variable.probability_dTVec  = probability_dTVec;
            end
            %% first do phase_I
            
            variable.arcF = false;
            
            bixelJApVecLastDose_sz              = effNumBix_lastDose_arcI.*(numVarMult+variable.totalNumOfShapes);
            bixelJApVecNextDose_sz              = effNumBix_nextDose_arcI.*(numVarMult+variable.totalNumOfShapes);
            variable.bixelJApVecLastDose_zeros  = zeros(bixelJApVecLastDose_sz,1);
            variable.bixelJApVecNextDose_zeros  = zeros(bixelJApVecNextDose_sz,1);
            
            variable.leftLeafPosI   = leftLeafPosI_arcI;
            variable.leftLeafPosF   = leftLeafPosF_arcI;
            variable.rightLeafPosI  = rightLeafPosI_arcI;
            variable.rightLeafPosF  = rightLeafPosF_arcI;
            
            variable.xPosLinearIndLeftLeafI     = xPosLinearIndLeftLeafI_arcI;
            variable.xPosLinearIndLeftLeafF     = xPosLinearIndLeftLeafF_arcI;
            variable.xPosLinearIndRightLeafI    = xPosLinearIndRightLeafI_arcI;
            variable.xPosLinearIndRightLeafF    = xPosLinearIndRightLeafF_arcI;
            
            variable.xPosIndLeftLeafI   = xPosIndLeftLeafI_arcI;
            variable.xPosIndLeftLeafF   = xPosIndLeftLeafF_arcI;
            variable.xPosIndRightLeafI  = xPosIndRightLeafI_arcI;
            variable.xPosIndRightLeafF  = xPosIndRightLeafF_arcI;
            
            variable.phase              = phase_I;
            variable.weight             = weight_I;
            variable.weightFactor       = static.weightFactor_I;
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                variable.DAOindex       = j+(phase_I-1)*numOfShapes;
                
                variable.jacobiScale    = jacobiScale_I;
            else
                
                variable.DAOindex_last      = DAOindex_last+(phase_I-1)*2;
                variable.DAOindex_next      = DAOindex_next+(phase_I-1)*2;
                
                variable.weight_last        = weight_last_I;
                variable.weight_next        = weight_next_I;
                variable.jacobiScale_last   = jacobiScale_last_I;
                variable.jacobiScale_next   = jacobiScale_next_I;
            end
            
            % calculate bixel weight and gradient
            tempResults = matRad_bixWeightAndGrad(static,variable);
            %[running.arcI,arcI_bixelJApVec_vec,arcI_bixelJApVec_i,arcI_bixelJApVec_j,arcI_bixelJApVec_offset,arcI_w] ...
            %    = matRad_bixWeightAndGradNEW(static,variable,running.arcI,arcI_bixelJApVec_vec,arcI_bixelJApVec_i,arcI_bixelJApVec_j,arcI_bixelJApVec_offset,arcI_w);
            
            % put tempResults into results
            % first do last dose beam
            results.arcI.lastDose.w{cellInd}                                                                                                                = tempResults.lastDose.w;
            results.arcI.lastDose.bixelJApVec_i{cellInd}((1:tempResults.lastDose.bixelJApVec_offset)+results.arcI.lastDose.bixelJApVec_offset{cellInd})     = full2UniqueLocalVar(tempResults.lastDose.bixelJApVec_i);
            results.arcI.lastDose.bixelJApVec_j{cellInd}((1:tempResults.lastDose.bixelJApVec_offset)+results.arcI.lastDose.bixelJApVec_offset{cellInd})     = tempResults.lastDose.bixelJApVec_j;
            results.arcI.lastDose.bixelJApVec_vec{cellInd}((1:tempResults.lastDose.bixelJApVec_offset)+results.arcI.lastDose.bixelJApVec_offset{cellInd})   = tempResults.lastDose.bixelJApVec_vec;
            results.arcI.lastDose.bixelJApVec_offset{cellInd}                                                                                               = results.arcI.lastDose.bixelJApVec_offset{cellInd}+tempResults.lastDose.bixelJApVec_offset;
            % now do next dose beam
            results.arcI.nextDose.w{cellInd}                                                                                                                = tempResults.nextDose.w;
            results.arcI.nextDose.bixelJApVec_i{cellInd}((1:tempResults.nextDose.bixelJApVec_offset)+results.arcI.nextDose.bixelJApVec_offset{cellInd})     = full2UniqueLocalVar(tempResults.nextDose.bixelJApVec_i);
            results.arcI.nextDose.bixelJApVec_j{cellInd}((1:tempResults.nextDose.bixelJApVec_offset)+results.arcI.nextDose.bixelJApVec_offset{cellInd})     = tempResults.nextDose.bixelJApVec_j;
            results.arcI.nextDose.bixelJApVec_vec{cellInd}((1:tempResults.nextDose.bixelJApVec_offset)+results.arcI.nextDose.bixelJApVec_offset{cellInd})   = tempResults.nextDose.bixelJApVec_vec;
            results.arcI.nextDose.bixelJApVec_offset{cellInd}                                                                                               = results.arcI.nextDose.bixelJApVec_offset{cellInd}+tempResults.nextDose.bixelJApVec_offset;
            
            % save shapeMap
            apertureInfo.beam(i).shape{phase_I}(j).shapeMap = apertureInfo.beam(i).shape{phase_I}(j).shapeMap+tempResults.shapeMap;
            
            % jacobian stuff
            if ~apertureInfo.propVMAT.beam(i).DAOBeam
                % if this is not an optimized beam, we need to add the
                % sumGradSq_weight to the previous and next optimized beams
                
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).sumGradSq_weight    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).sumGradSq_weight+tempResults.sumGradSq_weight_lastDAO;
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).sumGradSq_leaf      = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_lastDAO;
                
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).sumGradSq_weight    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).sumGradSq_weight+tempResults.sumGradSq_weight_nextDAO;
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).sumGradSq_leaf      = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_nextDAO;
            end
            
            % add the sumGradSq_leaf to the previous and next optimized
            % beams
            apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).sumGradSq_leaf = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_lastDAO;
            apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).sumGradSq_leaf = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_nextDAO;
            
            %% now do phase_F
            
            variable.arcF = true;
            
            bixelJApVecLastDose_sz              = effNumBix_lastDose_arcF.*(numVarMult+variable.totalNumOfShapes);
            bixelJApVecNextDose_sz              = effNumBix_nextDose_arcF.*(numVarMult+variable.totalNumOfShapes);
            variable.bixelJApVecLastDose_zeros  = zeros(bixelJApVecLastDose_sz,1);
            variable.bixelJApVecNextDose_zeros  = zeros(bixelJApVecNextDose_sz,1);
            
            variable.leftLeafPosI   = leftLeafPosI_arcF;
            variable.leftLeafPosF   = leftLeafPosF_arcF;
            variable.rightLeafPosI  = rightLeafPosI_arcF;
            variable.rightLeafPosF  = rightLeafPosF_arcF;
            
            variable.xPosLinearIndLeftLeafI     = xPosLinearIndLeftLeafI_arcF;
            variable.xPosLinearIndLeftLeafF     = xPosLinearIndLeftLeafF_arcF;
            variable.xPosLinearIndRightLeafI    = xPosLinearIndRightLeafI_arcF;
            variable.xPosLinearIndRightLeafF    = xPosLinearIndRightLeafF_arcF;
            
            variable.xPosIndLeftLeafI   = xPosIndLeftLeafI_arcF;
            variable.xPosIndLeftLeafF   = xPosIndLeftLeafF_arcF;
            variable.xPosIndRightLeafI  = xPosIndRightLeafI_arcF;
            variable.xPosIndRightLeafF  = xPosIndRightLeafF_arcF;
            
            variable.phase              = phase_F;
            variable.weight             = weight_F;
            variable.weightFactor       = static.weightFactor_F;
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                variable.DAOindex       = j+(phase_F-1)*numOfShapes;
                
                variable.jacobiScale    = jacobiScale_F;
            else
                
                variable.DAOindex_last      = DAOindex_last+(phase_F-1)*2;
                variable.DAOindex_next      = DAOindex_next+(phase_F-1)*2;
                
                variable.weight_last        = weight_last_F;
                variable.weight_next        = weight_next_F;
                variable.jacobiScale_last   = jacobiScale_last_F;
                variable.jacobiScale_next   = jacobiScale_next_F;
            end
            
            % calculate bixel weight and gradient
            tempResults = matRad_bixWeightAndGrad(static,variable);
            
            % put tempResults into results
            % first do last dose beam
            results.arcF.lastDose.w{cellInd}                                                                                                                = tempResults.lastDose.w;
            results.arcF.lastDose.bixelJApVec_i{cellInd}((1:tempResults.lastDose.bixelJApVec_offset)+results.arcF.lastDose.bixelJApVec_offset{cellInd})     = full2UniqueLocalVar(tempResults.lastDose.bixelJApVec_i);
            results.arcF.lastDose.bixelJApVec_j{cellInd}((1:tempResults.lastDose.bixelJApVec_offset)+results.arcF.lastDose.bixelJApVec_offset{cellInd})     = tempResults.lastDose.bixelJApVec_j;
            results.arcF.lastDose.bixelJApVec_vec{cellInd}((1:tempResults.lastDose.bixelJApVec_offset)+results.arcF.lastDose.bixelJApVec_offset{cellInd})   = tempResults.lastDose.bixelJApVec_vec;
            results.arcF.lastDose.bixelJApVec_offset{cellInd}                                                                                               = results.arcF.lastDose.bixelJApVec_offset{cellInd}+tempResults.lastDose.bixelJApVec_offset;
            % now do next dose beam
            results.arcF.nextDose.w{cellInd}                                                                                                                = tempResults.nextDose.w;
            results.arcF.nextDose.bixelJApVec_i{cellInd}((1:tempResults.nextDose.bixelJApVec_offset)+results.arcF.nextDose.bixelJApVec_offset{cellInd})     = full2UniqueLocalVar(tempResults.nextDose.bixelJApVec_i);
            results.arcF.nextDose.bixelJApVec_j{cellInd}((1:tempResults.nextDose.bixelJApVec_offset)+results.arcF.nextDose.bixelJApVec_offset{cellInd})     = tempResults.nextDose.bixelJApVec_j;
            results.arcF.nextDose.bixelJApVec_vec{cellInd}((1:tempResults.nextDose.bixelJApVec_offset)+results.arcF.nextDose.bixelJApVec_offset{cellInd})   = tempResults.nextDose.bixelJApVec_vec;
            results.arcF.nextDose.bixelJApVec_offset{cellInd}                                                                                               = results.arcF.nextDose.bixelJApVec_offset{cellInd}+tempResults.nextDose.bixelJApVec_offset;
            
            % save shapeMap
            apertureInfo.beam(i).shape{phase_F}(j).shapeMap     = apertureInfo.beam(i).shape{phase_F}(j).shapeMap+tempResults.shapeMap;
            
            % jacobian stuff
            if ~apertureInfo.propVMAT.beam(i).DAOBeam
                % if this is not an optimized beam, we need to add the
                % sumGradSq_weight to the previous and next optimized beams
                
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).sumGradSq_weight    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).sumGradSq_weight+tempResults.sumGradSq_weight_lastDAO;
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).sumGradSq_leaf      = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_lastDAO;
                
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).sumGradSq_weight    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).sumGradSq_weight+tempResults.sumGradSq_weight_nextDAO;
                apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).sumGradSq_leaf      = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_nextDAO;
            end
            
            % add the sumGradSq_leaf to the previous and next optimized
            % beams
            apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).sumGradSq_leaf = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_lastDAO;
            apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).sumGradSq_leaf = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).sumGradSq_leaf+tempResults.sumGradSq_leaf_nextDAO;
        end
        
        % cut out zeros, i.e. all elements after the offset
        results.arcI.lastDose.bixelJApVec_i{cellInd}((results.arcI.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcI.lastDose.bixelJApVec_j{cellInd}((results.arcI.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcI.lastDose.bixelJApVec_vec{cellInd}((results.arcI.lastDose.bixelJApVec_offset{cellInd}+1):end) = [];
        results.arcI.nextDose.bixelJApVec_i{cellInd}((results.arcI.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcI.nextDose.bixelJApVec_j{cellInd}((results.arcI.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcI.nextDose.bixelJApVec_vec{cellInd}((results.arcI.nextDose.bixelJApVec_offset{cellInd}+1):end) = [];
        results.arcF.lastDose.bixelJApVec_i{cellInd}((results.arcF.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcF.lastDose.bixelJApVec_j{cellInd}((results.arcF.lastDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcF.lastDose.bixelJApVec_vec{cellInd}((results.arcF.lastDose.bixelJApVec_offset{cellInd}+1):end) = [];
        results.arcF.nextDose.bixelJApVec_i{cellInd}((results.arcF.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcF.nextDose.bixelJApVec_j{cellInd}((results.arcF.nextDose.bixelJApVec_offset{cellInd}+1):end)   = [];
        results.arcF.nextDose.bixelJApVec_vec{cellInd}((results.arcF.nextDose.bixelJApVec_offset{cellInd}+1):end) = [];
    end
end

end