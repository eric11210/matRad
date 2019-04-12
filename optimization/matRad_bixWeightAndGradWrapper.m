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

bixelIndMap         = apertureInfo.beam(i).bixelIndMap-max(apertureInfo.beam(i).bixelIndMap(:))+apertureInfo.beam(i).numBixels;
bixelWidth          = apertureInfo.bixelWidth;
lim_l               = apertureInfo.beam(i).lim_l;
lim_r               = apertureInfo.beam(i).lim_r;
edges_l             = apertureInfo.beam(i).posOfCornerBixel(1)...
    + ((1:size(bixelIndMap,2))-1-1/2)*bixelWidth;
edges_r             = apertureInfo.beam(i).posOfCornerBixel(1)...
    + ((1:size(bixelIndMap,2))-1+1/2)*bixelWidth;
centres             = (edges_l+edges_r)/2;
widths              = edges_r-edges_l;
numRow              = apertureInfo.beam(i).numOfActiveLeafPairs;
numCol              = size(bixelIndMap,2);
numBix              = apertureInfo.beam(i).numBixels;
totalNumOfShapes    = apertureInfo.totalNumOfShapes;

% indices for various variables
tIx_Vec     = (1:totalNumOfShapes)+apertureInfo.beam(i).numUniqueVar-totalNumOfShapes;%(totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:totalNumOfShapes);

if ~apertureInfo.propVMAT.beam(i).DAOBeam
    time_last   = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time;
    time_next   = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time;
    time        = apertureInfo.beam(i).time;
    
    doseAngleBordersDiff        = apertureInfo.propVMAT.beam(i).doseAngleBordersDiff;
    doseAngleBordersDiff_last   = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff;
    doseAngleBordersDiff_next   = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff;
    timeFacCurr_last            = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr;
    timeFacCurr_next            = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr;
    fracFromLastDAO             = apertureInfo.propVMAT.beam(i).fracFromLastDAO;
    timeFracFromLastDAO         = apertureInfo.propVMAT.beam(i).timeFracFromLastDAO;
    timeFracFromNextDAO         = apertureInfo.propVMAT.beam(i).timeFracFromNextDAO;
    
    DAOindex_last   = 1;%apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
    DAOindex_next   = 2;%apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
    tIx_last        = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex+apertureInfo.beam(i).numUniqueVar-totalNumOfShapes;%(totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
    tIx_next        = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex+apertureInfo.beam(i).numUniqueVar-totalNumOfShapes;%(totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
    
    if apertureInfo.propVMAT.continuousAperture
        fracFromLastOpt = apertureInfo.propVMAT.beam(i).fracFromLastDAO;
        fracFromLastOptI = apertureInfo.propVMAT.beam(i).fracFromLastDAO_I*ones(numRow,1);
        fracFromLastOptF = apertureInfo.propVMAT.beam(i).fracFromLastDAO_F*ones(numRow,1);
        fracFromNextOptI = apertureInfo.propVMAT.beam(i).fracFromNextDAO_I*ones(numRow,1);
        fracFromNextOptF = apertureInfo.propVMAT.beam(i).fracFromNextDAO_F*ones(numRow,1);
    else
        fracFromLastOpt = apertureInfo.propVMAT.beam(i).fracFromLastDAO;
        fracFromLastOptI = apertureInfo.propVMAT.beam(i).fracFromLastDAO*ones(numRow,1);
        fracFromLastOptF = apertureInfo.propVMAT.beam(i).fracFromLastDAO*ones(numRow,1);
        fracFromNextOptI = (1-apertureInfo.propVMAT.beam(i).fracFromLastDAO)*ones(numRow,1);
        fracFromNextOptF = (1-apertureInfo.propVMAT.beam(i).fracFromLastDAO)*ones(numRow,1);
    end
end

%% prep static

static.lim_r            = lim_r;
static.edges_l          = edges_l;
static.edges_r          = edges_r;
static.centres          = centres;
static.widths           = widths;
static.bixelIndMap      = cast(bixelIndMap,'like',results.arcI.bixelJApVec_j{1});
static.numRow           = numRow;
static.numCol           = numCol;
static.numBix           = numBix;
static.DAOBeam          = apertureInfo.propVMAT.beam(i).DAOBeam;
static.totalNumOfBixels = apertureInfo.totalNumOfBixels;
static.bixIndVec        = 1:numBix;

if apertureInfo.propVMAT.beam(i).DAOBeam
    
    numVarMult              = 7;
else
    
    static.fracFromLastOpt  = fracFromLastOpt;
    static.fracFromLastOptI = fracFromLastOptI;
    static.fracFromLastOptF = fracFromLastOptF;
    static.fracFromNextOptI = fracFromNextOptI;
    static.fracFromNextOptF = fracFromNextOptF;
    
    static.time         = time;
    static.time_next    = time_next;
    static.time_last    = time_last;
    
    static.tIx_last         = tIx_last;
    static.tIx_next         = tIx_next;
    
    static.doseAngleBordersDiff  = doseAngleBordersDiff;
    static.doseAngleBordersDiff_last = doseAngleBordersDiff_last;
    static.doseAngleBordersDiff_next = doseAngleBordersDiff_next;
    static.timeFacCurr_last = timeFacCurr_last;
    static.timeFacCurr_next = timeFacCurr_next;
    static.fracFromLastDAO  = fracFromLastDAO;
    static.timeFracFromLastDAO = timeFracFromLastDAO;
    static.timeFracFromNextDAO = timeFracFromNextDAO;
    
    numVarMult              = 16;
end

% determine probabilities
[Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(apertureInfo.propVMAT.qij,apertureInfo.propVMAT.initProb,apertureInfo.beam(i).time,sum([apertureInfo.beam(1:(i-1)).time]));

for j = 1:numOfShapes
    
    for phase_I = 1:apertureInfo.numPhases
        
        transitions = apertureInfo.propVMAT.beam(i).transMask(phase_I,:);
        transitions(transitions == 0) = [];
        
        % loop over all (final) phases
        for phase_F = 1:apertureInfo.numPhases
            
            %% determine variable quantities
            
            weight_I        = apertureInfo.beam(i).shape{phase_I}(j).weight;
            weight_F        = apertureInfo.beam(i).shape{phase_F}(j).weight;
            weightFactor_I  = apertureInfo.propVMAT.beam(i).doseAngleBorderCentreDiff(1)./apertureInfo.propVMAT.beam(i).doseAngleBordersDiff;
            weightFactor_F  = apertureInfo.propVMAT.beam(i).doseAngleBorderCentreDiff(2)./apertureInfo.propVMAT.beam(i).doseAngleBordersDiff;
            leftLeafPos_I   = apertureInfo.beam(i).shape{phase_I}(j).leftLeafPos_I;
            leftLeafPos_F   = apertureInfo.beam(i).shape{phase_F}(j).leftLeafPos_F;
            rightLeafPos_I  = apertureInfo.beam(i).shape{phase_I}(j).rightLeafPos_I;
            rightLeafPos_F  = apertureInfo.beam(i).shape{phase_F}(j).rightLeafPos_F;
            
            %% sort out order, set up calculation of bixel weight and gradients
            
            % set the initial leaf positions to the minimum leaf positions
            % always, instead of the leaf positions at the actual beginning
            % of the arc
            % this simplifies the calculation
            % remember which one is actually I and F in leftMinInd
            [leftLeafPosI,leftMinInd] = min([leftLeafPos_I,leftLeafPos_F],[],2);
            leftLeafPosF = max([leftLeafPos_I,leftLeafPos_F],[],2);
            [rightLeafPosI,rightMinInd] = min([rightLeafPos_I,rightLeafPos_F],[],2);
            rightLeafPosF = max([rightLeafPos_I,rightLeafPos_F],[],2);
            
            leftLeafPosI = round2(leftLeafPosI,10);
            leftLeafPosF = round2(leftLeafPosF,10);
            rightLeafPosI = round2(rightLeafPosI,10);
            rightLeafPosF = round2(rightLeafPosF,10);
            
            leftLeafPosI(leftLeafPosI <= lim_l) = lim_l(leftLeafPosI <= lim_l);
            leftLeafPosF(leftLeafPosF <= lim_l) = lim_l(leftLeafPosF <= lim_l);
            rightLeafPosI(rightLeafPosI <= lim_l) = lim_l(rightLeafPosI <= lim_l);
            rightLeafPosF(rightLeafPosF <= lim_l) = lim_l(rightLeafPosF <= lim_l);
            leftLeafPosI(leftLeafPosI >= lim_r) = lim_r(leftLeafPosI >= lim_r);
            leftLeafPosF(leftLeafPosF >= lim_r) = lim_r(leftLeafPosF >= lim_r);
            rightLeafPosI(rightLeafPosI >= lim_r) = lim_r(rightLeafPosI >= lim_r);
            rightLeafPosF(rightLeafPosF >= lim_r) = lim_r(rightLeafPosF >= lim_r);
            
            % determine middle leaf positions
            leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
            rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
            
            % find bixel indices where leaves are located
            xPosIndLeftLeafI = min(floor((leftLeafPosI-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndLeftLeafM = min(floor((leftLeafPosM-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndLeftLeafF = max(ceil((leftLeafPosF-edges_r(1))./bixelWidth)+1,1);
            xPosIndRightLeafI = min(floor((rightLeafPosI-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndRightLeafM = min(floor((rightLeafPosM-edges_l(1))./bixelWidth)+1,numCol);
            xPosIndRightLeafF = max(ceil((rightLeafPosF-edges_r(1))./bixelWidth)+1,1);
            %
            xPosLinearIndLeftLeafI = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafI);
            xPosLinearIndLeftLeafM = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafM);
            xPosLinearIndLeftLeafF = sub2ind([numRow numCol],(1:numRow)',xPosIndLeftLeafF);
            xPosLinearIndRightLeafI = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafI);
            xPosLinearIndRightLeafM = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafM);
            xPosLinearIndRightLeafF = sub2ind([numRow numCol],(1:numRow)',xPosIndRightLeafF);
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                jacobiScale_I = apertureInfo.beam(i).shape{phase_I}(1).jacobiScale;
                jacobiScale_F = apertureInfo.beam(i).shape{phase_F}(1).jacobiScale;
                
                if apertureInfo.propVMAT.continuousAperture
                    vectorIx_LI   = numOfShapes.*apertureInfo.numPhases+2.*numRow.*(phase_I-1)+(1:numRow);%apertureInfo.beam(i).shape{phase_I}(j).vectorOffset(1) + ((1:n)-1);
                    vectorIx_LF   = numOfShapes.*apertureInfo.numPhases+2.*numRow.*(phase_F-1)+numRow+(1:numRow);%apertureInfo.beam(i).shape{phase_F}(j).vectorOffset(2) + ((1:n)-1);
                    vectorIx_RI   = vectorIx_LI+2.*numRow.*apertureInfo.numPhases;%vectorIx_LI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_RF   = vectorIx_LF+2.*numRow.*apertureInfo.numPhases;%vectorIx_LF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                else
                    vectorIx_LI   = apertureInfo.beam(i).shape{phase_I}(j).vectorOffset + ((1:numRow)-1);
                    vectorIx_LF   = apertureInfo.beam(i).shape{phase_F}(j).vectorOffset + ((1:numRow)-1);
                    vectorIx_RI   = vectorIx_LI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_RF   = vectorIx_LF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                end
                
                % change the vectorIx_xy elements to remember which
                % apertureVector elements the "new" I and F
                % if leftMinInd is 2, the I and F are switched
                tempL = vectorIx_LI;
                tempR = vectorIx_RI;
                vectorIx_LI(leftMinInd == 2)    = vectorIx_LF(leftMinInd == 2);
                vectorIx_LF(leftMinInd == 2)    = tempL(leftMinInd == 2);
                vectorIx_RI(rightMinInd == 2)   = vectorIx_RF(rightMinInd == 2);
                vectorIx_RF(rightMinInd == 2)   = tempR(rightMinInd == 2);
            else
                
                weight_last_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).weight;
                weight_last_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(j).weight;
                weight_next_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(j).weight;
                weight_next_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).weight;
                
                jacobiScale_last_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(1).jacobiScale;
                jacobiScale_last_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_F}(1).jacobiScale;
                jacobiScale_next_I = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_I}(1).jacobiScale;
                jacobiScale_next_F = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(1).jacobiScale;
                
                if apertureInfo.propVMAT.continuousAperture
                    vectorIx_LF_last  = 2.*apertureInfo.numPhases+2.*numRow.*(phase_I-1)+(1:numRow);%apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).vectorOffset(2) + ((1:numRow)-1);
                    vectorIx_LI_next  = 2.*apertureInfo.numPhases+2.*numRow.*(phase_F-1)+numRow+(1:numRow);%apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).vectorOffset(1) + ((1:numRow)-1);
                    vectorIx_RF_last  = vectorIx_LF_last+2.*numRow.*apertureInfo.numPhases;%vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_RI_next  = vectorIx_LI_next+2.*numRow.*apertureInfo.numPhases;%vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                else
                    vectorIx_LF_last  = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}(j).vectorOffset + ((1:numRow)-1);
                    vectorIx_LI_next  = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}(j).vectorOffset + ((1:numRow)-1);
                    vectorIx_RF_last  = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_RI_next  = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                end
                
                tempL = vectorIx_LF_last;
                tempR = vectorIx_RF_last;
                
                vectorIx_LF_last(leftMinInd == 2) = vectorIx_LI_next(leftMinInd == 2);
                vectorIx_LI_next(leftMinInd == 2) = tempL(leftMinInd == 2);
                
                vectorIx_RF_last(rightMinInd == 2) = vectorIx_RI_next(rightMinInd == 2);
                vectorIx_RI_next(rightMinInd == 2) = tempR(rightMinInd == 2);
            end
            
            % determine probabilities and derivatives
            probability         = Pi_T(phase_I).*Pij_transT(phase_I,phase_F);
            probability_dTVec   = sum(apertureInfo.propVMAT.jacobT(:,1:(i-1)),2).*Pi_T_dot(phase_I).*Pij_transT(phase_I,phase_F) + apertureInfo.propVMAT.jacobT(:,i).*Pi_T(phase_I).*Pij_transT_dot(phase_I,phase_F);
            
            % delete any derivatives with value less than eps
            delInd = abs(probability_dTVec) < eps;
            probability_dTVec(delInd)   = [];
            variable.tIx_Vec            = cast(tIx_Vec,'like',results.arcI.bixelJApVec_i{1});
            variable.tIx_Vec(delInd)    = [];
            variable.totalNumOfShapes   = numel(probability_dTVec);
            variable.bixelJApVec_sz     = numBix.*(numVarMult+variable.totalNumOfShapes);
            variable.numShapbixIndVec     = 1:(variable.totalNumOfShapes*numBix);
            
            %% first do phase_I
            
            variable.arcF = false;
            
            variable.leftLeafPosI   = leftLeafPosI;
            variable.leftLeafPosF   = leftLeafPosM;
            variable.rightLeafPosI  = rightLeafPosI;
            variable.rightLeafPosF  = rightLeafPosM;
            
            variable.xPosLinearIndLeftLeafI     = xPosLinearIndLeftLeafI;
            variable.xPosLinearIndLeftLeafF     = xPosLinearIndLeftLeafM;
            variable.xPosLinearIndRightLeafI    = xPosLinearIndRightLeafI;
            variable.xPosLinearIndRightLeafF    = xPosLinearIndRightLeafM;
            
            variable.xPosIndLeftLeafI   = xPosIndLeftLeafI;
            variable.xPosIndLeftLeafF   = xPosIndLeftLeafM;
            variable.xPosIndRightLeafI  = xPosIndRightLeafI;
            variable.xPosIndRightLeafF  = xPosIndRightLeafM;
            
            variable.phase              = phase_I;
            variable.probability        = probability;
            variable.probability_dTVec  = probability_dTVec;
            variable.weight             = weight_I;
            variable.weightFactor_I     = weightFactor_I;
            variable.weightFactor_F     = weightFactor_F;
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                variable.DAOindex       = j+(phase_I-1)*numOfShapes;
                
                variable.jacobiScale    = jacobiScale_I;
                variable.vectorIx_LI    = vectorIx_LI;
                variable.vectorIx_LF    = vectorIx_LF;
                variable.vectorIx_RI    = vectorIx_RI;
                variable.vectorIx_RF    = vectorIx_RF;
            else
                
                variable.DAOindex_last      = DAOindex_last+(phase_I-1)*2;
                variable.DAOindex_next      = DAOindex_next+(phase_I-1)*2;
                
                variable.vectorIx_LF_last   = vectorIx_LF_last;
                variable.vectorIx_LI_next   = vectorIx_LI_next;
                variable.vectorIx_RF_last   = vectorIx_RF_last;
                variable.vectorIx_RI_next   = vectorIx_RI_next;
                
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
            results.arcI.w{phase_I}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap)))         = results.arcI.w{phase_I}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap)))+tempResults.w;
            results.arcI.bixelJApVec_i{phase_I}((1:variable.bixelJApVec_sz)+results.arcI.bixelJApVec_offset{phase_I})   = tempResults.bixelJApVec_i;
            results.arcI.bixelJApVec_j{phase_I}((1:variable.bixelJApVec_sz)+results.arcI.bixelJApVec_offset{phase_I})   = tempResults.bixelJApVec_j;
            results.arcI.bixelJApVec_vec{phase_I}((1:variable.bixelJApVec_sz)+results.arcI.bixelJApVec_offset{phase_I}) = tempResults.bixelJApVec_vec;
            results.arcI.bixelJApVec_offset{phase_I}                                                                    = results.arcI.bixelJApVec_offset{phase_I}+variable.bixelJApVec_sz;
            
            apertureInfo.beam(i).shape{phase_I}(j).shapeMap     = apertureInfo.beam(i).shape{phase_I}(j).shapeMap+tempResults.shapeMap;
            apertureInfo.beam(i).shape{phase_I}(j).sumGradSq    = apertureInfo.beam(i).shape{phase_I}(j).sumGradSq+tempResults.sumGradSq;
            %% now do phase_F
            
            variable.arcF = true;
            
            variable.leftLeafPosI   = leftLeafPosM;
            variable.leftLeafPosF   = leftLeafPosF;
            variable.rightLeafPosI  = rightLeafPosM;
            variable.rightLeafPosF  = rightLeafPosF;
            
            variable.xPosLinearIndLeftLeafI     = xPosLinearIndLeftLeafM;
            variable.xPosLinearIndLeftLeafF     = xPosLinearIndLeftLeafF;
            variable.xPosLinearIndRightLeafI    = xPosLinearIndRightLeafM;
            variable.xPosLinearIndRightLeafF    = xPosLinearIndRightLeafF;
            
            variable.xPosIndLeftLeafI   = xPosIndLeftLeafM;
            variable.xPosIndLeftLeafF   = xPosIndLeftLeafF;
            variable.xPosIndRightLeafI  = xPosIndRightLeafM;
            variable.xPosIndRightLeafF  = xPosIndRightLeafF;
            
            variable.phase              = phase_F;
            variable.probability        = probability;
            variable.probability_dTVec  = probability_dTVec;
            variable.weight             = weight_F;
            variable.weightFactor_I     = weightFactor_I;
            variable.weightFactor_F     = weightFactor_F;
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                variable.DAOindex       = j+(phase_F-1)*numOfShapes;
                
                variable.jacobiScale    = jacobiScale_F;
                variable.vectorIx_LI    = vectorIx_LI;
                variable.vectorIx_LF    = vectorIx_LF;
                variable.vectorIx_RI    = vectorIx_RI;
                variable.vectorIx_RF    = vectorIx_RF;
            else
                
                variable.DAOindex_last      = DAOindex_last+(phase_F-1)*2;
                variable.DAOindex_next      = DAOindex_next+(phase_F-1)*2;
                
                variable.vectorIx_LF_last   = vectorIx_LF_last;
                variable.vectorIx_LI_next   = vectorIx_LI_next;
                variable.vectorIx_RF_last   = vectorIx_RF_last;
                variable.vectorIx_RI_next   = vectorIx_RI_next;
                
                variable.weight_last        = weight_last_F;
                variable.weight_next        = weight_next_F;
                variable.jacobiScale_last   = jacobiScale_last_F;
                variable.jacobiScale_next   = jacobiScale_next_F;
            end
            
            % calculate bixel weight and gradient
            tempResults = matRad_bixWeightAndGrad(static,variable);
            
            % put tempResults into results
            results.arcF.w{phase_F}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap)))         = results.arcF.w{phase_F}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap)))+tempResults.w;
            results.arcF.bixelJApVec_i{phase_F}((1:variable.bixelJApVec_sz)+results.arcF.bixelJApVec_offset{phase_F})   = tempResults.bixelJApVec_i;
            results.arcF.bixelJApVec_j{phase_F}((1:variable.bixelJApVec_sz)+results.arcF.bixelJApVec_offset{phase_F})   = tempResults.bixelJApVec_j;
            results.arcF.bixelJApVec_vec{phase_F}((1:variable.bixelJApVec_sz)+results.arcF.bixelJApVec_offset{phase_F}) = tempResults.bixelJApVec_vec;
            results.arcF.bixelJApVec_offset{phase_F}                                                                    = results.arcF.bixelJApVec_offset{phase_F}+variable.bixelJApVec_sz;
            
            apertureInfo.beam(i).shape{phase_F}(j).shapeMap     = apertureInfo.beam(i).shape{phase_F}(j).shapeMap+tempResults.shapeMap;
            apertureInfo.beam(i).shape{phase_F}(j).sumGradSq    = tempResults.sumGradSq;
            
        end
    end
end


end