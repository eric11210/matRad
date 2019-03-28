function [results,running] = matRad_bixWeightAndGradOLD(results,running,parameters)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the bixel weights from the aperture vector,
% and also the Jacobian matrix relating these two.
%
% call
%   [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sqrtSumGradSq,shapeMap,counters] = ...
%    matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j)
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

%% extract results and running

w               = results.w;
bixelJApVec_vec = results.bixelJApVec_vec;
bixelJApVec_i   = results.bixelJApVec_i;
bixelJApVec_j   = results.bixelJApVec_j;

bixelJApVec_offset  = running.bixelJApVec_offset;
sumGradSq           = running.sumGradSq;
shapeMap            = running.shapeMap;

%% extract variables from inputs

edges_l                 = parameters.edges_l;
edges_r                 = parameters.edges_r;
centres                 = parameters.centres;
numBix                  = parameters.numBix;

leftLeafPosI            = parameters.leftLeafPosI;
leftLeafPosF            = parameters.leftLeafPosF;
xPosLinearIndLeftLeafI  = parameters.xPosLinearIndLeftLeafI;
xPosLinearIndLeftLeafF  = parameters.xPosLinearIndLeftLeafF;

%% bixel weight calculation

%calculate fraction of fluence uncovered by left leaf
%initial computation
uncoveredByLeftLeaf = bsxfun(@minus,centres,leftLeafPosI)./repmat(leftLeafPosF-leftLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) + (leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafI)').*2);
uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) - (edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafF)').*2);
%round <0 to 0, >1 to 1
uncoveredByLeftLeaf(uncoveredByLeftLeaf < 0) = 0;
uncoveredByLeftLeaf(uncoveredByLeftLeaf > 1) = 1;

%calculate fraction of fluence covered by right leaf
%initial computation
coveredByRightLeaf = bsxfun(@minus,centres,rightLeafPosI)./repmat(rightLeafPosM-rightLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + (rightLeafPosI-edges_l(xPosIndRightLeafI)').^2./((rightLeafPosM-rightLeafPosI).*(widths(xPosIndRightLeafI)').*2);
coveredByRightLeaf(xPosLinearIndRightLeafM) = coveredByRightLeaf(xPosLinearIndRightLeafM) - (edges_r(xPosIndRightLeafM)'-rightLeafPosM).^2./((rightLeafPosM-rightLeafPosI).*(widths(xPosIndRightLeafM)').*2);
%round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = bsxfun(@minus,centres,leftLeafPosM)./(repmat(leftLeafPosM-leftLeafPosI,1,numBix)).^2;
dUl_dLM = bsxfun(@minus,leftLeafPosI,centres)./(repmat(leftLeafPosM-leftLeafPosI,1,numBix)).^2;

dCr_dRI = bsxfun(@minus,centres,rightLeafPosM)./(repmat(rightLeafPosM-rightLeafPosI,1,numBix)).^2;
dCr_dRM = bsxfun(@minus,rightLeafPosI,centres)./(repmat(rightLeafPosM-rightLeafPosI,1,numBix)).^2;

dUl_dLI(xPosLinearIndLeftLeafI) = dUl_dLI(xPosLinearIndLeftLeafI) + ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').*(2*leftLeafPosM-leftLeafPosI-edges_l(xPosIndLeftLeafI)'))./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLM(xPosLinearIndLeftLeafI) = dUl_dLM(xPosLinearIndLeftLeafI) - ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2)./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLI(xPosLinearIndLeftLeafM) = dUl_dLI(xPosLinearIndLeftLeafM) - ((edges_r(xPosIndLeftLeafM)'-leftLeafPosM).^2)./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafM)').*2);
dUl_dLM(xPosLinearIndLeftLeafM) = dUl_dLM(xPosLinearIndLeftLeafM) + ((edges_r(xPosIndLeftLeafM)'-leftLeafPosM).*(leftLeafPosM+edges_r(xPosIndLeftLeafM)'-2*leftLeafPosI))./((leftLeafPosM-leftLeafPosI).^2.*(widths(xPosIndLeftLeafM)').*2);

dCr_dRI(xPosLinearIndRightLeafI) = dCr_dRI(xPosLinearIndRightLeafI) + ((rightLeafPosI-edges_l(xPosIndRightLeafI)').*(2*rightLeafPosM-rightLeafPosI-edges_l(xPosIndRightLeafI)'))./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRM(xPosLinearIndRightLeafI) = dCr_dRM(xPosLinearIndRightLeafI) - ((rightLeafPosI-edges_l(xPosIndRightLeafI)').^2)./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRI(xPosLinearIndRightLeafM) = dCr_dRI(xPosLinearIndRightLeafM) - ((edges_r(xPosIndRightLeafM)'-rightLeafPosM).^2)./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafM)').*2);
dCr_dRM(xPosLinearIndRightLeafM) = dCr_dRM(xPosLinearIndRightLeafM) + ((edges_r(xPosIndRightLeafM)'-rightLeafPosM).*(rightLeafPosM+edges_r(xPosIndRightLeafM)'-2*rightLeafPosI))./((rightLeafPosM-rightLeafPosI).^2.*(widths(xPosIndRightLeafM)').*2);

for k = 1:n
    dUl_dLI(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLM(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLI(k,(xPosIndLeftLeafM(k)+1):numBix) = 0;
    dUl_dLM(k,(xPosIndLeftLeafM(k)+1):numBix) = 0;
    
    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafM(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary
        
        %19 July 2017 in journal
        dUl_dLI(k,xPosIndLeftLeafI(k)) = -1/(2*widths(xPosIndLeftLeafI(k))');
        dUl_dLM(k,xPosIndLeftLeafM(k)) = -1/(2*widths(xPosIndLeftLeafM(k))');
        if leftLeafPosM(k)-leftLeafPosI(k) <= eps(max(lim_r))
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafM(k)) = (edges_r(xPosIndLeftLeafM(k))-leftLeafPosM(k))./widths(xPosIndLeftLeafM(k));
        end
    end
    
    dCr_dRI(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRM(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRI(k,(xPosIndRightLeafM(k)+1):numBix) = 0;
    dCr_dRM(k,(xPosIndRightLeafM(k)+1):numBix) = 0;
    
    if xPosIndRightLeafI(k) >= xPosIndRightLeafM(k)
        dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
        dCr_dRM(k,xPosIndRightLeafM(k)) = -1/(2*widths(xPosIndRightLeafM(k))');
        if rightLeafPosM(k)-rightLeafPosI(k) <= eps(max(lim_r))
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafM(k)) = (edges_r(xPosIndRightLeafM(k))-rightLeafPosM(k))./widths(xPosIndRightLeafM(k));
        end
    end
end

% store information for Jacobi preconditioning
sumGradSq{phase_I} = sumGradSq{phase_I}+(weightFactor_I.*probability).^2.*mean([sum((dUl_dLI).^2,2); sum((dUl_dLM.*weightFactor_F).^2,2); sum((dUl_dLM.*weightFactor_I).^2,2); sum((dCr_dRI).^2,2); sum((dCr_dRM.*weightFactor_F).^2,2); sum((dCr_dRM.*weightFactor_I).^2,2)]);

%% save the bixel weights
%fluence is equal to fluence not covered by left leaf minus
%fluence covered by left leaf
shapeMap_nW = uncoveredByLeftLeaf-coveredByRightLeaf;
shapeMap_nW = round2(shapeMap_nW,15);
shapeMap_nW(isnan(shapeMap_nW)) = 0;

% find open bixels
%shapeMapIx = shapeMap > 0;
shapeMapIx = ~isnan(bixelIndMap);

currBixelIx = bixelIndMap(shapeMapIx);
w{phase_I}(currBixelIx) = w{phase_I}(currBixelIx) + shapeMap_nW(shapeMapIx).*weight_I.*weightFactor_I.*probability;
shapeMap{phase_I} = shapeMap{phase_I}+shapeMap_nW.*weight_I.*weightFactor_I.*probability;

%% save the gradients

bixelJApVec_offset_temp = bixelJApVec_offset{phase_I};
numSaveBixel = nnz(shapeMapIx);

if apertureInfo.propVMAT.beam(i).DAOBeam
    % indices
    vectorIxMat_LI = repmat(vectorIx_LI',1,numBix);
    vectorIxMat_LF = repmat(vectorIx_LF',1,numBix);
    vectorIxMat_RI = repmat(vectorIx_RI',1,numBix);
    vectorIxMat_RF = repmat(vectorIx_RF',1,numBix);
    
    % wrt weight
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_I;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex+(phase_I-1)*totalNumOfShapes;
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt initial left
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLI(shapeMapIx).*weight_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final left
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt initial right
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRI(shapeMapIx).*weight_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final right
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
else
    % indices
    vectorIxMat_LF_last = repmat(vectorIx_LF_last',1,numBix);
    vectorIxMat_LI_next = repmat(vectorIx_LI_next',1,numBix);
    vectorIxMat_RF_last = repmat(vectorIx_RF_last',1,numBix);
    vectorIxMat_RI_next = repmat(vectorIx_RI_next',1,numBix);
    
    % leaf interpolation fractions/weights
    fracFromLastOptI = repmat(fracFromLastOptI,1,numBix);
    fracFromLastOptF = repmat(fracFromLastOptF,1,numBix);
    fracFromNextOptI = repmat(fracFromNextOptI,1,numBix);
    fracFromNextOptF = repmat(fracFromNextOptF,1,numBix);
    
    % wrt last weight
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOpt*(time./time_last)*shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_last_I;
    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
    %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
    %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex_last+(phase_I-1)*totalNumOfShapes;
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt next weight
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = (1-fracFromLastOpt)*(time./time_next)*shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_next_I;
    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
    %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
    %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex_next+(phase_I-1)*totalNumOfShapes;
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    
    % updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
    %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
    
    % updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
    % updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
    
    % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt initial left (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptF(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final left (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptF(shapeMapIx).*dUl_dLM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt initial right (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptF(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final right (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptF(shapeMapIx).*dCr_dRM(shapeMapIx).*weight_I.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt last time
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = weightFactor_I.*probability.*(doseAngleBordersDiff.*timeFacCurr_last) ...
        .*(-fracFromLastDAO.*timeFracFromNextDAO.*(weight_last_I./doseAngleBordersDiff_next).*(time_next./time_last.^2) ...
        +(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next_I./doseAngleBordersDiff_last).*(1./time_next)) ...
        * shapeMap_nW(shapeMapIx);
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = tIx_last;
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt next time
    bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = weightFactor_I.*probability.*(doseAngleBordersDiff.*timeFacCurr_next) ...
        .*(fracFromLastDAO.*timeFracFromNextDAO.*(weight_last_I./doseAngleBordersDiff_next).*(1./time_last) ...
        -(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next_I./doseAngleBordersDiff_last).*(time_last./time_next.^2)) ...
        * shapeMap_nW(shapeMapIx);
    bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = tIx_next;
    bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel)) = bixelIndMap(shapeMapIx);
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
end

% wrt times (probability)
bixelJApVec_vec{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel*totalNumOfShapes)) = weight_I.*weightFactor_I.*shapeMap_nW(shapeMapIx)*probability_dTVec';
bixelJApVec_i{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel*totalNumOfShapes)) = ones(numSaveBixel,1)*tIx_Vec;
bixelJApVec_j{phase_I}(bixelJApVec_offset_temp+(1:numSaveBixel*totalNumOfShapes)) = bixelIndMap(shapeMapIx)*ones(1,totalNumOfShapes);
bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel*totalNumOfShapes;

% update counters
bixelJApVec_offset{phase_I} = bixelJApVec_offset_temp;

%% update results and running

results.w               = w;
results.bixelJApVec_vec = bixelJApVec_vec;
results.bixelJApVec_i   = bixelJApVec_i;
results.bixelJApVec_j   = bixelJApVec_j;

running.bixelJApVec_offset  = bixelJApVec_offset;
running.sumGradSq           = sumGradSq;
running.shapeMap            = shapeMap;

end