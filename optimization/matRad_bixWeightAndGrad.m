function [results,running] = matRad_bixWeightAndGrad(static,variable,results,running)

round2 = @(a,b) round(a*10^b)/10^b;

%% extract results and running

w                   = results.w;
bixelJApVec_vec     = results.bixelJApVec_vec ;
bixelJApVec_i       = results.bixelJApVec_i;
bixelJApVec_j       = results.bixelJApVec_j;
bixelJApVec_offset  = results.bixelJApVec_offset;

sumGradSq   = running.sumGradSq;
shapeMap    = running.shapeMap;

%% extract variables from input

lim_r               = static.lim_r;
edges_l             = static.edges_l;
edges_r             = static.edges_r;
centres             = static.centres;
widths              = static.widths;
n                   = static.n;
bixelIndMap         = static.bixelIndMap;
numBix              = static.numBix;
totalNumOfShapes    = static.totalNumOfShapes;
tIx_Vec             = static.tIx_Vec;

leftLeafPosI    = variable.leftLeafPosI;
leftLeafPosF    = variable.leftLeafPosF;
rightLeafPosI   = variable.rightLeafPosI;
rightLeafPosF   = variable.rightLeafPosF;

xPosLinearIndLeftLeafI  = variable.xPosLinearIndLeftLeafI;
xPosLinearIndLeftLeafF  = variable.xPosLinearIndLeftLeafF;
xPosLinearIndRightLeafI = variable.xPosLinearIndRightLeafI;
xPosLinearIndRightLeafF = variable.xPosLinearIndRightLeafF;

xPosIndLeftLeafI    = variable.xPosIndLeftLeafI;
xPosIndLeftLeafF    = variable.xPosIndLeftLeafF;
xPosIndRightLeafI   = variable.xPosIndRightLeafI;
xPosIndRightLeafF   = variable.xPosIndRightLeafF;

phase               = variable.phase;
probability         = variable.probability;
probability_dTVec   = variable.probability_dTVec;
weight              = variable.weight;
weightFactor_I      = variable.weightFactor_I;
weightFactor_F      = variable.weightFactor_F;


if static.DAOBeam
    
    DAOindex    = static.DAOindex;
    
    jacobiScale = variable.jacobiScale;
    vectorIx_LI = variable.vectorIx_LI;
    vectorIx_LF = variable.vectorIx_LF;
    vectorIx_RI = variable.vectorIx_RI;
    vectorIx_RF = variable.vectorIx_RF;
    
else
    
    vectorIx_LF_last = variable.vectorIx_LF_last;
    vectorIx_LI_next = variable.vectorIx_LI_next;
    vectorIx_RF_last = variable.vectorIx_RF_last;
    vectorIx_RI_next = variable.vectorIx_RI_next;
    
    fracFromLastOpt     = static.fracFromLastOpt;
    fracFromLastOptI    = static.fracFromLastOptI;
    fracFromLastOptF    = static.fracFromLastOptF;
    fracFromNextOptI    = static.fracFromNextOptI;
    fracFromNextOptF    = static.fracFromNextOptF;
    
    weight_last         = variable.weight_last;
    weight_next         = variable.weight_next;
    jacobiScale_last    = variable.jacobiScale_last;
    jacobiScale_next    = variable.jacobiScale_next;
    
    time        = static.time;
    time_next   = static.time_next;
    time_last   = static.time_last;
    
    DAOindex_last   = static.DAOindex_last;
    DAOindex_next   = static.DAOindex_next;
    tIx_last        = static.tIx_last;
    tIx_next        = static.tIx_next;
    
    doseAngleBordersDiff        = static.doseAngleBordersDiff;
    doseAngleBordersDiff_last   = static.doseAngleBordersDiff_last;
    doseAngleBordersDiff_next   = static.doseAngleBordersDiff_next;
    timeFacCurr_last            = static.timeFacCurr_last;
    timeFacCurr_next            = static.timeFacCurr_next;
    fracFromLastDAO             = static.fracFromLastDAO;
    timeFracFromLastDAO         = static.timeFracFromLastDAO;
    timeFracFromNextDAO         = static.timeFracFromNextDAO;
    
    
end

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
coveredByRightLeaf = bsxfun(@minus,centres,rightLeafPosI)./repmat(rightLeafPosF-rightLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + (rightLeafPosI-edges_l(xPosIndRightLeafI)').^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafI)').*2);
coveredByRightLeaf(xPosLinearIndRightLeafF) = coveredByRightLeaf(xPosLinearIndRightLeafF) - (edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafF)').*2);
%round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = bsxfun(@minus,centres,leftLeafPosF)./(repmat(leftLeafPosF-leftLeafPosI,1,numBix)).^2;
dUl_dLF = bsxfun(@minus,leftLeafPosI,centres)./(repmat(leftLeafPosF-leftLeafPosI,1,numBix)).^2;

dCr_dRI = bsxfun(@minus,centres,rightLeafPosF)./(repmat(rightLeafPosF-rightLeafPosI,1,numBix)).^2;
dCr_dRF = bsxfun(@minus,rightLeafPosI,centres)./(repmat(rightLeafPosF-rightLeafPosI,1,numBix)).^2;

dUl_dLI(xPosLinearIndLeftLeafI) = dUl_dLI(xPosLinearIndLeftLeafI) + ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').*(2*leftLeafPosF-leftLeafPosI-edges_l(xPosIndLeftLeafI)'))./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLF(xPosLinearIndLeftLeafI) = dUl_dLF(xPosLinearIndLeftLeafI) - ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2)./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLI(xPosLinearIndLeftLeafF) = dUl_dLI(xPosLinearIndLeftLeafF) - ((edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2)./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafF)').*2);
dUl_dLF(xPosLinearIndLeftLeafF) = dUl_dLF(xPosLinearIndLeftLeafF) + ((edges_r(xPosIndLeftLeafF)'-leftLeafPosF).*(leftLeafPosF+edges_r(xPosIndLeftLeafF)'-2*leftLeafPosI))./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafF)').*2);

dCr_dRI(xPosLinearIndRightLeafI) = dCr_dRI(xPosLinearIndRightLeafI) + ((rightLeafPosI-edges_l(xPosIndRightLeafI)').*(2*rightLeafPosF-rightLeafPosI-edges_l(xPosIndRightLeafI)'))./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRF(xPosLinearIndRightLeafI) = dCr_dRF(xPosLinearIndRightLeafI) - ((rightLeafPosI-edges_l(xPosIndRightLeafI)').^2)./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRI(xPosLinearIndRightLeafF) = dCr_dRI(xPosLinearIndRightLeafF) - ((edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2)./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafF)').*2);
dCr_dRF(xPosLinearIndRightLeafF) = dCr_dRF(xPosLinearIndRightLeafF) + ((edges_r(xPosIndRightLeafF)'-rightLeafPosF).*(rightLeafPosF+edges_r(xPosIndRightLeafF)'-2*rightLeafPosI))./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafF)').*2);

for k = 1:n
    dUl_dLI(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLF(k,1:(xPosIndLeftLeafI(k)-1)) = 0;
    dUl_dLI(k,(xPosIndLeftLeafF(k)+1):numBix) = 0;
    dUl_dLF(k,(xPosIndLeftLeafF(k)+1):numBix) = 0;
    
    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafF(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary
        
        %19 July 2017 in journal
        dUl_dLI(k,xPosIndLeftLeafI(k)) = -1/(2*widths(xPosIndLeftLeafI(k))');
        dUl_dLF(k,xPosIndLeftLeafF(k)) = -1/(2*widths(xPosIndLeftLeafF(k))');
        if leftLeafPosF(k)-leftLeafPosI(k) <= eps(max(lim_r))
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafF(k)) = (edges_r(xPosIndLeftLeafF(k))-leftLeafPosF(k))./widths(xPosIndLeftLeafF(k));
        end
    end
    
    dCr_dRI(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRF(k,1:(xPosIndRightLeafI(k)-1)) = 0;
    dCr_dRI(k,(xPosIndRightLeafF(k)+1):numBix) = 0;
    dCr_dRF(k,(xPosIndRightLeafF(k)+1):numBix) = 0;
    
    if xPosIndRightLeafI(k) >= xPosIndRightLeafF(k)
        dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
        dCr_dRF(k,xPosIndRightLeafF(k)) = -1/(2*widths(xPosIndRightLeafF(k))');
        if rightLeafPosF(k)-rightLeafPosI(k) <= eps(max(lim_r))
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafF(k)) = (edges_r(xPosIndRightLeafF(k))-rightLeafPosF(k))./widths(xPosIndRightLeafF(k));
        end
    end
end

%% swap I and F if calculating final phase stuff
if variable.arcF
    
    [weightFactor_I,weightFactor_F] = deal(weightFactor_F,weightFactor_I);
    
    [dUl_dLI,dUl_dLF] = deal(dUl_dLF,dUl_dLI);
    [dCr_dRI,dCr_dRF] = deal(dCr_dRF,dCr_dRI);
    
    if static.DAOBeam
        [vectorIx_LI,vectorIx_LF] = deal(vectorIx_LF,vectorIx_LI);
        [vectorIx_RI,vectorIx_RF] = deal(vectorIx_RF,vectorIx_RI);
        
    else
        [vectorIx_LF_last,vectorIx_LI_next] = deal(vectorIx_LI_next,vectorIx_LF_last);
        [vectorIx_RF_last,vectorIx_RI_next] = deal(vectorIx_RI_next,vectorIx_RF_last);
        
        [fracFromNextOptI,fracFromLastOptF] = deal(fracFromLastOptF,fracFromNextOptI);
        [fracFromNextOptF,fracFromLastOptI] = deal(fracFromLastOptI,fracFromNextOptF);
        
    end
end

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
w{phase}(currBixelIx) = w{phase}(currBixelIx) + shapeMap_nW(shapeMapIx).*weight.*weightFactor_I.*probability;
shapeMap{phase} = shapeMap{phase}+shapeMap_nW.*weight.*weightFactor_I.*probability;

% store information for Jacobi preconditioning
sumGradSq{phase} = sumGradSq{phase}+(weightFactor_I.*probability).^2.*mean([sum((dUl_dLI).^2,2); sum((dUl_dLF.*weightFactor_F).^2,2); sum((dUl_dLF.*weightFactor_I).^2,2); sum((dCr_dRI).^2,2); sum((dCr_dRF.*weightFactor_F).^2,2); sum((dCr_dRF.*weightFactor_I).^2,2)]);

%% save the gradients

bixelJApVec_offset_temp = bixelJApVec_offset{phase};
numSaveBixel = nnz(shapeMapIx);

if static.DAOBeam
    % indices
    vectorIxMat_LI = repmat(vectorIx_LI',1,numBix);
    vectorIxMat_LF = repmat(vectorIx_LF',1,numBix);
    vectorIxMat_RI = repmat(vectorIx_RI',1,numBix);
    vectorIxMat_RF = repmat(vectorIx_RF',1,numBix);
    
    % wrt weight
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex+(phase-1)*totalNumOfShapes;
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt initial left
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final left
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt initial right
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final right
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
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
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOpt*(time./time_last)*shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_last;
    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
    %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
    %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex_last+(phase-1)*totalNumOfShapes;
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt next weight
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = (1-fracFromLastOpt)*(time./time_next)*shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_next;
    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
    %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
    %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = DAOindex_next+(phase-1)*totalNumOfShapes;
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    
    % updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
    %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
    
    % updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
    % updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
    
    % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt initial left (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptI(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromLastOptF(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final left (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptI(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = fracFromNextOptF(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt initial right (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptI(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromLastOptF(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt final right (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % extra weightFactor_F for I->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptI(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = -fracFromNextOptF(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt last time
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = weightFactor_I.*probability.*(doseAngleBordersDiff.*timeFacCurr_last) ...
        .*(-fracFromLastDAO.*timeFracFromNextDAO.*(weight_last./doseAngleBordersDiff_next).*(time_next./time_last.^2) ...
        +(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next./doseAngleBordersDiff_last).*(1./time_next)) ...
        * shapeMap_nW(shapeMapIx);
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = tIx_last;
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
    
    % wrt next time
    bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = weightFactor_I.*probability.*(doseAngleBordersDiff.*timeFacCurr_next) ...
        .*(fracFromLastDAO.*timeFracFromNextDAO.*(weight_last./doseAngleBordersDiff_next).*(1./time_last) ...
        -(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next./doseAngleBordersDiff_last).*(time_last./time_next.^2)) ...
        * shapeMap_nW(shapeMapIx);
    bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = tIx_next;
    bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel)) = currBixelIx;
    bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel;
end

% wrt times (probability)
bixelJApVec_vec{phase}(bixelJApVec_offset_temp+(1:numSaveBixel*totalNumOfShapes)) = weight.*weightFactor_I.*shapeMap_nW(shapeMapIx)*probability_dTVec';
bixelJApVec_i{phase}(bixelJApVec_offset_temp+(1:numSaveBixel*totalNumOfShapes)) = ones(numSaveBixel,1)*tIx_Vec;
bixelJApVec_j{phase}(bixelJApVec_offset_temp+(1:numSaveBixel*totalNumOfShapes)) = currBixelIx*ones(1,totalNumOfShapes);
bixelJApVec_offset_temp = bixelJApVec_offset_temp+numSaveBixel*totalNumOfShapes;

% update counters
bixelJApVec_offset{phase} = bixelJApVec_offset_temp;

%% update results and running

results.w                   = w;
results.bixelJApVec_vec     = bixelJApVec_vec;
results.bixelJApVec_i       = bixelJApVec_i;
results.bixelJApVec_j       = bixelJApVec_j;
results.bixelJApVec_offset  = bixelJApVec_offset;

running.sumGradSq   = sumGradSq;
running.shapeMap    = shapeMap;

end