function results = matRad_bixWeightAndGrad(static,variable)

%round2 = @(a,b) round(a*10^b)/10^b;

%% extract variables from input

lim_r               = static.lim_r;
edges_l             = static.edges_l;
edges_r             = static.edges_r;
centres             = static.centres;
widths              = static.widths;
bixelIndMap         = static.bixelIndMap;
numRow              = static.numRow;
numCol              = static.numCol;
numBix              = static.numBix;
bixIndVec           = static.bixIndVec; % = 1:numBix

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

probability         = variable.probability;
probability_dTVec   = variable.probability_dTVec;
weight              = variable.weight;
weightFactor_I      = variable.weightFactor_I;
weightFactor_F      = variable.weightFactor_F;

totalNumOfShapes    = variable.totalNumOfShapes;
tIx_Vec             = variable.tIx_Vec;
numShapbixIndVec    = variable.numShapbixIndVec; % = 1:(totalNumOfShapes*numBix)


if static.DAOBeam
    
    DAOindex    = variable.DAOindex;
    
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
    
    DAOindex_last   = variable.DAOindex_last;
    DAOindex_next   = variable.DAOindex_next;
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

%% initialize results, extract and running

% REINITIALIZE THESE EVERY TIME THIS IS RUN
% SUM THEM INSIDE THE WRAPPER FUNCTION
bixelJApVec_vec     = zeros(variable.bixelJApVec_sz,1);
bixelJApVec_i       = zeros(variable.bixelJApVec_sz,1);
bixelJApVec_j       = zeros(variable.bixelJApVec_sz,1);

%% bixel weight calculation

% centres must be 1xn, leftLeafPosI must be mx1

%calculate fraction of fluence uncovered by left leaf
%initial computation
uncoveredByLeftLeaf = (centres-leftLeafPosI)./(leftLeafPosF-leftLeafPosI);% bsxfun(@minus,centres,leftLeafPosI)./repmat(leftLeafPosF-leftLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) + (leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafI)').*2);
uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) - (edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafF)').*2);
%round <0 to 0, >1 to 1
uncoveredByLeftLeaf(uncoveredByLeftLeaf < 0) = 0;
uncoveredByLeftLeaf(uncoveredByLeftLeaf > 1) = 1;

%calculate fraction of fluence covered by right leaf
%initial computation
coveredByRightLeaf = (centres-rightLeafPosI)./(rightLeafPosF-rightLeafPosI);%bsxfun(@minus,centres,rightLeafPosI)./repmat(rightLeafPosF-rightLeafPosI,1,numBix);
%correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + (rightLeafPosI-edges_l(xPosIndRightLeafI)').^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafI)').*2);
coveredByRightLeaf(xPosLinearIndRightLeafF) = coveredByRightLeaf(xPosLinearIndRightLeafF) - (edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafF)').*2);
%round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = (centres-leftLeafPosF)./(repmat(leftLeafPosF-leftLeafPosI,1,numCol)).^2;%bsxfun(@minus,centres,leftLeafPosF)./(repmat(leftLeafPosF-leftLeafPosI,1,numBix)).^2;
dUl_dLF = (leftLeafPosI-centres)./(repmat(leftLeafPosF-leftLeafPosI,1,numCol)).^2;%bsxfun(@minus,leftLeafPosI,centres)./(repmat(leftLeafPosF-leftLeafPosI,1,numBix)).^2;

dCr_dRI = (centres-rightLeafPosF)./(repmat(rightLeafPosF-rightLeafPosI,1,numCol)).^2;%bsxfun(@minus,centres,rightLeafPosF)./(repmat(rightLeafPosF-rightLeafPosI,1,numBix)).^2;
dCr_dRF = (rightLeafPosI-centres)./(repmat(rightLeafPosF-rightLeafPosI,1,numCol)).^2;%bsxfun(@minus,rightLeafPosI,centres)./(repmat(rightLeafPosF-rightLeafPosI,1,numBix)).^2;

dUl_dLI(xPosLinearIndLeftLeafI) = dUl_dLI(xPosLinearIndLeftLeafI) + ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').*(2*leftLeafPosF-leftLeafPosI-edges_l(xPosIndLeftLeafI)'))./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLF(xPosLinearIndLeftLeafI) = dUl_dLF(xPosLinearIndLeftLeafI) - ((leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2)./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafI)').*2);
dUl_dLI(xPosLinearIndLeftLeafF) = dUl_dLI(xPosLinearIndLeftLeafF) - ((edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2)./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafF)').*2);
dUl_dLF(xPosLinearIndLeftLeafF) = dUl_dLF(xPosLinearIndLeftLeafF) + ((edges_r(xPosIndLeftLeafF)'-leftLeafPosF).*(leftLeafPosF+edges_r(xPosIndLeftLeafF)'-2*leftLeafPosI))./((leftLeafPosF-leftLeafPosI).^2.*(widths(xPosIndLeftLeafF)').*2);

dCr_dRI(xPosLinearIndRightLeafI) = dCr_dRI(xPosLinearIndRightLeafI) + ((rightLeafPosI-edges_l(xPosIndRightLeafI)').*(2*rightLeafPosF-rightLeafPosI-edges_l(xPosIndRightLeafI)'))./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRF(xPosLinearIndRightLeafI) = dCr_dRF(xPosLinearIndRightLeafI) - ((rightLeafPosI-edges_l(xPosIndRightLeafI)').^2)./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafI)').*2);
dCr_dRI(xPosLinearIndRightLeafF) = dCr_dRI(xPosLinearIndRightLeafF) - ((edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2)./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafF)').*2);
dCr_dRF(xPosLinearIndRightLeafF) = dCr_dRF(xPosLinearIndRightLeafF) + ((edges_r(xPosIndRightLeafF)'-rightLeafPosF).*(rightLeafPosF+edges_r(xPosIndRightLeafF)'-2*rightLeafPosI))./((rightLeafPosF-rightLeafPosI).^2.*(widths(xPosIndRightLeafF)').*2);

dUl_dLI(uncoveredByLeftLeaf == 0 | uncoveredByLeftLeaf == 1) = 0;
dUl_dLF(uncoveredByLeftLeaf == 0 | uncoveredByLeftLeaf == 1) = 0;

dCr_dRI(coveredByRightLeaf == 0 | coveredByRightLeaf == 1) = 0;
dCr_dRF(coveredByRightLeaf == 0 | coveredByRightLeaf == 1) = 0;

for k = 1:numRow
    
    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafF(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary
        
        %19 July 2017 in journal
        dUl_dLI(k,xPosIndLeftLeafI(k)) = -1/(2*widths(xPosIndLeftLeafI(k))');
        dUl_dLF(k,xPosIndLeftLeafF(k)) = -1/(2*widths(xPosIndLeftLeafF(k))');
        if leftLeafPosF(k)-leftLeafPosI(k) <= eps(max(lim_r))
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafF(k)) = (edges_r(xPosIndLeftLeafF(k))-leftLeafPosF(k))./widths(xPosIndLeftLeafF(k));
            
            if xPosIndLeftLeafI(k) > xPosIndLeftLeafF(k)
                if xPosIndLeftLeafI(k) > 1
                    dUl_dLI(k,xPosIndLeftLeafI(k)-1) = -1/(2*widths(xPosIndLeftLeafI(k))');
                end
                if xPosIndLeftLeafF(k) < size(dUl_dLF,2)
                    dUl_dLF(k,xPosIndLeftLeafF(k)+1) = -1/(2*widths(xPosIndLeftLeafF(k))');
                end
            end
        end
    end
    
    if xPosIndRightLeafI(k) >= xPosIndRightLeafF(k)
        dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
        dCr_dRF(k,xPosIndRightLeafF(k)) = -1/(2*widths(xPosIndRightLeafF(k))');
        if rightLeafPosF(k)-rightLeafPosI(k) <= eps(max(lim_r))
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafF(k)) = (edges_r(xPosIndRightLeafF(k))-rightLeafPosF(k))./widths(xPosIndRightLeafF(k));
            
            if xPosIndRightLeafI(k) > xPosIndRightLeafF(k)
                if xPosIndRightLeafI(k) > 1
                    dCr_dRI(k,xPosIndRightLeafI(k)-1) = -1/(2*widths(xPosIndRightLeafI(k)-1)');
                end
                if xPosIndRightLeafF(k) < size(dCr_dRF,2)
                    dCr_dRF(k,xPosIndRightLeafF(k)+1) = -1/(2*widths(xPosIndRightLeafF(k)+1)');
                end
            end
        end
    end
end

%% swap I and F if calculating final phase stuff
if variable.arcF
    
    temp = weightFactor_I;
    weightFactor_I = weightFactor_F;
    weightFactor_F = temp;
    
    temp = dUl_dLI;
    dUl_dLI = dUl_dLF;
    dUl_dLF = temp;
    
    temp = dCr_dRI;
    dCr_dRI = dCr_dRF;
    dCr_dRF = temp;
    
    if static.DAOBeam
        
        temp = vectorIx_LI;
        vectorIx_LI = vectorIx_LF;
        vectorIx_LF = temp;
        
        temp = vectorIx_RI;
        vectorIx_RI = vectorIx_RF;
        vectorIx_RF = temp;
        
    else
        
        temp = vectorIx_LF_last;
        vectorIx_LF_last = vectorIx_LI_next;
        vectorIx_LI_next = temp;
        
        temp = vectorIx_RF_last;
        vectorIx_RF_last = vectorIx_RI_next;
        vectorIx_RI_next = temp;
        
        temp = fracFromNextOptI;
        fracFromNextOptI = fracFromLastOptF;
        fracFromLastOptF = temp;
        
        temp = fracFromNextOptF;
        fracFromNextOptF = fracFromLastOptI;
        fracFromLastOptI = temp;
    end
end

%% save the bixel weights
%fluence is equal to fluence not covered by left leaf minus
%fluence covered by left leaf
shapeMap_nW = uncoveredByLeftLeaf-coveredByRightLeaf;
%shapeMap_nW = round2(shapeMap_nW,15);
shapeMap_nW(isnan(shapeMap_nW)) = 0;

% find open bixels
%shapeMapIx = shapeMap > 0;
% shapeMapIx = ~isnan(bixelIndMap);
shapeMapIx = bixelIndMap ~= 0 & ~isnan(bixelIndMap);

currBixelIx = bixelIndMap(shapeMapIx);
w = shapeMap_nW(shapeMapIx).*weight.*weightFactor_I.*probability;
shapeMap = shapeMap_nW.*weight.*weightFactor_I.*probability;

% store information for Jacobi preconditioning
sumGradSq = (weightFactor_I.*probability).^2.*mean([sum((dUl_dLI).^2,2); sum((dUl_dLF.*weightFactor_F).^2,2); sum((dUl_dLF.*weightFactor_I).^2,2); sum((dCr_dRI).^2,2); sum((dCr_dRF.*weightFactor_F).^2,2); sum((dCr_dRF.*weightFactor_I).^2,2)]);

%% save the gradients

bixelJApVec_offset = 0;

if static.DAOBeam
    % indices
    vectorIxMat_LI = repmat(vectorIx_LI',1,numCol);
    vectorIxMat_LF = repmat(vectorIx_LF',1,numCol);
    vectorIxMat_RI = repmat(vectorIx_RI',1,numCol);
    vectorIxMat_RF = repmat(vectorIx_RF',1,numCol);
    
    % wrt weight
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = DAOindex; %DAOindex+(phase-1)*totalNumOfShapes;
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt initial left
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = dUl_dLI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LI(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % extra weightFactor_F for I->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LI(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt final left
    % extra weightFactor_I for F->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LF(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt initial right
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -dCr_dRI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RI(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % extra weightFactor_F for I->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RI(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt final right
    % extra weightFactor_I for F->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RF(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
else
    % indices
    vectorIxMat_LF_last = repmat(vectorIx_LF_last',1,numCol);
    vectorIxMat_LI_next = repmat(vectorIx_LI_next',1,numCol);
    vectorIxMat_RF_last = repmat(vectorIx_RF_last',1,numCol);
    vectorIxMat_RI_next = repmat(vectorIx_RI_next',1,numCol);
    
    % leaf interpolation fractions/weights
    fracFromLastOptI = repmat(fracFromLastOptI,1,numCol);
    fracFromLastOptF = repmat(fracFromLastOptF,1,numCol);
    fracFromNextOptI = repmat(fracFromNextOptI,1,numCol);
    fracFromNextOptF = repmat(fracFromNextOptF,1,numCol);
    
    % wrt last weight
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromLastOpt*(time./time_last)*shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_last;
    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
    %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
    %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = DAOindex_last;%DAOindex_last+(phase-1)*totalNumOfShapes;
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt next weight
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = (1-fracFromLastOpt)*(time./time_next)*shapeMap_nW(shapeMapIx).*weightFactor_I.*probability./jacobiScale_next;
    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
    %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(shapeMapIx) ...
    %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = DAOindex_next;%DAOindex_next+(phase-1)*totalNumOfShapes;
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    
    % updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
    %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
    
    % updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
    % updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
    
    % leftLeafPosM   = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM  = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt initial left (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromLastOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % extra weightFactor_F for I->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromLastOptI(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromLastOptF(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LF_last(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt final left (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromNextOptI(shapeMapIx).*dUl_dLI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % extra weightFactor_F for I->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromNextOptI(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = fracFromNextOptF(shapeMapIx).*dUl_dLF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_LI_next(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt initial right (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -fracFromLastOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % extra weightFactor_F for I->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -fracFromLastOptI(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -fracFromLastOptF(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RF_last(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt final right (optimization vector)
    % initial (interpolated arc)
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -fracFromNextOptI(shapeMapIx).*dCr_dRI(shapeMapIx).*weight.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % extra weightFactor_F for I->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -fracFromNextOptI(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_F.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    % final (interpolated arc)
    % extra weightFactor_I for F->M
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = -fracFromNextOptF(shapeMapIx).*dCr_dRF(shapeMapIx).*weight.*weightFactor_I.*weightFactor_I.*probability;
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = vectorIxMat_RI_next(shapeMapIx);
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt last time
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = weightFactor_I.*probability.*(doseAngleBordersDiff.*timeFacCurr_last) ...
        .*(-fracFromLastDAO.*timeFracFromNextDAO.*(weight_last./doseAngleBordersDiff_next).*(time_next./time_last.^2) ...
        +(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next./doseAngleBordersDiff_last).*(1./time_next)) ...
        * shapeMap_nW(shapeMapIx);
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = tIx_last;
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
    
    % wrt next time
    bixelJApVec_vec(bixelJApVec_offset+bixIndVec) = weightFactor_I.*probability.*(doseAngleBordersDiff.*timeFacCurr_next) ...
        .*(fracFromLastDAO.*timeFracFromNextDAO.*(weight_last./doseAngleBordersDiff_next).*(1./time_last) ...
        -(1-fracFromLastDAO).*timeFracFromLastDAO.*(weight_next./doseAngleBordersDiff_last).*(time_last./time_next.^2)) ...
        * shapeMap_nW(shapeMapIx);
    bixelJApVec_i(bixelJApVec_offset+bixIndVec) = tIx_next;
    bixelJApVec_j(bixelJApVec_offset+bixIndVec) = currBixelIx;
    bixelJApVec_offset = bixelJApVec_offset+numBix;
end

% wrt times (probability)
bixelJApVec_vec(bixelJApVec_offset+numShapbixIndVec) = reshape(weight.*weightFactor_I.*shapeMap_nW(shapeMapIx)*probability_dTVec',numBix*totalNumOfShapes,1);
bixelJApVec_i(bixelJApVec_offset+numShapbixIndVec) = repelem(tIx_Vec',numBix,1);%reshape(ones(numSaveBixel,1,'like',tIx_Vec)*tIx_Vec,numSaveBixel*totalNumOfShapes,1);
bixelJApVec_j(bixelJApVec_offset+numShapbixIndVec) = repmat(currBixelIx,totalNumOfShapes,1);%reshape(currBixelIx*ones(1,totalNumOfShapes,'currBixelIx'),numSaveBixel*totalNumOfShapes,1);
bixelJApVec_offset = bixelJApVec_offset+numBix*totalNumOfShapes;

%% update results and running

results.w                   = w;
results.bixelJApVec_vec     = bixelJApVec_vec;
results.bixelJApVec_i       = bixelJApVec_i;
results.bixelJApVec_j       = bixelJApVec_j;
results.bixelJApVec_offset  = bixelJApVec_offset;

results.sumGradSq   = sumGradSq;
results.shapeMap    = shapeMap;

end