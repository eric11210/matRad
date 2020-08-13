function results = matRad_bixWeightAndGrad(static,variable)

% called by matRad_bixWeightAndGradWrapper

%round2 = @(a,b) round(a*10^b)/10^b;

%% extract variables from input

lim_r               = static.lim_r;
edges_l             = static.edges_l;
edges_r             = static.edges_r;
centres             = static.centres;
widths              = static.widths;
lastBixelIndMap     = static.lastBixelIndMap;
nextBixelIndMap     = static.nextBixelIndMap;
numRow              = static.numRow;
numCol              = static.numCol;
lastNumBix          = static.lastNumBix;
nextNumBix          = static.nextNumBix;
lastBixIndVec       = static.lastBixIndVec; % = 1:numBix
nextBixIndVec       = static.nextBixIndVec; % = 1:numBix

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
if ~static.fixedGantrySpeed
    probability_dTVec   = variable.probability_dTVec;
    tIx_Vec             = variable.tIx_Vec;
end
weight              = variable.weight;
weightFactor        = variable.weightFactor;

totalNumOfShapes    = variable.totalNumOfShapes;
lastNumShapbixIndVec    = variable.lastNumShapbixIndVec; % = 1:(totalNumOfShapes*numBix)
nextNumShapbixIndVec    = variable.nextNumShapbixIndVec; % = 1:(totalNumOfShapes*numBix)


if static.DAOBeam
    
    DAOindex    = variable.DAOindex;
    jacobiScale = variable.jacobiScale;
else
    
    fracFromLastDAO_MU          = static.fracFromLastDAO_MU;
    fracFromNextDAO_MU          = static.fracFromNextDAO_MU;
    fracFromLastDAO_gantryRot   = static.fracFromLastDAO_gantryRot;
    fracFromNextDAO_gantryRot   = static.fracFromNextDAO_gantryRot;
    
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
    
    fluAngleBordersDiff         = static.fluAngleBordersDiff;
    fluAngleBordersDiff_last    = static.fluAngleBordersDiff_last;
    fluAngleBordersDiff_next    = static.fluAngleBordersDiff_next;
    timeFacCurr_last            = static.timeFacCurr_last;
    timeFacCurr_next            = static.timeFacCurr_next;
end

vectorIx_L_lastDAOI = variable.vectorIx_L_lastDAOI;
vectorIx_L_lastDAOF = variable.vectorIx_L_lastDAOF;
vectorIx_L_nextDAOI = variable.vectorIx_L_nextDAOI;
vectorIx_L_nextDAOF = variable.vectorIx_L_nextDAOF;

vectorIx_R_lastDAOI = variable.vectorIx_R_lastDAOI;
vectorIx_R_lastDAOF = variable.vectorIx_R_lastDAOF;
vectorIx_R_nextDAOI = variable.vectorIx_R_nextDAOI;
vectorIx_R_nextDAOF = variable.vectorIx_R_nextDAOF;

fracFromLastDAOI_leafI = static.fracFromLastDAOI_leafI;
fracFromLastDAOF_leafI = static.fracFromLastDAOF_leafI;
fracFromNextDAOI_leafF = static.fracFromNextDAOI_leafF;
fracFromNextDAOF_leafF = static.fracFromNextDAOF_leafF;

%% initialize results, extract and running

% REINITIALIZE THESE EVERY TIME THIS IS RUN
% SUM THEM INSIDE THE WRAPPER FUNCTION
bixelJApVecLastDose_vec     = variable.bixelJApVecLastDose_zeros;
bixelJApVecLastDose_i       = variable.bixelJApVecLastDose_zeros;
bixelJApVecLastDose_j       = variable.bixelJApVecLastDose_zeros;
bixelJApVecLastDose_offset  = 0;

bixelJApVecNextDose_vec     = variable.bixelJApVecNextDose_zeros;
bixelJApVecNextDose_i       = variable.bixelJApVecNextDose_zeros;
bixelJApVecNextDose_j       = variable.bixelJApVecNextDose_zeros;
bixelJApVecNextDose_offset  = 0;

%% bixel weight calculation

% centres must be 1xn, leftLeafPosI must be mx1

%calculate fraction of fluence uncovered by left leaf
%initial computation
uncoveredByLeftLeaf = (centres-leftLeafPosI)./(leftLeafPosF-leftLeafPosI);
%correct for overshoot in initial and final leaf positions
uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafI) + (leftLeafPosI-edges_l(xPosIndLeftLeafI)').^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafI)').*2);
uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) = uncoveredByLeftLeaf(xPosLinearIndLeftLeafF) - (edges_r(xPosIndLeftLeafF)'-leftLeafPosF).^2./((leftLeafPosF-leftLeafPosI).*(widths(xPosIndLeftLeafF)').*2);
%round <0 to 0, >1 to 1
uncoveredByLeftLeaf(uncoveredByLeftLeaf < 0) = 0;
uncoveredByLeftLeaf(uncoveredByLeftLeaf > 1) = 1;

%calculate fraction of fluence covered by right leaf
%initial computation
coveredByRightLeaf = (centres-rightLeafPosI)./(rightLeafPosF-rightLeafPosI);
%correct for overshoot in initial and final leaf positions
coveredByRightLeaf(xPosLinearIndRightLeafI) = coveredByRightLeaf(xPosLinearIndRightLeafI) + (rightLeafPosI-edges_l(xPosIndRightLeafI)').^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafI)').*2);
coveredByRightLeaf(xPosLinearIndRightLeafF) = coveredByRightLeaf(xPosLinearIndRightLeafF) - (edges_r(xPosIndRightLeafF)'-rightLeafPosF).^2./((rightLeafPosF-rightLeafPosI).*(widths(xPosIndRightLeafF)').*2);
%round <0 to 0, >1 to 1
coveredByRightLeaf(coveredByRightLeaf < 0) = 0;
coveredByRightLeaf(coveredByRightLeaf > 1) = 1;

%% gradient calculation

dUl_dLI = (centres-leftLeafPosF)./(repmat(leftLeafPosF-leftLeafPosI,1,numCol)).^2;
dUl_dLF = (leftLeafPosI-centres)./(repmat(leftLeafPosF-leftLeafPosI,1,numCol)).^2;

dCr_dRI = (centres-rightLeafPosF)./(repmat(rightLeafPosF-rightLeafPosI,1,numCol)).^2;
dCr_dRF = (rightLeafPosI-centres)./(repmat(rightLeafPosF-rightLeafPosI,1,numCol)).^2;

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
if ~variable.arcF
    fracFromFluI_leftLeafI = variable.fracFromFluI_leftLeafI_arcI;
    fracFromFluI_leftLeafF = variable.fracFromFluI_leftLeafF_arcI;
    fracFromFluF_leftLeafI = variable.fracFromFluF_leftLeafI_arcI;
    fracFromFluF_leftLeafF = variable.fracFromFluF_leftLeafF_arcI;
    
    fracFromFluI_rightLeafI = variable.fracFromFluI_rightLeafI_arcI;
    fracFromFluI_rightLeafF = variable.fracFromFluI_rightLeafF_arcI;
    fracFromFluF_rightLeafI = variable.fracFromFluF_rightLeafI_arcI;
    fracFromFluF_rightLeafF = variable.fracFromFluF_rightLeafF_arcI;
    
    fracIToLastDose     = static.fracToLastDose_arcI;
    fracIToNextDose     = static.fracToNextDose_arcI;
    
else
    
    fracFromFluI_leftLeafI = variable.fracFromFluI_leftLeafI_arcF;
    fracFromFluI_leftLeafF = variable.fracFromFluI_leftLeafF_arcF;
    fracFromFluF_leftLeafI = variable.fracFromFluF_leftLeafI_arcF;
    fracFromFluF_leftLeafF = variable.fracFromFluF_leftLeafF_arcF;
    
    fracFromFluI_rightLeafI = variable.fracFromFluI_rightLeafI_arcF;
    fracFromFluI_rightLeafF = variable.fracFromFluI_rightLeafF_arcF;
    fracFromFluF_rightLeafI = variable.fracFromFluF_rightLeafI_arcF;
    fracFromFluF_rightLeafF = variable.fracFromFluF_rightLeafF_arcF;
    
    fracIToLastDose     = static.fracToLastDose_arcF;
    fracIToNextDose     = static.fracToNextDose_arcF;
end

%% save the bixel weights
%fluence is equal to fluence not covered by left leaf minus
%fluence covered by left leaf
shapeMap_nW = uncoveredByLeftLeaf-coveredByRightLeaf;
%shapeMap_nW = round2(shapeMap_nW,15);
shapeMap_nW(isnan(shapeMap_nW) | shapeMap_nW < 0) = 0;

% find open bixels
%shapeMapIx = shapeMap > 0;
% shapeMapIx = ~isnan(bixelIndMap);
lastShapeMapIx  = lastBixelIndMap ~= 0 & ~isnan(lastBixelIndMap);
nextShapeMapIx  = nextBixelIndMap ~= 0 & ~isnan(nextBixelIndMap);
lastBixelIx     = lastBixelIndMap(lastShapeMapIx);
nextBixelIx     = nextBixelIndMap(nextShapeMapIx);

shapeMap    = shapeMap_nW.*weight.*weightFactor.*probability;
wLast       = fracIToLastDose.*shapeMap(lastShapeMapIx);
wNext       = fracIToNextDose.*shapeMap(nextShapeMapIx);

%% save the gradients

% repeat indices
vectorIxMat_L_lastDAOI = repmat(vectorIx_L_lastDAOI',1,numCol);
vectorIxMat_L_lastDAOF = repmat(vectorIx_L_lastDAOF',1,numCol);
vectorIxMat_L_nextDAOI = repmat(vectorIx_L_nextDAOI',1,numCol);
vectorIxMat_L_nextDAOF = repmat(vectorIx_L_nextDAOF',1,numCol);

vectorIxMat_R_lastDAOI = repmat(vectorIx_R_lastDAOI',1,numCol);
vectorIxMat_R_lastDAOF = repmat(vectorIx_R_lastDAOF',1,numCol);
vectorIxMat_R_nextDAOI = repmat(vectorIx_R_nextDAOI',1,numCol);
vectorIxMat_R_nextDAOF = repmat(vectorIx_R_nextDAOF',1,numCol);

% repeat fractions from fluence beams
fracMatFromFluI_leftLeafI = repmat(fracFromFluI_leftLeafI,1,numCol);
fracMatFromFluI_leftLeafF = repmat(fracFromFluI_leftLeafF,1,numCol);
fracMatFromFluF_leftLeafI = repmat(fracFromFluF_leftLeafI,1,numCol);
fracMatFromFluF_leftLeafF = repmat(fracFromFluF_leftLeafF,1,numCol);

fracMatFromFluI_rightLeafI = repmat(fracFromFluI_rightLeafI,1,numCol);
fracMatFromFluI_rightLeafF = repmat(fracFromFluI_rightLeafF,1,numCol);
fracMatFromFluF_rightLeafI = repmat(fracFromFluF_rightLeafI,1,numCol);
fracMatFromFluF_rightLeafF = repmat(fracFromFluF_rightLeafF,1,numCol);

for i = 1:2
    % repeat entire gradient calculation to give appropriate contributions
    % to last (i = 1) and next (i = 2) dose calculation beams
    
    % set curr factors, indices, etc., to last/next as appropriate
    if i == 1
        % we are doing the last dose beam
        fracIToCurrDose             = fracIToLastDose;
        currBixIndVec               = lastBixIndVec;
        currNumBix                  = lastNumBix;
        currBixelIx                 = lastBixelIx;
        currShapeMapIx              = lastShapeMapIx;
        currNumShapbixIndVec        = lastNumShapbixIndVec;
        bixelJApVecCurrDose_vec     = bixelJApVecLastDose_vec;
        bixelJApVecCurrDose_i       = bixelJApVecLastDose_i;
        bixelJApVecCurrDose_j       = bixelJApVecLastDose_j;
        bixelJApVecCurrDose_offset  = bixelJApVecLastDose_offset;
    else
        % we are doing the next dose beam
        fracIToCurrDose             = fracIToNextDose;
        currBixIndVec               = nextBixIndVec;
        currNumBix                  = nextNumBix;
        currBixelIx                 = nextBixelIx;
        currShapeMapIx              = nextShapeMapIx;
        currNumShapbixIndVec        = nextNumShapbixIndVec;
        bixelJApVecCurrDose_vec     = bixelJApVecNextDose_vec;
        bixelJApVecCurrDose_i       = bixelJApVecNextDose_i;
        bixelJApVecCurrDose_j       = bixelJApVecNextDose_j;
        bixelJApVecCurrDose_offset  = bixelJApVecNextDose_offset;
    end
    
    % if the current factor is 0, just skip it
    if fracIToCurrDose == 0
        continue
    end
    
    % vectorize matrices
    shapeMap_nW_vec                     = shapeMap_nW(currShapeMapIx);
    fracMatFromFluI_leftLeafI_vec       = fracMatFromFluI_leftLeafI(currShapeMapIx);
    fracMatFromFluI_leftLeafF_vec       = fracMatFromFluI_leftLeafF(currShapeMapIx);
    fracMatFromFluF_leftLeafI_vec       = fracMatFromFluF_leftLeafI(currShapeMapIx);
    fracMatFromFluF_leftLeafF_vec       = fracMatFromFluF_leftLeafF(currShapeMapIx);
    fracMatFromFluI_rightLeafI_vec      = fracMatFromFluI_rightLeafI(currShapeMapIx);
    fracMatFromFluI_rightLeafF_vec      = fracMatFromFluI_rightLeafF(currShapeMapIx);
    fracMatFromFluF_rightLeafI_vec      = fracMatFromFluF_rightLeafI(currShapeMapIx);
    fracMatFromFluF_rightLeafF_vec      = fracMatFromFluF_rightLeafF(currShapeMapIx);
    dUl_dLI_vec                         = dUl_dLI(currShapeMapIx);
    dUl_dLF_vec                         = dUl_dLF(currShapeMapIx);
    dCr_dRI_vec                         = dCr_dRI(currShapeMapIx);
    dCr_dRF_vec                         = dCr_dRF(currShapeMapIx);
    vectorIxMat_L_lastDAOI_vec          = vectorIxMat_L_lastDAOI(currShapeMapIx);
    vectorIxMat_L_lastDAOF_vec          = vectorIxMat_L_lastDAOF(currShapeMapIx);
    vectorIxMat_L_nextDAOI_vec          = vectorIxMat_L_nextDAOI(currShapeMapIx);
    vectorIxMat_L_nextDAOF_vec          = vectorIxMat_L_nextDAOF(currShapeMapIx);
    vectorIxMat_R_lastDAOI_vec          = vectorIxMat_R_lastDAOI(currShapeMapIx);
    vectorIxMat_R_lastDAOF_vec          = vectorIxMat_R_lastDAOF(currShapeMapIx);
    vectorIxMat_R_nextDAOI_vec          = vectorIxMat_R_nextDAOI(currShapeMapIx);
    vectorIxMat_R_nextDAOF_vec          = vectorIxMat_R_nextDAOF(currShapeMapIx);
    
    % the only gradients that differ between DAO and non-DAO beams are
    % those wrt weight (and time, through MU/weight rate interpolation)
    if static.DAOBeam
        
        % wrt weight
        bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose.*shapeMap_nW_vec.*weightFactor.*probability./jacobiScale;
        bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = DAOindex; %DAOindex+(phase-1)*totalNumOfShapes;
        bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
        bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
        
    else
        
        % wrt last weight
        bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*fracFromLastDAO_MU*(time./time_last)*shapeMap_nW_vec.*weightFactor.*probability./jacobiScale_last;
        bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = DAOindex_last;%DAOindex_last+(phase-1)*totalNumOfShapes;
        bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
        bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
        
        % wrt next weight
        bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*fracFromNextDAO_MU*(time./time_next)*shapeMap_nW_vec.*weightFactor.*probability./jacobiScale_next;
        bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = DAOindex_next;%DAOindex_next+(phase-1)*totalNumOfShapes;
        bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
        bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
        
        if ~static.fixedGantrySpeed
            % only do time variables for variable gantry speed
            
            % wrt last time
            bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*weightFactor.*probability.*(fluAngleBordersDiff.*timeFacCurr_last) ...
                .*(-fracFromLastDAO_MU.*fracFromNextDAO_gantryRot.*(weight_last./fluAngleBordersDiff_next).*(time_next./time_last.^2) ...
                +(fracFromNextDAO_MU).*fracFromLastDAO_gantryRot.*(weight_next./fluAngleBordersDiff_last).*(1./time_next)) ...
                * shapeMap_nW_vec;
            bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = tIx_last;
            bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
            bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
            
            % wrt next time
            bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*weightFactor.*probability.*(fluAngleBordersDiff.*timeFacCurr_next) ...
                .*(fracFromLastDAO_MU.*fracFromNextDAO_gantryRot.*(weight_last./fluAngleBordersDiff_next).*(1./time_last) ...
                -(fracFromNextDAO_MU).*fracFromLastDAO_gantryRot.*(weight_next./fluAngleBordersDiff_last).*(time_last./time_next.^2)) ...
                * shapeMap_nW_vec;
            bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = tIx_next;
            bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
            bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
        end
    end
    
    % leftLeafPos_I   = fracFromLastDAOI_leafI.*leftLeafPos_lastDAOI+fracFromLastDAOF_leafI.*leftLeafPos_lastDAOF;
    % rightLeafPos_I  = fracFromLastDAOI_leafI.*rightLeafPos_lastDAOI+fracFromLastDAOF_leafI.*rightLeafPos_lastDAOF;
    
    % leftLeafPos_F   = fracFromNextDAOI_leafF.*leftLeafPos_nextDAOI+fracFromNextDAOF_leafF.*leftLeafPos_nextDAOF;
    % rightLeafPos_F  = fracFromNextDAOI_leafF.*rightLeafPos_nextDAOI+fracFromNextDAOF_leafF.*rightLeafPos_nextDAOF;
    
    % leftLeafPosM    = weightFactor_F.*leftLeafPosI+weightFactor_I.*leftLeafPosF;
    % rightLeafPosM   = weightFactor_F.*rightLeafPosI+weightFactor_I.*rightLeafPosF;
    
    % wrt last initial left (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*fracFromLastDAOI_leafI.*(fracMatFromFluI_leftLeafI_vec.*dUl_dLI_vec+fracMatFromFluI_leftLeafF_vec.*dUl_dLF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_L_lastDAOI_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    % wrt last final left (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*fracFromLastDAOF_leafI.*(fracMatFromFluI_leftLeafI_vec.*dUl_dLI_vec+fracMatFromFluI_leftLeafF_vec.*dUl_dLF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_L_lastDAOF_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    % wrt next initial left (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*fracFromNextDAOI_leafF.*(fracMatFromFluF_leftLeafI_vec.*dUl_dLI_vec+fracMatFromFluF_leftLeafF_vec.*dUl_dLF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_L_nextDAOI_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    % wrt next final left (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = fracIToCurrDose*fracFromNextDAOF_leafF.*(fracMatFromFluF_leftLeafI_vec.*dUl_dLI_vec+fracMatFromFluF_leftLeafF_vec.*dUl_dLF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_L_nextDAOF_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    
    % wrt last initial right (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = -fracIToCurrDose*fracFromLastDAOI_leafI.*(fracMatFromFluI_rightLeafI_vec.*dCr_dRI_vec+fracMatFromFluI_rightLeafF_vec.*dCr_dRF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_R_lastDAOI_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    % wrt last final right (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = -fracIToCurrDose*fracFromLastDAOF_leafI.*(fracMatFromFluI_rightLeafI_vec.*dCr_dRI_vec+fracMatFromFluI_rightLeafF_vec.*dCr_dRF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_R_lastDAOF_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    % wrt next initial right (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = -fracIToCurrDose*fracFromNextDAOI_leafF.*(fracMatFromFluF_rightLeafI_vec.*dCr_dRI_vec+fracMatFromFluF_rightLeafF_vec.*dCr_dRF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_R_nextDAOI_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    % wrt next final right (optimization vector)
    bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currBixIndVec) = -fracIToCurrDose*fracFromNextDAOF_leafF.*(fracMatFromFluF_rightLeafI_vec.*dCr_dRI_vec+fracMatFromFluF_rightLeafF_vec.*dCr_dRF_vec).*weight.*weightFactor.*probability;
    bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currBixIndVec) = vectorIxMat_R_nextDAOF_vec;
    bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currBixIndVec) = currBixelIx;
    bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix;
    
    if ~static.fixedGantrySpeed
        % only do time variables for variable gantry speed
        
        % wrt times (probability)
        bixelJApVecCurrDose_vec(bixelJApVecCurrDose_offset+currNumShapbixIndVec) = reshape(fracIToCurrDose*weight.*weightFactor.*shapeMap_nW_vec*probability_dTVec',currNumBix*totalNumOfShapes,1);
        bixelJApVecCurrDose_i(bixelJApVecCurrDose_offset+currNumShapbixIndVec) = repelem(tIx_Vec',currNumBix,1);
        bixelJApVecCurrDose_j(bixelJApVecCurrDose_offset+currNumShapbixIndVec) = repmat(currBixelIx,totalNumOfShapes,1);
        bixelJApVecCurrDose_offset = bixelJApVecCurrDose_offset+currNumBix.*totalNumOfShapes;
    end
    
    % get the current bixelJApVecCurrDose stuff, put into last/next as
    % appropriate
    if i == 1
        % we are doing the last dose beam
        bixelJApVecLastDose_vec     = bixelJApVecCurrDose_vec;
        bixelJApVecLastDose_i       = bixelJApVecCurrDose_i;
        bixelJApVecLastDose_j       = bixelJApVecCurrDose_j;
        bixelJApVecLastDose_offset  = bixelJApVecCurrDose_offset;
    else
        % we are doing the next dose beam
        bixelJApVecNextDose_vec     = bixelJApVecCurrDose_vec;
        bixelJApVecNextDose_i       = bixelJApVecCurrDose_i;
        bixelJApVecNextDose_j       = bixelJApVecCurrDose_j;
        bixelJApVecNextDose_offset  = bixelJApVecCurrDose_offset;
    end
    
end

%% store information for Jacobi preconditioning

if ~static.DAOBeam
    % calculate the sum of the square of gradients wrt weight
    % include square of: probability and weight factor, etc.
    results.sumGradSq_weight_lastDAO = (weightFactor.*probability.*fracFromLastDAO_MU.*(time./time_last)).^2.*sum(shapeMap_nW(:).^2);
    results.sumGradSq_weight_nextDAO = (weightFactor.*probability.*fracFromNextDAO_MU.*(time./time_next)).^2.*sum(shapeMap_nW(:).^2);
end

% calculate the sum of the square of gradients wrt initial and final
% left and right leaf positions
% first sum over all like, then take mean
sumGadSq_nw_L_lastDAOI = fracFromLastDAOI_leafI.*sum([(dUl_dLI.*fracFromFluI_leftLeafI).^2 (dUl_dLF.*fracFromFluI_leftLeafF).^2],2);
sumGadSq_nw_L_lastDAOF = fracFromLastDAOF_leafI.*sum([(dUl_dLI.*fracFromFluI_leftLeafI).^2 (dUl_dLF.*fracFromFluI_leftLeafF).^2],2);
sumGadSq_nw_L_nextDAOI = fracFromNextDAOI_leafF.*sum([(dUl_dLI.*fracFromFluF_leftLeafI).^2 (dUl_dLF.*fracFromFluF_leftLeafF).^2],2);
sumGadSq_nw_L_nextDAOF = fracFromNextDAOF_leafF.*sum([(dUl_dLI.*fracFromFluF_leftLeafI).^2 (dUl_dLF.*fracFromFluF_leftLeafF).^2],2);
sumGadSq_nw_R_lastDAOI = fracFromLastDAOI_leafI.*sum([(dCr_dRI.*fracFromFluI_rightLeafI).^2 (dCr_dRF.*fracFromFluI_rightLeafF).^2],2);
sumGadSq_nw_R_lastDAOF = fracFromLastDAOF_leafI.*sum([(dCr_dRI.*fracFromFluI_rightLeafI).^2 (dCr_dRF.*fracFromFluI_rightLeafF).^2],2);
sumGadSq_nw_R_nextDAOI = fracFromNextDAOI_leafF.*sum([(dCr_dRI.*fracFromFluF_rightLeafI).^2 (dCr_dRF.*fracFromFluF_rightLeafF).^2],2);
sumGadSq_nw_R_nextDAOF = fracFromNextDAOF_leafF.*sum([(dCr_dRI.*fracFromFluF_rightLeafI).^2 (dCr_dRF.*fracFromFluF_rightLeafF).^2],2);

results.sumGradSq_leaf_lastDAO = (weight.*weightFactor.*probability).^2.*mean([sumGadSq_nw_L_lastDAOI; sumGadSq_nw_L_lastDAOF; sumGadSq_nw_R_lastDAOI; sumGadSq_nw_R_lastDAOF],1);
results.sumGradSq_leaf_nextDAO = (weight.*weightFactor.*probability).^2.*mean([sumGadSq_nw_L_nextDAOI; sumGadSq_nw_L_nextDAOF; sumGadSq_nw_R_nextDAOI; sumGadSq_nw_R_nextDAOF],1);

%% remove zeros

delInd                          = bixelJApVecLastDose_vec == 0;
bixelJApVecLastDose_i(delInd)   = [];
bixelJApVecLastDose_j(delInd)   = [];
bixelJApVecLastDose_vec(delInd) = [];
bixelJApVecLastDose_offset      = numel(bixelJApVecLastDose_i);

delInd                          = bixelJApVecNextDose_vec == 0;
bixelJApVecNextDose_i(delInd)   = [];
bixelJApVecNextDose_j(delInd)   = [];
bixelJApVecNextDose_vec(delInd) = [];
bixelJApVecNextDose_offset      = numel(bixelJApVecNextDose_i);

%% update results

results.lastDose.w                  = wLast;
results.lastDose.bixelJApVec_vec    = bixelJApVecLastDose_vec;
results.lastDose.bixelJApVec_i      = bixelJApVecLastDose_i;
results.lastDose.bixelJApVec_j      = bixelJApVecLastDose_j;
results.lastDose.bixelJApVec_offset = bixelJApVecLastDose_offset;

results.nextDose.w                  = wNext;
results.nextDose.bixelJApVec_vec    = bixelJApVecNextDose_vec;
results.nextDose.bixelJApVec_i      = bixelJApVecNextDose_i;
results.nextDose.bixelJApVec_j      = bixelJApVecNextDose_j;
results.nextDose.bixelJApVec_offset = bixelJApVecNextDose_offset;

results.shapeMap    = shapeMap;

end