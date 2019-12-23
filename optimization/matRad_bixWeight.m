function results = matRad_bixWeight(static,variable)

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
weight              = variable.weight;
weightFactor        = variable.weightFactor;

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

%% fix overshoot

for k = 1:numRow
    
    if xPosIndLeftLeafI(k) >= xPosIndLeftLeafF(k)
        % in discrete aperture, the xPosIndLeftLeafI is greater than
        % xPosIndLeftLeafM when leaf positions are at a bixel boundary
        
        %19 July 2017 in journal
        if leftLeafPosF(k)-leftLeafPosI(k) <= eps(max(lim_r))
            uncoveredByLeftLeaf(k,xPosIndLeftLeafI(k)) = (edges_r(xPosIndLeftLeafI(k))-leftLeafPosI(k))./widths(xPosIndLeftLeafI(k));
            uncoveredByLeftLeaf(k,xPosIndLeftLeafF(k)) = (edges_r(xPosIndLeftLeafF(k))-leftLeafPosF(k))./widths(xPosIndLeftLeafF(k));
        end
    end
    
    if xPosIndRightLeafI(k) >= xPosIndRightLeafF(k)
        if rightLeafPosF(k)-rightLeafPosI(k) <= eps(max(lim_r))
            coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
            coveredByRightLeaf(k,xPosIndRightLeafF(k)) = (edges_r(xPosIndRightLeafF(k))-rightLeafPosF(k))./widths(xPosIndRightLeafF(k));
        end
    end
end

%% swap I and F if calculating final phase stuff
if ~variable.arcF
    
    fracIToLastDose     = static.fracToLastDose_arcI;
    fracIToNextDose     = static.fracToNextDose_arcI;
    
else
    
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

shapeMap    = shapeMap_nW.*weight.*weightFactor.*probability;
wLast       = fracIToLastDose.*shapeMap(lastShapeMapIx);
wNext       = fracIToNextDose.*shapeMap(nextShapeMapIx);

%% update results

results.lastDose.w  = wLast;
results.nextDose.w  = wNext;
results.shapeMap    = shapeMap;


end