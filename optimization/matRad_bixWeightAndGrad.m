function [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,tempMap] = ...
    matRad_bixWeightAndGrad(apertureMotion,DAOBeam, ...
    lim_l,lim_r,edges_l,edges_r,centres,widths,n,numBix,bixelWidth, ...
    bixelIndMap,weight,jacobiScale, ...
    w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j, ...
    x1,x2,x3,x4,x5,x6,x7,x8)

round2 = @(a,b) round(a*10^b)/10^b;

if strcmp(apertureMotion,'continuous')
    
    [leftLeafPos_I,leftLeafPos_F,rightLeafPos_I,rightLeafPos_F] = deal(x1,x2,x3,x4);
    
    % set the initial leaf positions to the minimum leaf positions
    % always, instead of the leaf positions at the actual beginning
    % of the arc
    % this simplifies the calculation
    % remember which one is actually I and F in leftMinInd
    [leftLeafPosI,leftMinInd] = min([leftLeafPos_I,leftLeafPos_F],[],2);
    leftLeafPosF = max([leftLeafPos_I,leftLeafPos_F],[],2);
    [rightLeafPosI,rightMinInd] = min([rightLeafPos_I,rightLeafPos_F],[],2);
    rightLeafPosF = max([rightLeafPos_I,rightLeafPos_F],[],2);
    
    
    if DAOBeam
        [vectorIx_LI,vectorIx_LF,vectorIx_RI,vectorIx_RF] = deal(x5,x6,x7,x8);
        
        % change the vectorIx_xy elements to remember which
        % apertureVector elements the "new" I and F
        % if leftMinInd is 2, the I and F are switched
        tempL = vectorIx_LI;
        tempR = vectorIx_RI;
        vectorIx_LI(leftMinInd == 2) = vectorIx_LF(leftMinInd == 2);
        vectorIx_LF(leftMinInd == 2) = tempL(leftMinInd == 2);
        vectorIx_RI(rightMinInd == 2) = vectorIx_RF(rightMinInd == 2);
        vectorIx_RF(rightMinInd == 2) = tempR(rightMinInd == 2);
    else
        [vectorIx_LF_last,vectorIx_RF_last,vectorIx_LI_next,vectorIx_RI_next] = deal(x5,x6,x7,x8);
        
        tempL = vectorIx_LF_last;
        tempR = vectorIx_RF_last;
        
        vectorIx_LF_last(leftMinInd == 2) = vectorIx_LI_next(leftMinInd == 2);
        vectorIx_LI_next(leftMinInd == 2) = tempL(leftMinInd == 2);
        
        vectorIx_RF_last(rightMinInd == 2) = vectorIx_RI_next(rightMinInd == 2);
        vectorIx_RI_next(rightMinInd == 2) = tempR(rightMinInd == 2);
    end
    
    leftLeafPosI(leftLeafPosI <= lim_l) = lim_l(leftLeafPosI <= lim_l);
    leftLeafPosF(leftLeafPosF <= lim_l) = lim_l(leftLeafPosF <= lim_l);
    rightLeafPosI(rightLeafPosI <= lim_l) = lim_l(rightLeafPosI <= lim_l);
    rightLeafPosF(rightLeafPosF <= lim_l) = lim_l(rightLeafPosF <= lim_l);
    leftLeafPosI(leftLeafPosI >= lim_r) = lim_r(leftLeafPosI >= lim_r);
    leftLeafPosF(leftLeafPosF >= lim_r) = lim_r(leftLeafPosF >= lim_r);
    rightLeafPosI(rightLeafPosI >= lim_r) = lim_r(rightLeafPosI >= lim_r);
    rightLeafPosF(rightLeafPosF >= lim_r) = lim_r(rightLeafPosF >= lim_r);
    
    % find bixel indices where leaves are located
    xPosIndLeftLeafI = min(floor((leftLeafPosI-edges_l(1))./bixelWidth)+1,numBix);
    xPosIndLeftLeafF = max(ceil((leftLeafPosF-edges_r(1))./bixelWidth)+1,1);
    xPosIndRightLeafI = min(floor((rightLeafPosI-edges_l(1))./bixelWidth)+1,numBix);
    xPosIndRightLeafF = max(ceil((rightLeafPosF-edges_r(1))./bixelWidth)+1,1);
    %
    xPosLinearIndLeftLeafI = sub2ind([n numBix],(1:n)',xPosIndLeftLeafI);
    xPosLinearIndLeftLeafF = sub2ind([n numBix],(1:n)',xPosIndLeftLeafF);
    xPosLinearIndRightLeafI = sub2ind([n numBix],(1:n)',xPosIndRightLeafI);
    xPosLinearIndRightLeafF = sub2ind([n numBix],(1:n)',xPosIndRightLeafF);
    
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
    
    % gradients
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
        if xPosIndLeftLeafI(k) == xPosIndLeftLeafF(k)
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
        if xPosIndRightLeafI(k) == xPosIndRightLeafF(k)
            dCr_dRI(k,xPosIndRightLeafI(k)) = -1/(2*widths(xPosIndRightLeafI(k))');
            dCr_dRF(k,xPosIndRightLeafF(k)) = -1/(2*widths(xPosIndRightLeafF(k))');
            if rightLeafPosF(k)-rightLeafPosI(k) <= eps(max(lim_r))
                coveredByRightLeaf(k,xPosIndRightLeafI(k)) = (edges_r(xPosIndRightLeafI(k))-rightLeafPosI(k))./widths(xPosIndRightLeafI(k));
                coveredByRightLeaf(k,xPosIndRightLeafF(k)) = (edges_r(xPosIndRightLeafF(k))-rightLeafPosF(k))./widths(xPosIndRightLeafF(k));
            end
        end
    end
    
    
    % save the bixel weights
    %fluence is equal to fluence not covered by left leaf minus
    %fluence covered by left leaf
    tempMap = uncoveredByLeftLeaf-coveredByRightLeaf;
    tempMap = round2(tempMap,15);
    tempMap(isnan(tempMap)) = 0;
    
    % find open bixels
    %tempMapIx = tempMap > 0;
    tempMapIx = ~isnan(bixelIndMap);
    
    currBixelIx = bixelIndMap(tempMapIx);
    % look at bixelIndMap, probably have to add a factor to
    % correct this; either that or make w a cell array (lean
    % towards latter)
    % do currBixelIx+(phase-1)*dij.totalNumOfBixels?
    w{phase}(currBixelIx) = w{phase}(currBixelIx) + tempMap(tempMapIx)*weight;
    
    numSaveBixel = nnz(tempMapIx);
    
    if DAOBeam
        % indices
        vectorIx_LI = repmat(vectorIx_LI',1,numBix);
        vectorIx_LF = repmat(vectorIx_LF',1,numBix);
        vectorIx_RI = repmat(vectorIx_RI',1,numBix);
        vectorIx_RF = repmat(vectorIx_RF',1,numBix);
        
        % wrt weight
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = tempMap(tempMapIx)./jacobiScale;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = shapeInd;
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial left
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLI(tempMapIx)*weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLF(tempMapIx)*weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRI(tempMapIx)*weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRF(tempMapIx)*weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % store information for Jacobi preconditioning
        updatedInfo.beam(i).shape{phase}(j).sqrtSumGradSq = sqrt(mean([sum(dUl_dLI.^2,2); sum(dUl_dLF.^2,2); sum(dCr_dRI.^2,2); sum(dCr_dRF.^2,2)]));
        
        % increment shape index
        shapeInd = shapeInd +1;
    else
        % indices
        vectorIx_LF_last = repmat(vectorIx_LF_last',1,numBix);
        vectorIx_LI_next = repmat(vectorIx_LI_next',1,numBix);
        vectorIx_RF_last = repmat(vectorIx_RF_last',1,numBix);
        vectorIx_RI_next = repmat(vectorIx_RI_next',1,numBix);
        
        % leaf interpolation fractions/weights
        fracFromLastOptI = repmat(fracFromLastOptI,1,numBix);
        fracFromLastOptF = repmat(fracFromLastOptF,1,numBix);
        fracFromNextOptI = repmat(fracFromNextOptI,1,numBix);
        fracFromNextOptF = repmat(fracFromNextOptF,1,numBix);
        
        % wrt last weight
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOpt*(updatedInfo.beam(i).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time)*updatedInfo.beam(i).shape{phase}(j).shapeMap(tempMapIx) ...
            ./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).jacobiScale;
        %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
        %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(tempMapIx) ...
        %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOBeam])+(phase-1)*apertureInfo.totalNumOfShapes;
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt next weight
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = (1-fracFromLastOpt)*(updatedInfo.beam(i).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time)*updatedInfo.beam(i).shape{phase}(j).shapeMap(tempMapIx) ...
            ./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).jacobiScale;
        %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
        %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(tempMapIx) ...
        %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOBeam])+(phase-1)*apertureInfo.totalNumOfShapes;
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        
        %updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
        %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
        
        %updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
        %updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
        
        % wrt initial left (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOptI(tempMapIx).*dUl_dLI(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF_last(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOptF(tempMapIx).*dUl_dLF(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF_last(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final left (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromNextOptI(tempMapIx).*dUl_dLI(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI_next(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromNextOptF(tempMapIx).*dUl_dLF(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI_next(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt initial right (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromLastOptI(tempMapIx).*dCr_dRI(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF_last(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromLastOptF(tempMapIx).*dCr_dRF(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF_last(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt final right (optimization vector)
        % initial (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromNextOptI(tempMapIx).*dCr_dRI(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI_next(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        % final (interpolated arc)
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromNextOptF(tempMapIx).*dCr_dRF(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI_next(tempMapIx);
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt last time
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.propVMAT.beam(i).doseAngleBordersDiff.*updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr) ...
            .*(-updatedInfo.propVMAT.beam(i).fracFromLastDAO.*updatedInfo.propVMAT.beam(i).timeFracFromNextDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time.^2) ...
            +(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO).*updatedInfo.propVMAT.beam(i).timeFracFromLastDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(1./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time)) ...
            * updatedInfo.beam(i).shape{phase}(j).shapeMap(tempMapIx);
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOBeam])+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
        
        % wrt next time
        bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.propVMAT.beam(i).doseAngleBordersDiff.*updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr) ...
            .*(updatedInfo.propVMAT.beam(i).fracFromLastDAO.*updatedInfo.propVMAT.beam(i).timeFracFromNextDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(1./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time) ...
            -(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO).*updatedInfo.propVMAT.beam(i).timeFracFromLastDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time.^2)) ...
            * updatedInfo.beam(i).shape{phase}(j).shapeMap(tempMapIx);
        bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOBeam])+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;
        bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(tempMapIx);
        bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
    end
    
else
    
    [leftLeafPos,rightLeafPos] = deal(x1,x2);
    
end

end