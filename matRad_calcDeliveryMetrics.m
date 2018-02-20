function result = matRad_calcDeliveryMetrics(result,pln,stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad delivery metric calculation
%
% call
%   matRad_calcDeliveryMetrics(result,pln)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   pln:                matRad plan meta information struct
%
% output
%   All plans: total MU
%   VMAT plans: total time, leaf speed, MU rate, and gantry rotation speed
%   distributions
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

apertureInfo = result.apertureInfo;

fileName = apertureInfo.propVMAT.machineConstraintFile;
try
    load([pwd filesep fileName],'machine');
catch
    error(['Could not find the following machine file: ' fileName ]);
end


apertureInfo.planMU = 0;

l = 0;
if pln.propOpt.runVMAT
    %All of these are vectors
    %Each entry corresponds to a beam angle
    %Later, we will convert these to histograms, find max, mean, min, etc.
    gantryRot = zeros(1,size(pln.propStf.DAOGantryAngles,2)-1);
    MURate = gantryRot;
    times = gantryRot;
    angles = gantryRot;
    maxLeafSpeed = 0*pln.propStf.DAOGantryAngles;
    
    
    totTime = 0;
    
    for i = 1:size(apertureInfo.beam,2)
        apertureInfo.planMU = apertureInfo.planMU+apertureInfo.beam(i).MU;
        
        totTime = totTime+apertureInfo.beam(i).time;
        
        if apertureInfo.beam(i).numOfShapes
            l = l+1;
            
            gantryRot(l) = apertureInfo.beam(i).gantryRot;
            MURate(l) = apertureInfo.beam(i).MURate*60;
            times(l) = apertureInfo.beam(i).time;
            maxLeafSpeed(l) = apertureInfo.beam(i).maxLeafSpeed/10;
            angles(l) = apertureInfo.beam(i).gantryAngle;
        end
    end
    apertureInfo.time = totTime;
    
    
    %apertureInfoVec = apertureInfo.apertureVector;
    %{
    if pln.dynamic
        leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
        timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
        timeDoseBorderAngles = timeOptBorderAngles.*[apertureInfo.beam([apertureInfo.beam.optimizeBeam]).timeFacCurr]';
        
        leftLeafDiff = diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
        rightLeafDiff = diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
        
        leftLeafDiff = reshape(leftLeafDiff(repmat([apertureInfo.beam.optimizeBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
        rightLeafDiff = reshape(rightLeafDiff(repmat([apertureInfo.beam.optimizeBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
        
        lfspd = reshape([leftLeafDiff rightLeafDiff]./ ...
            repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
        %{
        v = [leftLeafDiff rightLeafDiff]./ ...
            repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2);
        v = reshape(v,[],91);
        dv = diff(v,1,2);
        dv = reshape(dv,[],1);
        %}
            
            optAngles = [apertureInfo.beam([apertureInfo.beam.optimizeBeam]).gantryAngle];
            optAnglesMat = reshape(repmat(optAngles,apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
    else
        leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
        optInd = [apertureInfo.beam.optimizeBeam];
        timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
        
        i = repelem(1:(apertureInfo.totalNumOfShapes-1),2);
        j = repelem(1:(apertureInfo.totalNumOfShapes),2);
        j(1) = [];
        j(end) = [];
        
        timeFac = [apertureInfo.beam(optInd).timeFac]';
        timeFac(timeFac == 0) = [];
        
        timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
        timeBNOptAngles = timeFacMatrix*timeOptBorderAngles;
        
        lfspd = reshape([abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
            abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
            repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
    end
    
    if pln.dynamic
        initBorders = unique([stf([stf.initializeBeam]).initAngleBorders]);
        forwardDir = 1-2*mod(1:(numel(initBorders)-1),2);
        numForward = zeros(numel(forwardDir),1);
        numBackward = zeros(numel(forwardDir),1);
        timeInInit = zeros(numel(forwardDir),1);
        
        plot(optAnglesMat,lfspd,'.')
        hold on
        counter = 1;
        for border = initBorders
            plot([border border],[-resultGUI.apertureInfo.machineConstraints.leafSpeed(2) resultGUI.apertureInfo.machineConstraints.leafSpeed(2)],'r-')
            
            if border < initBorders(end)
                curr_lfspd = lfspd(initBorders(counter) <= optAnglesMat & optAnglesMat <= initBorders(counter+1));
                
                numForward(counter) = nnz(curr_lfspd*forwardDir(counter) >= 0);
                numBackward(counter) = nnz(curr_lfspd*forwardDir(counter) < 0);
                timeInInit(counter) = sum(times(initBorders(counter) <= angles & angles <= initBorders(counter+1)));
                
                counter = counter+1;
            end
        end
        plot([min(initBorders)-5 max(initBorders)+5],[0 0],'k--')
        xlim([min(initBorders)-5 max(initBorders)+5])
        ylim([-resultGUI.apertureInfo.machineConstraints.leafSpeed(2)-5 resultGUI.apertureInfo.machineConstraints.leafSpeed(2)+5])
        xlabel('gantry angle (^\circ)')
        ylabel('leaf speed (cm/s)')
        
        figure
        plot(optAngles,gantryRot,'.')
        xlim([min(initBorders)-5 max(initBorders)+5])
        ylim([0 resultGUI.apertureInfo.machineConstraints.gantryRotationSpeed(2)+1])
        xlabel('gantry angle (^\circ)')
        ylabel('gantry rotation speed (^\circ/s)')
        
        figure
        plot(optAngles,MURate,'.')
        xlim([min(initBorders)-5 max(initBorders)+5])
        ylim([0 60*resultGUI.apertureInfo.machineConstraints.monitorUnitRate(2)+5])
        xlabel('gantry angle (^\circ)')
        ylabel('MU rate (MU/min)')
        
        apertureInfo.fracMaxMURate = sum(times(MURate > 60*machine.constraints.monitorUnitRate(2)*(1-1e-5)))./sum(times);
        apertureInfo.fracMinMURate = sum(times(MURate < 60*machine.constraints.monitorUnitRate(1)*(1+1e-5)))./sum(times);
        apertureInfo.fracMaxGantryRot = sum(times(gantryRot > machine.constraints.gantryRotationSpeed(2)*(1-1e-5)))./sum(times);
        apertureInfo.fracMaxLeafSpeed = sum(times(maxLeafSpeed > machine.constraints.leafSpeed(2)/10*(1-1e-5)))./sum(times);
        apertureInfo.fracHalfMaxLeafSpeed = sum(times(maxLeafSpeed > machine.constraints.leafSpeed(2)/10*(1-1e-5)/2))./sum(times);
        
        apertureInfo.fracForward = numForward./(numForward+numBackward);
        apertureInfo.fracBackward = 1-apertureInfo.fracForward;
        apertureInfo.totalFracForward = mean(apertureInfo.fracForward);
        %apertureInfo.totalFracForward = sum(apertureInfo.fracForward.*timeInInit)./sum(timeInInit);
        apertureInfo.totalFracBackward = 1-apertureInfo.totalFracForward;
    end
    %}
    
end

result.apertureInfo = apertureInfo;

