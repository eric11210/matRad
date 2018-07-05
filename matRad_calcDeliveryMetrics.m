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

apertureInfo.planMU = 0;
%apertureInfo.planArea = 0;
%apertureInfo.planModulation = 0;

%apertureMU = nan(1000,1);
%apertureArea = nan(1000,1);
%l = 1;
if pln.propOpt.runVMAT
    
    for i = 1:numel(apertureInfo.beam)
        apertureInfo.planMU = apertureInfo.planMU+apertureInfo.beam(i).shape{1}(1).MU;
    end
    
    %initInd = find([apertureInfo.beam.initializeBeam]);
    %{
    for i = initInd
        apertureInfo.beam(i).beamMU = 0;
        apertureInfo.beam(i).beamArea = 0;
        apertureInfo.beam(i).unionShapeMap = 0*apertureInfo.beam(i).shape(1).shapeMap;
        
        for j = stf(i).beamChildrenIndex'
            apertureInfo.beam(j).shape(1).MU = apertureInfo.weightToMU*apertureInfo.beam(j).shape(1).weight;
            apertureInfo.beam(j).shape(1).apertureArea = (apertureInfo.bixelWidth/10)^2*sum(apertureInfo.beam(j).shape(1).shapeMap(:));
            
            apertureInfo.beam(i).beamMU = apertureInfo.beam(i).beamMU+apertureInfo.beam(j).shape(1).MU;
            apertureInfo.beam(i).beamArea = apertureInfo.beam(i).beamArea+apertureInfo.beam(j).shape(1).MU*apertureInfo.beam(j).shape(1).apertureArea;
            
            apertureInfo.beam(i).unionShapeMap = max(apertureInfo.beam(i).unionShapeMap,apertureInfo.beam(j).shape(1).shapeMap);
            
            apertureMU(l) = apertureInfo.beam(j).shape(1).MU;
            apertureArea(l) = apertureInfo.beam(j).shape(1).apertureArea;
            l = l+1;
        end
        apertureInfo.beam(i).beamArea = apertureInfo.beam(i).beamArea./apertureInfo.beam(i).beamMU;
        apertureInfo.beam(i).unionArea = (apertureInfo.bixelWidth/10)^2*sum(apertureInfo.beam(i).unionShapeMap(:));
        apertureInfo.beam(i).beamModulation = 1-apertureInfo.beam(i).beamArea./apertureInfo.beam(i).unionArea;
        apertureInfo.beam(i).beamK = apertureInfo.beam(i).numOfShapes*(1-apertureInfo.beam(i).beamModulation);
        
        apertureInfo.planMU = apertureInfo.planMU+apertureInfo.beam(i).beamMU;
        apertureInfo.planArea = apertureInfo.planArea+apertureInfo.beam(i).beamArea*apertureInfo.beam(i).beamMU;
        apertureInfo.planModulation = apertureInfo.planModulation+apertureInfo.beam(i).beamModulation*apertureInfo.beam(i).beamMU;
    end
    %}
end

%apertureInfo.planArea = apertureInfo.planArea./apertureInfo.planMU;
%apertureInfo.planModulation = apertureInfo.planModulation./apertureInfo.planMU;


%beamMU = [apertureInfo.beam(:).beamMU]';
%beamArea = [apertureInfo.beam(:).beamArea]';
%beamModulation = [apertureInfo.beam(:).beamModulation]';
%beamK = [apertureInfo.beam(:).beamK]';


%apertureMU(isnan(apertureMU)) = [];
%apertureArea(isnan(apertureArea)) = [];

%{
figure
histogram2(beamMU,beamArea,5,'DisplayStyle','tile')
xlabel('Monitor units')
ylabel('Beam area (cm^2)')
colorbar

figure
plot(apertureMU,apertureArea,'.')
xlabel('Monitor units')
ylabel('Aperture area (cm^2)')

figure
histogram2(beamMU,beamModulation,5,'DisplayStyle','tile')
xlabel('Monitor units')
ylabel('Beam modulation')
colorbar

figure
histogram2(beamArea,beamModulation,5,'DisplayStyle','tile')
xlabel('Area (cm^2)')
ylabel('Beam modulation')
colorbar

figure
plot(beamK,'.')
xlabel('Beam number')
ylabel('k = N*(1-BM)')

%}
l = 0;
if pln.propOpt.runVMAT
    
    fileName = apertureInfo.propVMAT.machineConstraintFile;
    try
        load(fileName,'machine');
    catch
        error(['Could not find the following machine file: ' fileName ]);
    end
    
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
        
        totTime = totTime+apertureInfo.beam(i).time; %time until next optimized beam
        
        if apertureInfo.beam(i).numOfShapes %only optimized beams have their time in the data struct
            l = l+1;
            gantryRot(l) = apertureInfo.beam(i).gantryRot;
            MURate(l) = apertureInfo.beam(i).shape{1}(1).MURate*60;
            times(l) = apertureInfo.beam(i).time;
            maxLeafSpeed(l) = apertureInfo.beam(i).maxLeafSpeed/10;
            angles(l) = apertureInfo.beam(i).gantryAngle;
        end
    end
    apertureInfo.time = totTime;
    
    
    apertureInfoVec = apertureInfo.apertureVector;
    if pln.propOpt.VMAToptions.continuousAperture
        leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
        timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
        timeDoseBorderAngles = timeOptBorderAngles.*[apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFacCurr]';
        
        leftLeafDiff = diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
        rightLeafDiff = diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
        
        leftLeafDiff = reshape(leftLeafDiff(repmat([apertureInfo.propVMAT.beam.DAOBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
        rightLeafDiff = reshape(rightLeafDiff(repmat([apertureInfo.propVMAT.beam.DAOBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
        
        lfspd = reshape([leftLeafDiff rightLeafDiff]./ ...
            repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
        %{
        v = [leftLeafDiff rightLeafDiff]./ ...
            repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2);
        v = reshape(v,[],91);
        dv = diff(v,1,2);
        dv = reshape(dv,[],1);
        %}
            
            optAngles = [apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).gantryAngle];
            optAnglesMat = reshape(repmat(optAngles,apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
    else
        leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
        optInd = [apertureInfo.propVMAT.beam.DAOBeam];
        timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
        
        i = repelem(1:(apertureInfo.totalNumOfShapes-1),2);
        j = repelem(1:(apertureInfo.totalNumOfShapes),2);
        j(1) = [];
        j(end) = [];
        
        timeFac = [apertureInfo.propVMAT.beam(optInd).timeFac]';
        timeFac(timeFac == 0) = [];
        
        timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
        timeBNOptAngles = timeFacMatrix*timeOptBorderAngles;
        
        lfspd = reshape([abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
            abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
            repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
    end
    
    if pln.propOpt.VMAToptions.continuousAperture
        
        FMOBorders = zeros(1,2*numel(pln.propStf.FMOGantryAngles));
        counter = 1;
        for i = 1:numel(pln.propStf.gantryAngles)
            if stf(i).propVMAT.FMOBeam
                FMOBorders(counter) = stf(i).propVMAT.FMOAngleBorders(1);
                FMOBorders(counter+1) = stf(i).propVMAT.FMOAngleBorders(2);
                counter = counter+2;
            else
                continue
            end
        end
        FMOBorders = unique(FMOBorders);
        forwardDir = 1-2*mod(1:(numel(FMOBorders)-1),2);
        numForward = zeros(numel(forwardDir),1);
        numBackward = zeros(numel(forwardDir),1);
        timeInInit = zeros(numel(forwardDir),1);
        
        plot(optAnglesMat,lfspd,'.')
        hold on
        counter = 1;
        for border = FMOBorders
            plot([border border],[-machine.constraints.leafSpeed(2) machine.constraints.leafSpeed(2)],'r-')
            
            if border < FMOBorders(end)
                curr_lfspd = lfspd(FMOBorders(counter) <= optAnglesMat & optAnglesMat <= FMOBorders(counter+1));
                
                numForward(counter) = nnz(curr_lfspd*forwardDir(counter) >= 0);
                numBackward(counter) = nnz(curr_lfspd*forwardDir(counter) < 0);
                timeInInit(counter) = sum(times(FMOBorders(counter) <= angles & angles <= FMOBorders(counter+1)));
                
                counter = counter+1;
            end
        end
        figure
        plot([min(FMOBorders)-5 max(FMOBorders)+5],[0 0],'k--')
        xlim([min(FMOBorders)-5 max(FMOBorders)+5])
        ylim([-machine.constraints.leafSpeed(2)-5 machine.constraints.leafSpeed(2)+5])
        xlabel('gantry angle (^\circ)')
        ylabel('leaf speed (cm/s)')
        
        figure
        plot(optAngles,gantryRot,'.')
        xlim([min(FMOBorders)-5 max(FMOBorders)+5])
        ylim([0 machine.constraints.gantryRotationSpeed(2)+1])
        xlabel('gantry angle (^\circ)')
        ylabel('gantry rotation speed (^\circ/s)')
        
        figure
        plot(optAngles,MURate,'.')
        xlim([min(FMOBorders)-5 max(FMOBorders)+5])
        ylim([0 60*machine.constraints.monitorUnitRate(2)+5])
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
    
end

result.apertureInfo = apertureInfo;

