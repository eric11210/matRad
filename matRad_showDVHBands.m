function matRad_showDVHBands(dvh,pdvh,percentileIndices,cst,pln,lineStyleIndicator)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dvh visualizaion
% 
% call
%   matRad_showDVH(dvh,cst,pln,lineStyleIndicator)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   pln:                matRad pln struct
%   lineStyleIndicator: integer (1,2,3,4) to indicate the current linestyle
%                       (hint: use different lineStyles to overlay
%                       different dvhs)
%
% output
%   graphical display of DVH   
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('lineStyleIndicator','var') || isempty(lineStyleIndicator)
    lineStyleIndicator = 1;
end

% create new figure and set default line style indicator if not explictly
% specified
hold on;

numOfVois = size(cst,1);
        
%% print the dvh
colorMx    = colorcube;
colorMx    = colorMx(1:floor(64/numOfVois):64,:);

lineStyles = {'-',':','--','-.'};

maxDVHvol  = 0;
maxDVHdose = 0;

for i = 1:numOfVois
    if cst{i,5}.Visible
        % cut off at the first zero value where there is no more signal
        % behind
        ix      = max([1 find(dvh(i).volumePoints>0,1,'last')]);
        dvhOfMean = [dvh(i).doseGrid(1:ix).*pln.numOfFractions; dvh(i).volumePoints(1:ix)];
        
        plot(dvhOfMean(1,:),dvhOfMean(2,:),'LineWidth',4,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2})
        
        ix = pdvh(i).volumePoints(percentileIndices(1),:) > 0;
        lowerDVH = [pdvh(i).doseGrid(ix); pdvh(i).volumePoints(percentileIndices(1),ix)];
        
        ix = pdvh(i).volumePoints(percentileIndices(2),:) > 0;
        medianDVH = [pdvh(i).doseGrid(ix); pdvh(i).volumePoints(percentileIndices(2),ix)];
        
        ix = pdvh(i).volumePoints(percentileIndices(3),:) > 0;
        upperDVH = [pdvh(i).doseGrid(ix); pdvh(i).volumePoints(percentileIndices(3),ix)];
        
        plot(lowerDVH(1,:),lowerDVH(2,:),'LineWidth',4,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator+1},'DisplayName',cst{i,2},'HandleVisibility','off')
        
        plot(medianDVH(1,:),medianDVH(2,:),'LineWidth',4,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator+2},'DisplayName',cst{i,2},'HandleVisibility','off')
        
        plot(upperDVH(1,:),upperDVH(2,:),'LineWidth',4,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator+1},'DisplayName',cst{i,2},'HandleVisibility','off')
        
        maxDVHvol  = max(maxDVHvol,max(upperDVH(2,:)));
        maxDVHdose = max(maxDVHdose,max(upperDVH(1,:)));
    end
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff

ylim([0 1.1*maxDVHvol]);
xlim([0 1.1*maxDVHdose]);

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if strcmp(pln.propOpt.bioOptimization,'none')
     xlabel('Dose [Gy]','FontSize',fontSizeValue);
else
     xlabel('RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
end
