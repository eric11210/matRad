function [fVar,gVar] = matRad_varObjAndGradFunc(apertureInfo,dij,cst)

global matRad_global_fVar;
global matRad_global_fVarGrad;
global matRad_global_recalcVar;

if isempty(matRad_global_recalcVar)
    recalcVar = true;
else
    recalcVar = matRad_global_recalcVar;
end

if recalcVar
    
    % calculate variance and gradient
    [dVarSum,dVarSumGrad] = matRad_doseVarianceSum(apertureInfo,dij);
    
    % initialize fVar, fVarGrad
    fVar = 0;
    gVar = zeros(1,size(apertureInfo.apertureVector,1));
    
    % compute objective function for every VOI.
    for  i = 1:size(cst,1)
        
        % Only take target VOI.
        if ~isempty(cst{i,4}{1}) && isequal(cst{i,3},'TARGET')
            
            % loop over the number of constraints for the current VOI
            for j = 1:numel(cst{i,6})
                
                % only perform objective computations for square deviation
                % objectives
                if strcmp(cst{i,6}(j).type,'square deviation')
                    
                    dVar_i      = numel(cst{i,4}{1}).*dVarSum./numel(dij.targetVox);
                    dVarGrad_i  = numel(cst{i,4}{1}).*dVarSumGrad./numel(dij.targetVox);
                    
                    fVar    = fVar + (cst{i,6}(j).penalty./dij.numOfFractions).*(dVar_i./numel(cst{i,4}{1}));
                    gVar    = gVar + (cst{i,6}(j).penalty./dij.numOfFractions).*(dVarGrad_i./numel(cst{i,4}{1}));
                    
                end
            end
        end
    end
    
    % transpose gradient
    gVar = gVar';
    
    % update global variables
    matRad_global_fVar      = fVar;
    matRad_global_fVarGrad  = gVar;
    matRad_global_recalcVar = false;
    
else
    
    % get objective function, gradient from variables
    fVar    = matRad_global_fVar;
    gVar    = matRad_global_fVarGrad;
end

end