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
    [dVar,dVarGrad] = matRad_doseVariance(apertureInfo,dij);
    
    % initialize fVar, fVarGrad
    fVar = 0;
    gVar = zeros(1,size(apertureInfo.apertureVector,1));
    
    % compute objective function for every VOI.
    for  i = 1:size(cst,1)
        
        % Only take OAR or target VOI.
        if ~isempty(cst{i,4}{1}) && isequal(cst{i,3},'OAR')
            
            % loop over the number of constraints for the current VOI
            for j = 1:numel(cst{i,6})
                
                % only perform objective computations for objectives
                if isempty(strfind(cst{i,6}(j).type,'constraint'))
                    
                    dVar_i      = dVar(cst{i,4}{1});
                    dVarGrad_i  = dVarGrad(cst{i,4}{1},:);
                    
                    fVar    = fVar + (cst{i,6}(j).penalty./dij.numOfFractions).*(sum(dVar_i,1)./numel(dVar_i));
                    gVar    = gVar + (cst{i,6}(j).penalty./dij.numOfFractions).*(sum(dVarGrad_i,1)./numel(dVarGrad_i));
                    
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