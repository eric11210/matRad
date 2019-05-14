function f = matRad_daoObjFunc(apertureInfoVec,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: objective function for direct aperture optimization
%
% call
%   f = matRad_daoObjFunc(apertureInfoVect,apertureInfo,dij,cst,options)  
%
% input
%   apertureInfoVect:   aperture info in form of vector
%   apertureInfo:       aperture info struct
%   dij:                matRad dij struct as generated by bixel-based dose calculation
%   cst:                matRad cst struct
%   options:            option struct defining the type of optimization
%
% output
%   f: objective function value
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
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

% read in the global apertureInfo and apertureVector variables
global matRad_global_apertureInfo;
% update apertureInfo from the global variable
apertureInfo = matRad_global_apertureInfo;

% also read in the global recalc variance variable
global matRad_global_recalcVar;

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
    matRad_global_apertureInfo = apertureInfo;
    
    % recalculate the variance
    matRad_global_recalcVar = true;
end

% bixel based objective function calculation
f = matRad_objFuncWrapper(apertureInfo.bixelWeights,dij,cst,options);

if apertureInfo.varOpt
    % add in variance term
    [fVar,~]    = matRad_varObjAndGradFunc(apertureInfo,dij,cst);
    f           = f+fVar;
end
