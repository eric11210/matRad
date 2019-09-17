function [ct,cst] = matRad_importXCATfromBin(importOptions)

if importOptions.repPhase && ~contains(importOptions.fnameXcatRoot,'rep')
    importOptions.fnameXcatRoot = [importOptions.fnameXcatRoot '_rep'];
end

%% Convert binaries from XCAT to DICOM
dirDICOM = fullfile(fileparts(mfilename('fullpath')),'DICOM',filesep);

xcatLog = matRad_xcatBin2DICOM(importOptions);


%% Convert DICOM to ct structure
fprintf('matRad: Converting DICOM to ct struct ... \n');

% determine the number of frames we're going to generate
if importOptions.vmcDef
    % if we're using the vmc++ deformation code, we only need one
    numFrames = 1;
else
    % otherwise we need all of them
    numFrames = xcatLog.numFrames;
end

for frame = 1:numFrames
    
    fprintf('Frame %d of %d.\n',frame,numFrames);
    
    ctList = xcatLog.ctList(:,frame);
    
    ctTemp = matRad_importDicomCt(ctList,xcatLog.resolution,false);
    
    if ~exist('ct','var')
        ct = ctTemp;
        ct.cube = cell(numFrames,1);
        ct.cubeHU = cell(numFrames,1);
    end
    
    ct.cube{frame} = ctTemp.cube{1};
    ct.cubeHU{frame} = ctTemp.cubeHU{1};
end

% but use the correct number of frames here
ct.numOfCtScen = importOptions.numPhases;

%% Import vector files
fprintf('matRad: Importing motion vectors from XCAT files ... \n');

ct = matRad_parseMVF(ct,xcatLog,importOptions);

fprintf('\nDone!\n');

%% pad ct and interpolate mvf with zeros from vmc++
fprintf('matRad: Padding ct and interpolating motion vectors ... \n');

ct = matRad_padCtInterpMvf(ct,xcatLog,importOptions);

%% Bin ct frames based on tumour motion data from XCAT
fprintf('matRad: Binning frames into phases ... \n');

ct = matRad_binFrames2Phases(ct,importOptions);
% add field in importOptions if we don't want to bin it, i.e. frame = bin

%% Import DICOM structure set
fprintf('matRad: Importing DICOM structure set ... \n');

%fnameRTSS = fullfile(dirDICOM,sprintf('%s_RTSS.dcm',importOptions.fnameXcatRoot));
fnameRTSS = fullfile(dirDICOM,'XCAT_RTSS.dcm');
structures = matRad_importDicomRtss(fnameRTSS,ct.dicomInfo);

for i = 1:numel(structures)
    
    fprintf('Structure %d of %d.\n',i,numel(structures));
    structures(i).indices = matRad_convRtssContours2Indices(structures(i),ct);
end
fprintf('Done!\n');

fprintf('matRad: Adding tumour volumes to structure set... \n');
[structures,ct] = matRad_tumourContourWrapper(structures,ct,importOptions);
fprintf('Done!\n');

cst = matRad_createCst(structures);

%% Save file

saveName = importOptions.fnameXcatRoot;

if ~importOptions.keepAllFrames && importOptions.averageCT && ~importOptions.vmcDef
    saveName = [saveName '_avCT'];
end

if ~importOptions.keepAllFrames && importOptions.averageMVF
    saveName = [saveName '_avMVF'];
end

% rename vmcDef to ITV?
if ~importOptions.vmcDef
    saveName = [saveName '_ITV'];
end

save(saveName,'ct','cst','importOptions','-v7.3');

end

