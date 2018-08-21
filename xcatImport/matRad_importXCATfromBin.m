function [ct,cst] = matRad_importXCATfromBin(importOptions)

if importOptions.repPhase && ~strfind(importOptions.fnameXcatRoot,'rep')
    importOptions.fnameXcatRoot = [importOptions.fnameXcatRoot '_rep'];
end

%% Convert binaries from XCAT to DICOM
dirDICOM = fullfile(fileparts(mfilename('fullpath')),'DICOM',filesep);

xcatLog = matRad_xcatBin2DICOM(importOptions);


%% Convert DICOM to ct structure
fprintf('matRad: Converting DICOM to ct struct ... \n');

for frame = 1:xcatLog.numFrames
    
    fprintf('Phase %d of %d.\n',frame,xcatLog.numFrames);
    
    ctList = xcatLog.ctList(:,frame);
    
    ctTemp = matRad_importDicomCt(ctList,xcatLog.resolution,false);
    
    if ~exist('ct','var')
        ct = ctTemp;
        ct.cube = cell(xcatLog.numFrames,1);
        ct.cubeHU = cell(xcatLog.numFrames,1);
    end
    
    ct.cube{frame} = ctTemp.cube{1};
    ct.cubeHU{frame} = ctTemp.cubeHU{1};
end

ct.numOfCtScen = xcatLog.numFrames;

%% Import vector files
fprintf('matRad: Importing motion vectors from XCAT files ... \n');

ct = matRad_parseMVF(ct,xcatLog,importOptions);

fprintf('Done!\n');

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

if ~importOptions.keepAllFrames && importOptions.averageCT
    saveName = [saveName '_avCT'];
end

if ~importOptions.keepAllFrames && importOptions.averageMVF
    saveName = [saveName '_avMVF'];
end

save(saveName,'ct','cst','importOptions','-v7.3');

end

