function [ct,cst] = matRad_importXCATfromBin(fnameXcatRoot)

%% Convert binaries from XCAT to DICOM
dirDICOM = fullfile(fileparts(mfilename('fullpath')),'DICOM',filesep);

xcatLog = matRad_xcatBin2DICOM(fnameXcatRoot);


%% Convert DICOM to ct structure
fprintf('matRad: Converting DICOM to ct struct ... \n');

for phase = 1:xcatLog.numPhases
    
    fprintf('Phase %d of %d.\n',phase,xcatLog.numPhases);
    
    ctList = xcatLog.ctList(:,phase);
    
    ctTemp = matRad_importDicomCt(ctList,xcatLog.resolution,false);
    
    if ~exist('ct','var')
        ct = ctTemp;
    end
    
    ct.cube{phase} = ctTemp.cube{1};
    ct.cubeHU{phase} = ctTemp.cubeHU{1};
end

ct.numOfCtScen = xcatLog.numPhases;

%delete(fullfile(dirDICOM,'*Slice*.dcm'));

%% Import vector files
fprintf('matRad: Importing motion vectors from XCAT files ... \n');

ct = matRad_parseMVF(ct,xcatLog,fnameXcatRoot);

fprintf('Done!\n');

%% Import DICOM structure set
fprintf('matRad: Importing DICOM structure set ... \n');

fnameRTSS = fullfile(dirDICOM,sprintf('%s_RTSS.dcm',fnameXcatRoot));
structures = matRad_importDicomRtss(fnameRTSS,ct.dicomInfo);

for i = 1:numel(structures)
    
    fprintf('Structure %d of %d.\n',i,numel(structures));
    structures(i).indices = matRad_convRtssContours2Indices(structures(i),ct);
end
fprintf('Done!\n');
cst = matRad_createCst(structures);


save(fnameXcatRoot,'ct','cst','-v7.3');

end

