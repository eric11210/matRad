function xcatLog = matRad_xcatBin2DICOM(fnameXcatRoot)

dirXCAT = fullfile(fileparts(mfilename('fullpath')),'XCAT',filesep);
dirDICOM = fullfile(fileparts(mfilename('fullpath')),'DICOM',filesep);

fnameXcatLog = fullfile(dirXCAT,sprintf('%s_log',fnameXcatRoot));
xcatLog = matRad_xcatReadLog(fnameXcatLog);

dicomInfo = matRad_dicomHeader(xcatLog);
dicomInfo.PatientName.FamilyName = fnameXcatRoot;
dicomInfo.PatientName.GivenName = fnameXcatRoot;

numVox = xcatLog.dim.x*xcatLog.dim.y*xcatLog.dim.z;

ctList_all = ls(fullfile(dirDICOM,sprintf('%s_Phase*_Slice*.dcm',fnameXcatRoot)));

xcatLog.ctList = cell(xcatLog.dim.z,xcatLog.numPhases);


fprintf('matRad: Converting XCAT binaries to DICOM ... \n');

for phase = 1:xcatLog.numPhases
    
    if size(ctList_all,1) ~= xcatLog.numPhases*xcatLog.dim.z
        fnameXcatBin = fullfile(dirXCAT,sprintf('%s_atn_%d.bin',fnameXcatRoot,phase));
        
        fid = fopen(fnameXcatBin);
        phant = fread(fid,numVox,'single');
        phant = reshape(phant,[xcatLog.dim.x xcatLog.dim.y xcatLog.dim.z]);
        phant = permute(phant,[2 1 3]);
        fclose(fid);
        
        phant = 1000*(phant-xcatLog.muWater)./xcatLog.muWater;
        
        dicomInfo.RescaleSlope = (max(phant(:))-min(phant(:)))./(dicomInfo.LargestImagePixelValue-dicomInfo.SmallestImagePixelValue);
        dicomInfo.RescaleIntercept = min(phant(:))-dicomInfo.SmallestImagePixelValue*dicomInfo.RescaleSlope;
    end
    
    for slice = 1:xcatLog.dim.z
        
        fnameDICOM = fullfile(dirDICOM,sprintf('%s_Phase%02d_Slice%03d.dcm',fnameXcatRoot,phase,slice));
        xcatLog.ctList{slice,phase} = fnameDICOM;
        
        if size(ctList_all,1) ~= xcatLog.numPhases*xcatLog.dim.z
            
            % x is R-L, y is A-P, z is I-S
            % original CT, contours:    start @ z = 1165 mm, res = [0.85
            % 0.85 1] mm
            %dicomInfo.ImagePositionPatient = [-217; -217; slice-161];
            % new CT:                   start @ z = 1164 mm, res = [3 3 3]
            % mm
            dicomInfo.ImagePositionPatient = [-217-0.85/2+dicomInfo.PixelSpacing(1)/2; -217-0.85/2+dicomInfo.PixelSpacing(2)/2; dicomInfo.SliceThickness*(slice-1)-161];
            dicomInfo.SliceLocation = dicomInfo.SliceThickness*(slice-1)-161;
            dicomInfo.SeriesNumber = phase;
            dicomInfo.InstanceNumber = slice;
            dicomInfo.SeriesDescription = sprintf('%.1f%%',10*(phase-1));
            dicomInfo.SeriesInstanceUID = dicomuid;
            dicomInfo.ImageComments = sprintf('%.1f%%',10*(phase-1));
            
            image = phant(:,:,slice);
            
            %image=imrotate(image,-90);
            %image=fliplr(image);
            
            image = (image-dicomInfo.RescaleIntercept)./dicomInfo.RescaleSlope;
            dicomwrite(uint16(image),fnameDICOM,dicomInfo,'CreateMode','Copy');
            
            matRad_progress((phase-1)*xcatLog.dim.z+slice,xcatLog.numPhases*xcatLog.dim.z);
        end
    end
end

fprintf('Done!\n');


end

