function xcatLog = matRad_xcatBin2DICOM(importOptions)

dirXCAT = fullfile(fileparts(mfilename('fullpath')),'XCAT',filesep);
dirDICOM = fullfile(fileparts(mfilename('fullpath')),'DICOM',filesep);

fnameXcatLog = fullfile(dirXCAT,sprintf('%s_log',importOptions.fnameXcatRoot));
xcatLog = matRad_xcatReadLog(fnameXcatLog);

dicomInfo = matRad_dicomHeader(xcatLog);
dicomInfo.PatientName.FamilyName = importOptions.fnameXcatRoot;
dicomInfo.PatientName.GivenName = importOptions.fnameXcatRoot;

numVox = xcatLog.dim.x*xcatLog.dim.y*xcatLog.dim.z;

ctList_all = ls(fullfile(dirDICOM,sprintf('%s_Phase*_Slice*.dcm',importOptions.fnameXcatRoot)));

xcatLog.ctList = cell(xcatLog.dim.z,xcatLog.numFrames);


% first determine the padding in the x- and y-direction
% we need to loop over all phases to ensure that they all get to be the
% same size in the end
zeroIndPreX = xcatLog.dim.x;
zeroIndPostX = 0;
zeroIndPreY = xcatLog.dim.y;
zeroIndPostY = 0;

fprintf('matRad: Determing padding in CT images ... \n');
for phase = 1:xcatLog.numFrames
        
        if importOptions.massConserve
            fnameXcatBin = fullfile(dirXCAT,sprintf('%s_atnMC_%d.bin',importOptions.fnameXcatRoot,phase));
        else
            fnameXcatBin = fullfile(dirXCAT,sprintf('%s_atn_%d.bin',importOptions.fnameXcatRoot,phase));
        end
        
        fid = fopen(fnameXcatBin);
        phantBody = fread(fid,numVox,'single');
        phantBody = reshape(phantBody,[xcatLog.dim.x xcatLog.dim.y xcatLog.dim.z]);
        phantBody = permute(phantBody,[2 1 3]);
        fclose(fid);
        
        if importOptions.addTumour
            
            if importOptions.massConserve
                fnameXcatBin = fullfile(dirXCAT,sprintf('%s_tumour_atnMC_%d.bin',importOptions.fnameXcatRoot,phase));
            else
                fnameXcatBin = fullfile(dirXCAT,sprintf('%s_tumour_atn_%d.bin',importOptions.fnameXcatRoot,phase));
            end
            
            fid = fopen(fnameXcatBin);
            phantTumour = fread(fid,numVox,'single');
            phantTumour = reshape(phantTumour,[xcatLog.dim.x xcatLog.dim.y xcatLog.dim.z]);
            phantTumour = permute(phantTumour,[2 1 3]);
            fclose(fid);
            
            phant = phantBody+phantTumour;
        else
            
            phant = phantBody;
        end
        
        % find padding in x- and y-direction
        sumPhantZ = sum(phant,3);
        sumPhantXZ = sum(sumPhantZ,2);
        sumPhantYZ = sum(sumPhantZ,1);
        zeroIndPreX = min(find(sumPhantYZ ~= 0,1,'first')-1,zeroIndPreX);
        zeroIndPostX = max(find(sumPhantYZ ~= 0,1,'last')+1,zeroIndPostX);
        zeroIndPreY = min(find(sumPhantXZ ~= 0,1,'first')-1,zeroIndPreY);
        zeroIndPostY = max(find(sumPhantXZ ~= 0,1,'last')+1,zeroIndPostY);
    
    matRad_progress(phase,xcatLog.numFrames);
end

% save the indices in the xcatLog
xcatLog.zeroIndPreX = zeroIndPreX;
xcatLog.zeroIndPostX = zeroIndPostX;
xcatLog.zeroIndPreY = zeroIndPreY;
xcatLog.zeroIndPostY = zeroIndPostY;

%
dicomInfo.Rows = zeroIndPostY-zeroIndPreY-1;
dicomInfo.Columns = zeroIndPostX-zeroIndPreX-1;
dicomInfo.Width = zeroIndPostX-zeroIndPreX-1;
dicomInfo.Height = zeroIndPostY-zeroIndPreY-1;

fprintf('matRad: Converting XCAT binaries to DICOM ... \n');
for phase = 1:xcatLog.numFrames
    
    if size(ctList_all,1) ~= xcatLog.numFrames*xcatLog.dim.z
        
        if importOptions.massConserve
            fnameXcatBin = fullfile(dirXCAT,sprintf('%s_atnMC_%d.bin',importOptions.fnameXcatRoot,phase));
        else
            fnameXcatBin = fullfile(dirXCAT,sprintf('%s_atn_%d.bin',importOptions.fnameXcatRoot,phase));
        end
        
        fid = fopen(fnameXcatBin);
        phantBody = fread(fid,numVox,'single');
        phantBody = reshape(phantBody,[xcatLog.dim.x xcatLog.dim.y xcatLog.dim.z]);
        phantBody = permute(phantBody,[2 1 3]);
        fclose(fid);
        
        if importOptions.addTumour
            
            if importOptions.massConserve
                fnameXcatBin = fullfile(dirXCAT,sprintf('%s_tumour_atnMC_%d.bin',importOptions.fnameXcatRoot,phase));
            else
                fnameXcatBin = fullfile(dirXCAT,sprintf('%s_tumour_atn_%d.bin',importOptions.fnameXcatRoot,phase));
            end
            
            fid = fopen(fnameXcatBin);
            phantTumour = fread(fid,numVox,'single');
            phantTumour = reshape(phantTumour,[xcatLog.dim.x xcatLog.dim.y xcatLog.dim.z]);
            phantTumour = permute(phantTumour,[2 1 3]);
            fclose(fid);
            
            phant = phantBody+phantTumour;
        else
            
            phant = phantBody;
        end
        
        % remove extra padding
        phant(:,zeroIndPostX:end,:) = [];
        phant(:,1:zeroIndPreX,:)    = [];
        phant(zeroIndPostY:end,:,:) = [];
        phant(1:zeroIndPreY,:,:) = [];
        
        % convert density to HU
        phant = 1000*(phant-xcatLog.muWater)./xcatLog.muWater;
        
        % use entire range of int
        dicomInfo.RescaleSlope = (max(phant(:))-min(phant(:)))./(dicomInfo.LargestImagePixelValue-dicomInfo.SmallestImagePixelValue);
        dicomInfo.RescaleIntercept = min(phant(:))-dicomInfo.SmallestImagePixelValue*dicomInfo.RescaleSlope;
    end
    
    for slice = 1:xcatLog.dim.z
        
        fnameDICOM = fullfile(dirDICOM,sprintf('%s_Phase%02d_Slice%03d.dcm',importOptions.fnameXcatRoot,phase,slice));
        xcatLog.ctList{slice,phase} = fnameDICOM;
        if size(ctList_all,1) ~= xcatLog.numFrames*xcatLog.dim.z
            % x is R-L, y is A-P, z is I-S
            % original CT, contours:    start @ z = 1165 mm, res = [0.85
            % 0.85 1] mm
            %dicomInfo.ImagePositionPatient = [-217; -217; slice-161];
            % new CT:                   start @ z = 1164 mm, res = [3 3 3]
            % mm
            % also account for removal of the padding in the CT image
            dicomInfo.ImagePositionPatient = [-217-0.85/2+dicomInfo.PixelSpacing(1)*(1/2+zeroIndPreX); -217-0.85/2+dicomInfo.PixelSpacing(2)*(1/2+zeroIndPreY); dicomInfo.SliceThickness*(slice-1)-161];
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
            
            matRad_progress((phase-1)*xcatLog.dim.z+slice,xcatLog.numFrames*xcatLog.dim.z);
        end
    end
end

%
xcatLog.dim.x = dicomInfo.Width;
xcatLog.dim.y = dicomInfo.Height;

fprintf('Done!\n');



end

