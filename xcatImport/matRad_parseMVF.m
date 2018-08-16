function ct = matRad_parseMVF(ct,xcatLog,fnameXcatRoot)

dirXCAT = fullfile(fileparts(mfilename('fullpath')),'XCAT',filesep);

fnameXcatPar = fullfile(dirXCAT,sprintf('%s.par',fnameXcatRoot));
[tumourPosInit,tumourRadius] = matRad_xcatReadPar(fnameXcatPar);
tumourPosInit = tumourPosInit-[xcatLog.zeroIndPreX xcatLog.zeroIndPreY 0];

tumourPos = zeros(xcatLog.numFrames+1,3);
tumourPos(1,:) = tumourPosInit;
tumourPos(xcatLog.numFrames+1,:) = tumourPos(1,:);
t = xcatLog.deltaT.*((1:(xcatLog.numFrames+1))-1);

ct.motionVecX = cell(xcatLog.numFrames,1);
ct.motionVecY = cell(xcatLog.numFrames,1);
ct.motionVecZ = cell(xcatLog.numFrames,1);

for phase = 1:xcatLog.numFrames
    
    fprintf('\nPhase %d of %d.\n',phase,xcatLog.numFrames);
    
    %initialize motion vectors
    ct.motionVecX{phase} = zeros(ct.cubeDim);
    ct.motionVecY{phase} = zeros(ct.cubeDim);
    ct.motionVecZ{phase} = zeros(ct.cubeDim);
    
    %assume no motion for unspecified voxels
    for i = 1:ct.cubeDim(1)
        ct.motionVecY{phase}(i,:,:) = i;
    end
    for i = 1:ct.cubeDim(2)
        ct.motionVecX{phase}(:,i,:) = i;
    end
    for i = 1:ct.cubeDim(3)
        ct.motionVecZ{phase}(:,:,i) = i;
    end
    
    if phase == 1
        continue
    end
    
    resRatioX = ct.resolution.x./xcatLog.resolution.x;
    resRatioY = ct.resolution.y./xcatLog.resolution.y;
    resRatioZ = ct.resolution.z./xcatLog.resolution.z;
    
    if ~( round(resRatioX) == resRatioX && round(resRatioY) == resRatioY && round(resRatioZ) == resRatioZ )
        error('Non-integral ratios of resolutions.  Inexact motion vector mapping may occur.');
    end
    
    fnameXcatMVF = fullfile(dirXCAT,sprintf('%s_vec_frame1_to_frame%d.txt',fnameXcatRoot,phase));
    
    fid = fopen(fnameXcatMVF,'r');
    
    %first line is empty
    tline = fgetl(fid);
    %second line contains frame number
    tline = fgetl(fid);
    %third and subsequent lines contain vector information
    tline = fgetl(fid);
    
    i = 1;
    
    while ischar(tline)
        
        %{
    if ~strcmp(tline(1:5),'known')
        error('Unknown vector encountered');
    end
        %}
        
        spacesInd = strfind(tline,' ');
        
        % each line contains the voxel coordinate (in the original resolution)
        % of each transformed voxel, and also the voxel sub-coordinates in the
        % transformed frame
        
        % also correct for the removal of the padding
        xCoord_vox = 1+str2double(tline((spacesInd(2)+1):(spacesInd(3)-1)))./resRatioX-xcatLog.zeroIndPreX;
        yCoord_vox = 1+str2double(tline((spacesInd(3)+1):(spacesInd(4)-1)))./resRatioY-xcatLog.zeroIndPreY;
        zCoord_vox = 1+str2double(tline((spacesInd(4)+1):(spacesInd(5)-1)))./resRatioZ;
        
        if round(xCoord_vox) == xCoord_vox && round(yCoord_vox) == yCoord_vox && round(zCoord_vox) == zCoord_vox ...
                && 1 <= xCoord_vox && xCoord_vox <= ct.cubeDim(2) && 1<= yCoord_vox && yCoord_vox <= ct.cubeDim(1)
            %only store motion vectors for voxels which map exactly to a voxel
            %in the interpolated resolution (note that the voxel in the new
            %frame will probably not map exactly to a voxel)
            xCoordNewFrame_vox = 1+str2double(tline((spacesInd(7)+1):(spacesInd(8)-1)))./resRatioX-xcatLog.zeroIndPreX;
            yCoordNewFrame_vox = 1+str2double(tline((spacesInd(8)+1):(spacesInd(9)-1)))./resRatioY-xcatLog.zeroIndPreY;
            zCoordNewFrame_vox = 1+str2double(tline((spacesInd(9)+1):end))./resRatioZ;
            
            ct.motionVecX{phase}(yCoord_vox,xCoord_vox,zCoord_vox) = xCoordNewFrame_vox;
            ct.motionVecY{phase}(yCoord_vox,xCoord_vox,zCoord_vox) = yCoordNewFrame_vox;
            ct.motionVecZ{phase}(yCoord_vox,xCoord_vox,zCoord_vox) = zCoordNewFrame_vox;
            
            matRad_progress(i,numel(ct.motionVecX{phase}));
            i = i+1;
            tline = fgetl(fid);
        else
            tline = fgetl(fid);
            continue
        end
        
    end
    
    tumourPos(phase,1) = interp3(ct.motionVecX{phase},tumourPos(1,1),tumourPos(1,2),tumourPos(1,3));
    tumourPos(phase,2) = interp3(ct.motionVecY{phase},tumourPos(1,1),tumourPos(1,2),tumourPos(1,3));
    tumourPos(phase,3) = interp3(ct.motionVecZ{phase},tumourPos(1,1),tumourPos(1,2),tumourPos(1,3));
    
    fclose(fid);
    
end

ct.tumourMotion.coordsVox   = tumourPos;
ct.tumourMotion.t           = t;
ct.tumourMotion.radius      = tumourRadius;

end