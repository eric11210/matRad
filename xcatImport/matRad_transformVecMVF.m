function transPos = matRad_transformVecMVF(initPos,phase,xcatLog,fnameXcatRoot)

dirXCAT = fullfile(fileparts(mfilename('fullpath')),'XCAT',filesep);

fnameXcatMVF = fullfile(dirXCAT,sprintf('%s_vec_frame1_to_frame%d.txt',fnameXcatRoot,phase));

fid = fopen(fnameXcatMVF,'r');

%first line is empty
tline = fgetl(fid);
%second line contains frame number
tline = fgetl(fid);
%third and subsequent lines contain vector information
tline = fgetl(fid);

%initialize motion vectors
motionVecX = zeros(xcatLog.dim.y,xcatLog.dim.x,xcatLog.dim.z);
motionVecY = zeros(xcatLog.dim.y,xcatLog.dim.x,xcatLog.dim.z);
motionVecZ = zeros(xcatLog.dim.y,xcatLog.dim.x,xcatLog.dim.z);

%assume no motion for unspecified voxels
for i = 1:xcatLog.dim.y
    motionVecY(i,:,:) = i;
end
for i = 1:xcatLog.dim.x
    motionVecX(:,i,:) = i;
end
for i = 1:xcatLog.dim.z
    motionVecZ(:,:,i) = i;
end

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
    xCoord_vox = 1+str2double(tline((spacesInd(2)+1):(spacesInd(3)-1)));
    yCoord_vox = 1+str2double(tline((spacesInd(3)+1):(spacesInd(4)-1)));
    zCoord_vox = 1+str2double(tline((spacesInd(4)+1):(spacesInd(5)-1)));
    
    if round(xCoord_vox) == xCoord_vox && round(yCoord_vox) == yCoord_vox && round(zCoord_vox) == zCoord_vox
        %only store motion vectors for voxels which map exactly to a voxel
        %in the interpolated resolution (note that the voxel in the new
        %frame will probably not map exactly to a voxel)
        xCoordNewFrame_vox = 1+str2double(tline((spacesInd(7)+1):(spacesInd(8)-1)));
        yCoordNewFrame_vox = 1+str2double(tline((spacesInd(8)+1):(spacesInd(9)-1)));
        zCoordNewFrame_vox = 1+str2double(tline((spacesInd(9)+1):end));
        
        motionVecX(yCoord_vox,xCoord_vox,zCoord_vox) = xCoordNewFrame_vox;
        motionVecY(yCoord_vox,xCoord_vox,zCoord_vox) = yCoordNewFrame_vox;
        motionVecZ(yCoord_vox,xCoord_vox,zCoord_vox) = zCoordNewFrame_vox;
        
        matRad_progress(i,numel(motionVecX));
        i = i+1;
        tline = fgetl(fid);
    else
        tline = fgetl(fid);
        continue
    end
    
end

fclose(fid);

transPos = zeros(1,3);

transPos(1) = interp3(motionVecX,initPos(1),initPos(2),initPos(3));
transPos(2) = interp3(motionVecY,initPos(1),initPos(2),initPos(3));
transPos(3) = interp3(motionVecZ,initPos(1),initPos(2),initPos(3));


end