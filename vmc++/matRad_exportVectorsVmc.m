function matRad_exportVectorsVmc(ct,frame,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad text vectors export for vmc++
% 
% call
%   matRad_exportVectorsVmc(ct,filename)
%
% input
%   ct:             matRad ct struct
%   filename:       path where vectorsFile is created
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the vectors file is a text file which contains the x-, y-, and
% z-deformations of the vertices of each voxel (in cm)

% the ordering of vertices is the same as it is for the CT


%% interpolation coordinates

% determine the coordinates at the voxel centres
iVoxel = 1:ct.cubeDim(2);
jVoxel = 1:ct.cubeDim(1);
kVoxel = 1:ct.cubeDim(3);
% create meshgrid
[iVoxel_grid,jVoxel_grid,kVoxel_grid] = meshgrid(iVoxel,jVoxel,kVoxel);

% determine the coordinates at the voxel vertices
iVertex = 0.5:1:(ct.cubeDim(2)+0.5);
jVertex = 0.5:1:(ct.cubeDim(1)+0.5);
kVertex = 0.5:1:(ct.cubeDim(3)+0.5);
% create meshgrid
[iVertex_grid,jVertex_grid,kVertex_grid] = meshgrid(iVertex,jVertex,kVertex);

% determine total number of vertices
numVertices = numel(iVertex_grid);

%% get deformation at voxel centres and vertices

% determine the deformation (in cm) at the voxel centres
defXVoxel = (ct.motionVecX{frame}-iVoxel_grid).*ct.resolution.x/10;
defYVoxel = (ct.motionVecY{frame}-jVoxel_grid).*ct.resolution.y/10;
defZVoxel = (ct.motionVecZ{frame}-kVoxel_grid).*ct.resolution.z/10;

% determine the deformation at the voxel vertices
defXVertex = interp3(iVoxel_grid,jVoxel_grid,kVoxel_grid,defXVoxel,iVertex_grid,jVertex_grid,kVertex_grid,'linear',0);
defYVertex = interp3(iVoxel_grid,jVoxel_grid,kVoxel_grid,defYVoxel,iVertex_grid,jVertex_grid,kVertex_grid,'linear',0);
defZVertex = interp3(iVoxel_grid,jVoxel_grid,kVoxel_grid,defZVoxel,iVertex_grid,jVertex_grid,kVertex_grid,'linear',0);

%% writing prep

% get deformations into the correct format/order
defXVertex_writeFormat = permute(defXVertex,[2 1 3]);
defXVertex_writeFormat = reshape(defXVertex_writeFormat,[],1);
defYVertex_writeFormat = permute(defYVertex,[2 1 3]);
defYVertex_writeFormat = reshape(defYVertex_writeFormat,[],1);
defZVertex_writeFormat = permute(defZVertex,[2 1 3]);
defZVertex_writeFormat = reshape(defZVertex_writeFormat,[],1);

%% writing

% open file
fid = fopen(filename,'w');

% first write the number of vertices
%fprintf(fid,'%d\n',numVertices);
fwrite(fid,numVertices,'int32');

% now loop through and write the deformations for each vertex
for i = 1:numVertices
    %fprintf(fid,'%6.6f, %6.6f, %6.6f\n',defXVertex_writeFormat(i),defYVertex_writeFormat(i),defZVertex_writeFormat(i));
    fwrite(fid,defXVertex_writeFormat(i),'float32');
    fwrite(fid,defYVertex_writeFormat(i),'float32');
    fwrite(fid,defZVertex_writeFormat(i),'float32');
end

% close file
fclose(fid);

end
