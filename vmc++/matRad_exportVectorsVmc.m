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

%% check for collapsed voxels

collapsedVoxel = false(size(iVoxel_grid));

diffX = diff(ct.motionVecX{frame},1,2);
diffY = diff(ct.motionVecY{frame},1,1);
diffZ = diff(ct.motionVecZ{frame},1,3);

maskX = diffX <= 0;
maskY = diffY <= 0;
maskZ = diffZ <= 0;

collapsedVoxel(:,1:(ct.cubeDim(2)-1),:)     = maskX | collapsedVoxel(:,1:(ct.cubeDim(2)-1),:);
collapsedVoxel(:,2:ct.cubeDim(2),:)         = maskX | collapsedVoxel(:,2:ct.cubeDim(2),:);
collapsedVoxel(1:(ct.cubeDim(1)-1),:,:)     = maskY | collapsedVoxel(1:(ct.cubeDim(1)-1),:,:);
collapsedVoxel(2:ct.cubeDim(1),:,:)         = maskY | collapsedVoxel(2:ct.cubeDim(1),:,:);
collapsedVoxel(:,:,1:(ct.cubeDim(3)-1))     = maskZ | collapsedVoxel(:,:,1:(ct.cubeDim(3)-1));
collapsedVoxel(:,:,2:ct.cubeDim(3))         = maskZ | collapsedVoxel(:,:,2:ct.cubeDim(3));

% Emily's way
collapsedVoxelE = false(size(iVoxel_grid));

for i = 2:(ct.cubeDim(2)-1)
    for j = 2:(ct.cubeDim(1)-1)
        for k = 2:(ct.cubeDim(3)-1)
            
            jacobian = zeros(3,3);
            
            jacobian(1,1) = 1+(defXVoxel(j,i+1,k)-defXVoxel(j,i-1,k))./(2.*ct.resolution.x./10);
            jacobian(1,2) = (defYVoxel(j,i+1,k)-defYVoxel(j,i-1,k))./(2.*ct.resolution.y./10);
            jacobian(1,3) = (defZVoxel(j,i+1,k)-defZVoxel(j,i-1,k))./(2.*ct.resolution.z./10);
            
            jacobian(2,1) = (defXVoxel(j+1,i,k)-defXVoxel(j-1,i,k))./(2.*ct.resolution.x./10);
            jacobian(2,2) = 1+(defYVoxel(j+1,i,k)-defYVoxel(j-1,i,k))./(2.*ct.resolution.y./10);
            jacobian(2,3) = (defZVoxel(j+1,i,k)-defZVoxel(j-1,i,k))./(2.*ct.resolution.z./10);
            
            jacobian(3,1) = (defXVoxel(j,i,k+1)-defXVoxel(j,i,k-1))./(2.*ct.resolution.x./10);
            jacobian(3,2) = (defYVoxel(j,i,k+1)-defYVoxel(j,i,k-1))./(2.*ct.resolution.y./10);
            jacobian(3,3) = 1+(defZVoxel(j,i,k+1)-defZVoxel(j,i,k-1))./(2.*ct.resolution.z./10);
            
            value = det(jacobian);
            
            if value <= 0
                collapsedVoxelE(j,i,k) = true;
            end
            
        end
    end
end

collapsedVertex = false(size(iVertex_grid));

for i = 2:(numel(iVertex)-1)
    for j = 2:(numel(jVertex)-1)
        for k = 2:(numel(kVertex)-1)
            
            jacobian = zeros(3,3);
            
            jacobian(1,1) = 1+(defXVertex(j,i+1,k)-defXVertex(j,i-1,k))./(2.*ct.resolution.x./10);
            jacobian(1,2) = (defYVertex(j,i+1,k)-defYVertex(j,i-1,k))./(2.*ct.resolution.y./10);
            jacobian(1,3) = (defZVertex(j,i+1,k)-defZVertex(j,i-1,k))./(2.*ct.resolution.z./10);
            
            jacobian(2,1) = (defXVertex(j+1,i,k)-defXVertex(j-1,i,k))./(2.*ct.resolution.x./10);
            jacobian(2,2) = 1+(defYVertex(j+1,i,k)-defYVertex(j-1,i,k))./(2.*ct.resolution.y./10);
            jacobian(2,3) = (defZVertex(j+1,i,k)-defZVertex(j-1,i,k))./(2.*ct.resolution.z./10);
            
            jacobian(3,1) = (defXVertex(j,i,k+1)-defXVertex(j,i,k-1))./(2.*ct.resolution.x./10);
            jacobian(3,2) = (defYVertex(j,i,k+1)-defYVertex(j,i,k-1))./(2.*ct.resolution.y./10);
            jacobian(3,3) = 1+(defZVertex(j,i,k+1)-defZVertex(j,i,k-1))./(2.*ct.resolution.z./10);
            
            value = det(jacobian);
            
            if value <= 0
                collapsedVertex(j,i,k) = true;
            end
            
        end
    end
end

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
fprintf(fid,'%d\n',numVertices);

% now loop through and write the deformations for each vertex
for i = 1:numVertices
    fprintf(fid,'%6.6f, %6.6f, %6.6f\n',defXVertex_writeFormat(i),defYVertex_writeFormat(i),defZVertex_writeFormat(i));
end

% close file
fclose(fid);

end
