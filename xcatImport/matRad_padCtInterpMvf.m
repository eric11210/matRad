function ct = matRad_padCtInterpMvf(ct,xcatLog,importOptions)

% we want to pad ct for two reasons (both related to def_tetra in vmc++):
% 1) we need a 2-voxel thick pad on all edges for def_tetra to work
% 2) we also cannot have any voxels/tetrahedrons collapse

% at the same time, interpolate all voxels with 0 density (air) so that
% voxels/tetrahedrons don't collapse

% we also wish to remove instances of voxel collapse that happen in XCAT

%% get padding information

% find max extent of deformed CT
xIndPre     = 1;
yIndPre     = 1;
zIndPre     = 1;
xIndPost    = ct.cubeDim(2);
yIndPost    = ct.cubeDim(1);
zIndPost    = ct.cubeDim(3);

% loop over all phases
for frame = 1:xcatLog.numFrames
    
    % for each frame, find max extent of deformed CT
    xIndPre_frame = min(ct.motionVecX{frame}(:));
    yIndPre_frame = min(ct.motionVecY{frame}(:));
    zIndPre_frame = min(ct.motionVecZ{frame}(:));
    
    xIndPost_frame = max(ct.motionVecX{frame}(:));
    yIndPost_frame = max(ct.motionVecY{frame}(:));
    zIndPost_frame = max(ct.motionVecZ{frame}(:));
    
    % update max
    xIndPre = min([xIndPre xIndPre_frame]);
    yIndPre = min([yIndPre yIndPre_frame]);
    zIndPre = min([zIndPre zIndPre_frame]);
    
    xIndPost = max([xIndPost xIndPost_frame]);
    yIndPost = max([yIndPost yIndPost_frame]);
    zIndPost = max([zIndPost zIndPost_frame]);
    
end

% now round down for the pre, up for the post
% also subtract/add two, since we need a 2-voxel thick layer
xIndPre = floor(xIndPre)-2;
yIndPre = floor(yIndPre)-2;
zIndPre = floor(zIndPre)-2;

xIndPost = ceil(xIndPost)+2;
yIndPost = ceil(yIndPost)+2;
zIndPost = ceil(zIndPost)+2;

% now determine by how much we need to pad each direction
xPadPre = 1-xIndPre;
yPadPre = 1-yIndPre;
zPadPre = 1-zIndPre;

xPadPost = xIndPost-ct.cubeDim(2);
yPadPost = yIndPost-ct.cubeDim(1);
zPadPost = zIndPost-ct.cubeDim(3);

%% do the padding and interpolation

% first fix all the metadata
ct.cubeDim(2) = ct.cubeDim(2)+xPadPre+xPadPost;
ct.cubeDim(1) = ct.cubeDim(1)+yPadPre+yPadPost;
ct.cubeDim(3) = ct.cubeDim(3)+zPadPre+zPadPost;

xPosPre = ct.x(1)-ct.resolution.x.*xPadPre;
yPosPre = ct.y(1)-ct.resolution.y.*yPadPre;
zPosPre = ct.z(1)-ct.resolution.z.*zPadPre;

xPosPost = ct.x(end)+ct.resolution.x.*xPadPost;
yPosPost = ct.y(end)+ct.resolution.y.*yPadPost;
zPosPost = ct.z(end)+ct.resolution.z.*zPadPost;

ct.x = xPosPre:ct.resolution.x:xPosPost;
ct.y = yPosPre:ct.resolution.y:yPosPost;
ct.z = zPosPre:ct.resolution.z:zPosPost;

ct.tumourMotion.coordsVox(:,1) = ct.tumourMotion.coordsVox(:,1)+xPadPre;
ct.tumourMotion.coordsVox(:,2) = ct.tumourMotion.coordsVox(:,2)+yPadPre;
ct.tumourMotion.coordsVox(:,3) = ct.tumourMotion.coordsVox(:,3)+zPadPre;

% prepare for the padding by defining base motionVecs (i.e. what they would
% look like with no padding)
[baseMotionVecX,baseMotionVecY,baseMotionVecZ] = meshgrid(1:ct.cubeDim(2),1:ct.cubeDim(1),1:ct.cubeDim(3));

% determine all voxels not part of the 2-voxel thick pad
notPad = false(ct.cubeDim);
notPad(3:(ct.cubeDim(1)-2),3:(ct.cubeDim(2)-2),3:(ct.cubeDim(3)-2)) = true;

if importOptions.vmcDef
    % if we're doing vmc deformation code, we only need/have the first
    % phase
    ct.cube{1} = padarray(ct.cube{1},[yPadPre xPadPre zPadPre],0,'pre');
    ct.cube{1} = padarray(ct.cube{1},[yPadPost xPadPost zPadPost],0,'post');
    
    ct.cubeHU{1} = padarray(ct.cubeHU{1},[yPadPre xPadPre zPadPre],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'pre');
    ct.cubeHU{1} = padarray(ct.cubeHU{1},[yPadPost xPadPost zPadPost],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'post');
end

% now pad the motionVecs
for frame = 1:xcatLog.numFrames
    
    fprintf('Frame %d of %d... ',frame,xcatLog.numFrames);
    
    if ~importOptions.vmcDef
        ct.cube{frame} = padarray(ct.cube{frame},[yPadPre xPadPre zPadPre],0,'pre');
        ct.cube{frame} = padarray(ct.cube{frame},[yPadPost xPadPost zPadPost],0,'post');
        
        ct.cubeHU{frame} = padarray(ct.cubeHU{frame},[yPadPre xPadPre zPadPre],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'pre');
        ct.cubeHU{frame} = padarray(ct.cubeHU{frame},[yPadPost xPadPost zPadPost],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'post');
    end
    
    ct.motionVecX{frame} = padarray(ct.motionVecX{frame},[yPadPre xPadPre zPadPre],0,'pre');
    ct.motionVecX{frame} = padarray(ct.motionVecX{frame},[yPadPost xPadPost zPadPost],0,'post');
    
    ct.motionVecY{frame} = padarray(ct.motionVecY{frame},[yPadPre xPadPre zPadPre],0,'pre');
    ct.motionVecY{frame} = padarray(ct.motionVecY{frame},[yPadPost xPadPost zPadPost],0,'post');
    
    ct.motionVecZ{frame} = padarray(ct.motionVecZ{frame},[yPadPre xPadPre zPadPre],0,'pre');
    ct.motionVecZ{frame} = padarray(ct.motionVecZ{frame},[yPadPost xPadPost zPadPost],0,'post');
    
    % shift all motionVecs to the new coordinate system
    ct.motionVecX{frame} = ct.motionVecX{frame}+xPadPre;
    ct.motionVecY{frame} = ct.motionVecY{frame}+yPadPre;
    ct.motionVecZ{frame} = ct.motionVecZ{frame}+zPadPre;
    
    % give correct values for the zeros we just introduced
    newVoxels = ct.motionVecX{frame} == xPadPre; % motionVecX is only xPadPre where we just padded it
    
    ct.motionVecX{frame}(newVoxels) = baseMotionVecX(newVoxels);
    ct.motionVecY{frame}(newVoxels) = baseMotionVecY(newVoxels);
    ct.motionVecZ{frame}(newVoxels) = baseMotionVecZ(newVoxels);
    
    
    %% now onto interpolation
    
    % now determine which voxels need to be interpolated
    % these are all voxels which have zero density EXCEPT for the
    % 2-voxel thick pad
    
    % start out with no voxels
    doInterp = false(ct.cubeDim);
    
    % now set to true all voxels which have a zero density
    doInterp(ct.cube{1} == 0) = true;
    
    % now add the skin, two voxels deep
    doInterp = addSkin(doInterp,2);
    
    % now set to false all voxels which have a non-zero deformation
    % i.e., their motionVec does not equal the baseMotionVecX
    doInterp(ct.motionVecX{frame} ~= baseMotionVecX) = false;
    doInterp(ct.motionVecY{frame} ~= baseMotionVecY) = false;
    doInterp(ct.motionVecZ{frame} ~= baseMotionVecZ) = false;
    
    % now take out all voxels in the pad
    doInterp(~notPad) = false;
    
    % interpolate motion vectors at the doInterp points
    [motionVecX_interp, motionVecY_interp, motionVecZ_interp] = interpSomeVoxels(ct.motionVecX{frame}, ct.motionVecY{frame}, ct.motionVecZ{frame}, baseMotionVecX, baseMotionVecY, baseMotionVecZ, doInterp);
    
    % put interpolated vectors into ct struct
    ct.motionVecX{frame} = motionVecX_interp;
    ct.motionVecY{frame} = motionVecY_interp;
    ct.motionVecZ{frame} = motionVecZ_interp;
    
    % fix the 2-voxel thick pad
    ct.motionVecX{frame}(~notPad) = baseMotionVecX(~notPad);
    ct.motionVecY{frame}(~notPad) = baseMotionVecY(~notPad);
    ct.motionVecZ{frame}(~notPad) = baseMotionVecZ(~notPad);
    
    %% correct collapsed voxels
    
    % keep old motionVecs before fixing collapse
    ct.origMotionVecX{frame} = ct.motionVecX{frame};
    ct.origMotionVecY{frame} = ct.motionVecY{frame};
    ct.origMotionVecZ{frame} = ct.motionVecZ{frame};
    
    % find collapsed voxels using adjacency test
    collapsedVoxel = findCollapsedVoxels(ct.motionVecX{frame}, ct.motionVecY{frame}, ct.motionVecZ{frame}, 'adjacent');
    
    % now that we know which voxels have collapsed, we can interpolate the
    % deformation vectors from those which have not collapsed
    
    % interpolate motion vectors at the doInterp points
    [motionVecX_interp, motionVecY_interp, motionVecZ_interp] = interpSomeVoxels(ct.motionVecX{frame}, ct.motionVecY{frame}, ct.motionVecZ{frame}, baseMotionVecX, baseMotionVecY, baseMotionVecZ, collapsedVoxel);
    
    % put interpolated vectors into ct struct
    ct.motionVecX{frame} = motionVecX_interp;
    ct.motionVecY{frame} = motionVecY_interp;
    ct.motionVecZ{frame} = motionVecZ_interp;
    
    % check if any collapsed voxels persist using tetrahedral volume test
    collapsedVoxel = findCollapsedVoxels(ct.motionVecX{frame}, ct.motionVecY{frame}, ct.motionVecZ{frame}, 'tetra', ct.resolution);
    
    % initialize completion flag
    voxelsAreCollapsed =  any(collapsedVoxel(:));
    
    % initialize loop variables
    iteration = 0;
    motionVecX_smooth = ct.motionVecX{frame};
    motionVecY_smooth = ct.motionVecY{frame};
    motionVecZ_smooth = ct.motionVecZ{frame};
    
    % remove any further instances of voxel collapse by smoothing out the
    % vectors
    % do while loop, increasing range each time
    while voxelsAreCollapsed
        
        % update iteration, range
        iteration = iteration+1;
        width = iteration*2+1;
        
        % do average filter
        % add in the pad to fix issues with convolution (MATLAB will 0-pad
        % anyway)
        
        % subtract off coordinates first
        
        kernel = ones(width,width,width)./(width^3);
        motionVecX_smooth = convn(padarray(ct.motionVecX{frame},[iteration iteration iteration],'replicate','both'),kernel,'same');
        motionVecY_smooth = convn(padarray(ct.motionVecY{frame},[iteration iteration iteration],'replicate','both'),kernel,'same');
        motionVecZ_smooth = convn(padarray(ct.motionVecZ{frame},[iteration iteration iteration],'replicate','both'),kernel,'same');
        
        % remove the extra padding
        motionVecX_smooth = motionVecX_smooth((iteration+1):(end-iteration),(iteration+1):(end-iteration),(iteration+1):(end-iteration));
        motionVecY_smooth = motionVecY_smooth((iteration+1):(end-iteration),(iteration+1):(end-iteration),(iteration+1):(end-iteration));
        motionVecZ_smooth = motionVecZ_smooth((iteration+1):(end-iteration),(iteration+1):(end-iteration),(iteration+1):(end-iteration));
        
        % make sure the pad is not changed from the 0-motion
        motionVecX_smooth(~notPad) = ct.motionVecX{frame}(~notPad);
        motionVecY_smooth(~notPad) = ct.motionVecY{frame}(~notPad);
        motionVecZ_smooth(~notPad) = ct.motionVecZ{frame}(~notPad);
        
        % check if any collapsed voxels persist using Jacobian test
        collapsedVoxel = findCollapsedVoxels(motionVecX_smooth, motionVecY_smooth, motionVecZ_smooth, 'tetra', ct.resolution);
        
        % update completion flag
        voxelsAreCollapsed =  any(collapsedVoxel(:));
        
    end
    
    % update motionVecs in ct
    ct.motionVecX{frame} = motionVecX_smooth;
    ct.motionVecY{frame} = motionVecY_smooth;
    ct.motionVecZ{frame} = motionVecZ_smooth;
    
    fprintf('No collapse after %d iterations.\n',iteration);
end

end

function skinMask = addSkin(mask,skinThickness)

thinEdge    = false(size(mask));
thickEdge   = false(size(mask));

for k = 1:size(mask,3)
    
    thinEdge(:,:,k) = edge(mask(:,:,k));
    
    for i = 1:size(mask,2)
        for j = 1:size(mask,1)
            
            if thinEdge(j,i,k)
                
                j0 = max([j-skinThickness 1]);
                j1 = min([j+skinThickness size(mask,1)]);
                i0 = max([i-skinThickness 1]);
                i1 = min([i+skinThickness size(mask,2)]);
                k0 = max([k-skinThickness 1]);
                k1 = min([k+skinThickness size(mask,3)]);
                
                thickEdge(j0:j1,i0:i1,k0:k1) = true;
            end
            
        end
    end
    
end

skinMask = thickEdge | mask;

end

function collapsedVoxel = findCollapsedVoxels(motionVecX, motionVecY, motionVecZ, test, resolution)

% find cube dimensions
cubeDim = size(motionVecX);

% define collapsed voxel Boolean array
collapsedVoxel = false(cubeDim);

switch test
    case 'adjacent'
        % voxels are collapsed if the motion vec of a particular component is
        % not strictly monotonically increasing in that particular direction
        
        % first find the cumulative max for each component in the same
        % direction
        cumMaxX = cummax(motionVecX,2);
        cumMaxY = cummax(motionVecY,1);
        cumMaxZ = cummax(motionVecZ,3);
        
        % the cumulative max should be strictly monotonically increasing. if it
        % is not, then there *might* be voxel collapse (this will generate
        % false positives but it should catch all instances of voxel collapse)
        % the cummax is also useful because it accounts for voxel collapse
        % which is "> 1 voxel deep". i.e., it can detect problems where a voxel
        % might move past several voxels (which do not collapse between
        % themselves)
        
        diffX = diff(cumMaxX,1,2);
        diffY = diff(cumMaxY,1,1);
        diffZ = diff(cumMaxZ,1,3);
        
        maskX = diffX <= 0;
        maskY = diffY <= 0;
        maskZ = diffZ <= 0;
        
        collapsedVoxel(:,1:(cubeDim(2)-1),:)     = maskX | collapsedVoxel(:,1:(cubeDim(2)-1),:);
        collapsedVoxel(:,2:cubeDim(2),:)         = maskX | collapsedVoxel(:,2:cubeDim(2),:);
        collapsedVoxel(1:(cubeDim(1)-1),:,:)     = maskY | collapsedVoxel(1:(cubeDim(1)-1),:,:);
        collapsedVoxel(2:cubeDim(1),:,:)         = maskY | collapsedVoxel(2:cubeDim(1),:,:);
        collapsedVoxel(:,:,1:(cubeDim(3)-1))     = maskZ | collapsedVoxel(:,:,1:(cubeDim(3)-1));
        collapsedVoxel(:,:,2:cubeDim(3))         = maskZ | collapsedVoxel(:,:,2:cubeDim(3));
        
    case 'jacobian'
        
        % determine the coordinates at the voxel centres
        iVoxel = 1:cubeDim(2);
        jVoxel = 1:cubeDim(1);
        kVoxel = 1:cubeDim(3);
        % create meshgrid
        [iVoxel_grid,jVoxel_grid,kVoxel_grid] = meshgrid(iVoxel,jVoxel,kVoxel);
        
        % determine the deformation (in cm) at the voxel centres
        defXVoxel = (motionVecX-iVoxel_grid).*resolution.x/10;
        defYVoxel = (motionVecY-jVoxel_grid).*resolution.y/10;
        defZVoxel = (motionVecZ-kVoxel_grid).*resolution.z/10;
        
        for i = 2:(cubeDim(2)-1)
            for j = 2:(cubeDim(1)-1)
                for k = 2:(cubeDim(3)-1)
                    
                    jacobian = zeros(3,3);
                    
                    jacobian(1,1) = 1+(defXVoxel(j,i+1,k)-defXVoxel(j,i-1,k))./(2.*resolution.x./10);
                    jacobian(1,2) = (defYVoxel(j,i+1,k)-defYVoxel(j,i-1,k))./(2.*resolution.y./10);
                    jacobian(1,3) = (defZVoxel(j,i+1,k)-defZVoxel(j,i-1,k))./(2.*resolution.z./10);
                    
                    jacobian(2,1) = (defXVoxel(j+1,i,k)-defXVoxel(j-1,i,k))./(2.*resolution.x./10);
                    jacobian(2,2) = 1+(defYVoxel(j+1,i,k)-defYVoxel(j-1,i,k))./(2.*resolution.y./10);
                    jacobian(2,3) = (defZVoxel(j+1,i,k)-defZVoxel(j-1,i,k))./(2.*resolution.z./10);
                    
                    jacobian(3,1) = (defXVoxel(j,i,k+1)-defXVoxel(j,i,k-1))./(2.*resolution.x./10);
                    jacobian(3,2) = (defYVoxel(j,i,k+1)-defYVoxel(j,i,k-1))./(2.*resolution.y./10);
                    jacobian(3,3) = 1+(defZVoxel(j,i,k+1)-defZVoxel(j,i,k-1))./(2.*resolution.z./10);
                    
                    value = det(jacobian);
                    
                    if value <= 0
                        collapsedVoxel(j,i,k) = true;
                    end
                    
                end
            end
        end
        
    case 'tetra'
        
        % determine the coordinates at the voxel centres
        iVoxel = 1:cubeDim(2);
        jVoxel = 1:cubeDim(1);
        kVoxel = 1:cubeDim(3);
        % create meshgrid
        [iVoxel_grid,jVoxel_grid,kVoxel_grid] = meshgrid(iVoxel,jVoxel,kVoxel);
        
        % determine the coordinates at the voxel vertices
        iVertex = 0.5:1:(cubeDim(2)+0.5);
        jVertex = 0.5:1:(cubeDim(1)+0.5);
        kVertex = 0.5:1:(cubeDim(3)+0.5);
        % create meshgrid
        [iVertex_grid,jVertex_grid,kVertex_grid] = meshgrid(iVertex,jVertex,kVertex);
        
        % determine the deformation (in cm) at the voxel centres
        defXVoxel = (motionVecX-iVoxel_grid).*resolution.x/10;
        defYVoxel = (motionVecY-jVoxel_grid).*resolution.y/10;
        defZVoxel = (motionVecZ-kVoxel_grid).*resolution.z/10;
        
        % determine the deformation at the voxel vertices
        defXVertex = interp3(iVoxel_grid,jVoxel_grid,kVoxel_grid,defXVoxel,iVertex_grid,jVertex_grid,kVertex_grid,'linear',0);
        defYVertex = interp3(iVoxel_grid,jVoxel_grid,kVoxel_grid,defYVoxel,iVertex_grid,jVertex_grid,kVertex_grid,'linear',0);
        defZVertex = interp3(iVoxel_grid,jVoxel_grid,kVoxel_grid,defZVoxel,iVertex_grid,jVertex_grid,kVertex_grid,'linear',0);
        
        % determine the motionVec at the voxel vertices
        motionVecX_vertex = defXVertex.*10./resolution.x + iVertex_grid;
        motionVecY_vertex = defYVertex.*10./resolution.x + jVertex_grid;
        motionVecZ_vertex = defZVertex.*10./resolution.x + kVertex_grid;
        
        % tetrahedron offsets
        iTetra = [  0 0 1 0;
                    0 0 1 1;
                    0 0 0 1;
                    0 1 1 1;
                    0 1 1 0;
                    1 1 1 0];
        
        jTetra = [  0 1 1 1;
                    0 1 1 1;
                    0 0 1 1;
                    0 0 1 1;
                    0 1 0 0;
                    0 0 1 0];
        
        kTetra = [  0 0 0 1;
                    0 1 0 1;
                    0 1 1 1;
                    0 0 1 0;
                    0 1 0 1;
                    1 0 1 1];
        
        for i = 1:cubeDim(2)
            for j = 1:cubeDim(1)
                for k = 1:cubeDim(3)
                    
                    for tetra = 1:6
                        
                        % get tetrahedron vertex deformation
                        vert0Def = [motionVecX_vertex(j+jTetra(tetra,1),i+iTetra(tetra,1),k+kTetra(tetra,1)) motionVecY_vertex(j+jTetra(tetra,1),i+iTetra(tetra,1),k+kTetra(tetra,1)) motionVecZ_vertex(j+jTetra(tetra,1),i+iTetra(tetra,1),k+kTetra(tetra,1))];
                        vert1Def = [motionVecX_vertex(j+jTetra(tetra,2),i+iTetra(tetra,2),k+kTetra(tetra,2)) motionVecY_vertex(j+jTetra(tetra,2),i+iTetra(tetra,2),k+kTetra(tetra,2)) motionVecZ_vertex(j+jTetra(tetra,2),i+iTetra(tetra,2),k+kTetra(tetra,2))];
                        vert2Def = [motionVecX_vertex(j+jTetra(tetra,3),i+iTetra(tetra,3),k+kTetra(tetra,3)) motionVecY_vertex(j+jTetra(tetra,3),i+iTetra(tetra,3),k+kTetra(tetra,3)) motionVecZ_vertex(j+jTetra(tetra,3),i+iTetra(tetra,3),k+kTetra(tetra,3))];
                        vert3Def = [motionVecX_vertex(j+jTetra(tetra,4),i+iTetra(tetra,4),k+kTetra(tetra,4)) motionVecY_vertex(j+jTetra(tetra,4),i+iTetra(tetra,4),k+kTetra(tetra,4)) motionVecZ_vertex(j+jTetra(tetra,4),i+iTetra(tetra,4),k+kTetra(tetra,4))];
                        
                        % get volume of tetrahedron
                        volume = tetraVol(vert0Def,vert1Def,vert2Def,vert3Def);
                        
                        % if the volume is less than or equal to 0, then the voxel is collapsed
                        if volume <= 0
                            collapsedVoxel(j,i,k) = true;
                        end
                    end
                end
            end
        end
        
end

end


function volume = tetraVol(vert0Def,vert1Def,vert2Def,vert3Def)

v1 = vert1Def-vert0Def;
v2 = vert2Def-vert0Def;
v3 = vert3Def-vert0Def;

volume = -det([v1; v2; v3])/6;

end

function [motionVecX_interp, motionVecY_interp, motionVecZ_interp] = interpSomeVoxels(motionVecX, motionVecY, motionVecZ, X, Y, Z, doInterp)

% set interpolated motion vecs to be the old ones
motionVecX_interp = motionVecX;
motionVecY_interp = motionVecY;
motionVecZ_interp = motionVecZ;

% use mask to delete voxels not being interpolated
X_masked = X;
Y_masked = Y;
Z_masked = Z;

motionVecX_masked = motionVecX;
motionVecY_masked = motionVecY;
motionVecZ_masked = motionVecZ;

X_masked(doInterp) = [];
Y_masked(doInterp) = [];
Z_masked(doInterp) = [];

motionVecX_masked(doInterp) = [];
motionVecY_masked(doInterp) = [];
motionVecZ_masked(doInterp) = [];

X_masked = X_masked(:);
Y_masked = Y_masked(:);
Z_masked = Z_masked(:);

motionVecX_masked = motionVecX_masked(:);
motionVecY_masked = motionVecY_masked(:);
motionVecZ_masked = motionVecZ_masked(:);

% set the grid vectors to be very far apart so they don't interfere with
% each other
X_far = X.*10^6;
Y_far = Y.*10^6;
Z_far = Z.*10^6;

X_masked_far = X_masked.*10^6;
Y_masked_far = Y_masked.*10^6;
Z_masked_far = Z_masked.*10^6;

% do interpolation
motionVecX_interp(doInterp) = griddata(X_masked,Y_masked_far,Z_masked_far,motionVecX_masked,X(doInterp),Y_far(doInterp),Z_far(doInterp));
motionVecY_interp(doInterp) = griddata(X_masked_far,Y_masked,Z_masked_far,motionVecY_masked,X_far(doInterp),Y(doInterp),Z_far(doInterp));
motionVecZ_interp(doInterp) = griddata(X_masked_far,Y_masked_far,Z_masked,motionVecZ_masked,X_far(doInterp),Y_far(doInterp),Z(doInterp));

end