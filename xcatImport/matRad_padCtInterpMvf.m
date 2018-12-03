function ct = matRad_padCtInterpMvf(ct,xcatLog)

% we want to pad ct for two reasons (both related to def_tetra in vmc++):
% 1) we need a 2-voxel thick pad on all edges for def_tetra to work
% 2) we also cannot have any voxels/tetrahedrons collapse

% at the same time, interpolate all voxels with 0 density (air) so that
% voxels/tetrahedrons don't collapse


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

%%%%% TEMPORARY
    ct.cube{1} = padarray(ct.cube{1},[yPadPre xPadPre zPadPre],0,'pre');
    ct.cube{1} = padarray(ct.cube{1},[yPadPost xPadPost zPadPost],0,'post');
    
    ct.cubeHU{1} = padarray(ct.cubeHU{1},[yPadPre xPadPre zPadPre],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'pre');
    ct.cubeHU{1} = padarray(ct.cubeHU{1},[yPadPost xPadPost zPadPost],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'post');


% now pad the motionVecs
for frame = 1:xcatLog.numFrames
    
    frame = 5;
    
    fprintf('\nFrame %d of %d.\n',frame,xcatLog.numFrames);
    
    ct.cube{frame} = padarray(ct.cube{frame},[yPadPre xPadPre zPadPre],0,'pre');
    ct.cube{frame} = padarray(ct.cube{frame},[yPadPost xPadPost zPadPost],0,'post');
    
    ct.cubeHU{frame} = padarray(ct.cubeHU{frame},[yPadPre xPadPre zPadPre],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'pre');
    ct.cubeHU{frame} = padarray(ct.cubeHU{frame},[yPadPost xPadPost zPadPost],interp1(ct.hlut(:,2),ct.hlut(:,1),0),'post');
    
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
    
    % use mask to delete voxels not being interpolated
    baseMotionVecX_masked = baseMotionVecX;
    baseMotionVecY_masked = baseMotionVecY;
    baseMotionVecZ_masked = baseMotionVecZ;
    
    motionVecX_masked = ct.motionVecX{frame};
    motionVecY_masked = ct.motionVecY{frame};
    motionVecZ_masked = ct.motionVecZ{frame};
    
    baseMotionVecX_masked(doInterp) = [];
    baseMotionVecY_masked(doInterp) = [];
    baseMotionVecZ_masked(doInterp) = [];
    
    motionVecX_masked(doInterp) = [];
    motionVecY_masked(doInterp) = [];
    motionVecZ_masked(doInterp) = [];
    
    % create the scattered interpolants
    motionVecXInterpolant = scatteredInterpolant(baseMotionVecX_masked',baseMotionVecY_masked',baseMotionVecZ_masked',motionVecX_masked');
    motionVecYInterpolant = scatteredInterpolant(baseMotionVecX_masked',baseMotionVecY_masked',baseMotionVecZ_masked',motionVecY_masked');
    motionVecZInterpolant = scatteredInterpolant(baseMotionVecX_masked',baseMotionVecY_masked',baseMotionVecZ_masked',motionVecZ_masked');
    
    % finally, use the interpolants to update the values in the motion vecs
    ct.motionVecX{frame}(doInterp) = motionVecXInterpolant(baseMotionVecX(doInterp),baseMotionVecY(doInterp),baseMotionVecZ(doInterp));
    ct.motionVecY{frame}(doInterp) = motionVecYInterpolant(baseMotionVecX(doInterp),baseMotionVecY(doInterp),baseMotionVecZ(doInterp));
    ct.motionVecZ{frame}(doInterp) = motionVecZInterpolant(baseMotionVecX(doInterp),baseMotionVecY(doInterp),baseMotionVecZ(doInterp));
    
    % now do a mean filter on the whole thing
    %ct.motionVecX{frame} = convn(ct.motionVecX{frame},ones(3,3,3)./3^3,'same');
    %ct.motionVecY{frame} = convn(ct.motionVecY{frame},ones(3,3,3)./3^3,'same');
    %ct.motionVecZ{frame} = convn(ct.motionVecZ{frame},ones(3,3,3)./3^3,'same');
    
    % now do a median filter on the whole thing
    %ct.motionVecX{frame} = medfilt3(ct.motionVecX{frame});
    %ct.motionVecY{frame} = medfilt3(ct.motionVecY{frame});
    %ct.motionVecZ{frame} = medfilt3(ct.motionVecZ{frame});
    
    % fix the 2-voxel thick pad
    ct.motionVecX{frame}(~notPad) = baseMotionVecX(~notPad);
    ct.motionVecY{frame}(~notPad) = baseMotionVecY(~notPad);
    ct.motionVecZ{frame}(~notPad) = baseMotionVecZ(~notPad);
    
    %% check for collapsed voxels
    collapsedVoxel = false(ct.cubeDim);
    
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