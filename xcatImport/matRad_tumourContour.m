function indV = matRad_tumourContour(ct,radius,tumourMotion)

tumourCloudxLine = (-radius:(ct.resolution.x/10):radius);
tumourCloudyLine = (-radius:(ct.resolution.y/10):radius);
tumourCloudzLine = (-radius:(ct.resolution.z/10):radius);

[tumourCloudxGrid,tumourCloudyGrid,tumourCloudzGrid] = meshgrid(tumourCloudxLine,tumourCloudyLine,tumourCloudzLine);

tumourCloudMask = sqrt(tumourCloudxGrid.^2+tumourCloudyGrid.^2+tumourCloudzGrid.^2) <= radius;

tumourCloudx = tumourCloudxGrid(tumourCloudMask)+interp1(ct.x,ct.tumourMotion.coordsVox(1,1));
tumourCloudy = tumourCloudyGrid(tumourCloudMask)+interp1(ct.y,ct.tumourMotion.coordsVox(1,2));
tumourCloudz = tumourCloudzGrid(tumourCloudMask)+interp1(ct.z,ct.tumourMotion.coordsVox(1,3));

tumourCloudx_coord = interp1(ct.x,1:ct.cubeDim(2),tumourCloudx);
tumourCloudy_coord = interp1(ct.y,1:ct.cubeDim(1),tumourCloudy);
tumourCloudz_coord = interp1(ct.z,1:ct.cubeDim(3),tumourCloudz);

if tumourMotion
    
    tumourCloudx_coord_p1 = tumourCloudx_coord;
    tumourCloudy_coord_p1 = tumourCloudy_coord;
    tumourCloudz_coord_p1 = tumourCloudz_coord;
    
    tumourCloudx_coord = zeros(numel(tumourCloudx_coord_p1),ct.tumourMotion.numFrames);
    tumourCloudy_coord = zeros(numel(tumourCloudx_coord_p1),ct.tumourMotion.numFrames);
    tumourCloudz_coord = zeros(numel(tumourCloudx_coord_p1),ct.tumourMotion.numFrames);
    
    for frame = 1:ct.tumourMotion.numFrames
        
        tumourCloudx_coord(:,frame) = interp3(ct.motionVecX{frame},tumourCloudx_coord_p1,tumourCloudy_coord_p1,tumourCloudz_coord_p1);
        tumourCloudy_coord(:,frame) = interp3(ct.motionVecY{frame},tumourCloudx_coord_p1,tumourCloudy_coord_p1,tumourCloudz_coord_p1);
        tumourCloudz_coord(:,frame) = interp3(ct.motionVecZ{frame},tumourCloudx_coord_p1,tumourCloudy_coord_p1,tumourCloudz_coord_p1);
    end
end

tumourCloudx_coord = round(tumourCloudx_coord(:));
tumourCloudy_coord = round(tumourCloudy_coord(:));
tumourCloudz_coord = round(tumourCloudz_coord(:));

indV = unique(sub2ind(ct.cubeDim,tumourCloudy_coord,tumourCloudx_coord,tumourCloudz_coord));

end