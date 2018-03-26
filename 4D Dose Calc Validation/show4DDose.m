w = zeros(dij.totalNumOfBixels,1);
w(90) = 1;

physicalDose = cell(dij.numOfScenarios,1);

for i = 1:dij.numOfScenarios
    
    physicalDose{i} = dij.physicalDose{i}*w;
    physicalDose{i} = reshape(physicalDose{i},dij.dimensions);
    
    imCt = imrotate(squeeze(ct.cube{i}(:,56,:)),90);
    imDose = imrotate(squeeze(physicalDose{i}(:,56,:)),90);
    figure
    imagesc(imCt)
    hold on
    super = imagesc(imDose);
    set(super,'AlphaData',0.5);
    
    title(sprintf('Phase %d of %d',i,dij.numOfScenarios));
    xlabel('Increasing posterior voxel')
    ylabel('Increasing inferior voxel')
    
    savefig(sprintf('Phase %d of %d',i,dij.numOfScenarios));
    pause(1)
end