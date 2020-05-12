w = zeros(dij.totalNumOfBixels,1);
w(110) = 1;
%w(66) = 1;

physicalDose = cell(dij.numOfScenarios,1);

for i = 1:dij.numPhases
    
    physicalDose{i} = dij.physicalDose{i}*w;
    physicalDose{i} = reshape(physicalDose{i},dij.dimensions);
    
    imCt    = imrotate(squeeze(ct.cubeHU{1}(:,53,:)),90);
    imDose  = imrotate(squeeze(physicalDose{i}(:,53,:)),90);
    alpha   = 1.*(imDose > 0);
    
    x = [ct.y(1) ct.y(end)];
    y = [ct.z(1) ct.z(end)];
    
    % create figure
    figure;
    % create two axes
    ax1 = axes;
    ax2 = axes;
    % put the images in the axes
    imagesc(ax1,x,y,imCt)
    imagesc(ax2,x,y,imDose,'AlphaData',alpha)
    % hide the top axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    % give each one its own colormap
    colormap(ax1,'gray')
    colormap(ax2,'default')
    % set aspect ratio
    ax1.DataAspectRatio = [1 1 1];
    ax2.DataAspectRatio = [1 1 1];
    
    % then add colorbars and get everything lined up
    cb1 = colorbar(ax1,'WestOutside');
    cb2 = colorbar(ax2,'EastOutside');
    set([ax1 ax2],'Position',mean([ax1.Position; ax2.Position],1));
    
    caxis(ax1,[-1000 600])
    caxis(ax2,[0 0.6])
    
    %title(sprintf('Phase %d of %d',i,dij.numPhases));
    xlabel(ax1,'posterior position / simm')
    ylabel(ax1,'inferior position / simm')
    
    %set(gcf, 'Colormap', gray);
    
    %savefig(sprintf('Phase %d of %d',i,dij.numOfScenarios));
    pause(1)
    
end