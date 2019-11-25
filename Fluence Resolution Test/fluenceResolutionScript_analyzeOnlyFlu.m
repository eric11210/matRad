% setup
maxFluGantryAngleSpacingS = [4 2 1 0.5 0.25];

fluence = zeros(size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,2), numel(maxFluGantryAngleSpacingS));
weight = zeros(size(maxFluGantryAngleSpacingS));

i = 1;
for maxFluGantryAngleSpacing = maxFluGantryAngleSpacingS
    
    % determine name
    fname = sprintf('Max fluence gantry angle spacing = %.3f.mat',maxFluGantryAngleSpacing);
    
    % load in recalculated results
    load(fname,'resultGUI','recalc');
    
    for i_flu = 1:numel(recalc.apertureInfo.beam)
        for phase_flu = 1:recalc.apertureInfo.numPhases
            fluence(:,:,i) = fluence(:,:,i)+recalc.apertureInfo.beam(i_flu).shape{phase_flu}(1).shapeMap;
            weight(i) = weight(i)+recalc.apertureInfo.beam(i_flu).shape{phase_flu}(1).weight;
        end
    end
    
    i = i+1;
end

save('Conventional Fluence Res Results','fluence','weight')
