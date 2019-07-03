% setup
numPhases_vec = 1:5;

fluence = zeros(size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape{1}(1).shapeMap,2), numel(numPhases_vec));
weight = zeros(size(numPhases_vec));

i = 1;
for numPhases = numPhases_vec
    
    % determine name
    fname = sprintf('%d phases.mat',numPhases);
    
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

save('Phase Res Results','fluence','weight')
