%% display PDVHs and save

maxFluGantryAngleSpacingS = pln.propOpt.VMAToptions.fluGantryAngleSpacing./(1:round(pln.propOpt.prop4D.motionModel.deltaT_sample/0.04));

for maxFluGantryAngleSpacing = maxFluGantryAngleSpacingS
    
    % load
    fname = sprintf('convFrac max flu gantry angle spacing = %.4f.mat',maxFluGantryAngleSpacing);
    load(fname)
    
    % display
    matRad_showDVHBands(recalc.resultGUI.dvh_mean_MC,recalc.resultGUI.pdvh_MC,[1 3 5],cst,pln);
    title(sprintf('PDVH for \\Delta\\theta_{flu} = %f',maxFluGantryAngleSpacing))
    set(gcf, 'Position', get(0, 'Screensize'));
    
    %save
    savefig(sprintf('PDVH for flu gantry angle spacing = %f.fig',maxFluGantryAngleSpacing))
end
