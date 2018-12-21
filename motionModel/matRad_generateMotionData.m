function data = matRad_generateMotionData(dataOptions)

dt = 1./dataOptions.samplingFreq;

data.t_cut = (0:dt:(dataOptions.totTime-dt))';

switch dataOptions.function
    case 'sawtooth'
        
        data.x_cut = (dataOptions.amplitude./2).*sawtooth((2.*pi./dataOptions.period).*data.t_cut);
    case 'square'
        
        data.x_cut = (dataOptions.amplitude./2).*square((2.*pi./dataOptions.period).*data.t_cut);
end


end