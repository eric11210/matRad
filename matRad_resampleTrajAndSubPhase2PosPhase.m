function l = matRad_resampleTrajAndSubPhase2PosPhase(l_deltaT,t,motionModel)
% l_deltaT and t_deltaT are spaced at resolution deltaT
% t is at some other resolution, we want to interpolate l to this
% resolution, and convert from subphase to posPhase

% first convert subphase to posSubPhase
l_deltaT = motionModel.indices.subPhase2PosSubPhase(l_deltaT);

t_deltaT = motionModel.deltaT_sample.*(0:(numel(l_deltaT)-1))';

% now perform linear interpolation
l_linear = round(interp1(t_deltaT,l_deltaT,t));

% now try sinc interpolation for fun
[T,T_deltaT] = ndgrid(t,t_deltaT);
l_sinc = round(sinc((T-T_deltaT)./motionModel.deltaT_sample)*l_deltaT);

% now convert posSubPhase to posPhase
l = motionModel.indices.posSubPhase2PosPhase(l_linear);

end

