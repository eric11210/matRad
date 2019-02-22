function lineVelocity = matRad_FSMgetFitVelocity(x,deltaT)
% fit a single line to some data, return the velocity

% generate times vector
t = (0:(numel(x)-1))'.*deltaT;

% fit line to data
p = polyfit(t,x,1);

% get velocity
lineVelocity = p(1);

end