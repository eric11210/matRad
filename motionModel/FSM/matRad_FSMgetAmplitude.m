function lambda = matRad_FSMgetAmplitude(x1,x2,deltaT)
% fit two different lines to data, return mean amplitude of the lines

% generate times vectors
t1 = (0:(numel(x1)-1))'.*deltaT;
t2 = (0:(numel(x2)-1))'.*deltaT + t1(end);

% fit lines to data
p1 = polyfit(t1,x1,1);
p2 = polyfit(t2,x2,1);

% get x-values at beginning and end of lines
x10 = polyval(p1,t1(1));
x11 = polyval(p1,t1(end));
x20 = polyval(p2,t2(1));
x21 = polyval(p2,t2(end));

% calculate amplitudes
lambda1 = norm(x11-x10);
lambda2 = norm(x21-x20);

% calculate mean amplitude
lambda = (lambda1 + lambda2)./2;

end