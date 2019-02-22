function xFit = matRad_FSMpartialSplineFunction(param,t)
% free spline function

% get parameters
tBreak  = param(1);
m1      = param(2);
if numel(param) == 3
    b1      = 0;
    m2      = param(3);
else
    b1      = param(3);
    m2      = param(4);
end
b2      = (m1-m2).*tBreak+b1;

% calculate values fitted spline function
xFit = zeros(size(t));
xFit(t <  tBreak)   = m1.*t(t <  tBreak) + b1;
xFit(t >= tBreak)   = m2.*t(t >= tBreak) + b2;

end