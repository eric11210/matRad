function [tKnotNew,xKnotNew] = matRad_FSMpartialSplineFit(x1,x2,t1,t2,xKnot,tKnot)
% do a partial spline fit, return index of break and position of new knot

%% first do two separate fits
if nargin < 5
    % set xKnot to 0
    xKnot = 0;
    % set tKnot to 0
    tKnot = 0;
end

% offset x,t by the knot positiong
% this way we can set the y-intercept of the first line to be 0
x1 = x1-xKnot;
x2 = x2-xKnot;
t1 = t1-tKnot;
t2 = t2-tKnot;

% fit lines to data
p1 = polyfit(t1,x1,1);
p2 = polyfit(t2,x2,1);

%% now do full fit

% generate times vector
t = [t1; t2];

% generate positions vector
x = [x1; x2];

% initialize parameters
initParam(1) = t2(1);
initParam(2) = p1(1);
lbT = t(1)+0.05.*(t(end)-t(1));%min([mean(t) t(end-round(K/2))]);
ubT = t(end)-0.05.*(t(end)-t(1));%max([mean(t) t(round(K/2))]);
if nargin == 5
    initParam(3) = p2(1);
    
    % bounds
    lb = [lbT -inf -inf]';
    ub = [ubT inf inf]';
else
    initParam(3) = p1(2);
    initParam(4) = p2(1);
    
    % bounds
    lb = [lbT -inf -inf -inf]';
    ub = [ubT inf inf inf]';
end

% no display
options.Display = 'none';

% do fit
fitParam = lsqcurvefit(@matRad_FSMpartialSplineFunction,initParam,t,x,lb,ub,options);

% get the new knot time position, including the offset
xKnotNew = xKnot+matRad_FSMpartialSplineFunction(fitParam,fitParam(1));
tKnotNew = tKnot+fitParam(1);


end