function tBreak = matRad_FSMfreeSplineFit(x,t,K,m1,m2)
% do a free spline fit, time of break

% initialize parameters
initParam       = zeros(4,1);
if K <= numel(t)
    initParam(1) = t(K);%t(end-K+1);
else
    initParam(1) = mean(t);
end
initParam(2)    = m1;
initParam(3)    = x(1);
initParam(4)    = m2;

% bounds
lbT = t(1)+0.05.*(t(end)-t(1));%min([mean(t) t(end-round(K/2))]);
ubT = t(end)-0.05.*(t(end)-t(1));%max([mean(t) t(round(K/2))]);
lb = [lbT -inf -inf -inf]';
ub = [ubT inf inf inf]';

% no display
options.Display = 'none';

% do fit
fitParam = lsqcurvefit(@matRad_FSMfreeSplineFunction,initParam,t,x,lb,ub,options);

% get the break time
tBreak = fitParam(1);

end