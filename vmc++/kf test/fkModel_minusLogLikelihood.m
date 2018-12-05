function mLL = fkModel_minusLogLikelihood(parameters)

global DELTA;

sigma   = parameters(1);
alpha1  = parameters(2);
DELTA1  = parameters(3);
alpha2  = parameters(4);
DELTA2  = parameters(5);

mLL = -sum(log(alpha1*normpdf(DELTA,DELTA1,sigma)+alpha2*normpdf(DELTA,DELTA2,sigma)+(1-alpha1-alpha2)*normpdf(DELTA,0,sigma)));

end

