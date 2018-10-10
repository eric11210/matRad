function mLL = minusLogLikelihood(x)

global D;

sigma = x(1);
alpha1 = x(2);
delta1 = x(3);
alpha2 = x(4);
delta2 = x(5);

mLL = -sum(log(alpha1*normpdf(D,delta1,sigma)+alpha2*normpdf(D,delta2,sigma)+(1-alpha1-alpha2)*normpdf(D,0,sigma)));

end

