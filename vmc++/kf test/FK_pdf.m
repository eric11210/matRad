function pdf = FK_pdf(delta,parameters)

sigma   = parameters(1);
alpha1  = parameters(2);
delta1  = parameters(3);
alpha2  = parameters(4);
delta2  = parameters(5);

pdf = alpha1*normpdf(delta,delta1,sigma)+alpha2*normpdf(delta,delta2,sigma)+(1-alpha1-alpha2)*normpdf(delta,0,sigma);

end

