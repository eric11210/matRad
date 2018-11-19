function mLL = fkModel_minusLogLikelihood_sigDef(parameters)

global delta;

sigma   = 0.2;
alpha1  = parameters(1);
delta1  = parameters(2);
alpha2  = parameters(3);
delta2  = parameters(4);
%{
if alpha1 > 0
    mLL_comp1   =   log(alpha1)             - log(sigma*sqrt(2*pi))     - (delta-delta1).^2 ./(2.*sigma.^2);
else
    mLL_comp1   =   0;
end

if alpha2 > 0
    mLL_comp2   =   log(alpha2)             - log(sigma*sqrt(2*pi))     - (delta-delta2).^2 ./(2.*sigma.^2);
else
    mLL_comp2   =   0;
end

if 1-alpha1-alpha2 > 0
    mLL_comp3   =   log(1-alpha1-alpha2)    - log(sigma*sqrt(2*pi))     - (delta).^2        ./(2.*sigma.^2);
else
    mLL_comp3   =   0;
end
%}

%mLL = -sum(mLL_comp1+mLL_comp2+mLL_comp3);

mLL = -sum(log(alpha1*normpdf(delta,delta1,sigma)+alpha2*normpdf(delta,delta2,sigma)+(1-alpha1-alpha2)*normpdf(delta,0,sigma)));

end

