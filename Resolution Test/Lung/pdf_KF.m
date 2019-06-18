function pdf = pdf_KF(delta,sigma,alpha1,delta1,alpha2,delta2)

pdf = alpha1*normpdf(delta,delta1,sigma)+alpha2*normpdf(delta,delta2,sigma)+(1-alpha1-alpha2)*normpdf(delta,0,sigma);

end

