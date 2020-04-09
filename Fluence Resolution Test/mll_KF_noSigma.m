function mLL = mll_KF_noSigma(x)

global D;

sigma   = 1;
alpha1  = x(1);
delta1  = x(2);
alpha2  = x(3);
delta2  = x(4);

pdf             = pdf_KF(D,sigma,alpha1,delta1,alpha2,delta2);
pdf(pdf == 0)   = eps;

mLL = -sum(log(pdf));

end

