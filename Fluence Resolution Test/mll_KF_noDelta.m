function mLL = mll_KF_noDelta(x)

global D;

sigma = x(1);
alpha1 = 0;
delta1 = 0;
alpha2 = 0;
delta2 = 0;

mLL = -sum(log(pdf_KF(D,sigma,alpha1,delta1,alpha2,delta2)));


end

