function mLL = mll_KF_full(x)

global D;

sigma   = x(1);
alpha1  = x(2);
delta1  = x(3);
alpha2  = x(4);
delta2  = x(5);

pdf             = pdf_KF(D,sigma,alpha1,delta1,alpha2,delta2);
pdf(pdf == 0)   = [];

mLL = -sum(log(pdf));

end

