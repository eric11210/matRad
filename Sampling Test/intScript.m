%% parameters

% mean and variance of the dose from MC
% dVar might move around, but dMean should be fixed
dMean   = [0.0001 0.001 0.01 0.1:0.1:2 3:1:10 100 1000];
dSig    = 0.1;%[0.0001 0.001 0.01 0.1:0.1:2];

% threshold dose
dTh     = 1;

syms t;

%% set up arrays
mySampVar   = zeros(numel(dMean),numel(dSig));
matMean     = zeros(numel(dMean),numel(dSig));
matVar      = zeros(numel(dMean),numel(dSig));
matSampVar  = zeros(numel(dMean),numel(dSig));


%% calculations

for i = 1:numel(dMean)
    for j = 1:numel(dSig)
        
        % my own calc
        erfArg = (dTh-dMean(i))./(dSig(j).*sqrt(2));
        
        erfp_coeff  = (dMean(i).*dTh - dMean(i).^2);
        erfp_term   = (1+erf(erfArg))./2;
        
        erfm_coeff  = dSig(j).^2;
        erfm_term   = (1-erf(erfArg))./2;
        
        gauss_coeff = dMean(i).*dSig(j).^2;
        gauss_term  = normpdf(dTh,dMean(i),dSig(j));
        % clear NaNs from bixelDoseError = 0
        gauss_term(isnan(gauss_term)) = 0;
        
        mySampVar(i,j) = erfp_coeff.*erfp_term + erfm_coeff.*erfm_term + gauss_coeff.*gauss_term;
        
        
        % MATLAB's calc
        f = (1./sqrt(2*pi*dSig(j)^2))*exp(-(t-dMean(i)).^2./(2*dSig(j)^2));
        
        matMean(i,j)        = int(t*f,t,[-inf inf]);
        matSampVar(i,j)     = dTh.*int(t*f,t,[-inf dTh]) + int(t^2*f,t,[dTh inf])-matMean(i,j).^2;
        
        
    end
end