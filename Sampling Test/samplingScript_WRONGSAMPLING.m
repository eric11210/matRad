%%%% NOTE: there is a discrepancy here which samples all negative values of
%%%% dMC, even ones which have an absolute value greater than dTh. I could
%%%% fix this, but it would require more work than is worth it, since
%%%% negative doses would never appear in MC anyway.

%% parameters

% mean and variance of the dose from MC
% dVar might move around, but dMean should be fixed
dMean   = [0.0001 0.001 0.01 0.1:0.1:2 3:1:10 100 1000];
dSig    = [0.0001 0.001 0.01 0.1:0.1:2];

% threshold dose
dTh     = 1;

% number of trials
nTrials = 1000000;

%% set up arrays
dSampled_cVar       = zeros(numel(dMean),numel(dSig));
dSampled_var        = zeros(numel(dMean),numel(dSig));
dSampled_expectSVar = zeros(numel(dMean),numel(dSig));
bias_sVar           = zeros(numel(dMean),numel(dSig));
dSampled_expectD    = zeros(numel(dMean),numel(dSig));


%% simulation

for i = 1:numel(dMean)
    for j = 1:numel(dSig)
        
        % simulate MC
        dMC = normrnd(dMean(i),dSig(j),[nTrials 1]);
        
        % sampling
        r = rand([nTrials 1]);
        probToKeep = abs(dMC./dTh);
        
        dSampled = dMC;
        dSampled(dMC < dTh & r < probToKeep) = sign(dMC(dMC < dTh & r < probToKeep)).*dTh;
        dSampled(dMC < dTh & r > probToKeep) = 0;
        
        % sample variance
        erfArg  = (dTh-dMC)./(dSig(j).*sqrt(2));
        erfMArg = dMC./(dSig(j).*sqrt(2));
        
        erfp_coeff  = -dMC.^2;
        erfp_term   = (1+erf(erfArg))./2;
        
        erfm_coeff  = dSig(j).^2;
        erfm_term   = (1-erf(erfArg))./2;
        
        erfM_coeff  = dTh.*dMC;
        erfM_term   = (2.*erf(erfMArg)+erf(erfArg)-1)./2;
        
        gauss_coeff = dMC.*dSig(j).^2;
        gauss_term  = normpdf(dTh,dMC,dSig(j));
        % clear NaNs from bixelDoseError = 0
        gauss_term(isnan(gauss_term)) = 0;
        
        gaussM_coeff = 2.*dTh.*dSig(j).^2;
        gaussM_term  = normpdf(0,dMC,dSig(j));
        % clear NaNs from bixelDoseError = 0
        gaussM_term(isnan(gaussM_term)) = 0;
        
        dSampled_sVar = erfp_coeff.*erfp_term + erfm_coeff.*erfm_term + erfM_coeff.*erfM_term + gauss_coeff.*gauss_term + gaussM_coeff.*gaussM_term;
        
        %% calculated stuff
        
        % variance
        erfArg = (dTh-dMean(i))./(dSig(j).*sqrt(2));
        erfMArg = dMean(i)./(dSig(j).*sqrt(2));
        
        erfp_coeff  = -dMean(i).^2;
        erfp_term   = (1+erf(erfArg))./2;
        
        erfm_coeff  = dSig(j).^2;
        erfm_term   = (1-erf(erfArg))./2;
        
        erfM_coeff  = dTh.*dMean(i);
        erfM_term   = (2.*erf(erfMArg)+erf(erfArg)-1)./2;
        
        gauss_coeff = dMean(i).*dSig(j).^2;
        gauss_term  = normpdf(dTh,dMean(i),dSig(j));
        % clear NaNs from bixelDoseError = 0
        gauss_term(isnan(gauss_term)) = 0;
        
        gaussM_coeff = 2.*dTh.*dSig(j).^2;
        gaussM_term  = normpdf(0,dMean(i),dSig(j));
        % clear NaNs from bixelDoseError = 0
        gaussM_term(isnan(gaussM_term)) = 0;
        
        dSampled_var(i,j) = erfp_coeff.*erfp_term + erfm_coeff.*erfm_term + erfM_coeff.*erfM_term + gauss_coeff.*gauss_term + gaussM_coeff.*gaussM_term;
        
        
        % expectation value of dose
        dSampled_expectD(i,j) = mean(dSampled);
        
        % calculated variance of dose
        dSampled_cVar(i,j) = var(dSampled);
        
        % expectation value of sample variance
        dSampled_expectSVar(i,j) = mean(dSampled_sVar);
        
        % bias of sample variance
        bias_sVar(i,j) = dSampled_expectSVar(i,j)-dSampled_var(i,j);
        
    end
end

%% limiting behaviour
var_ldMean = repmat(sqrt(2/pi).*dSig.*dTh,[numel(dMean),1]);
var_hdMean = repmat(dSig.^2,[numel(dMean),1]);

%{
%% make plots

% dMean << dTh
plot(dMean,dSampled_var)
hold
plot(dMean,var_ldMean,'k')
xlim([0 1])
legend(cellstr(num2str(dSig')),'location','best')
title('Behaviour of variance when $\overline{D} << D_{th}$','interpreter','latex')
xlabel('$\overline{D}$','interpreter','latex')
ylabel('$\sigma_D^2$','interpreter','latex')

% dMean >> dTh
%}
