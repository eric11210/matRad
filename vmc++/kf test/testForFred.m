%% filenames

fname_1802_9373_j0 = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_simpleSource_1802_9373j0.dos');

fname_1802_9373_j1 = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_simpleSource_1802_9373j1.dos');

%% setup

Nx = 101;
Ny = 101;
Nz = 150;

%% stuff for KF model
global DELTA;
parameters_init = [1    0   0       0   0]';
lb              = [0    0   -inf    0   -inf]';
ub              = [inf  1   inf     1   inf]';
A               = [0    1   0       1   0];
b               = 1;

%% read simpleSource seed1

fid = fopen(fname_1802_9373_j0,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_1802_9373_j0       = fread(fid, no_regions, 'float32');
doseError_1802_9373_j0  = fread(fid, no_regions, 'float32');

dose_1802_9373_j0 = reshape(dose_1802_9373_j0,[Nx Ny Nz]);
dose_1802_9373_j0 = permute(dose_1802_9373_j0,[2 1 3]);

doseError_1802_9373_j0 = reshape(doseError_1802_9373_j0,[Nx Ny Nz]);
doseError_1802_9373_j0 = permute(doseError_1802_9373_j0,[2 1 3]);

fclose(fid);

clear fid

%% read simpleSource seed1

fid = fopen(fname_1802_9373_j1,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_1802_9373_j1       = fread(fid, no_regions, 'float32');
doseError_1802_9373_j1  = fread(fid, no_regions, 'float32');

dose_1802_9373_j1 = reshape(dose_1802_9373_j1,[Nx Ny Nz]);
dose_1802_9373_j1 = permute(dose_1802_9373_j1,[2 1 3]);

doseError_1802_9373_j1 = reshape(doseError_1802_9373_j1,[Nx Ny Nz]);
doseError_1802_9373_j1 = permute(doseError_1802_9373_j1,[2 1 3]);

fid = fclose(fid);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff       = dose_1802_9373_j0-dose_1802_9373_j1;
doseDiffError  = sqrt(doseError_1802_9373_j0.^2+doseError_1802_9373_j1.^2);

delta = doseDiff./doseDiffError;

delta = delta(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose = max([max(dose_1802_9373_j0(:)) max(dose_1802_9373_j1(:))]);

deleteInd = dose_1802_9373_j0 < 0.5.*maxDose | dose_1802_9373_j1 < 0.5.*maxDose | doseDiffError == 0;

delta_HD               = delta;
delta_HD(deleteInd)    = [];

DELTA = delta_HD;

parameters = fmincon(@fkModel_minusLogLikelihood,parameters_init,A,b,[],[],lb,ub);

sigma  = parameters(1);
alpha1 = parameters(2);
delta1 = parameters(3);
alpha2 = parameters(4);
delta2 = parameters(5);

% make figure
xBound  = ceil(max(abs(delta_HD)));
x       = -xBound:0.01:xBound;
pdf     = FK_pdf(x,parameters);

figure
hold on
plot(x,pdf,'r-')
histogram(delta_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('\\sigma = %f',sigma))

