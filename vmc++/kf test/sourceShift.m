%% filenames

fname_sourceShift_base = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_sourceShift_base.dos');

fname_sourceShift_shift = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_sourceShift_shift.dos');

%% setup

Nx = 101;
Ny = 101;
Nz = 150;

%% stuff for KF model
global delta;
parameters_init = [0.25   -0.2       0.25   0.2]';
lb              = [0   -inf    0   -inf]';
ub              = [1   inf     1   inf]';
A               = [1   0       1   0];
b               = 1;

%% read sourceShift_base

fid = fopen(fname_sourceShift_base,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_sourceShift_base       = fread(fid, no_regions, 'float32');
doseError_sourceShift_base  = fread(fid, no_regions, 'float32');

dose_sourceShift_base = reshape(dose_sourceShift_base,[Nx Ny Nz]);
dose_sourceShift_base = permute(dose_sourceShift_base,[2 1 3]);

doseError_sourceShift_base = reshape(doseError_sourceShift_base,[Nx Ny Nz]);
doseError_sourceShift_base = permute(doseError_sourceShift_base,[2 1 3]);

fclose(fid);

clear fid

%% read sourceShift_shift

fid = fopen(fname_sourceShift_shift,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_sourceShift_shift       = fread(fid, no_regions, 'float32');
doseError_sourceShift_shift  = fread(fid, no_regions, 'float32');

dose_sourceShift_shift = reshape(dose_sourceShift_shift,[Nx Ny Nz]);
dose_sourceShift_shift = permute(dose_sourceShift_shift,[2 1 3]);

doseError_sourceShift_shift = reshape(doseError_sourceShift_shift,[Nx Ny Nz]);
doseError_sourceShift_shift = permute(doseError_sourceShift_shift,[2 1 3]);

fid = fclose(fid);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff_sourceShift       = dose_sourceShift_base-dose_sourceShift_shift;
doseDiffError_sourceShift  = sqrt(doseError_sourceShift_base.^2+doseError_sourceShift_shift.^2);

delta_sourceShift = doseDiff_sourceShift./doseDiffError_sourceShift;

delta_sourceShift = delta_sourceShift(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose_sourceShift = max([max(dose_sourceShift_base(:)) max(dose_sourceShift_shift(:))]);

deleteInd = dose_sourceShift_base < 0.5.*maxDose_sourceShift | dose_sourceShift_shift < 0.5.*maxDose_sourceShift | doseDiffError_sourceShift == 0;

delta_sourceShift_HD               = delta_sourceShift;
delta_sourceShift_HD(deleteInd)    = [];

delta = delta_sourceShift_HD;

parameters_sourceShift = fmincon(@fkModel_minusLogLikelihood_sigDef,parameters_init,A,b,[],[],lb,ub);

sigma_sourceShift  = 0.2;
alpha1_sourceShift = parameters_sourceShift(1);
delta1_sourceShift = parameters_sourceShift(2);
alpha2_sourceShift = parameters_sourceShift(3);
delta2_sourceShift = parameters_sourceShift(4);

% make figure
xBound = ceil(max(abs(delta_sourceShift_HD)));
x                   = -xBound:0.01:xBound;
pdf_sourceShift    = FK_pdf(x,[sigma_sourceShift; parameters_sourceShift]);

figure
hold on
plot(x,pdf_sourceShift,'r-')
histogram(delta_sourceShift_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('source shift, \\sigma = %f',sigma_sourceShift))
