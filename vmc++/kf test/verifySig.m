%% filenames

fname_simpleSource_seed1 = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_simpleSource_seed1.dos');

fname_simpleSource_seed2 = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_simpleSource_seed2.dos');

fname_phspSource_seed1 = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_phspSource_seed1.dos');

fname_phspSource_seed2 = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_phspSource_seed2.dos');

%% setup

Nx = 101;
Ny = 101;
Nz = 150;

%% stuff for KF model
global delta;
parameters_init = [1    0   0       0   0]';
lb              = [0    0   -inf    0   -inf]';
ub              = [inf  1   inf     1   inf]';
A               = [0    1   0       1   0];
b               = 1;

%% read simpleSource seed1

fid = fopen(fname_simpleSource_seed1,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_simpleSource_seed1       = fread(fid, no_regions, 'float32');
doseError_simpleSource_seed1  = fread(fid, no_regions, 'float32');

dose_simpleSource_seed1 = reshape(dose_simpleSource_seed1,[Nx Ny Nz]);
dose_simpleSource_seed1 = permute(dose_simpleSource_seed1,[2 1 3]);

doseError_simpleSource_seed1 = reshape(doseError_simpleSource_seed1,[Nx Ny Nz]);
doseError_simpleSource_seed1 = permute(doseError_simpleSource_seed1,[2 1 3]);

fclose(fid);

clear fid

%% read simpleSource seed1

fid = fopen(fname_simpleSource_seed2,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_simpleSource_seed2       = fread(fid, no_regions, 'float32');
doseError_simpleSource_seed2  = fread(fid, no_regions, 'float32');

dose_simpleSource_seed2 = reshape(dose_simpleSource_seed2,[Nx Ny Nz]);
dose_simpleSource_seed2 = permute(dose_simpleSource_seed2,[2 1 3]);

doseError_simpleSource_seed2 = reshape(doseError_simpleSource_seed2,[Nx Ny Nz]);
doseError_simpleSource_seed2 = permute(doseError_simpleSource_seed2,[2 1 3]);

fid = fclose(fid);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff_simpleSource       = dose_simpleSource_seed1-dose_simpleSource_seed2;
doseDiffError_simpleSource  = sqrt(doseError_simpleSource_seed1.^2+doseError_simpleSource_seed2.^2);

delta_simpleSource = doseDiff_simpleSource./doseDiffError_simpleSource;

delta_simpleSource = delta_simpleSource(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose_simpleSource = max([max(dose_simpleSource_seed1(:)) max(dose_simpleSource_seed2(:))]);

deleteInd = dose_simpleSource_seed1 < 0.5.*maxDose_simpleSource | dose_simpleSource_seed2 < 0.5.*maxDose_simpleSource | doseDiffError_simpleSource == 0;

delta_simpleSource_HD               = delta_simpleSource;
delta_simpleSource_HD(deleteInd)    = [];

delta = delta_simpleSource_HD;

parameters_simpleSource = fmincon(@fkModel_minusLogLikelihood,parameters_init,A,b,[],[],lb,ub);

sigma_simpleSource  = parameters_simpleSource(1);
alpha1_simpleSource = parameters_simpleSource(2);
delta1_simpleSource = parameters_simpleSource(3);
alpha2_simpleSource = parameters_simpleSource(4);
delta2_simpleSource = parameters_simpleSource(5);

% make figure
x                   = -2:0.01:2;
pdf_simpleSource    = FK_pdf(x,parameters_simpleSource);

figure
hold on
plot(x,pdf_simpleSource,'r-')
histogram(delta_simpleSource_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('simple source, \\sigma = %f',sigma_simpleSource))

%% read pshpSource seed1

fid = fopen(fname_phspSource_seed1,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_phspSource_seed1       = fread(fid, no_regions, 'float32');
doseError_phspSource_seed1  = fread(fid, no_regions, 'float32');

dose_phspSource_seed1 = reshape(dose_phspSource_seed1,[Nx Ny Nz]);
dose_phspSource_seed1 = permute(dose_phspSource_seed1,[2 1 3]);

doseError_phspSource_seed1 = reshape(doseError_phspSource_seed1,[Nx Ny Nz]);
doseError_phspSource_seed1 = permute(doseError_phspSource_seed1,[2 1 3]);

fclose(fid);

clear fid

%% read pshpSource seed1

fid = fopen(fname_phspSource_seed2,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_phspSource_seed2       = fread(fid, no_regions, 'float32');
doseError_phspSource_seed2  = fread(fid, no_regions, 'float32');

dose_phspSource_seed2 = reshape(dose_phspSource_seed2,[Nx Ny Nz]);
dose_phspSource_seed2 = permute(dose_phspSource_seed2,[2 1 3]);

doseError_phspSource_seed2 = reshape(doseError_phspSource_seed2,[Nx Ny Nz]);
doseError_phspSource_seed2 = permute(doseError_phspSource_seed2,[2 1 3]);

fid = fclose(fid);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff_phspSource       = dose_phspSource_seed1-dose_phspSource_seed2;
doseDiffError_phspSource  = sqrt(doseError_phspSource_seed1.^2+doseError_phspSource_seed2.^2);

delta_phspSource = doseDiff_phspSource./doseDiffError_phspSource;

delta_phspSource = delta_phspSource(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose_phspSource = max([max(dose_phspSource_seed1(:)) max(dose_phspSource_seed2(:))]);

deleteInd = dose_phspSource_seed1 < 0.5.*maxDose_phspSource | dose_phspSource_seed2 < 0.5.*maxDose_phspSource | doseDiffError_phspSource == 0;

delta_phspSource_HD               = delta_phspSource;
delta_phspSource_HD(deleteInd)    = [];

delta = delta_phspSource_HD;

parameters_phspSource = fmincon(@fkModel_minusLogLikelihood,parameters_init,A,b,[],[],lb,ub);

sigma_phspSource  = parameters_phspSource(1);
alpha1_phspSource = parameters_phspSource(2);
delta1_phspSource = parameters_phspSource(3);
alpha2_phspSource = parameters_phspSource(4);
delta2_phspSource = parameters_phspSource(5);

% make figure
x                   = -2:0.01:2;
pdf_phspSource      = FK_pdf(x,parameters_phspSource);

figure
hold on
plot(x,pdf_phspSource,'r-')
histogram(delta_phspSource_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('phsp source, \\sigma = %f',sigma_phspSource))

