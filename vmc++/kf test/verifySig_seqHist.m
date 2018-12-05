%% filenames

fname_simpleSource_1million = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_simpleSource_1million.dos');
fname_simpleSource_2million = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_simpleSource_2million.dos');

fname_phspSource_1million = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_phspSource_1million.dos');
fname_phspSource_2million = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_phspSource_2million.dos');

fname_dupedPhspSource_1million = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_dupedPhspSource_1million.dos');
fname_dupedPhspSource_2million = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_dupedPhspSource_2million.dos');

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

%%

%%% SIMPLE SOURCE %%%

%% read 1 million histories

fid = fopen(fname_simpleSource_1million,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_simpleSource_1million       = fread(fid, no_regions, 'float32');
doseError_simpleSource_1million  = fread(fid, no_regions, 'float32');

dose_simpleSource_1million = reshape(dose_simpleSource_1million,[Nx Ny Nz]);
dose_simpleSource_1million = permute(dose_simpleSource_1million,[2 1 3]);

doseError_simpleSource_1million = reshape(doseError_simpleSource_1million,[Nx Ny Nz]);
doseError_simpleSource_1million = permute(doseError_simpleSource_1million,[2 1 3]);

fclose(fid);

clear fid

%% read 2 million histories

fid = fopen(fname_simpleSource_2million,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_simpleSource_2million       = fread(fid, no_regions, 'float32');
doseError_simpleSource_2million  = fread(fid, no_regions, 'float32');

dose_simpleSource_2million = reshape(dose_simpleSource_2million,[Nx Ny Nz]);
dose_simpleSource_2million = permute(dose_simpleSource_2million,[2 1 3]);

doseError_simpleSource_2million = reshape(doseError_simpleSource_2million,[Nx Ny Nz]);
doseError_simpleSource_2million = permute(doseError_simpleSource_2million,[2 1 3]);

fclose(fid);

%% calculate 1st, 2nd million

dose_simpleSource_1stMillion         = dose_simpleSource_1million(:);
doseError_simpleSource_1stMillion    = doseError_simpleSource_1million(:);

dose_simpleSource_2ndMillion         = 2.*dose_simpleSource_2million(:)-dose_simpleSource_1million(:);
doseError_simpleSource_2ndMillion    = doseError_simpleSource_1million(:);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff_simpleSource       = dose_simpleSource_2ndMillion-dose_simpleSource_1stMillion;
doseDiffError_simpleSource  = sqrt(doseError_simpleSource_2ndMillion.^2+doseError_simpleSource_1stMillion.^2);

delta_simpleSource = doseDiff_simpleSource./doseDiffError_simpleSource;

delta_simpleSource = delta_simpleSource(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose_simpleSource = max([max(dose_simpleSource_1million(:)) max(dose_simpleSource_2million(:))]);

deleteInd = dose_simpleSource_1million(:) < 0.5.*maxDose_simpleSource | dose_simpleSource_2million(:) < 0.5.*maxDose_simpleSource | doseDiffError_simpleSource == 0;

delta_simpleSource_HD               = delta_simpleSource;
delta_simpleSource_HD(deleteInd)    = [];

DELTA = delta_simpleSource_HD;

parameters = fmincon(@fkModel_minusLogLikelihood,parameters_init,A,b,[],[],lb,ub);

sigma_simpleSource  = parameters(1);
alpha1_simpleSource = parameters(2);
delta1_simpleSource = parameters(3);
alpha2_simpleSource = parameters(4);
delta2_simpleSource = parameters(5);

% make figure
xBound = ceil(max(abs(delta_simpleSource_HD)));
x      = -xBound:0.01:xBound;
pdf    = FK_pdf(x,parameters);

figure
hold on
plot(x,pdf,'r-')
histogram(delta_simpleSource_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('simple source, \\sigma = %f',sigma_simpleSource))

%%

%%% PHSP SOURCE %%%

%% read 1 million histories

fid = fopen(fname_phspSource_1million,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_phspSource_1million       = fread(fid, no_regions, 'float32');
doseError_phspSource_1million  = fread(fid, no_regions, 'float32');

dose_phspSource_1million = reshape(dose_phspSource_1million,[Nx Ny Nz]);
dose_phspSource_1million = permute(dose_phspSource_1million,[2 1 3]);

doseError_phspSource_1million = reshape(doseError_phspSource_1million,[Nx Ny Nz]);
doseError_phspSource_1million = permute(doseError_phspSource_1million,[2 1 3]);

fclose(fid);

clear fid

%% read 2 million histories

fid = fopen(fname_phspSource_2million,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_phspSource_2million       = fread(fid, no_regions, 'float32');
doseError_phspSource_2million  = fread(fid, no_regions, 'float32');

dose_phspSource_2million = reshape(dose_phspSource_2million,[Nx Ny Nz]);
dose_phspSource_2million = permute(dose_phspSource_2million,[2 1 3]);

doseError_phspSource_2million = reshape(doseError_phspSource_2million,[Nx Ny Nz]);
doseError_phspSource_2million = permute(doseError_phspSource_2million,[2 1 3]);

fclose(fid);

%% calculate 1st, 2nd million

dose_phspSource_1stMillion         = dose_phspSource_1million(:);
doseError_phspSource_1stMillion    = doseError_phspSource_1million(:);

dose_phspSource_2ndMillion         = 2.*dose_phspSource_2million(:)-dose_phspSource_1million(:);
doseError_phspSource_2ndMillion    = doseError_phspSource_1million(:);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff_phspSource       = dose_phspSource_2ndMillion-dose_phspSource_1stMillion;
doseDiffError_phspSource  = sqrt(doseError_phspSource_2ndMillion.^2+doseError_phspSource_1stMillion.^2);

delta_phspSource = doseDiff_phspSource./doseDiffError_phspSource;

delta_phspSource = delta_phspSource(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose_phspSource = max([max(dose_phspSource_1million(:)) max(dose_phspSource_2million(:))]);

deleteInd = dose_phspSource_1million(:) < 0.5.*maxDose_phspSource | dose_phspSource_2million(:) < 0.5.*maxDose_phspSource | doseDiffError_phspSource == 0;

delta_phspSource_HD               = delta_phspSource;
delta_phspSource_HD(deleteInd)    = [];

DELTA = delta_phspSource_HD;

parameters = fmincon(@fkModel_minusLogLikelihood,parameters_init,A,b,[],[],lb,ub);

sigma_phspSource  = parameters(1);
alpha1_phspSource = parameters(2);
delta1_phspSource = parameters(3);
alpha2_phspSource = parameters(4);
delta2_phspSource = parameters(5);

% make figure
xBound = ceil(max(abs(delta_phspSource_HD)));
x      = -xBound:0.01:xBound;
pdf    = FK_pdf(x,parameters);

figure
hold on
plot(x,pdf,'r-')
histogram(delta_phspSource_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('phsp source, \\sigma = %f',sigma_phspSource))

%%

%%% DUPED PHSP SOURCE %%%

%% read 1 million histories

fid = fopen(fname_dupedPhspSource_1million,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_dupedPhspSource_1million       = fread(fid, no_regions, 'float32');
doseError_dupedPhspSource_1million  = fread(fid, no_regions, 'float32');

dose_dupedPhspSource_1million = reshape(dose_dupedPhspSource_1million,[Nx Ny Nz]);
dose_dupedPhspSource_1million = permute(dose_dupedPhspSource_1million,[2 1 3]);

doseError_dupedPhspSource_1million = reshape(doseError_dupedPhspSource_1million,[Nx Ny Nz]);
doseError_dupedPhspSource_1million = permute(doseError_dupedPhspSource_1million,[2 1 3]);

fclose(fid);

clear fid

%% read 2 million histories

fid = fopen(fname_dupedPhspSource_2million,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_dupedPhspSource_2million       = fread(fid, no_regions, 'float32');
doseError_dupedPhspSource_2million  = fread(fid, no_regions, 'float32');

dose_dupedPhspSource_2million = reshape(dose_dupedPhspSource_2million,[Nx Ny Nz]);
dose_dupedPhspSource_2million = permute(dose_dupedPhspSource_2million,[2 1 3]);

doseError_dupedPhspSource_2million = reshape(doseError_dupedPhspSource_2million,[Nx Ny Nz]);
doseError_dupedPhspSource_2million = permute(doseError_dupedPhspSource_2million,[2 1 3]);

fclose(fid);

%% calculate 1st, 2nd million

dose_dupedPhspSource_1stMillion         = dose_dupedPhspSource_1million(:);
doseError_dupedPhspSource_1stMillion    = doseError_dupedPhspSource_1million(:);

dose_dupedPhspSource_2ndMillion         = 2.*dose_dupedPhspSource_2million(:)-dose_dupedPhspSource_1million(:);
doseError_dupedPhspSource_2ndMillion    = doseError_dupedPhspSource_1million(:);

%% calculate delta, difference divided by combined statistical uncertainty

doseDiff_dupedPhspSource       = dose_dupedPhspSource_2ndMillion-dose_dupedPhspSource_1stMillion;
doseDiffError_dupedPhspSource  = sqrt(doseError_dupedPhspSource_2ndMillion.^2+doseError_dupedPhspSource_1stMillion.^2);

delta_dupedPhspSource = doseDiff_dupedPhspSource./doseDiffError_dupedPhspSource;

delta_dupedPhspSource = delta_dupedPhspSource(:);

%% fit delta to KF model

% first look at voxels with D > 50% Dmax

maxDose_dupedPhspSource = max([max(dose_dupedPhspSource_1million(:)) max(dose_dupedPhspSource_2million(:))]);

deleteInd = dose_dupedPhspSource_1million(:) < 0.5.*maxDose_dupedPhspSource | dose_dupedPhspSource_2million(:) < 0.5.*maxDose_dupedPhspSource | doseDiffError_dupedPhspSource == 0;

delta_dupedPhspSource_HD               = delta_dupedPhspSource;
delta_dupedPhspSource_HD(deleteInd)    = [];

DELTA = delta_dupedPhspSource_HD;

parameters = fmincon(@fkModel_minusLogLikelihood,parameters_init,A,b,[],[],lb,ub);

sigma_dupedPhspSource  = parameters(1);
alpha1_dupedPhspSource = parameters(2);
delta1_dupedPhspSource = parameters(3);
alpha2_dupedPhspSource = parameters(4);
delta2_dupedPhspSource = parameters(5);

% make figure
xBound = ceil(max(abs(delta_dupedPhspSource_HD)));
x      = -xBound:0.01:xBound;
pdf    = FK_pdf(x,parameters);

figure
hold on
plot(x,pdf,'r-')
histogram(delta_dupedPhspSource_HD,100,'Normalization','pdf');
xlabel('\Delta')
ylabel('pdf')
title(sprintf('duped phsp source, \\sigma = %f',sigma_dupedPhspSource))

