%% filenames

fname_defDose = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\','testDeformation_allDirections_def.dos');

fname_dose = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\','testDeformation_allDirections.dos');

%% setup

Nx = 41;
Ny = 41;
Nz = 42;

%% read defDose

fid_defDose = fopen(fname_defDose,'r');

Header      = fread(fid_defDose,1,'int32');
no_regions  = Header(1);

defDose       = fread(fid_defDose, no_regions, 'float32');
defDoseError  = fread(fid_defDose, no_regions, 'float32');

defDose = reshape(defDose,[Nx Ny Nz]);
defDose = permute(defDose,[2 1 3]);

defDoseError = reshape(defDoseError,[Nx Ny Nz]);
defDoseError = permute(defDoseError,[2 1 3]);

fid_defDose = fclose(fid_defDose);

%% read dose

fid_dose = fopen(fname_dose,'r');

Header      = fread(fid_dose,1,'int32');
no_regions  = Header(1);

dose       = fread(fid_dose, no_regions, 'float32');
doseError  = fread(fid_dose, no_regions, 'float32');

dose = reshape(dose,[Nx Ny Nz]);
dose = permute(dose,[2 1 3]);

doseError = reshape(doseError,[Nx Ny Nz]);
doseError = permute(doseError,[2 1 3]);

fid_dose = fclose(fid_dose);

%% differences
doseDiff = dose-defDose;
doseDiffError = sqrt(doseError.^2+defDoseError.^2);

relDiff = doseDiff./doseDiffError;
relDiff(doseDiff == 0) = 0;
relDiff = relDiff(:);

mu = mean(relDiff);
sig1 = nnz(abs(relDiff) <= 1)./numel(relDiff);
sig2 = nnz(abs(relDiff) <= 2)./numel(relDiff);
sig3 = nnz(abs(relDiff) <= 3)./numel(relDiff);
