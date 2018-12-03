%% filenames

fname_vmc = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','noVR.dos');

fname_DOSXYZnrc = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','noVR_DOSXYZnrc.dos');

%% setup

Nx = 41;
Ny = 41;
Nz = 60;

%% read vmc

fid = fopen(fname_vmc,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_vmc       = fread(fid, no_regions, 'float32');
doseError_vmc  = fread(fid, no_regions, 'float32');

dose_vmc = reshape(dose_vmc,[Nx Ny Nz]);
dose_vmc = permute(dose_vmc,[2 1 3]);

doseError_vmc = reshape(doseError_vmc,[Nx Ny Nz]);
doseError_vmc = permute(doseError_vmc,[2 1 3]);

fclose(fid);

relError_vmc = doseError_vmc./dose_vmc;

clear fid

%% read DOSXYZnrc

fid = fopen(fname_DOSXYZnrc,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

dose_DOSXYZnrc       = fread(fid, no_regions, 'float32');
doseError_DOSXYZnrc  = fread(fid, no_regions, 'float32');

dose_DOSXYZnrc = reshape(dose_DOSXYZnrc,[Nx Ny Nz]);
dose_DOSXYZnrc = permute(dose_DOSXYZnrc,[2 1 3]);

doseError_DOSXYZnrc = reshape(doseError_DOSXYZnrc,[Nx Ny Nz]);
doseError_DOSXYZnrc = permute(doseError_DOSXYZnrc,[2 1 3]);

fid = fclose(fid);

relError_DOSXYZnrc  = doseError_DOSXYZnrc./dose_DOSXYZnrc;

