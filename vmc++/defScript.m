fname_defCT = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\phantoms','testDeformation_def.ct');
fname_defvec = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\vectors','testDeformation_def.vectors');
fname_defDose = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\','testDeformation_def.dos');

fname_CT = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\phantoms','testDeformation.ct');
fname_dose = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\','testDeformation.dos');

%% defCt
fid_defCT = fopen(fname_defCT,'r');

Nx = fread(fid_defCT,1,'int32');
Ny = fread(fid_defCT,1,'int32');
Nz = fread(fid_defCT,1,'int32');

X = zeros(Nx+1,1);
Y = zeros(Ny+1,1);
Z = zeros(Nz+1,1);

for i = 1:(Nx+1)
    X(i) = fread(fid_defCT,1,'float32');
end

for i = 1:(Ny+1)
    Y(i) = fread(fid_defCT,1,'float32');
end

for i = 1:(Nz+1)
    Z(i) = fread(fid_defCT,1,'float32');
end

N = Nx*Ny*Nz;

dens = fread(fid_defCT,N,'float32');

dens = reshape(dens,[Ny Nx Nz]);
dens = permute(dens,[2 1 3]);

fid_defCT = fclose(fid_defCT);

%% vectors
fid_defvec = fopen(fname_defvec,'r');

%Nvec = (Nx+1)*(Ny+1)*(Nz+1);

% get number of elements
tline = fgetl(fid_defvec);
Nvec_f = str2double(tline);


dX = zeros(Nvec_f,1);
dY = zeros(Nvec_f,1);
dZ = zeros(Nvec_f,1);

for i = 1:Nvec_f
    
    % get next line
    tline = fgetl(fid_defvec);
    
    % find comma indices
    commaInd = strfind(tline,',');
    
    % determine deformations
    dx_str = tline(1:(commaInd(1)-1));
    dy_str = tline((commaInd(1)+2):(commaInd(2)-1));
    dz_str = tline((commaInd(2)+2):end);
    
    % enter deformations
    dX(i) = str2double(dx_str);
    dY(i) = str2double(dy_str);
    dZ(i) = str2double(dz_str);
    
end

dX = reshape(dX,[Ny+1 Nx+1 Nz+1]);
dX = permute(dX,[2 1 3]);

dY = reshape(dY,[Ny+1 Nx+1 Nz+1]);
dY = permute(dY,[2 1 3]);

dZ = reshape(dZ,[Ny+1 Nx+1 Nz+1]);
dZ = permute(dZ,[2 1 3]);

fid_defvec = fclose(fid_defvec);

%% ct
fid_CT = fopen(fname_CT,'w');

dZVec = squeeze(min(min(dZ,[],1),[],2));
ZnoDef = Z+dZVec;
densnoDef = 2*dens;
densnoDef = permute(densnoDef,[2 1 3]);
densnoDef = reshape(densnoDef,[],1);

fwrite(fid_CT,[Ny Nx Nz],'int32');

fwrite(fid_CT,X,'float32');
fwrite(fid_CT,Y,'float32');
fwrite(fid_CT,ZnoDef,'float32');

fwrite(fid_CT,densnoDef,'float32');

fid_CT = fclose(fid_CT);

%% defDose
fid_defDose = fopen(fname_defDose,'r');

Header      = fread(fid_defDose,1,'int32');
no_regions  = Header(1);

defDose       = fread(fid_defDose, no_regions, 'float32');
defDoseError  = fread(fid_defDose, no_regions, 'float32');

defDose = reshape(defDose,[Ny Nx Nz]);
defDose = permute(defDose,[2 1 3]);

defDoseError = reshape(defDoseError,[Ny Nx Nz]);
defDoseError = permute(defDoseError,[2 1 3]);

fid_defDose = fclose(fid_defDose);

%% dose
fid_dose = fopen(fname_dose,'r');

Header      = fread(fid_dose,1,'int32');
no_regions  = Header(1);

dose       = fread(fid_dose, no_regions, 'float32');
doseError  = fread(fid_dose, no_regions, 'float32');

dose = reshape(dose,[Ny Nx Nz]);
dose = permute(dose,[2 1 3]);

doseError = reshape(doseError,[Ny Nx Nz]);
doseError = permute(doseError,[2 1 3]);

fid_dose = fclose(fid_dose);

%% differences
diff = dose-defDose;
diffError = sqrt(doseError.^2+defDoseError.^2);

relDiff = diff./diffError;
relDiff(diff == 0) = 0;
relDiff = relDiff(:);

mu = mean(relDiff);
sig1 = nnz(abs(relDiff) <= 1)./numel(relDiff);
sig2 = nnz(abs(relDiff) <= 2)./numel(relDiff);
sig3 = nnz(abs(relDiff) <= 3)./numel(relDiff);

