fname_CT = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\phantoms','test_def.ct');
fname_vec = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\vectors','test_def.vectors');

%% ct
fid_CT = fopen(fname_CT,'r');

Nx = fread(fid_CT,1,'int32');
Ny = fread(fid_CT,1,'int32');
Nz = fread(fid_CT,1,'int32');

X = zeros(Nx+1,1);
Y = zeros(Ny+1,1);
Z = zeros(Nz+1,1);

for i = 1:(Nx+1)
    X(i) = fread(fid_CT,1,'float32');
end

for i = 1:(Ny+1)
    Y(i) = fread(fid_CT,1,'float32');
end

for i = 1:(Nz+1)
    Z(i) = fread(fid_CT,1,'float32');
end

N = Nx*Ny*Nz;

dens = fread(fid_CT,N,'float32');

dens = reshape(dens,[Ny Nx Nz]);
dens = permute(dens,[2 1 3]);

fid_CT = fclose(fid_CT);

%% vectors
fid_vec = fopen(fname_vec,'r');

%Nvec = (Nx+1)*(Ny+1)*(Nz+1);

% get number of elements
tline = fgetl(fid_vec);
Nvec_f = str2double(tline);


dX = zeros(Nvec_f,1);
dY = zeros(Nvec_f,1);
dZ = zeros(Nvec_f,1);

for i = 1:Nvec_f
    
    % get next line
    tline = fgetl(fid_vec);
    
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


fid_vec = fclose(fid_vec);