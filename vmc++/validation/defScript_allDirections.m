%% filenames

fname_defCT = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\phantoms','testDeformation_allDirections_def.ct');
fname_defvec = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\vectors','testDeformation_allDirections_def.vectors');

fname_CT = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\phantoms','testDeformation_allDirections.ct');

%% setup for deformed CT

Nx = 41;
Ny = 41;
Nz = 42;

dX = 0.5;
dY = 0.5;
dZ = 0.5;

halfLengthX = Nx*dX/2;
halfLengthY = Ny*dY/2;
halfLengthZ = Nz*dZ/2;

X = (-halfLengthX:dX:halfLengthX)';
Y = (-halfLengthY:dY:halfLengthY)';
Z = (-halfLengthZ:dZ:halfLengthZ)';

X_centre = cumsum(diff(X))-diff(X)/2+X(1);
Y_centre = cumsum(diff(Y))-diff(Y)/2+Y(1);
Z_centre = cumsum(diff(Z))-diff(Z)/2+Z(1);

dens = zeros(Ny,Nx,Nz);
dens(3:(Ny-2),3:(Nx-2),3:(Nz-2)) = 1;

dens_writeFormat = permute(dens,[2 1 3]);
dens_writeFormat = reshape(dens_writeFormat,[],1);

%% setup for deformed vectors

xCompFac = 1.25;
yCompFac = 1.05;
zCompFac = 2;

xShift = 0*dX;
yShift = 0*dY;
zShift = 0*dZ;

dX_noDef = dX./xCompFac;
dY_noDef = dY./yCompFac;
dZ_noDef = dZ./zCompFac;

halfLengthX_noDef = Nx*dX_noDef/2;
halfLengthY_noDef = Ny*dY_noDef/2;
halfLengthZ_noDef = Nz*dZ_noDef/2;

X_noDef = (-halfLengthX_noDef:dX_noDef:halfLengthX_noDef)'+xShift;
Y_noDef = (-halfLengthY_noDef:dY_noDef:halfLengthY_noDef)'+yShift;
Z_noDef = (-halfLengthZ_noDef:dZ_noDef:halfLengthZ_noDef)'+zShift;

X_noDef(1:2) = X(1:2);
X_noDef(Nx:(Nx+1)) = X(Nx:(Nx+1));
Y_noDef(1:2) = Y(1:2);
Y_noDef(Ny:(Ny+1)) = Y(Ny:(Ny+1));
Z_noDef(1:2) = Z(1:2);
Z_noDef(Nz:(Nz+1)) = Z(Nz:(Nz+1));

[xGrid,yGrid,zGrid] = meshgrid(X,Y,Z);

[xGrid_noDef,yGrid_noDef,zGrid_noDef] = meshgrid(X_noDef,Y_noDef,Z_noDef);

defX = xGrid_noDef-xGrid;
defY = yGrid_noDef-yGrid;
defZ = zGrid_noDef-zGrid;

zeroMask = zeros(Ny+1,Nx+1,Nz+1);
zeroMask(3:(Ny-1),3:(Nx-1),3:(Nz-1)) = 1;
zeroMask = ~zeroMask;

defX(zeroMask) = 0;
defY(zeroMask) = 0;
defZ(zeroMask) = 0;

defX_writeFormat = permute(defX,[2 1 3]);
defX_writeFormat = reshape(defX_writeFormat,[],1);
defY_writeFormat = permute(defY,[2 1 3]);
defY_writeFormat = reshape(defY_writeFormat,[],1);
defZ_writeFormat = permute(defZ,[2 1 3]);
defZ_writeFormat = reshape(defZ_writeFormat,[],1);

numPts = numel(defX);


%% setup for undeformed CT

X_noDef_centre = cumsum(diff(X_noDef))-diff(X_noDef)/2+X_noDef(1);
Y_noDef_centre = cumsum(diff(Y_noDef))-diff(Y_noDef)/2+Y_noDef(1);
Z_noDef_centre = cumsum(diff(Z_noDef))-diff(Z_noDef)/2+Z_noDef(1);

dens_noDef = xCompFac.*yCompFac.*zCompFac.*dens;
dens_noDef_writeFormat = permute(dens_noDef,[2 1 3]);
dens_noDef_writeFormat = reshape(dens_noDef_writeFormat,[],1);

%% write defCt

fid = fopen(fname_defCT,'w');

fwrite(fid,[Ny Nx Nz],'int32');

fwrite(fid,X,'float32');
fwrite(fid,Y,'float32');
fwrite(fid,Z,'float32');

fwrite(fid,dens_writeFormat,'float32');

fclose(fid);

%% write defVec

fid = fopen(fname_defvec,'w');

fprintf(fid,'%d\n',numPts);

for i = 1:numPts
    fprintf(fid,'%6.6f, %6.6f, %6.6f\n',defX_writeFormat(i),defY_writeFormat(i),defZ_writeFormat(i));
end

fclose(fid);

%% write CT

fid = fopen(fname_CT,'w');

fwrite(fid,[Ny Nx Nz],'int32');

fwrite(fid,X_noDef,'float32');
fwrite(fid,Y_noDef,'float32');
fwrite(fid,Z_noDef,'float32');

fwrite(fid,dens_noDef_writeFormat,'float32');

fclose(fid);
