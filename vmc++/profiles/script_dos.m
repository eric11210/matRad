% binary
fname_full = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run',[fname '.dos']);
fid = fopen(fname_full,'r');

Header      = fread(fid,1,'int32');
no_regions  = Header(1);

bixelDose       = fread(fid, no_regions, 'float32');
bixelDoseError  = fread(fid, no_regions, 'float32');

d = reshape(bixelDose,[351 351 351]);
d = permute(d,[2 1 3]);

dError = reshape(bixelDoseError,[351 351 351]);
dError = permute(dError,[2 1 3]);

% read in reversed
% reverse x-y