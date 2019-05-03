fname_defvec_txt = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\vectors','testDeformation_def.vectors');
fname_defvec_bin = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\vectors','testDeformation_def_bin.vectors');

fid_defvec_txt = fopen(fname_defvec_txt,'r');
fid_defvec_bin = fopen(fname_defvec_bin,'w');

% get number of elements
tline = fgetl(fid_defvec_txt);
Nvec_f = str2double(tline);

% write number of elements
fwrite(fid_defvec_bin,Nvec_f,'int32');

for i = 1:Nvec_f
    
    % get next line
    tline = fgetl(fid_defvec_txt);
    
    % find comma indices
    commaInd = strfind(tline,',');
    
    % determine deformations
    dx_str = tline(1:(commaInd(1)-1));
    dy_str = tline((commaInd(1)+2):(commaInd(2)-1));
    dz_str = tline((commaInd(2)+2):end);
    
    % get deformations
    dX = str2double(dx_str);
    dY = str2double(dy_str);
    dZ = str2double(dz_str);
    
    % write deformations
    fwrite(fid_defvec_bin,dX,'float32');
    fwrite(fid_defvec_bin,dY,'float32');
    fwrite(fid_defvec_bin,dZ,'float32');
    
end

% close files
fid_defvec_txt = fclose(fid_defvec_txt);
fid_defvec_bin = fclose(fid_defvec_bin);