% z-profile

fname_full = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run',[fname '.prof0']);
fid = fopen(fname_full);

% first four lines are header
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

i = 1;
z = zeros(101,1);
profile_z = zeros(101,1);
while ischar(tline) && ~isempty(tline)
    
    temp = str2num(tline);
    
    z(i) = temp(1);
    profile_z(i) = temp(4);
    i = i+1;
    
    tline = fgetl(fid);
end

profile_z = profile_z./max(profile_z);

fclose(fid);

% x-profile

fname_full = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run',[fname '.prof1']);
fid = fopen(fname_full);

% first four lines are header
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

i = 1;
x = zeros(101,1);
profile_x = zeros(101,1);
while ischar(tline) && ~isempty(tline)
    
    temp = str2num(tline);
    
    x(i) = temp(1);
    profile_x(i) = temp(4);
    i = i+1;
    
    tline = fgetl(fid);
end

profile_x = profile_x./max(profile_x);

fclose(fid);

% y-profile

fname_full = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run',[fname '.prof2']);
fid = fopen(fname_full);

% first four lines are header
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

i = 1;
y = zeros(101,1);
profile_y = zeros(101,1);
while ischar(tline) && ~isempty(tline)
    
    temp = str2num(tline);
    
    y(i) = temp(1);
    profile_y(i) = temp(4);
    i = i+1;
    
    tline = fgetl(fid);
end

profile_y = profile_y./max(profile_y);

fclose(fid);

figure
plot(z,profile_z)
xlabel('z / cm')
ylabel('rel dose')
savefig([fname ' z_profile'])

figure
plot(x,profile_x)
xlabel('x / cm')
ylabel('rel dose')
savefig([fname ' x_profile'])

figure
plot(y,profile_y)
xlabel('y / cm')
ylabel('rel dose')
savefig([fname ' y_profile'])




