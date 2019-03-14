function data = matRad_readMotionData(fileInfo)

p   = fileInfo.p;
f   = fileInfo.f;
m   = fileInfo.m;

if fileInfo.processed
    fname_wild = sprintf('p*_f%d_m%d_proc2_pca.csv',f,m);
    path = 'ck_processed';
else
    fname_wild = sprintf('p*_f%d_m%d.csv',f,m);
    path = 'ck_original';
end

matchedFiles = dir(fullfile(fileparts(mfilename('fullpath')),path,fname_wild));

if size(matchedFiles,1) > 1
    error('More than one patient with this fraction and marker number.');
elseif size(matchedFiles,1) == 0
    error('No patient with this fraction and marker number.');
end

contents = dlmread(fullfile(fileparts(mfilename('fullpath')),path,matchedFiles.name));

t = contents(:,1);
x = contents(:,2);
%x = -x; % flip x-direction for this patient

t_start = 2503; % 2151 start at inhale for f = 42
t_tot   = 60.*10;
t_end   = t_start+t_tot;

ind_start   = round(interp1(t,1:numel(t),t_start));
ind_end     = round(interp1(t,1:numel(t),t_end));

ind_start   = max([ind_start 1]);
ind_end     = min([ind_end numel(t)]);

data.t_cut = t(ind_start:ind_end);
data.x_cut = x(ind_start:ind_end);

% trim up to first exhale, use min/peakfinder?

end

