function data = matRad_readMotionData(fileInfo)

p   = fileInfo.p;
f   = fileInfo.f;
m   = fileInfo.m;

if fileInfo.processed
    fname = sprintf('p%d_f%d_m%d_proc2_pca.csv',p,f,m);
    path = 'ck_processed';
else
    fname = sprintf('p%d_f%d_m%d.csv',p,f,m);
    path = 'ck_original';
end

contents = dlmread(fullfile(fileparts(mfilename('fullpath')),path,fname));

t = contents(:,1);
x = contents(:,2);

t_start = t(57135);
t_tot   = 60.*5;
t_end   = t_start+t_tot;

ind_start   = round(interp1(t,1:numel(t),t_start));
ind_end     = round(interp1(t,1:numel(t),t_end));

ind_start   = max([ind_start 1]);
ind_end     = min([ind_end numel(t)]);

data.t_cut = t(ind_start:ind_end);
data.x_cut = x(ind_start:ind_end);

end

