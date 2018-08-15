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

data.t_cut = t(57135:end);
data.x_cut = x(57135:end);

end

