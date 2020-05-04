function data = matRad_readMotionData(fileInfo,percExtTarg)

%% tested properties

% fraction number
fTested         = [19 21];
% marker number
mTested         = [1 3];
% do we need to flip the vertical axis? the positive direction is exhale
flipTested      = [false true];
% this should start in the middle of an inhaling (negative slope) cycle
t_startTested   = [2503 1998.5];

%% get appropriate properties

f = fileInfo.f;
m = fileInfo.m;

ind = fTested == f & mTested == m;

flip    = flipTested(ind);
t_start = t_startTested(ind);

t_tot = fileInfo.t_tot;

t_end = t_start+t_tot;

%% actual function

% read training data
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

% check if we need to flip x-direction for this patient
if flip
    x = -x;
end

ind_start   = round(interp1(t,1:numel(t),t_start));
ind_end     = round(interp1(t,1:numel(t),t_end));

ind_start  = max([ind_start 1]);
ind_end    = min([ind_end numel(t)]);

t_cut = t(ind_start:ind_end);
x_cut = x(ind_start:ind_end);

t_cut = t_cut-t_cut(1);
x_cut = scaleData(x_cut,percExtTarg);

data.t_cut = t_cut;
data.x_cut = x_cut;

end

function xScaled = scaleData(x,percExtTarg)

% we want all of the data scaled between 0 and 1, with the exception of the
% percentage percExtTarg at both extrema

% find cutoffs
cutOff_top      = max(x);
cutOff_bottom   = min(x);

stepSize = (cutOff_top-cutOff_bottom)./1e4;

percAbove = 0;
percBelow = 0;

while percAbove < percExtTarg || percBelow < percExtTarg
    
    if percAbove < percExtTarg
        cutOff_top = cutOff_top-stepSize;
    end
    
    if percBelow < percExtTarg
        cutOff_bottom = cutOff_bottom+stepSize;
    end
    
    percAbove = 100.*nnz(x > cutOff_top)./numel(x);
    percBelow = 100.*nnz(x < cutOff_bottom)./numel(x);
    
end

% scale the data
m = 1.6888./(cutOff_top-cutOff_bottom);
b = -cutOff_bottom.*m;

xScaled = m.*x+b;

end

