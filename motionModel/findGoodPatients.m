m = 1;

fStart  = 22;
fEnd    = 121;
fVec    = fStart:fEnd;

for f = fVec
    
    fname_wild = sprintf('p*_f%d_m%d_proc2_pca.csv',f,m);
    path = 'ck_processed';
    matchedFiles = dir(fullfile(fileparts(mfilename('fullpath')),path,fname_wild));
    
    if size(matchedFiles,1) > 1
        error('More than one patient with this fraction and marker number.');
    elseif size(matchedFiles,1) == 0
        error('No patient with this fraction and marker number.');
    end
    
    contents = dlmread(fullfile(fileparts(mfilename('fullpath')),path,matchedFiles.name));
    
    t = contents(:,1);
    x = contents(:,2);
    x = -x; % flip x-direction for this patient
    
    plot(t,x)
    title(sprintf('fraction %d',f));
    waitforbuttonpress;
    
end