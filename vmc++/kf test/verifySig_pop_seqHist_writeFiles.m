%% setup

nHist1 = 21;
nHist2 = 50;

%% get base file contents

baseFileName = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\kf test\','kfMethod_getSig_1million.vmc');

fid = fopen(baseFileName,'r+');

i = 1;
tline = fgetl(fid);
baseFileContents{i} = tline;

while ischar(tline)
    % loop through file
    
    if contains(tline,'ncase')
        % find the line that has ncase
        
        ncaseLine = i;
    end
    
    i = i+1;
    tline = fgetl(fid);
    baseFileContents{i} = tline;
    
end

fclose(fid);


%% copy base file to new file, then modify new file

for nHist = nHist1:nHist2
    
    newFileName = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\run\',sprintf('kfMethod_getSig_%dmillion.vmc',nHist));
    
    newFileContents = baseFileContents;
    
    newFileContents{ncaseLine} = sprintf('ncase = %d',nHist.*10^6);
    
    fid = fopen(newFileName,'w');
    
    for i = 1:numel(newFileContents)
        
        if newFileContents{i} ~= -1
            fprintf(fid,'%s\n',newFileContents{i});
        end
    end
    
    fclose(fid);
    
end

%% write batch file

batchFileName = fullfile('C:\Users\eric\Documents\GitHub\matRad\vmc++\','kfMethod_getSig_allmillion.bat');

verboseString = '';

parallelProcesses = cell(1,nHist2-nHist1+1);
i = 1;
for nHist = nHist1:nHist2
    parallelProcesses{1,i} = sprintf('start "" 9>"%%lock%%%d" %s .\\bin\\vmc_Windows.exe -i kfMethod_getSig_%dmillion',i,verboseString,nHist);
    i = i+1;
end

batchFile = {...
    ['@echo off'],...
    ['setlocal'],...
    ['set "lock=%temp%\wait%random%.lock"'],...
    [''],...
    parallelProcesses{:},...
    [''],...
    [':Wait for all processes to finish (wait until lock files are no longer locked)'],...
    ['1>nul 2>nul ping /n 2 ::1'],...
    ['for %%N in (',strjoin(arrayfun(@(x) num2str(x),(1:(nHist2-nHist1+1)),'UniformOutput',false),' '),') do ('],...
    ['  (call ) 9>"%lock%%%N" || goto :Wait'],...
    [') 2>nul'],...
    [''],...
    ['del "%lock%*"'],...
    ['']...
    %,['echo Done - ready to continue processing']
    };


fid = fopen(batchFileName,'wt');
for i = 1:length(batchFile)
  fprintf(fid,'%s\n',batchFile{i});
end
fclose(fid);
