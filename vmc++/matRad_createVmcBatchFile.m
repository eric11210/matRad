function matRad_createVmcBatchFile(parallelSimulations,filepath,inOutName_base,verboseLevel)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad batchfile creation
%
% call
%   matRad_createVmcBatchFile(parallelSimulations,filepath,verboseLevel)
%
% input
%   parallelSimulations: no of parallel simulations
%
%   filepath:            path where batchfile is created (this has to be the 
%                        path of the vmc++ folder)
%
%   inOutName_base:      base filename for the input/output files
%
%   verboseLevel:        optional. number specifying the amount of output 
%                        printed to the command prompt
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set verbose level
if nargin < 3
    verboseString = '/B'; % to not open terminals by default
else
    if verboseLevel > 1 % open terminals is verboseLevel > 1
        verboseString = '';
    else
        verboseString = '/B';
    end        
end

if ispc % parallelization only possible on windows systems

    %parallelProcesses = cell(1,2*parallelSimulations);
    parallelProcesses = cell(1,parallelSimulations);
    for i = 1:parallelSimulations
        parallelProcesses{1,i} = sprintf('start "" 9>"%%lock%%%d" %s .\\bin\\vmc_Windows.exe -i %s_%d',i,verboseString,inOutName_base,i);
        %parallelProcesses{1,2*i} = 'timeout /T 2';
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
        ['for %%N in (',strjoin(arrayfun(@(x) num2str(x),(1:parallelSimulations),'UniformOutput',false),' '),') do ('],...
        ['  (call ) 9>"%lock%%%N" || goto :Wait'],...
        [') 2>nul'],...
        [''],...
        ['del "%lock%*"'],...
        ['']...
        ['exit']...
        %['copy NUL EmptyFile.txt']...
        %,['echo Done - ready to continue processing']
        };
    
elseif isunix
    
    batchFile = {sprintf('./bin/vmc_Linux.exe %s_1',inOutName_base)};

end

% write batch file
fid = fopen(filepath,'wt');
for i = 1 : length(batchFile)
  fprintf(fid,'%s\n',batchFile{i});
end
fclose(fid);

if isunix
   system(['chmod a+x ' filepath]);
end

end
