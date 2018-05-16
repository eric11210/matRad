function [tumourPos,t] = matRad_getTumourMotion(fnameXcatRoot)

dirXCAT = fullfile(fileparts(mfilename('fullpath')),'XCAT',filesep);

fnameXcatLog = fullfile(dirXCAT,sprintf('%s_log',fnameXcatRoot));
xcatLog = matRad_xcatReadLog(fnameXcatLog);

fnameXcatPar = fullfile(dirXCAT,sprintf('%s.par',fnameXcatRoot));
tumourPosInit = matRad_xcatReadPar(fnameXcatPar);

tumourPos = zeros(xcatLog.numPhases,3);
tumourPos(1,:) = tumourPosInit;

for phase = 2:xcatLog.numPhases
    
    fprintf('\nPhase %d of %d.\n',phase,xcatLog.numPhases);
    
    tumourPos(phase,:) = matRad_transformVecMVF(tumourPosInit,phase,xcatLog,fnameXcatRoot);
end
% add extra phase at the end to make it periodic
tumourPos(phase+1,:) = tumourPos(1,:);

t = xcatLog.deltaT.*((1:(xcatLog.numPhases+1))-1);


end