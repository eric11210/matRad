function [tPos,tRad] = matRad_xcatReadPar(fnameLog)

fid = fopen(fnameLog,'r');

tline = fgetl(fid);

tPos = zeros(1,3);

while(ischar(tline))
    
    if contains(tline,'x_location')
        
        equalInd    = find(tline == '=');
        hashInd     = find(tline == '#');
        tPos(1)     = str2double(tline((equalInd(1)+1):(hashInd(1)-1)));
    
    elseif contains(tline,'y_location')
        
        equalInd    = find(tline == '=');
        hashInd     = find(tline == '#');
        tPos(2)     = str2double(tline((equalInd(1)+1):(hashInd(1)-1)));
        
    elseif contains(tline,'z_location')
        
        equalInd    = find(tline == '=');
        hashInd     = find(tline == '#');
        tPos(3)     = str2double(tline((equalInd(1)+1):(hashInd(1)-1)));
        
    elseif contains(tline,'lesn_diameter')
        
        equalInd    = find(tline == '=');
        hashInd     = find(tline == '#');
        tRad        = str2double(tline((equalInd(1)+1):(hashInd(1)-1)))/2;
        
    end
    
    tline = fgetl(fid);
end

fclose(fid);

tPos = tPos+1;

end