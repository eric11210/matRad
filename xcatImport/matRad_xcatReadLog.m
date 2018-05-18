function xcatLog = matRad_xcatReadLog(fnameLog)

fid = fopen(fnameLog,'r');

tline = fgetl(fid);
tline = fgetl(fid);

muPixel = false;

while(ischar(tline))
    
    if ~isempty(strfind(tline,'Total Output Frames'))
        
        equalInd = find(tline == '=');
        numPhases = str2double(tline((equalInd+1):end));
        
    elseif ~isempty(strfind(tline,'Time Per Frame'))
        
        equalInd = find(tline == '=');
        sInd = find(tline == 's');
        deltaT = str2double(tline((equalInd(3)+1):(sInd(1)-1)));
    
    elseif ~isempty(strfind(tline,'pixel width')) && isempty(strfind(tline,'voxel volume'))
        
        equalInd = find(tline == '=');
        parenthInd = find(tline == '(');
        resolution.x = str2double(tline((equalInd+1):(parenthInd-1)))*10; % mm
        resolution.y = str2double(tline((equalInd+1):(parenthInd-1)))*10; % mm
        
    elseif ~isempty(strfind(tline,'slice width')) && isempty(strfind(tline,'voxel volume'))
        
        equalInd = find(tline == '=');
        parenthInd = find(tline == '(');
        resolution.z = str2double(tline((equalInd+1):(parenthInd-1)))*10; % mm
        
    elseif ~isempty(strfind(tline,'array_size'))
        
        equalInd = find(tline == '=');
        dim.x = str2double(tline((equalInd+1):end));
        dim.y = str2double(tline((equalInd+1):end));
        
    elseif ~isempty(strfind(tline,'starting slice number'))
        
        equalInd = find(tline == '=');
        startingSlice = str2double(tline((equalInd+1):end));
        
    elseif ~isempty(strfind(tline,'ending slice number'))
        
        equalInd = find(tline == '=');
        endingSlice = str2double(tline((equalInd+1):end));
        
    elseif ~isempty(strfind(tline,'Linear Attenuation Coefficients (1/pixel):'))
        
        % next lines are the mu's in pixel units
        muPixel = true;
        
    end
    
    if muPixel && ~isempty(strfind(tline,'Body (water)'))
        
        equalInd = find(tline == '=');
        muWater = str2double(tline((equalInd+1):end));
        
    end
    
    
    tline = fgetl(fid);
end

dim.z = endingSlice-startingSlice+1;


xcatLog.dim = dim;
xcatLog.resolution = resolution;
xcatLog.muWater = muWater;
xcatLog.numFrames = numPhases;
xcatLog.deltaT = deltaT;

fclose(fid);

end

