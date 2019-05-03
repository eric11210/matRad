function xcatLog = matRad_xcatReadLog(fnameLog)

fid = fopen(fnameLog,'r');

tline = fgetl(fid);
tline = fgetl(fid);

muPixel = false;

while(ischar(tline))
    
    if contains(tline,'Total Output Frames')
        
        equalInd = find(tline == '=');
        numPhases = str2double(tline((equalInd+1):end));
        
    elseif contains(tline,'Time Per Frame')
        
        equalInd = find(tline == '=');
        sInd = find(tline == 's');
        deltaT = str2double(tline((equalInd(3)+1):(sInd(1)-1)));
    
    elseif contains(tline,'pixel width') && ~contains(tline,'voxel volume')
        
        equalInd = find(tline == '=');
        parenthInd = find(tline == '(');
        resolution.x = str2double(tline((equalInd+1):(parenthInd-1)))*10; % mm
        resolution.y = str2double(tline((equalInd+1):(parenthInd-1)))*10; % mm
        
    elseif contains(tline,'slice width') && ~contains(tline,'voxel volume')
        
        equalInd = find(tline == '=');
        parenthInd = find(tline == '(');
        resolution.z = str2double(tline((equalInd+1):(parenthInd-1)))*10; % mm
        
    elseif contains(tline,'array_size')
        
        equalInd = find(tline == '=');
        dim.x = str2double(tline((equalInd+1):end));
        dim.y = str2double(tline((equalInd+1):end));
        
    elseif contains(tline,'starting slice number')
        
        equalInd = find(tline == '=');
        startingSlice = str2double(tline((equalInd+1):end));
        
    elseif contains(tline,'ending slice number')
        
        equalInd = find(tline == '=');
        endingSlice = str2double(tline((equalInd+1):end));
        
    elseif contains(tline,'Linear Attenuation Coefficients (1/pixel):')
        
        % next lines are the mu's in pixel units
        muPixel = true;
        
    end
    
    if muPixel && contains(tline,'Body (water)')
        
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

