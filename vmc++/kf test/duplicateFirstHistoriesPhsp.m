function duplicateFirstHistoriesPhsp(numHist,fname_full,fname_dupe)

if nargin < 2
    fname_full = 'C:\Users\eric\Documents\GitHub\matRad\vmc++\run\phsp\5cmx5cm_SSD50cm.egsphsp1';
    
    fname_dupe = sprintf('C:\\Users\\eric\\Documents\\GitHub\\matRad\\vmc++\\run\\phsp\\5cmx5cm_SSD50cm_%dx2.egsphsp1',numHist);
end

% open full phsp
fid_full = fopen(fname_full,'r');

% get header of full phsp
[fid_full, header_full] = getHeader(fid_full);
mode = char(header_full.MODE_RW(5));

% initialize header, will change it later
header_dupe             = header_full;
header_dupe.NPPHSP      = 0;
header_dupe.NPHOTPHSP   = 0;
header_dupe.EKMAXPHSP   = 0;
header_dupe.EKMINPHSPE  = 1000;

% write duped header
writeHeader(fname_dupe,header_dupe,0);

% initialize buffer for all records to be duped
buffer_dupe = cell(numHist,1);

% now loop until we have found numHist photons
notFinished = true;
while notFinished
    
    % extract record
    [fid_full, record] = getRecord(fid_full,mode);
    
    % is this the first particle scored from a new primary history?
    if record.E < 0
        
        firstParticle = true;
    end
    
    % modify max/min energies, increment number of photons
    % must determine particle type using LATCH
    LATCH = de2bi(record.LATCH,32);
    if LATCH(30:31) == [0 0]
        % photon
        header_dupe.EKMAXPHSP   = max(header_dupe.EKMAXPHSP,abs(record.E));
        header_dupe.NPHOTPHSP   = header_dupe.NPHOTPHSP+1;
        
    elseif LATCH(30:31) == [0 1]
        % electron
        header_dupe.EKMINPHSPE  = min(header_dupe.EKMINPHSPE,abs(record.E)-0.511);
        header_dupe.EKMAXPHSP   = max(header_dupe.EKMAXPHSP,abs(record.E)-0.511);
        
        % IF ELECTRON, CONTINUE
        continue
        
    elseif LATCH(30:31) == [1 0]
        % positron
        header_dupe.EKMINPHSPE  = min(header_dupe.EKMINPHSPE,abs(record.E)-0.511);
        header_dupe.EKMAXPHSP   = max(header_dupe.EKMAXPHSP,abs(record.E)-0.511);
        
        % IF POSITRON, CONTINUE
        continue
        
    elseif LATCH(30:31) == [1 1]
        error('Electron and positron???')
        
    end
    
    % we should only be looking at photons now
    
    % increment number of particles
    header_dupe.NPPHSP = header_dupe.NPPHSP+1;
    
    if firstParticle
        % this is the first from a new primary history
        % make the energy negative
        record.E = -abs(record.E);
        % make sure next particle is not negative
        firstParticle = false;
    else
        % these particles should already have positive energy
        if record.E < 0
            % SHOULDN'T HAPPEN
            warning('NEGATIVE ENERGY')
        end
    end
    
    % put record in buffer
    buffer_dupe{header_dupe.NPPHSP} = record;
    
    % we are finished when header_dupe.NPPHSP == numHist
    notFinished = header_dupe.NPPHSP ~= numHist;
    
    % display progress
    if mod(header_dupe.NPPHSP,max(1,round(numHist/200))) == 0
        matRad_progress(header_dupe.NPPHSP/max(1,round(numHist/200)),...
            floor(numHist/max(1,round(numHist/200))));
    end
    
    
end

% close full phsp
fclose(fid_full);

% now write buffer twice
writeBuffer(fname_dupe,buffer_dupe,mode);
writeBuffer(fname_dupe,buffer_dupe,mode);

% double number of particles in duped phsp
header_dupe.NPPHSP      = 2.*header_dupe.NPPHSP;
header_dupe.NPHOTPHSP   = 2.*header_dupe.NPHOTPHSP;

% clean up duped phsp header
if header_dupe.EKMINPHSPE == 1000
    header_dupe.EKMINPHSPE = 0;
end

writeHeader(fname_dupe,header_dupe,1);

end

% read/write functions

function [fid, header] = getHeader(fid)

header.MODE_RW      = fread(fid,5,'uint8');
header.NPPHSP       = fread(fid,1,'int32');
header.NPHOTPHSP    = fread(fid,1,'int32');
header.EKMAXPHSP    = fread(fid,1,'float32');
header.EKMINPHSPE   = fread(fid,1,'float32');
header.NINCPHSP     = fread(fid,1,'float32');
header.garbage      = fread(fid,3,'int8');

end

function [fid, record] = getRecord(fid,mode)

record.LATCH    = fread(fid,1,'uint32');
record.E        = fread(fid,1,'float32');
record.X        = fread(fid,1,'float32');
record.Y        = fread(fid,1,'float32');
record.U        = fread(fid,1,'float32');
record.V        = fread(fid,1,'float32');
record.WT       = fread(fid,1,'float32');

if mode == 2
    record.ZLAST    = fread(fid,1,'float32');
end

end

function writeHeader(fname,header,rewrite)

if rewrite
    fid = fopen(fname,'r+');
else
    fid = fopen(fname,'w');
end

fseek(fid,0,'bof');

fwrite(fid,header.MODE_RW,'uint8');
fwrite(fid,header.NPPHSP,'int32');
fwrite(fid,header.NPHOTPHSP,'int32');
fwrite(fid,header.EKMAXPHSP,'float32');
fwrite(fid,header.EKMINPHSPE,'float32');
fwrite(fid,header.NINCPHSP,'float32');
fwrite(fid,header.garbage,'int8');

fclose(fid);

end

function writeRecord(fid,record,mode)

fwrite(fid,record.LATCH,'uint32');
fwrite(fid,record.E,'float32');
fwrite(fid,record.X,'float32');
fwrite(fid,record.Y,'float32');
fwrite(fid,record.U,'float32');
fwrite(fid,record.V,'float32');
fwrite(fid,record.WT,'float32');

if mode == 2
    fwrite(fid,record.ZLAST,'float32');
end

end

function writeBuffer(fname,buffer,mode)

fid = fopen(fname,'r+');
fseek(fid,0,'eof');

numRecords = sum(~cellfun(@isempty,buffer));

for i = 1:numRecords
    writeRecord(fid,buffer{i},mode);
end

fclose(fid);

end