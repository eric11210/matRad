function dij = matRad_loadDij(fnameBase)

%% load metadata first

% get dij path
dijPath = 'E:\matRad\dij';

% construct full fnameBase
fnameBase = fullfile(dijPath,fnameBase);

% save remainder of dij struct
load(fnameBase,'dij');

%% now save actual dij matrix

% get fname
fname = sprintf('%s.dij',fnameBase);

% open dij file for writing
fid = fopen(fname,'r');

% allocate physicalDose
physicalDose = cell(1,dij.numPhases);

% loop through phases
for phase = 1:dij.numPhases
    
    % read total number of non-zero dij elements
    physicalDose_phase_nnz = fread(fid,1,'uint32');
    
    % prepare physicalDose sparse arrays
    physicalDose_phase_i    = zeros(physicalDose_phase_nnz,1);
    physicalDose_phase_j    = zeros(physicalDose_phase_nnz,1);
    physicalDose_phase_s    = zeros(physicalDose_phase_nnz,1);
    
    % set offset to 0
    offset = 0;
    
    % loop through bixels
    for bixel = 1:dij.totalNumOfBixels
        
        % read number of non-zero elements for this bixel
        physicalDose_phase_bixel_numberNonZero = fread(fid,1,'uint32');
        
        % read dij non-zero elements
        physicalDose_phase_bixel_nonZero = fread(fid,physicalDose_phase_bixel_numberNonZero,'double');
        
        % read non-zero voxel indices
        voxelIndices_phase_bixel_nonZero = fread(fid,physicalDose_phase_bixel_numberNonZero,'uint32');
        
        % determine indices for physicalDose sparse arrays
        ind = offset+(1:physicalDose_phase_bixel_numberNonZero);
        
        % update offset
        offset = offset+physicalDose_phase_bixel_numberNonZero;
        
        % enter physicalDose sparse arrays
        physicalDose_phase_i(ind)   = voxelIndices_phase_bixel_nonZero;
        physicalDose_phase_j(ind)   = bixel;
        physicalDose_phase_s(ind)   = physicalDose_phase_bixel_nonZero;
    end
    
    % put dij matrix for phase in physicalDose
    physicalDose{phase} = sparse(physicalDose_phase_i,physicalDose_phase_j,physicalDose_phase_s,dij.numOfVoxels,dij.totalNumOfBixels);
    
end

% close dij file
fclose(fid);

% put physicalDose in dij struct
dij.physicalDose = physicalDose;

end

