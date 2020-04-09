function matRad_saveDij(dij,fnameBase)

%% save metadata first

% get dij path
dijPath = 'E:\matRad\dij';

% construct full fnameBase
fnameBase = fullfile(dijPath,fnameBase);

% get physicalDose
physicalDose = dij.physicalDose;

% delete physicalDose field from dij
dij = rmfield(dij,'physicalDose');

% save remainder of dij struct
save(fnameBase,'dij');

%% now save actual dij matrix

% get fname
fname = sprintf('%s.dij',fnameBase);

% open dij file for writing
fid = fopen(fname,'w');

% get voxel indices
voxelIndices = 1:dij.numOfVoxels;

% loop through phases
for phase = 1:dij.numPhases
    
    % get dij matrix for phase
    physicalDose_phase = physicalDose{phase};
    
    % get total number of non-zero dij elements
    physicalDose_phase_nnz = nnz(physicalDose_phase);
    
    % save total number of non-zero dij elements
    fwrite(fid,physicalDose_phase_nnz,'uint32');
    
    % loop through bixels
    for bixel = 1:dij.totalNumOfBixels
        
        % get dij column for particular phase and bixel
        physicalDose_phase_bixel = physicalDose_phase(:,bixel);
        
        % find non-zedro dij elements
        nonZero = physicalDose_phase_bixel > 0;
        
        % get non-zero dij elements
        physicalDose_phase_bixel_nonZero = full(physicalDose_phase_bixel(nonZero));
        
        % get non-zero voxel indices
        voxelIndices_phase_bixel_nonZero = voxelIndices(nonZero);
        
        % get number of non-zero dij elements
        physicalDose_phase_bixel_numberNonZero = numel(physicalDose_phase_bixel_nonZero);
        
        % write number of non-zero elements for this bixel
        fwrite(fid,physicalDose_phase_bixel_numberNonZero,'uint32');
        
        % write dij non-zero elements
        fwrite(fid,physicalDose_phase_bixel_nonZero,'double');
        
        % write voxel non-zero indices
        fwrite(fid,voxelIndices_phase_bixel_nonZero,'uint32');
    end
end

% close dij file
fclose(fid);

end

