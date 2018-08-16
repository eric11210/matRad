function ct = matRad_binFrames2Phases(ct,importOptions)

if importOptions.binFrames2Phases
    % get data from ct
    tumourPos = ct.tumourMotion.coordsVox;
    t = ct.tumourMotion.t;
    
    % extract principal component
    [~,components] = pca(tumourPos);
    x_XCAT = components(:,1);
    
    % get boundaries on bins
    lBoundsMax = max(x_XCAT);
    lBoundsMin = min(x_XCAT);
    nBins = importOptions.numPhases/2;
    lBounds = linspace(lBoundsMax,lBoundsMin,nBins+1);
    
    % preliminary binning of data
    l_XCAT = zeros(size(x_XCAT));
    for bin = 1:nBins
        l_XCAT(lBounds(bin+1) <= x_XCAT & x_XCAT <= lBounds(bin)) = bin;
    end
    
    % include inhale/exhale information:
    % breathing in is the beginning of the cycle, so it should have
    % lower phase number
    % breathing out is the end of the cycle, so it should have higher phase
    % number
    
    % define breathing in/out by peaks: points after a min and before a max are
    % exhale, points after a max and before a min are inhale.
    %[~,ind_maxPeaks] = findpeaks(x_XCAT,'MinPeakProminence',0.5);
    [~,ind_minPeaks] = findpeaks(-x_XCAT,'MinPeakProminence',0.5);
    
    l_XCAT(ind_minPeaks:end) = importOptions.numPhases+1-l_XCAT(ind_minPeaks:end);
    
    % find time at which mean x value of each phase is reached
    x_l = zeros(importOptions.numPhases,1);
    t_x_l = zeros(importOptions.numPhases,1);
    t_x_l_XCAT = zeros(importOptions.numPhases,1);
    ind_x_l_XCAT = zeros(importOptions.numPhases,1);
    for phase = 1:importOptions.numPhases
        x_l(phase) = mean(x_XCAT(l_XCAT == phase));
        t_x_l(phase) = interp1(x_XCAT(l_XCAT == phase),t(l_XCAT == phase),x_l(phase));
        t_x_l_XCAT(phase) = interp1(t,t,t_x_l(phase),'nearest','extrap');
        ind_x_l_XCAT(phase) = find(t == t_x_l_XCAT(phase));
    end
    
    
    
    
    ct.tumourMotion.frames2Phases = l_XCAT(1:(end-1));
    ct.tumourMotion.nFramesPerPhase = accumarray(ct.tumourMotion.frames2Phases,1);
    %ct.tumourMotion.nFramesPerPhase = repelem(ct.tumourMotion.nFramesPerPhase,ct.tumourMotion.nFramesPerPhase);
    ct.tumourMotion.numFrames = numel(l_XCAT)-1;
    ct.tumourMotion.numPhases = importOptions.numPhases;
    
    
    ct.cube_new = cell(importOptions.numPhases,1);
    ct.cubeHU_new = cell(importOptions.numPhases,1);
    ct.motionVecX_new = cell(importOptions.numPhases,1);
    ct.motionVecY_new = cell(importOptions.numPhases,1);
    ct.motionVecZ_new = cell(importOptions.numPhases,1);
    
    for phase = 1:importOptions.numPhases
        % create new cubes for each phase, that are the frame snapshot
        % at the representative time t_x_l
        ct.cube_new{phase} = ct.cube{ind_x_l_XCAT(phase)};
        ct.cubeHU_new{phase} = ct.cubeHU{ind_x_l_XCAT(phase)};
        
    end
    
    if ~importOptions.keepAllFrames
        if importOptions.averageCT || importOptions.averageMVF
            % overwrite new cubes if we want to do an average instead
            for frame = 1:ct.tumourMotion.numFrames
                phase = ct.tumourMotion.frames2Phases(frame);
                
                if importOptions.averageCT
                    if find(ct.tumourMotion.frames2Phases == phase,1,'first') == frame
                        ct.cube_new{phase} = zeros(ct.cubeDim);
                        ct.cubeHU_new{phase} = zeros(ct.cubeDim);
                    end
                    
                    ct.cube_new{phase} = ct.cube_new{phase}+ct.cube{frame}./ct.tumourMotion.nFramesPerPhase(phase);
                    ct.cubeHU_new{phase} = ct.cubeHU_new{phase}+ct.cubeHU{frame}./ct.tumourMotion.nFramesPerPhase(phase);
                end
                
                if importOptions.averageMVF
                    
                    if find(ct.tumourMotion.frames2Phases == phase,1,'first') == frame
                        ct.motionVecX_new{phase} = zeros(ct.cubeDim);
                        ct.motionVecY_new{phase} = zeros(ct.cubeDim);
                        ct.motionVecZ_new{phase} = zeros(ct.cubeDim);
                    end
                    
                    ct.motionVecX_new{phase} = ct.motionVecX_new{phase}+ct.motionVecX{frame}./ct.tumourMotion.nFramesPerPhase(phase);
                    ct.motionVecY_new{phase} = ct.motionVecY_new{phase}+ct.motionVecY{frame}./ct.tumourMotion.nFramesPerPhase(phase);
                    ct.motionVecZ_new{phase} = ct.motionVecZ_new{phase}+ct.motionVecZ{frame}./ct.tumourMotion.nFramesPerPhase(phase);
                end
                
            end
        end
        
        if ~importOptions.averageCT
            fprintf('\nMake sure to re-run XCAT option = 0 at the new times to get the correct CTs!\n');
        end
        if ~importOptions.averageMVF
            fprintf('\nMake sure to re-run XCAT option = 4 at the new times to get the correct MVFs!\n');
        end
        
        % overwrite old ct info with new
        ct.cube         = ct.cube_new;
        ct.cubeHU       = ct.cubeHU_new;
        ct.motionVecX   = ct.motionVecX_new;
        ct.motionVecY   = ct.motionVecY_new;
        ct.motionVecZ   = ct.motionVecZ_new;
        
    end
    
    % now delete _new fields
    ct = rmfield(ct,{'cube_new' 'cubeHU_new' 'motionVecX_new' 'motionVecY_new' 'motionVecZ_new'});
    
    ct.tumourMotion.tPhase = t_x_l;
    ct.tumourMotion.indPhase = ind_x_l_XCAT;
    ct.tumourMotion.tPhaseRounded = t_x_l_XCAT;
    
    ct.tumourMotion.XCATPar.out_period = 2.*ct.tumourMotion.tPhase;
    %ct.tumourMotion.XCATPar.hrt_start_ph_index = ct.tumourMotion.tPhase(1)./1; % the 1 is the hrt_period
    %ct.tumourMotion.XCATPar.resp_start_ph_index = ct.tumourMotion.tPhase(1)./5; % the 1 is the resp_period
else
    
    ct.tumourMotion.frames2Phases = (1:importOptions.numPhases)';
    ct.tumourMotion.nFramesPerPhase = accumarray(ct.tumourMotion.frames2Phases,1);
    ct.tumourMotion.numFrames = importOptions.numPhases;
    ct.tumourMotion.numPhases = importOptions.numPhases;
end

end