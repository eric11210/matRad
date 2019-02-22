function [fmTStart,fmIndStart,fmS,xKnot] = matRad_FSMlocalPruning(fmTStart,fmIndStart,fmS,xKnot,velocityRangeS,dV)

% update beginning times, indices, statesof the most
% recent lines, starting with the most recent one
% also knots
fmTStart    = circshift(fmTStart,-1);
fmIndStart  = circshift(fmIndStart,-1);
fmS         = circshift(fmS,-1);
xKnot       = circshift(xKnot,-1);


% check if current state is irregular; we might be able to go right into
% regular breathing if IRR was a mistake
if fmS(1) == 4 && (velocityRangeS(1,1) < dV && dV < velocityRangeS(1,2) || velocityRangeS(2,1) < dV && dV < velocityRangeS(2,2))
    
    %% enter regular breathing
    
    % now find out which state it is
    
    % first inhale state
    if velocityRangeS(1,1) < dV && dV < velocityRangeS(1,2)
        
        % set next state
        fSNext = 1;
    end
    
    % then exhale state
    if velocityRangeS(2,1) < dV && dV < velocityRangeS(2,2)
        
        % set next state
        fSNext = 2;
    end
    
    fmS(1) = fSNext;
    
    % now that we've rearranged the states, check again to see if
    % states are the same
    if fmS(1) == fmS(2) && fmIndStart(2) > 1
        
        % update beginning times, indices, statesof the most
        % recent lines, starting with the most recent one
        % also knots
        fmTStart    = circshift(fmTStart,-1);
        fmIndStart  = circshift(fmIndStart,-1);
        fmS         = circshift(fmS,-1);
        xKnot       = circshift(xKnot,-1);
        
    end
    
end


end

