function [Pij_deltaTSample, Pi_deltaTSample] = generateProbMatrix(l_sample,nSubPhases)


if true
    
    % count number of transitions from one phase to the next
    zij_deltaTSample = zeros(nSubPhases,nSubPhases);
    
    for i = 1:numel(l_sample)
        
        if i == 1
            currPhase = l_sample(i);
            nextPhase = l_sample(i+1);
        elseif i == numel(l_sample)
            %prevPhase = p_sample(i-1);
            currPhase = l_sample(i);
        else
            %prevPhase = p_sample(i-1);
            currPhase = l_sample(i);
            nextPhase = l_sample(i+1);
        end
        
        if i ~= numel(l_sample)
            zij_deltaTSample(currPhase,nextPhase) = zij_deltaTSample(currPhase,nextPhase)+1;
        end
    end
    
    
    % construction of probability transition matrix
    zi_deltaTSample = repmat(sum(zij_deltaTSample,2),1,nSubPhases);
    Pij_deltaTSample = zij_deltaTSample./zi_deltaTSample;
    Pij_deltaTSample(zi_deltaTSample == 0) = 0;
    
    %{
% estimates on variance, look into this
uij_deltaTSample = (zij_deltaTSample+1);
ui_deltaTSample = repmat(sum(uij_deltaTSample,2),1,nSubPhases+2);

varPij_deltaTSample = uij_deltaTSample.*(ui_deltaTSample-uij_deltaTSample)./(ui_deltaTSample.^2.*(ui_deltaTSample+1));

relVarPij_deltaTSample = varPij_deltaTSample./Pij_deltaTSample;
relVarPij_deltaTSample(Pij_deltaTSample == 0) = 0;
    %}
    % construction of probability vector (time-homogeneous)
    Pi_deltaTSample = zeros(nSubPhases,1);
    
    for i = 1:nSubPhases
        Pi_deltaTSample(i) = nnz(l_sample == i)./numel(l_sample);
    end
    
else
    
    % count number of transitions from one phase to the next
    zij_deltaTSample = zeros(nSubPhases+2,nSubPhases+2);
    
    for i = 1:numel(l_sample)
        
        if i == 1
            currPhase = l_sample(i);
            nextPhase = l_sample(i+1);
        elseif i == numel(l_sample)
            %prevPhase = p_sample(i-1);
            currPhase = l_sample(i);
        else
            %prevPhase = p_sample(i-1);
            currPhase = l_sample(i);
            nextPhase = l_sample(i+1);
        end
        
        if i ~= numel(l_sample)
            zij_deltaTSample(currPhase,nextPhase) = zij_deltaTSample(currPhase,nextPhase)+1;
        end
    end
    
    
    % construction of probability transition matrix
    zi_deltaTSample = repmat(sum(zij_deltaTSample,2),1,nSubPhases+2);
    Pij_deltaTSample = zij_deltaTSample./zi_deltaTSample;
    Pij_deltaTSample(zi_deltaTSample == 0) = 0;
    
    %{
% estimates on variance, look into this
uij_deltaTSample = (zij_deltaTSample+1);
ui_deltaTSample = repmat(sum(uij_deltaTSample,2),1,nSubPhases+2);

varPij_deltaTSample = uij_deltaTSample.*(ui_deltaTSample-uij_deltaTSample)./(ui_deltaTSample.^2.*(ui_deltaTSample+1));

relVarPij_deltaTSample = varPij_deltaTSample./Pij_deltaTSample;
relVarPij_deltaTSample(Pij_deltaTSample == 0) = 0;
    %}
    % construction of probability vector (time-homogeneous)
    Pi_deltaTSample = zeros(nSubPhases+2,1);
    
    for i = 1:nSubPhases+2
        Pi_deltaTSample(i) = nnz(l_sample == i)./numel(l_sample);
    end
    
end

end