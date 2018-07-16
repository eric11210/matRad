function [Pij_mean,Pij_mean_std,Pij_std,Pij_std_std,Pi_mean,Pi_std] = stdProbMatrix(Pij_sample,nSubPhases,nHistories,nSteps_sample,l_init)

Pij_simulatedArray = zeros(nSubPhases,nSubPhases,nHistories);
Pi_simulatedArray = zeros(nSubPhases,nHistories);

fprintf('Monte Carlo estimation.\n');
for history = 1:nHistories
    
    % simulate a sequence of the Markov chain with Monte Carlo
    l_simulated = MCMC(Pij_sample,nSteps_sample,l_init);
    
    % extra probability matrices from the simulated sequence
    [Pij_simulated, Pi_simulated] = generateProbMatrix(l_simulated,nSubPhases);
    
    % insert in arrays
    Pij_simulatedArray(:,:,history) = Pij_simulated;
    Pi_simulatedArray(:,history) = Pi_simulated;
    
    matRad_progress(history,nHistories);
end

% calculate mean and standard deviation of simulated matrices
Pij_mean = mean(Pij_simulatedArray,3);
Pij_std = std(Pij_simulatedArray,0,3);
Pij_mean_std = Pij_std./sqrt(nHistories);
Pij_mu4 = sum((Pij_simulatedArray-repmat(Pij_mean,1,1,nHistories)).^4,3)./(nHistories-1);
Pij_std_std = sqrt((Pij_mu4-(nHistories-3).*(Pij_std.^4)./(nHistories-1))./nHistories);


Pi_mean = mean(Pi_simulatedArray,2);
Pi_std = std(Pi_simulatedArray,0,2);


end