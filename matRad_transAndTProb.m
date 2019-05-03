function [Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(qij,initProb,transT,T)

% calculate transition probability matrix and derivative
Pij_transT      = expm(qij.*transT);
Pij_transT_dot  = qij*Pij_transT;

% calculate probability to arrive at i at T and derivative
Pi_T        = initProb*expm(qij.*T);
Pi_T_dot    = Pi_T*qij;

% [subPhase2PosPhase_gridJ, subPhase2PosPhase_gridI] = meshgrid(subPhase2PosPhase);

% norm = repmat(accumarray(subPhase2PosPhase_gridI(:),Pij(:)),[1 nPosPhases];

% PIJ = accumarray([subPhase2PosPhase_gridI(:) subPhase2PosPhase_gridJ(:)],Pij(:))./norm;

% PIJ(isnan(PIJ)) = 0;

% PI = accumarray(subPhase2PosPhase,Pi);

end