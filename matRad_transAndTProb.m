function [Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(qij,initProb,transT,T)

% calculate transition probability matrix and derivative
Pij_transT  = expm(qij.*transT);
Pij_transT_dot = qij*Pij_transT;

% calculate probability to arrive at i at T and derivative
Pi_T    = initProb*expm(qij.*T);
Pi_T_dot   = Pi_T*qij;


end