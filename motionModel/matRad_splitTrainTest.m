function [FS_sample,t_sample,x_sample,v_sample,l_sample,indFirstCycle_train,indLastCycle_train,indFirstCycle_test,indLastCycle_test] = matRad_splitTrainTest(FS_sample,t_sample,x_sample,v_sample,l_sample,trainRatio)
% find the starting and stopping indices corresponding to full breathing
% cycles

% also split data into training and testing

% both training and testing data should have only full cycles

% default training ratio is 0.5
if nargin < 6
    trainRatio = 0.5;
end

%% first determine overall starting and stopping indices and cut

% determine the overall starting index
% this is defined by the first time the patient inhales AFTER an EOE
indFirstEOE     = find(FS_sample == 3,1,'first');
indFirstCycle   = find(FS_sample(indFirstEOE:end) == 1,1,'first')+indFirstEOE-1;

% determine the overall stopping index
% this is defined by the last time the patient EOEs AFTER an inhale
indLastInh      = find(FS_sample == 1,1,'last');
indLastCycle    = find(FS_sample(1:indLastInh) == 3,1,'last');

% now cut
t_sample((indLastCycle+1):end)  = [];
x_sample((indLastCycle+1):end)  = [];
v_sample((indLastCycle+1):end)  = [];
l_sample((indLastCycle+1):end)  = [];
FS_sample((indLastCycle+1):end) = [];

t_sample(1:(indFirstCycle-1))   = [];
x_sample(1:(indFirstCycle-1))   = [];
v_sample(1:(indFirstCycle-1))   = [];
l_sample(1:(indFirstCycle-1))   = [];
FS_sample(1:(indFirstCycle-1))  = [];

%% now split training data

% determine number of points in the training and testing data
numTrainPts = round(numel(FS_sample).*trainRatio);
%numTestPts  = numel(FS_sample)-numTrainPts;

% the first part of the data will be used for training, the second for
% testing
% split it now
FS_sample_train = FS_sample(1:numTrainPts);
FS_sample_test  = FS_sample((numTrainPts+1):end);

% find the stopping index for the training data
indLastInh_train    = find(FS_sample_train == 1,1,'last');
indLastCycle_train  = find(FS_sample_train(1:indLastInh_train) == 3,1,'last');

% now find the starting index for the testing data
indFirstEOE_test    = find(FS_sample_test == 3,1,'first');
indFirstCycle_test  = find(FS_sample_test(indFirstEOE_test:end) == 1,1,'first')+indFirstEOE_test+numTrainPts-1;

% the starting index for the training data is just 1
indFirstCycle_train = 1;

% likewise, the stopping index for the testing data is just the number of
% elements in the data
indLastCycle_test = numel(FS_sample);

end

