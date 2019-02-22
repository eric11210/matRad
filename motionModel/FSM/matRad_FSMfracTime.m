function FS_sample = matRad_FSMfracTime(FS_sample,deltaT_sample,options)

%% find maximum time spent in each state

% create vector containing maximum time spent in each state
maxTime = zeros(4,1);

% initialize time spent in current state
timeInState = 1;

% loop through entire trace
for i = 1:(numel(FS_sample)-1)
    
    % determine current and next states
    currState = FS_sample(i);
    nextState = FS_sample(i+1);
    
    % determine if two states are equal
    if currState == nextState
        % they are equal
        
        % increment currTime by 1
        timeInState = timeInState+1;
    else
        % they are not equal
        
        % replace maxTime of current state by new max
        maxTime(currState) = max([maxTime(currState) timeInState]);
        
        % reset currTime to 1
        timeInState = 1;
    end
    
end

%% split states into fractional times

% determine bounds on regular states
% all IRR are put into the same time fraction
fracTimeBounds = zeros(3,options.nTimeFracs+1);

for state = 1:3
    
    fracTimeBounds(state,:) = linspace(0,maxTime(state),options.nTimeFracs+1)+1;
end
fracTimeBoundsL = fracTimeBounds(:,1:(end-1));
fracTimeBoundsU = fracTimeBounds(:,2:end);

% initialize time spent in current state
timeInState = 1;

% initialize timeFrac vector
timeFrac_sample = zeros(size(FS_sample));

% loop through entire trace
for i = 1:numel(FS_sample)
    
    % determine current and next states
    currState = FS_sample(i);
    if i == numel(FS_sample)
        nextState = currState;
    else
        nextState = FS_sample(i+1);
    end
    
    % first determine if state is IRR
    if currState == 4
        % state is IRR
        
        % IRR states get a timeFrac of 1
        timeFrac_sample(i) = 1;
        
    else
        % state is not IRR, proceed as normal
        
        % determine which time frac this is
        timeFracInd = find(fracTimeBoundsL(currState,:) <= timeInState & timeInState < fracTimeBoundsU(currState,:));
        
        % put this timeFracInd in the trace
        timeFrac_sample(i) = timeFracInd;
        
        % determine if two states are equal
        if currState == nextState
            % they are equal
            
            % increment currTime by 1
            timeInState = timeInState+1;
        else
            % they are not equal
            
            % reset currTime to 1
            timeInState = 1;
        end
        
    end
    
end

% modify FS_sample using the fractional times
FS_sample = (FS_sample-1).*options.nTimeFracs+timeFrac_sample;


end

