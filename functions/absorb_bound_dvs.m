function [P_left, P_right] = absorb_bound_dvs(P_left, P_right, bound, isRight, crossStep)
%% =======================================================================
% Absorbing bound on simulated accumulators
%
% Inputs:
%   P_left, P_right - time steps x trials, left and right accumulator
%   bound           - bound at each time step
%   isRight         - trials x 1, choice (true = right)
%   crossStep       - trials x 1, time step of the first bound crossing
%                     (NaN = no crossing)
%
% Winning accumulator (chosen side): stays at the bound from its crossing
% on. Losing accumulator: keeps accumulating until it reaches the bound
% itself, then stays at the bound (otherwise until stimulus offset).
% Trials without a crossing are not changed.
%% =======================================================================

[nSteps, nTrials] = size(P_left);

bound     = bound(:);
timeStep  = (1:nSteps)';
isRight   = logical(isRight(:))';
crossed   = isfinite(crossStep(:))';

crossStep = crossStep(:)';
crossStep(~crossed) = nSteps + 1;

% Winning and losing accumulator
winning = P_left;   winning(:,isRight) = P_right(:,isRight);
losing  = P_right;  losing(:,isRight)  = P_left(:,isRight);

% Winning accumulator: bound value at the crossing, from the crossing on
boundAtCrossing = nan(1,nTrials);
boundAtCrossing(crossed) = bound(crossStep(crossed));

afterCrossing = timeStep >= crossStep;
boundMatrix   = repmat(boundAtCrossing, nSteps, 1);
winning(afterCrossing) = boundMatrix(afterCrossing);

% Losing accumulator: first time step after the winner's crossing at
% which it reaches the bound
reachesBound = (timeStep > crossStep) & (losing >= bound);
[losingReaches, losingStep] = max(reachesBound, [], 1);
losingReaches = logical(losingReaches) & crossed;

boundAtLosing = nan(1,nTrials);
boundAtLosing(losingReaches) = bound(losingStep(losingReaches));

afterLosing = losingReaches & (timeStep >= losingStep);
boundMatrix = repmat(boundAtLosing, nSteps, 1);
losing(afterLosing) = boundMatrix(afterLosing);

% Back to left / right accumulator
P_left  = winning;  P_left(:,isRight)  = losing(:,isRight);
P_right = losing;   P_right(:,isRight) = winning(:,isRight);

end
