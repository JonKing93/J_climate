function[sigScaled] = scaleSignals(signals, eigVals)
%% Scales signals to the standard deviation of the standardized data.
%
% [sigScaled] = scaleSignals(signals, eigVals)
%
% ----- Inputs -----
%
% signals: a set of signals from a PCA analysis. If signals is a 3D stack
%   of matrices each dim1 x dim2 matrix correpsonds to one set of signals.
%
% eigVals: A matrix containing the associated eigenvalues for each set of 
%   Each column contains the eigenvalues for one set of signals. The
%   eigenvalues SHOULD NOT be normalized.
%
%
% ----- Outputs -----
%
% sigScaled: The set of scaled signals. These signals may be plotted
%   directly against the standardized dataset.

[nsignals , nstack] = setup(signals, eigVals);

% Preallocate signals
sigScaled = NaN(size(signals));

% Take square root of eigVals
eigVals = sqrt(eigVals);

% for each matrix of signals...
for j = 1:nstack
    for k = 1:nsignals % For each signal
        % Scale the signal by the corresponding eigenvalue
        sigScaled(:,k,j) = signals(:,k,j) ./ sqrt(eigVals(k,j));
    end
end

end

%%%%% Helper Functions %%%%%
function[nsignals, nstack] = setup(signals, eigVals)

% Ensure signals is 3D
if ndims(signals) > 3
    error('signals must be 3D');
end

% Ensure eigVals is a matrix
if ~ismatrix(eigVals)
    error('eigVals must be a matrix');
end

% Ensure there are no NaNs
if NaNcheck(signals)
    error('signals cannot contain NaN');
end
if NaNcheck(eigVals)
    error('eigVals cannot contain NaN');
end

% Ensure eigVals are positive
if any(any( eigVals < 0))
    error('The eigenvalues should all be positive');
end

% Ensure the signals and eigenvalues align properly
[~, nsignals, nstack] = size(signals);
[neig, neigSets] = size(eigVals);
if nsignals ~= neig || nstack ~= neigSets
    error('signals and eigVals do not have matching dimensions');
end
end
