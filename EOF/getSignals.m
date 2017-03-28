function[signals] = getSignals(Datax0, modes)
%% Gets the signals from a standardized dataset and its eigenvectors.
%
% [signals] = getSignals(Datax0, modes)
%
% ----- Inputs -----
%
% Datax0: A standardized dataset. A 2D matrix with columns that represent
%       individual time series.
%
% modes: A set of the modes of the standardized dataset. Each mode is along one column
%
%
% ----- Outputs -----
%
% signals: The signals that arise from redistributing the dataset over the modes.
%       If all modes are used, no information is lost and the signals may 
%       be used to completely reconstruct the original data.


% Check for errors
errorCheck(Datax0, modes);

% Get the signals
signals = Datax0 * modes;

end

%%%%% Helper functions %%%%%
function[] = errorCheck(Datax0, eigVecs)

% Ensure Datax0 is 2D
if ~ismatrix(Datax0)
    error('Datax0 must be a 2D matrix');
end

% Ensure eigVecs is 2D
if ~ismatrix(eigVecs)
    error('eigVecs must be a 2D matrix');
end

% Check for NaNs
if NaNcheck(Datax0)
    error('Datax0 cannot contain NaNs');
end

if NaNcheck(eigVecs)
    error('eigVecs cannot contain NaNs');
end

% Ensure matrix dimensions align
[~,x] = size(Datax0);
[y,~] = size(eigVecs);
if x ~= y
    error('Data and eigVecs have mismatched sizes');
end

end