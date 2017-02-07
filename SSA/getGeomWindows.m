function[W] = getGeomWindows(minVal, maxVal, power, npoints, a, keepModes)
%% Calculates window sizes using a power series for MS_SSA given a fixed value "a" and desired number of modes.
%
% ----- INPUTS -----
%
% minVal: The minimum desired window size. The power series will be
%   calculated beginning with this value
%
% maxVal: An integer, the maximum desired window size
%   If the string 'max' is entered, maxVal will be set to npoints
%
% power: The geometric factor for the power series
%
% npoints: The number of points in the data time series
%
% a: The fixed value of "a" in MS_SSA (See MS_SSA documentation)
%
% keepModes: The desired modes to be saved for MS_SSA.
%   (e.g. keepModes = [1 2 55] will cause modes 1, 2, and 55 to be saved.)
%
% ----- OUTPUTS -----
%
% W: A vector of window sizes that will contain the desired number of 

% Run a quick error check
[maxVal] = setup(minVal, maxVal, power, npoints, a, keepModes);

% Calculate the power series
W = minVal;
k = W(end);

while k*power < maxVal
    W = [W; W(end)*power]; %#ok<AGROW>
    k = W(end);
end

% Calculate the possible sizes of the SSA partition "M" given the fixed "a"
% and calculated power series.
M = W ./ a;

% Get the indices of windows in W for which the desired number of saved
% modes exceeds the size of the partition window "M"
removeDex = M < max(keepModes);

% Remove these values from the power series
W(removeDex) = [];

% If no windows meet the criteria, print a warning
if isempty(W)
    warning('No window sizes meet the criteria');
end


end

% ----- Helper functions -----
function[maxVal] = setup(minVal, maxVal, power, npoints, a, keepModes)

% Ensure minVal is a positive integer
if minVal <=0 || mod(minVal,1)~=0
    error('minVal must be a positive integer');
end

% Check if maxVal contains a flagged string
if ischar(maxVal)
    switch maxVal
        case 'max'
            maxVal = npoints;
        otherwise
            error('Unrecognized value of maxVal');
    end
end

% Ensure maxVal does not exceed the number of points
if maxVal > npoints
    error('maxVal must be smaller than npoints');
end

% Ensure power is positive
if power <= 0
    error('power must be positive');
end

% Ensure a is positive
if a <=0
    error('''a'' must be positive');
end

% Ensure keepModes is a vector
[x,y] = size(keepModes);
if ~ismatrix(keepModes) || (x>1 && y>1)
    error('keepModes must be a vector');
end
nmodes = max([x,y]);

% Ensure keepModes only contains positive integers
for j = 1:nmodes
    if keepModes(j)<=0 || mod(keepModes(j),1)~=0
        error('keepModes must only contain positive integers');
    end
end

end