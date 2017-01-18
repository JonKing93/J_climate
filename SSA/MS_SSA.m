function[winEigvals, winEigvecs, winDatam0] = MS_SSA(Data, a, W, algorithm, keepModes, varargin)
%% Runs a multiscale SSA Analysis on a set of time series.
%
% ----- Inputs -----
%
% Data: A 2D matrix, each column is a time series
% 
% a: The integer to which the ratio W / M is fixed. (suggested, a=3)
%
% W: The desired window sizes (often a geometric series)
%
% algorithm: The desired algorithm to use for the SSA trajectories
%   'VG': Vautard-Ghil
%   'BK': Broomhead-King
%
% keepModes: A vector of integers indexing the SSA modes to keep
%   (e.g. keepVecs = [1 2 55] will save SSA EOF modes 1, 2 and 55
%   
%   'all': will store every eigenmode of every analysis
%       (For large datasets, this may exceed system memory)
%
% Beta: (Optional Argument) A vector containing scaling factors to decrease
%   computation time. For a given window size W, windowed samples of the
%   time series X(t) will begin at each point Beta*t
%
%   Example 1: Beta = 1, W = 60, npoints = 100
%       The windows will be selected at each time point giving the vectors
%       1:60, 2:61, 3:62, ... 40:99, and 41:100.
%
%   Example 2: Beta = 5, W = 60, npoints = 100
%       The window will be selected every 5 time points giving the vectors
%       1:60, 6:65, 11:70, ... 36:95, and 41:100
%
%       Total Computation time decreases by a factor of 5 from Example 1
%
%   IMPORTANT: For best results, Beta must be much smaller than the window
%       size (Beta <<< W).
%
% ----- Outputs -----
%
% winEigvals: A cell vector, each cell contains the eigenvalues of each
%   time series for the SSA EOFs computed for each window size.
%
%   If only 1 window size is input, returns a matrix.
%
% winEigvecs: A cell vector, each cell contains the eigenvectors of each
%   time series for the SSA EOFs computed for each window size.
%
%   If only 1 window size is input, returns a matrix.
%

% Error check and get sizes
[M, nSSA, nwins, nseries, nvecs, nEigs, keepModes] = setup(Data, a, W, keepModes, varargin{:});

winEigvals = cell(nwins,1); % The number of eigs will change as partition
winEigvecs = cell(nwins,1); % size changes, thus store them in cells.

% For each window size in the analysis...
for k = 1:nwins

    % Preallocate
    winData = NaN( W(k), nSSA(k), nseries);
    eigvals = NaN( nvecs, nSSA(k), nseries);
    eigvecs = NaN( M(k), nvecs, nSSA(k), nseries);
    
    % Sample each time series for each window of length W
    for i = 1:nSSA(k)
        winData(:, i, :) = Data(i:W(k)+i-1,:);
    end
    
    % For each time series...
    for j = 1:nseries
        
        % Run an SSA on the windowed samples, obtain eigs
        [S, V, winDatam0] = simpleSSA(winData(:,:,j), M(k), algorithm, nEigs);
        
        % Save the desired number of eigs
        eigvals(:,:,j) = S(keepModes, :);
        eigvecs(:,:,:,j) = V(:,keepModes,:);
    end
    
    % Store the differently sized sets of eigs in a cell
    winEigvals{k} = eigvals;
    winEigvecs{k} = eigvecs;
end

% Return a matrix if only 1 window is given. (Friendly for users
% unaccustomed to MATLAB cells.)
if nwins == 1
    winEigvals = winEigvals{:};
    winEigvecs = winEigvecs{:};
end
    
end

% ----- Helper functions -----
function[M, nSSA, nwins, nseries, nvecs, nEigs, keepModes] = ...
    setup(Data, a, W, keepModes, varargin)

% Check that Data is 2D
if ~ismatrix(Data)
    error('Data matrix must be 2D');
end

% Check that W is a simple vector
[x,y] = size(W);
if ~ismatrix(W) || (x>1 && y>1) 
    error('W must be a single vector of window sizes');
end
nwins = max([x,y]);

% Check that the window size does not exceed data length
[npoints, nseries] = size(Data);
if max(W) > npoints
    error('Maximum window size exceeds data length');
end

% Calculate M for each window
M = W ./ a;

% Round M if not an integer, notify user
if sum( mod(M,1)~=0) ~=0
    M = round(M);
    fprintf('Rounding SSA partition (M) values to nearest integer');
end

% Calculate the numer of SSA analyses for each window
nSSA = npoints-W+1;

% Get the maximum index of an eig for the svds
if strcmp(keepModes, 'all')
    nEigs = M;
    keepModes = 1:M;
else
    nEigs = max(keepModes);
end

% Check that keepModes is a vector
[x,y] = size(keepModes);
if ~ismatrix(keepModes) || (x>1 && y>1)
    error('keepModes must be a vector');
end

% Get the number of saved vectors
nvecs = max([x,y]);

% Check that the desired saved modes do not exceed the maximum mode
if nEigs > min(M)
    error('Desired saved mode %.0f exceeds number of modes for window size %.0f',...
        max(keepModes), W( M == min(M)));
end

end

