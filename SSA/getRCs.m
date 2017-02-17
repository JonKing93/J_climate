function[RCs, signals] = getRCs(eigvecs, traj, algorithm)
%% Calculates the RCs for a SSA output
% 
% [RCs, signals] = getRCs(eigvecs, traj, algorithm)
%
% ----- Inputs -----
%
% eigvecs: The eigenvectors for a set of trajectory matrices
%
% traj: A set of trajectory matrices. The RCs of the dim1 x dim2 matrices
% will be calculated.
%
% algorithm: The algorithm used to calculate the trajectory matrices
%   'BK': Broomhead - King
%   'VG': Vautard - Ghil
%
% ----- Outputs -----
% 
% RCs: The reconstructed components for a set of trajectories and
% eigenvectors. Each column of the output is a particular RC
%
% signals: The signals calculated and used to construct the RCs. Each
% column corresponds to the signal for a particular mode of the trajectory

% Error check, get some initial sizes
[M, N, ncols, nseries] = setup(eigvecs, traj, algorithm);

% Preallocate
signals = NaN(M, ncols, nseries);
RCs = NaN(M,N, nseries);

% For each series, calculate signals and RCs
for j = 1:nseries
    % Get signals
    signals(:,:,j) = eigvecs(:,:,j)' * traj(:,:,j);
    
    
    % Build the the RC matrix
    for k=1:M
        notNAN = ~isnan(signals(k,:,j));
        RCs(k,:,j) = conv( eigvecs(:,k,j), signals(k, notNAN, j) );
    end

    % Apply convolution weights (i.e. normalize the RCs)
    for k=1:N
        % Case 1: incomplete overlap of functions because we are at the front
        % of the vector
        if k <= M-1
            RCs(:,k,j) = RCs(:,k,j) ./ k;
        
        % Case 2: Complete overlap of the two functions
        elseif k >= M && k <= N-M+1
            RCs(:,k,j) = RCs(:,k,j) ./ M;
        
        % Case 3: Incomplete overlap because we are at the end of the vector
        else
            RCs(:,k,j) = RCs(:,k,j) ./ (N-k+1);
        end
    end
end

% Return with individual RCs in columns
RCs = permute(RCs, [2,1,3]);
signals = permute(signals, [2,1,3]);

end

% ----- Helper functions -----
function[M, N, ncols, nseries] = setup(eigvecs, traj, algorithm)

% Ensure input matrices are 3D
if ndims(eigvecs) > 3 || ndims(traj) > 3
    error('Too many dimensions in input matrices');
end

% Calculate some sizes
[M, ncols, nseries] = size(traj);

% Get the number of points
switch algorithm
    case 'VG'
        N = ncols - M + 1;
    case 'BK'
        N = ncols + M - 1;
    otherwise
        error('Unrecognized algorithm');
end

end






