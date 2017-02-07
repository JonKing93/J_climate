function[traj, C] = buildTrajandCov(Data_m0, M, algorithm)
%% Builds a trajectory matrix using either the BK (Broomhead King) or VG (Vautard Ghil) algorithm,
% and constructs the associated covariance matrix quickly by using the
% inherent symmetry of each method.
%
% [traj, C] = buildTrajandCov(Data_m0, M, algorithm)
%
% ----- Inputs -----
%
% Data_m0: A set of time series with the means removed. Each column is a
%   time series.
%
% M: The size of window used to build SSA trajectories
%
% algorithm:
%   'BK': Broomhead King
%   'VG': Vautard Ghil 
%
% ----- Outputs -----
% 
% traj: The trajectory matrix.
%
% C: The covariance matrix associated with each method

% Run an error check and get pre-allocation sizes
[N, nseries] = errCheck(Data_m0, M);

switch algorithm
    % In the VG algorithm, the window slides off the ends of the time
    % series.
    case 'VG'

        % Set the matrix size
        ncols = N+M-1;
        
        % Preallocate
        traj = NaN(M, ncols, nseries);
        C = NaN(M,M,nseries);
        
        for j = 1:nseries  % Get trajectory and covariance of each series
                 
            % Fill matrix
            for k = 1:M-1    % Before the window completely covers the series
                traj(M-k+1:M, k, j) = Data_m0(1:k, j);
            end
        
            for k = M: ncols-M   % When there is complete overlap
                traj(:,k,j) = Data_m0(k-M+1:k, j);
            end
        
            for k = ncols-M+1:ncols    % When the window has fallen off
                traj(1:ncols-k+1,k,j) = Data_m0(k-M+1:N,j);
            end
        
            % Calculate covariance matrix. This will be a symmetric
            % Toeplitz matrix, which is used to speed computation.
            for k = 1:M
                % Get the first row
                notNAN = ( ~isnan(traj(1,:,j)) & ~isnan(traj(k,:,j)) );
                C(1,k,j) = sum( traj(1,notNAN,j) .* traj(k,notNAN,j)) ./ sum(notNAN);
            end
            
            % Get the final covariance matrix using Toeplitz symmetry
            C(:,:,j) = toeplitz( C(1,:,j));
        
        end
        
        
 
    % In the BK algorithm, the window stops when it reaches either end of
    % the time series.   
    case 'BK'
        
        ncols = N-M+1;
        
        % Preallocate
        traj = NaN(M, ncols, nseries);
        C = NaN(M,M,nseries);
        
        for j = 1:nseries
        
            % Fill matrix
            for k = 1:ncols  % For each value in the window M
                traj(:,k,j) = Data_m0(k:M+k-1,j);
            end
        
            % Calculate covariance matrix
            C(:,:,j) = cov(traj(:,:,j)');
        end
        
        
        
    otherwise
        error('Unrecognized algorithm');
end
end

% ----- Helper Functions -----
function[N, nseries] = errCheck(Data_m0, M)

% Check that Data_m0 is a single time series 
if ~ismatrix(Data_m0)
    error('Data_m0 must be a 2D collection of time series');
end

% Get size of Data
[N, nseries] = size(Data_m0);

% Check that the window size does not exceed the time series size
if M > N
    error('Window size exceeds time series length');
end

end




 