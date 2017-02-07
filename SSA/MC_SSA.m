function[surrEig] = MC_SSA(Data_m0, eigvecs, MC, noise, window, algorithm)
%% Runs a Monte Carlo singular spectrum analysis. Generates surrogate eigenvalues.
%
% [surrEig] = MC_SSA(Data_m0, eigvecs, MC, noise, window, algorithm)
%
% ----- Inputs -----
% 
% Data_m0: A set of time series with the means removed. Each column is a
%   time series.
%
% eigvecs: The eigenvectors for each time series. Each dim1 x dim2 matrix
% corresponds to the eigenvectors of the corresponding time series. 
%
% MC: The Monte Carlo number for the MC_SSA
% 
% noise: The type of noise with which to construct surrogate trajectories
%   'white': White Gaussian noise
%   'red': Red lag-1 autocorrelated noise
%
% window: The size of window used to build SSA trajectories
%
% algorithm: The algorithm used to build SSA trajectories
%   'BK': Broomhead-King
%   'VG': Vautard-Ghil
%
% ----- Outputs -----
%
% surrEig: The matrix of surrogate eigenvectors
% 

% Error check and get some initial sizes
[npoints, nseries] = errCheck(Data_m0, MC);

% Preallocate surrogate matrices
surr = NaN(npoints, MC, nseries);
surrC = NaN(window, window, MC, nseries);
surrProj = NaN(window, window, MC, nseries);
surrEig = NaN(MC, window, nseries);

% Calculate the appropriate autocorrelation coefficient for the noise type
switch noise
    % White noise
    case 'white'
        ar1 = zeros(nseries, 1);
        
    % Red autocorrelated lag-1 noise    
    case 'red'
        ar1 = corr( Data_m0(1:end-1,:), Data_m0(2:end,:) );
        ar1 = diag(ar1);
        
    otherwise
        error('Unrecognized noise type');
end

% Get the variances of the time series
standev = std(Data_m0);

% Initialize the series
surr(1,:,:) = 0;

% For each time series...
for j = 1:nseries
    
    % Build the surrogate series
    for k = 2:npoints
        surr(k,:,j) = ( ar1(j) .* surr(k-1,:,j) ) + randn(1,MC);
    end
    
    % Remove mean and variance
    surr(:,:,j) = zscore(surr(:,:,j));
    
    % Scale the surrogate series to the standard deviation of the data
    % series
    surr(:,:,j) = surr(:,:,j) * standev(j);
    
    
    % For each surrogate series...
    for k = 1:MC
        
        % Get the trajectory and associated covariance matrix
        [~, surrC(:,:,k,j)] = buildTrajandCov(surr(:,k,j), window, algorithm);
        
        % Project onto the data eigenvector basis
        surrProj(:,:,k,j) = eigvecs(:,:,j)' * surrC(:,:,k,j) * eigvecs(:,:,j);
        
        % Extract the diagonal elements
        surrEig(k,:,j) = diag( surrProj(:,:,k,j) );
    end
    
    % Sort the projection diagonals for each time series
    surrEig(:,:,j) = sort( surrEig(:,:,j) );
end


end


% ----- Helper functions -----
function[npoints, nseries] = errCheck(Data_m0, MC)

% Ensure Data_m0 is 2D
if ~ismatrix(Data_m0)
    error('MC_SSA is for 2D sets of time series');
end

% Ensure MC is positive
if MC < 1
    error('Monte Carlo number must be positive');
end

% Get some sizes
[npoints, nseries] = size(Data_m0);
end