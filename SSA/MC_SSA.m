function[surrVals, iterSigVals, iterTrueConf] = MC_SSA(ts_m0, singVecs, MC, noise, window, algorithm, pval, varargin)
%% Runs a Monte Carlo singular spectrum analysis. Generates surrogate eigenvalues.
%
% [surrVals, iterSigVals, iterTrueConf] = MC_SSA(ts_m0, singVecs, MC, noise, window, algorithm, pval)
% Conducts an MC-SSA significance test on a set of singular vectors.
% Returns the set of surrogate singular values as well as surrogate
% eigenvalues at the significance level for each iteration.
%
% [...] = MC_SSA(..., showProgress)
% A flag to show the current Monte Carlo iteration number.
%
% [surrVals] = MC_SSA(..., convergeTest)
% A flag to block the saving of convergence information. May improve
% runtime for large analyses, but causes a loss of information.
%
%
% ----- Inputs -----
% 
% ts_m0: A time series with mean of 0.
%
% singVecs: The singular vectors for the time series.
%
% MC: The number of Monte Carlo iterations in the MC-SSA
% 
% noise: The type of noise with which to construct surrogate trajectories
%       'white': White Gaussian noise
%       'red': Red lag-1 autocorrelated noise with added white noise
%
% window: The embedding dimension used to construct trajectories
%
% algorithm: The algorithm used to build SSA trajectories
%       'BK': Broomhead-King
%       'VG': Vautard-Ghil
%
% pval: The significance level to record for convergence testing.
%
% showProgress: A flag to display the current Monte Carlo iteration number
%       'noProgress' (Default): Do not display number
%       'showProgress': Display number
%
% convergeTest: A flag to toggle significance testing.
%       'convergeTest' (Default): Record convergence information
%       'noConvergeTest': Block the storage of convergence information
%
%
% ----- Outputs -----
%
% surrVals: The matrix of surrogate singular values

% Error check
[showProgress, convergeTest] = parseInputs(varargin{:});
errCheck(ts_m0, MC, pval);

% Preallocate surrogate matrices
surrC = NaN(window, window, MC);
surrProj = NaN(window, window, MC);
surrVals = NaN(MC, window);
if convergeTest
    iterSigVals = NaN(MC, window);
    iterTrueConf = NaN(MC,1);
else
    iterSigVals = [];
    iterTrueConf = [];
end

% Build the surrogate series
surr = randNoiseSeries(noise, ts_m0, MC);

% For each surrogate series...
for k = 1:MC
    
    % Show progress if desired
    if showProgress
        fprintf('Monte Carlo Iteration: %i / %i\r\n', k, MC);
    end

    % Get the trajectory and associated covariance matrix
    [~, surrC(:,:,k)] = getTandC(surr(:,k), window, algorithm);

    % Project onto the data eigenvector basis
    surrProj(:,:,k) = singVecs' * surrC(:,:,k) * singVecs;

    % Extract the diagonal elements
    surrVals(k,:) = diag( surrProj(:,:,k) );
    
    % Sort surrVals every iteration if testing convergence
    if convergeTest
        % Sort the new set of surrogate values
        surrVals = sort(surrVals);
        
        % Calculate the confidence level threshold
        thresh = ceil(k* (1-pval));
        iterTrueConf(k) = thresh/k;
        
        % Get the set of values on the confidence interval
        iterSigVals(k,:) = iterSigVals(thresh,:);
    end        
end

% Sort the projection diagonals if no convergence tests
if ~convergeTest
    surrVals = sort(surrVals);
end

end


% ----- Helper functions -----
function[showProgress, convergeTest] = parseInputs(varargin)
inArgs = varargin;
showProgress = false;
convergeTest = true;

if ~isempty(inArgs)
    for k = 1:length(inArgs)
        arg = inArgs{k};
        if strcmpi(arg, 'showProgress')
            showProgress = true;
        elseif strcmpi(arg, 'noProgress')
            showProgress = false;
        elseif strcmpi(arg, 'noConvergeTest')
            convergeTest = false;
        elseif strcmpi(arg, 'convergeTest')
            % Do nothing
        else
            error('Unrecognized Input');
        end
    end
end
end
            

function[] = errCheck(ts_m0, MC, pval)

if pval<=0 || pval>=1
    error('pval must be on the interval (0,1)');
end

% Ensure Data_m0 is 2D
if ~isvector(ts_m0)
    error('ts_m0 must be a vector');
end

% Ensure MC is positive
if MC < 1
    error('Monte Carlo number must be positive');
end
end