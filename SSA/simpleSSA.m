function[eigvals, eigvecs, Data_m0, traj, C] = simpleSSA(Data, M, varargin)
%% Runs a single-spectrum analysis.
%
% [eigvals, eigvecs, Data_m0, traj, C] = simpleSSA(Data, M)
% performs singular spectrum analysis on each time series in a set of data
% series.
%
% [...] = simpleSSA(..., algorithm)
% performs singular spectrum analysis using a trajectory matrix constructed
% via a particular algorithm.
%
% [...] = simpleSSA(..., nEigs)
% will perform singular spectrum analysis only a specfied subset of
% eigenvectors for the dataset.
%
% ----- Inputs -----
%
% Data: A matrix of time series. Each column of corresponds to a particular
%   time series.
%
% M: The window size used to construct the trajectory matrix
%
% % *** Optional Inputs ***
%
% algorithm: a flag for the method used to construct the trajectory matrix
%     'VG': (Default) Vautard-Ghil
%     'BK': Broomhead-King
%
% nEigs: An integer specifying the number of eigenvalues and eigenvectors
%   to obtain. Will use a svds rather than svd analysis. 
%   (GENERALLY NOT RECOMMENDED).
%
% ----- Outputs -----
%
% eigvals: The singular values of each SSA of each time series
%
% eigvecs: The singular vectors of each SSA of each time series 
%
% Data_m0: a set of time series with means removed
%
% traj: A set of trajectory matrices for the time series. Each dim1 x dim2
%   matrix is the trajectory matrix for one time series.
%
% C: The covariance matrix for each trajectory matrix. Each dim1 x dim2
%   matrix is the covariance matrix for a particular time series.
%
% ----- Additional Reading -----

%% Error check
[algorithm, svdtype] = setup(Data, M, varargin);


%% Remove mean from data
Data_m0 = detrend(Data, 'constant');


%% Get trajectories and covariance
[traj, C] = buildTrajandCov(Data_m0, M, algorithm);


%% Run an svd of each series or svds as appropriate
[eigvals, eigvecs] = quickSVD(C, svdtype{:});

end


% ----- Helper functions -----
function[algorithm, svdtype] = setup(Data, M, inArgs)

% Ensure that Data is not 3D or greater
if ~ismatrix(Data)
    error('Data must be a matrix');
end

% Ensure M is a positive, integer, scalar no longer than the data series
if ~isscalar(M)
    error('M must be a scalar');
elseif M<1  
    error('M must be positive');
elseif mod(M,1) ~= 0
    error('M must be an integer');
elseif M > size(Data,1)
    error('M cannot exceed the length of the data series');
end
if length(M) ~= 1
    error('M must be singular');
end

% Set the defaults
algorithm = 'VG';
svdtype = {'svd'};

% Parse the inputs
if ~isempty(inArgs)
    nIn = length(inArgs);
    if nIn > 2
        error('Too many input arguments')
    else
        algset = false;
        eigset = false;
        for k = 1:nIn
            if ~algset && ischar(inArgs{k})
                algorithm = inArgs{k};
                algset = true;
            elseif isscalar(inArgs{k}) && ~eigset
                svdtype = [{'svds'}, {inArgs{k}}];
                eigset = true;
            else
                error('Unrecognized input');
            end
        end
    end
end
    

end