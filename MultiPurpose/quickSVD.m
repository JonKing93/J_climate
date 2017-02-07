function[eigvals, eigvecs] = quickSVD(C, varargin)
%% Run singular value decomposition on a matrix. Returns eigenvalues, eigenvectors with majority positive loadings.
%
% [eigvals, eigvecs] = quickSVD(C, decompSpecs)
%
% ----- Inputs -----
%
% C: A matrix or 3D stack of matrices. If C is a 3D stack of matrices, an
%       svd will be performed on each dim1 x dim2 matrix.
%
% *** Optional Inputs ***
%
% decompSpecs: Additional specifications for the svd
%
%       Decomposition type: svd or svds (these may have differing runtimes)
%           'svd': (default) Uses the svd function for decomposition
%           'svds': Uses the svds function for decomposition
%             
%       Decomposition size: The number of eigenvectors found (may affect runtime)
%           'all': (default) Performs the full decomposition
%           'econ': Performs the economy size decomposition
%           nEigs: An integer specifying the number of leading eigenvectors to find.
%                     (neigs is only recommended for use with the 'svds' switch)
%                     This may save memory for large analyses.
%
%
% ----- Outputs -----
%
% eigvals: A 2D matrix of eigenvalues. Each column corresponds to the
%       values for a particular dim1 x dim2 matrix in C.
%
% eigvecs: A 3D matrix of eigenvectors. Each dim1 x dim2 matrix corresponds
%       to a particular matrix in C.
%
%
% ----- Background Readings -----
%
% IN PROGRESS...
%
%


[nCols, numMat, decompType, decompArg] = setup(C, varargin{:});

% Preallocate
eigvals = NaN(nCols, numMat);
eigvecs = NaN(nCols, nCols, numMat);

% Run the SVD on each matrix
for k = 1:numMat
    
    switch decompType
        case 'svd'
            if strcmp(decompArg,'all')
                [~,S,V] = svd(C(:,:,k));           
            else
                [~,S,V] = svd(C(:,:,k), 0);
            end
            % Reduce dataset if desired
            if isscalar(decompArg)
                S = S(1:nCols,:);
                V = V(1:nCols,:);
            end
            
        case 'svds'
            if isscalar(decompArg)
                [~,S,V] = svds(C(:,:,k), decompArg);
            else
                [~,S,V] = svds(C(:,:,k));
            end
    end
                    
    % Pull the eigenvalues off of the diagonal
    eigvals(:,k) = diag(S);
    
    % Set the eigenvectors so that the majority of loadings are positive
    eigvecs(:,:,k) = posColSign(V);
        
end
end
    
    
%%%%% Helper Functions %%%%%

function[nCols, numMat, decompType, decompArg] = setup(C,varargin)
%% Error checking and determination of input specifications

% Ensure C is 3D 
if ndims(C) > 3
    error('quickSVD is for 3D matrices only');
end

% Ensure no NaNs are in C
if NaNcheck(C)
    error('C may not contain NaN entries');
end

% Get the size of C
[~, nCols, numMat] = size(C);

% Set defaults
decompType = 'svd';
decompArg = 'all';

% Done if there were no additional specifications
if isempty(varargin)
    return;
end

% Set additional specifications
for j = 1:length(varargin)
    spec = varargin{j};
    
    if strcmp(spec, 'svd') % Do nothing (default)
        
    elseif strcmp(spec, 'svds')
        decompType = 'svds';
        
    elseif strcmp(spec, 'all') % Do nothing (default)
        
    elseif strcmp(spec, 'econ')
        decompArg = 'econ';
        
    elseif isscalar(spec) % for nEigs, will not need as large of a matrix size
        decompArg = spec;
        if spec > nCols   % Throw error if spec is too high
            error('The desired number of eigs is greater than the possible number of eigs');
        end
        nCols = spec;
        
    else
        error('Unrecognized input argument');
    end
end

% Check that svds is not set to the full decomposition
if strcmp(decompType, 'svds') && strcmp(decompArg, 'all')
    error('svds can only perform the economy decomposition');
end

% Send warning if using svd with nEigs
if strcmp(decompType, 'svd') && isscalar(decompArg)
    warning('Calculating a limited number of eigs works better with svds rather than svd');
end


end
