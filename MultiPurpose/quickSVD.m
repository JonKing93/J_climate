function[eigvals, eigvecs] = quickSVD(C, varargin)
%% Run singular value decomposition on a matrix. Returns eigenvalues, eigenvectors with majority positive loadings.
%
% [eigvals, eigvecs] = quickSVD(C)
%   performs a full SVD on a matrix or set of matrices.
%
% [eigvals, eigvecs] = quickSVD(C, svdType, svdSize)
% allows the user to choose between svd and svds decompositions, and select
% the size of each type of decomposition.
%
% ----- Inputs -----
%
% C: A matrix or 3D stack of matrices. If C is a 3D stack of matrices, an
%       svd will be performed on each dim1 x dim2 matrix.
%
% *** Optional Inputs ***
%
% svdType: A flag for the type of decomposition
%       'svd': (Default) Uses the svd decomposition
%       'svds': Uses the function for the svds decomposition.
%
% svdSize: Specify the number of singular values found.
%       For svd: 'all': (Default) Performs the full decomposition
%                'econ': Performs the economy sized decomposition
%                nEigs: an integer specifying the number of singular values
%                       to find (not recommended for svd)
%
%       For svds: 'econ': (Default) Performs the economy sized decomposition
%                 nEigs: An integer specifying the number of singular values to find
%
% ----- Outputs -----
%
% eigvals: A 2D matrix of eigenvalues. Each column corresponds to the
%       values for a particular dim1 x dim2 matrix in C.
%
% eigvecs: A 3D matrix of eigenvectors. Each dim1 x dim2 matrix corresponds
%       to a particular matrix in C.
%
% ---
% Jonathan King, 2017

% Error check, parse inputs, get number of matrices
[numMat, decompType, decompArg] = setup(C, varargin{:});

% Run the SVD on each matrix
for k = 1:numMat
    
    switch decompType
        case 'svd'
            if strcmp(decompArg,'all')
                [~,S,V] = svd(C(:,:,k));           
            else
                [~,S,V] = svd(C(:,:,k), 0);
            end
            
        case 'svds'
            if isscalar(decompArg)
                [~,S,V] = svds(C(:,:,k), decompArg);
            else
                [~,S,V] = svds(C(:,:,k));
            end 
    end
   
    % Preallocate if this is the first run
    if k==1
        eigvals = NaN( length(diag(S)), numMat);
        eigvecs = NaN( [size(V), numMat]);
    end
                  
    % Pull the eigenvalues off of the diagonal
    eigvals(:,k) = diag(S);
    
    % Set the eigenvectors so that the majority of loadings are positive
    eigvecs(:,:,k) = posColSign(V);
        
end
end
    
    
%%%%% Helper Functions %%%%%

function[nMat, svdType, nEigs] = setup(C,varargin)
%% Error checking and determination of input specifications

% Ensure C is 3D 
if ndims(C) > 3
    error('quickSVD is for 3D matrices only');
end

% Ensure no NaNs are in C
if NaNcheck(C)
    error('C may not contain NaN entries');
end

% Get the number of matrix stacks
nMat = size(C,3);

% Set defaults
svdType = 'svd';
nEigs = 'all';

% Check the additional specifications
if ~isempty(varargin)
    if length(varargin) > 2
        error('Too many inputs');
    else
        % SVD type
        if strcmp( varargin{1}, 'svd')
            % Do nothing, this is the default
        elseif strcmp( varargin{1}, 'svds')
            svdType = 'svds';
            nEigs = 'econ';
        else
            error('Unrecognized svd type');
        end
        % nEigs
        if length(varargin) == 1
            % Do nothing
        elseif strcmpi( varargin{2}, 'all') && strcmp(svdType, 'svd')
            % Do nothing
        elseif strcmpi( varargin{2}, 'econ')
            nEigs = 'econ';
        elseif isscalar(varargin{2}) && strcmp(svdType, 'svd')
            error('Use the svds switch to find a specific number of singular values');
        elseif strcmpi(varargin{2}, 'all')
            error('Use the svd switch to perform a full decomposition');
        elseif isscalar(varargin{2})
            nEigs = varargin{2};
        else
            error('Unrecognized input');
        end
    end
end

end
