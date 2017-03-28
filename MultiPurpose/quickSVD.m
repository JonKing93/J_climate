function[eigvals, eigvecs] = quickSVD(C, varargin)
%% Run singular value decomposition on a matrix. Returns eigenvalues, eigenvectors with majority positive loadings.
%
% [eigvals, eigvecs] = quickSVD(C)
%   performs a full SVD on a matrix or set of matrices.
%
% [...] = quickSVD(C, 'svds', 'econ')
% Uses an economy sized svds decomposition rather than the full svd.
%
% [...] = quickSVD(C, 'svds', nEigs)
% Uses the svds decomposition and determines the first nEigs eigenvalues.
%
% [...] = quickSVD(C, 'svd')
% Performs the normal svd.
%
% ----- Inputs -----
%
% C: A matrix. May not contain NaN.
%
%
% ----- Outputs -----
%
% eigvals: A vector containing the eigenvalues of C.
%
% eigvecs: A matrix containing the eigenvectors of C. Each column is one
%       eigenvector.

% Parse Inputs
[svdFunc, svdsArg] = parseInputs( varargin(:));


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
function[svdFunc, svdsArg] = parseInputs( varargin)
inArgs = varargin;

% Set the defaults
svdFunc = 'svd';
svdsArg = NaN;

if ~isempty(inArgs)
   for k = 1:length(inArgs)
       arg = inArgs{k};
       
       if strcmpi(arg,'svd')
           % Do nothing
       elseif strcmpi(arg, 'svds')
            if length(inArgs) >= k+1 && ( isscalar(inArgs{k+1}) || strcmpi(inArgs{k+1},'econ') )
                svdFunc = 'svds';
                svdsArg = inArgs{k+1};
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end
        else
            error('Unrecognized Input');
       end
   end
end
end
          


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
