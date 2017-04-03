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
[svdFunc, svdsArg] = parseInputs( varargin{:});

% Error check
errCheck(C);

% Run the SVD on each matrix
switch svdFunc
    case 'svd'
        [~,S,V] = svd(C);

    case 'svds'
        [~,S,V] = svds(C, svdsArg); 
end

% Pull the eigenvalues off of the diagonal
eigvals = diag(S);

% Set the eigenvectors so that the majority of elements are positive
eigvecs = posColSign(V);
        
end
    
    
%%%%% Helper Functions %%%%%
function[svdFunc, svdsArg] = parseInputs( varargin)
inArgs = varargin;

% Set the defaults
svdFunc = 'svd';
svdsArg = NaN;

if ~isempty(inArgs)
    isSvdsArg = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
       
        if isSvdsArg
            if isscalar(arg) || strcmpi(arg,'econ')
                svdFunc = 'svds';
                svdsArg = arg;
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end
        elseif strcmpi(arg,'svd')
            % Do nothing
        elseif strcmpi(arg, 'svds') 
            if length(inArgs)>=k+1
                isSvdsArg = true;
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end            
        else
            error('Unrecognized Input');
       end
   end
end
end
          
function[] = errCheck(C)
% C is a matrix
if ~ismatrix(C)
    error('C must be a matrix');
% No NaN
elseif hasNaN(C)
    error('C may not contain NaN');
end
end

