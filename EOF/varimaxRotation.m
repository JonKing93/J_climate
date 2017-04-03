function[rotModes, rotEigvals, rotExpVar, rotMatrix] = varimaxRotation(scaModes, eigVals)
%% Performs a VARIMAX rotation on a set of scaled eigenvectors and eigenvalues.
%
% [rotModes, rotMatrix] = varimaxRotation(scaledEigvecs) returns a set fo
%   eigenvectors rotated using the varimax criterion and the corresponding
%   rotation matrix
%
% [rotModes, rotMatrix, rotEigvals] = varimaxRotation(scaledEigvecs, eigVals)
%   also returns the rotated eigenvalues when the initial values are given
%   as input.
%
% ----- Inputs -----
%
% scaModes: A set of modes scaled for varimax rotation by the square root
% of associated eigenvalues. Each column contains 1 mode.
%   matrix contains a set of eigenvectors
%
% *** Optiontal Inputs ***
%
% eigVals: The eigenvalues / loadings associated with a set of eigenvector
%   modes. Each column corresponds to one set of eigenvectors. If the
%   eigenvalues are given, varimaxRotation 
%
%
% ----- Outputs -----
%
% rotEigvecs: The rotated eigenvectors.
%
% rotEigvals: The rotated eigenvalues.
%
% rotExpVar: The explained variance of the rotated eigenvalues.
%
% rotMatrix: The matrix used to perform the eigenvector rotation.

% Error checking, determine if eigenvalues were given
if nargin == 1
    [nSets, haveEigvals] = setup(scaModes);
else
    [nSets, haveEigvals] = setup(scaModes, eigVals);
end

% Perform the rotation
[rotModes, rotMatrix] = rotatefactors(scaModes);

% Calculate the rotated eigenvalues
if haveEigvals
    rotEigvals = NaN(size(eigVals));
    rotExpVar = NaN(size(eigVals));
    for k = 1:nSets
        rotEigvals(:,k) = diag(  diag(eigVals(:,k)) * rotMatrix  );
        rotExpVar(:,k) = 100 * ( rotEigvals(:,k) ./ sum(eigVals(:,k)) );
    end
end

end

%%%%% Helper Functions %%%%%
function[nSets, haveEigvals] = setup(scaledEigvecs, eigVals)
% Check if eigenvalue checking is needed
if nargin == 1
    haveEigvals = false;
else
    haveEigvals = true;
end

% Check for proper number of dimensions
if ndims(scaledEigvecs) > 3
    error('scaledEigvecs must be 3D');
end
if haveEigvals && ~ismatrix(eigVals)
    error('eigVals must be a matrix');
end

% Check for NaN
if hasNaN(scaledEigvecs)
    error('scaledEigVecs cannot contain NaN');
end
if haveEigvals && NaNcheck(eigVals)
    error('eigVals cannot contain NaN');
end

% Ensure the eigenvalues are all positive
if haveEigvals && any(any( eigVals < 0))
    error('The eigenvalues must be all positive');
end

% Get the size of eigenvectors
[~, nVecs, nSets] = size(scaledEigvecs);

% Ensure sizes align
if haveEigvals
    [neigs, nValSets] = size(eigVals);
    if nValSets ~= nSets || neigs ~= nVecs
        error('scaledEigvecs and eigVals do not have corresponding sizes');
    end
end

end

