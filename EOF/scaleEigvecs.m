function[scaledVecs] = scaleEigvecs(eigVecs, eigVals)
%% Scales the eigenvectors by the square root of the eigenvalues.
% A necessary step before EOF rotation.
%
% [scaledVecs] = scaleEigvecs(eigVecs, eigVals)
%
% ----- Inputs -----
%
% eigVecs: A set of eigenvectors for Varimax rotation. Each column is an
%   eigenvector. If eigVecs is 3D, each dim1 x dim2 matrix is a set of
%   eigenvectors.
%
% eigVals: The eigenvalues associated with each eigenvector. Each column of
%   eigVals corresponds to a set of eigenvectors.
%
%
% ----- Outputs -----
%
% scaledVecs: The scaled eigenvectors.

[npoints,nVecs,nSets] = setup(eigVecs, eigVals);

scaledVecs = NaN(npoints, nVecs, nSets);
for k = 1:nSets
    scaledVecs(:,:,k) = eigVecs(:,:,k) .* ...
        repmat( sqrt(eigVals(:,k))', [npoints, 1] );
end

end

%%%%% Helper Functions %%%%%
function[npoints, nVecs, nSets] = setup(eigVecs, eigVals)

% Check eigVecs is 3D
if ndims(eigVecs)>3
    error('eigVecs must be 3D');
end
% Check eigVals is a matrix
if ~ismatrix(eigVals)
    error('eigVals must be a matrix');
end

% Ensure there are no NaNs
if NaNcheck(eigVecs)
    error('eigVecs cannot contain NaN');
end
if NaNcheck(eigVals)
    error('eigVals cannot contain NaN');
end

% Check that the dimensions of the eigvecs and eigvals align
[npoints, nVecs, nSets] = size(eigVecs);
[nEigvals, nValsets] = size(eigVals);

if nSets ~= nValsets || nEigvals ~= nVecs
    error('The size of eigVals and eigVecs do not correctly align');
end
end
    