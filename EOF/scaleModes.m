function[scaledVecs] = scaleModes(modes, eigVals)
%% Scales the EOF modes by the square root of the eigenvalues.
% A necessary step before EOF rotation.
%
% [scaledModes] = scaleModes(modes, eigVals)
%
%
% ----- Inputs -----
%
% modes: A set of EOF modes for Varimax rotation. Each column is a mode
%
% eigVals: A vector containing the eigenvalues associated with each mode.
%
%
% ----- Outputs -----
%
% scaledModes: The scaled eigenvectors.

[npoints,nVecs,nSets] = setup(modes, eigVals);

scaledVecs = NaN(npoints, nVecs, nSets);
for k = 1:nSets
    scaledVecs(:,:,k) = modes(:,:,k) .* ...
        repmat( sqrt(eigVals(:,k))', [npoints, 1] );
end

end

%%%%% Helper Functions %%%%%
function[npoints, nVecs, nSets] = setup(modes, eigVals)

% Check eigVecs is a matrix
if ~ismatrix(modes)
    error('eigVecs must be a matrix');
end

% Check eigVals is a vector
if ~isvector(eigVals)
    error('eigVals must be a vector');
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
    