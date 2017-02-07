W = 243;
M = floor(W/11);
nPoints = 11775;
nChunks = floor(nPoints / W);
Data = Srm0;

chunkArray = NaN(W, nChunks);
for k = 1:nChunks
    chunkArray(:,k) = Sr( k*W-(W-1): k*W);
end

chunkArray = detrend(chunkArray,'constant');

[eigvals, eigvecs] = simpleSSA(chunkArray, M, 'VG');

tic
surrEigvals = MC_SSA(chunkArray, eigvecs, 1000, 'red', M, 'VG');
toc

[sigEigdex, highsurrtail, ~] = sigTestMCSSA(0.95, eigvals, surrEigvals);

[maxFreq, maxPeriod] = maxFreqPeriod(