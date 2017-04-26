function[s] = fieldcorr_Analysis(ts, field)
%% Computes the correlation of a time series with a field and tests the 
% significance of the correlations.
%
% [s] = fieldcorr_Analysis(ts, field)
% Calculates the correlation of a time series and an n-dimensional field.
% Checks that 
%
% [s] = fieldcorr_Analysis(ts, field, MC, noiseType, p)
%
% [s] = fieldcorr_Analysis(..., 'lead', leads, iTS, iField)
%
% [s] = fieldcorr_Analysis(..., 'noSpatial')
%
% [s] = fieldcorr_Analysis(..., 'noConvergeTest')
%
% [s] = fieldcorr_Analysis(..., 'fieldDim', d)
%
% [s] = fieldcorr_Analysis(..., 'corrArgs', {corrArgs})
%
% [s] = fieldcorr_Analysis(..., 'fdr', ...........? )
%



% Perform the field correlation at the desired lags

