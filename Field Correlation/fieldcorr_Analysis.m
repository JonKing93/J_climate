function[s] = fieldcorr_Analysis(ts, field)
%% Computes the correlation of a time series with a field and tests the 
% significance of the correlations. Includes options to consider the
% effects of multiple comparisons, spatial inter-correlation, lag-1
% autocorrelation, and false detection rate.
%
% [s] = fieldcorr_Analysis(ts, field)
% Calculates the correlation of a time series and an n-dimensional field.
% Checks that 
%
% [s] = fieldcorr_Analysis(ts, field, MC, noiseType, p)
%
% [s] = fieldcorr_Analysis(..., 'leadTS', tsLeads, iTS, iField)
% Shifts the time series by the desired number of leads/lags, and
% performs the field correlation and significance tests on each lead/lag.
%
% [s] = fieldcorr_Analysis(..., 'leadField', fieldLeads, iTS, iField)
% Shifts the field by the desired number of leads/lags, and performs the
% field correlation and significance tests on each lead/lag. May be used
% simultaneously with the 'leadTS' flag.
%
% [s] = fieldcorr_Analysis(..., 'noSpatial')
% Blocks the spatial significance test.
%
% [s] = fieldcorr_Analysis(..., 'noConvergeTest')
% Blocks the recording of Monte Carlo convergence data. May improve runtime
% for very large analyses.
%
% [s] = fieldcorr_Analysis(..., 'fieldDim', d)
% Performs the correlations along the specified field dimension. By
% default, correlation are computed along the first field dimension.
%
% [s] = fieldcorr_Analysis(..., 'corrArgs', {corrArgs})
% Performs the analysis using alternative methods of correlation. See the
% MATLAB help page on "corr" for details.
%
% [s] = fieldcorr_Analysis(..., 'fdr', ...........? )
%



% Perform the field correlation at the desired lags

