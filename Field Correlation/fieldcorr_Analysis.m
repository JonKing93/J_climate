function[s] = fieldcorr_Analysis(ts, field, varargin)
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
% [s] = fieldcorr_Analysis(..., 'fdr', q, fdrType)
% Restricts significant correlations using a false detection rate procedure.
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
%
% ----- Inputs -----
%
% ts: a time series. This is a single vector.
%
% field: A field. The time series at individual points in the field should
%       have the same length as ts. Field may be n-dimensional. By default,
%       fieldcorr assumes that the time series runs parallel to the first dimension.
%
% d: A scalar specifying the dimension of the field on which to perform the
%       correlation.
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
% q: The rate of false discovery. This is the percent of null hypotheses
%       that are falsely rejected and attributed as significant. (q = 0.05
%       is a commonly used value.) q must be on the interval (0,1)
%
% fdrType: A flag for the desired fdr procedure
%       'BH': The Benjamini-Hochberg procedure. Guaranteed to control the 
%               false detection rate when all data are independent or
%               positively correlated. Moderately conservative.
%       'BY': The Benjamini-Yekutieli procedure. Guaranteed for any type of
%               data dependency. A very conservative test.
%       'BKY': Benjamini-Krieger-Yekutieli. Guaranteed only for independent
%               datasets. Uses a two-stage process to provide a less overly
%               conservative approach than BH and BY.






