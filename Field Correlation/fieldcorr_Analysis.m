function[s] = fieldcorr_Analysis(ts, field, varargin)
%% Computes the correlation of a time series with a field and tests the 
% significance of the correlations. Includes options to consider the
% effects of multiple comparisons, spatial inter-correlation, lag-1
% autocorrelation, and false detection rate.
%
%
% [s] = fieldcorr_Analysis(ts, field, MC, noiseType, p)
% Calculates the correlation coefficients and associated p-values between
% a time series and n-dimensional field. Considers the effect of multiple
% hypothesis testing and uses a Monte Carlo test to consider spatial
% interdependence to determine significant correlations.
%
% [s] = fieldcorr_Analysis(ts, field, 'noSpatial')
% Blocks the spatial interdependence test.
%
% [s] = fieldcorr_Analysis(..., 'fdr', q, fdrType)
% Further restricts significant correlations using a false detection rate procedure.
%
% [s] = fieldcorr_Analysis(..., 'lagTS', tsLags, iTS, iField)
% Performs all tests at specified lags / leads of the time series.
%
% [s] = fieldcorr_Analysis(..., 'lagField', fieldLags, iTS, iField)
% Performs all tests at specified leads / lags of the field.
%
% [s] = fieldcorr_Analysis(..., 'lagTS', tsLags, 'lagField', fieldLags, iTS, iField)
% Performs all tests on leads and lags of both the time series and field.
%
% [s] = fieldcorr_Analysis(..., 'noConvergeTest')
% Blocks the recording of Monte Carlo convergence data. May improve runtime
% for very large analyses.
%
% [s] = fieldcorr_Analysis(..., 'fieldDim', d)
% Performs the correlations along the specified field dimension. By
% default, correlation coefficiens are computed along the first field dimension.
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
%
% tsLags: The desired lags (negative integers) and leads (positive integers)
%   at which to compute correlation coefficients. (e.g. lags = [-3 0 10]
%   will calculate the correlation with a lag of 3 time indices, no lags, and
%   lead of 10 indices for the time series.)
%
% fLags: The desired lags / leads for the field.
%
% iTS: a set of time indices for time series. See the assignLagIndices function.
%
% iField: A set of time indices for the field. These must overlap
%    with the tsIndices. See the assignLagIndices function.
%
% d: A scalar specifying the dimension of the field on which to perform the
%       correlation.
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
%
% ----- References -----
%
% See Contents folder
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu


% Parse Inputs
[runSpatial, MC, noiseType, p, q, fdrType, noLags, lagArgs, nlags, testConverge, dim, corrArgs] = ...
    parseInputs(varargin{:});

% Declare the initial structure
s = struct();

% Make the field 2D along the dimension of interest
[field, dSize, dOrder] = dimNTodim2(field, dim);

% Run the field correlation at all lags
if isempty(lagArgs)
    [s.corrmaps, s.pmaps] = fieldcorr(ts, field, 'corrArgs', corrArgs{:});
else 
    [s.corrmaps, s.pmaps, s.N] = lagFieldcorr(ts, iTS, field, iField, lagArgs{:}, 'corrArgs', corrArgs);
end

% Global (at all lags) significance test for multiple comparisons
nTests = numel(s.pmaps(~isnan(s.pmaps)));
nPass = numel(  s.pmaps(~isnan(s.pmaps)) < p);
s.isSigGlobal = finiteTestsAreSig(nTests, p, nPass);

% Local (for each individual lag) significance test for multiple comparison
s.isSig = NaN(nlags,1);
for k = 1:nlags
    nPass = numel(s.pmaps
    s.isSig(k) = finiteTestsAreSig(





end

% ----- Helper Functions -----

function[runSpatial, MC, noiseType, p, q, fdrType, tsLags, fieldLags, iTS, ...
    iField, testConverge, dim, corrArgs] = parseInputs(varargin)

inArgs = varargin;

% Set defaults
runSpatial = NaN;
MC = NaN;
noiseType = NaN;
p = NaN;
q = NaN;
fdrType = NaN;
tsLags = [];
fieldLags = [];
iTS = NaN;
iField = NaN;
testConverge = true;
dim = 1;
corrArgs = {'type', 'Pearson'};
bothLags = false;

if ~isempty(inArgs)
   
    % Set flag switches
    isnoisetype = false;
    isp = false;
    isq = false;
    isfdrtype = false;
    istslags = false;
    isfieldlags = false;
    isits = false;
    isifield = false;
    isdim = false;
    iscorrargs = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        % True switches
        if isnoisetype
            noiseType = arg;
            isnoisetype = false;
            isp = true;
        elseif isp
            p = arg;
            isp = false;
        elseif isq
            q = arg;
            isq = false;
            isfdrtype = true;
        elseif isfdrtype
            fdrType = arg;
            isfdrtype = false;
        elseif istslags
            tsLags = arg;
            istslags = false;
            if ~bothLags
                isits = true;
            end
        elseif isfieldlags
            fieldLags = arg;
            isfieldlags = false;
            isits = true;
        elseif isits
            iTS = arg;
            isits = false;
            isifield = true;
        elseif isifield
            iField = arg;
            isifield = false;
        elseif isdim
            dim = arg;
            isdim = false;
        elseif iscorrargs
            corrArgs = arg;
            iscorrargs = false;
               
        % Switch triggers
        elseif k==1 && ~strcmpi(arg, 'noSpatial')
            if length(inArgs)>= k+2
                isnoisetype = true;
                runSpatial = true;
            else 
                error('noiseType and p are not specified');
            end
            
        elseif k==1 && strcmpi(arg, 'noSpatial')
            runSpatial = false;
            
        elseif strcmpi(arg, 'fdr')
            if length(inArgs)>=k+2
                isq = true;
            else
                error('q and fdrType are not specified');
            end
            
        elseif strcmpi(arg, 'lagTS')
            if length(inArgs)>=k+3
               if strcmpi(inArgs{k+2},'lagField')
                   bothLags = true;
                   istslags = true;
               else
                   istslags = true;
               end
            else
                error('tsLags, iTS and iField are not all specified');
            end
            
        elseif strcmpi(arg, 'lagField')
            if length(inArgs)>=k+3
                isfieldlags = true;                
            else
                error('fieldLags, iTS and iField are not all specified');
            end
    
        elseif strcmpi(arg, 'noConvergeTest')
            testConverge = false;
            
        elseif strcmpi(arg, 'fieldDim')
            isdim = true;
            
        elseif strcmpi(arg, 'corrArgs')
            iscorrargs = true;
            
        else
            error('Unrecognized Input');
        end
    end
else
    error('Insufficient Inputs');
end

end

