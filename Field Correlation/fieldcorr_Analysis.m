function[s] = fieldcorr_Analysis(ts, field, p, varargin)
%% Computes the correlation of a time series with a field and tests the 
% significance of the correlations. Includes options to consider the
% effects of multiple comparisons, spatial inter-correlation, lag-1
% autocorrelation, and false detection rate.
%
%
% [s] = fieldcorr_Analysis(ts, field, p, MC, noiseType)
% Calculates the correlation coefficients and associated p-values between
% a time series and n-dimensional field. Considers the effect of multiple
% hypothesis testing and uses a Monte Carlo test to consider spatial
% interdependence to determine significant correlations.
%
% [s] = fieldcorr_Analysis(ts, field, p, 'noSpatial')
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
% [s] = fieldcorr_Analysis(..., 'lagArgs', {lagArgs})
% Runs the lagged field correlation using the 'fixedN' or 'restrictBounds'
% arguments. See the help secton of "lagFieldcorr" for details.
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
% p: A vector containing significance levels to be tested. Each p value 
%       must be on the interval (0,1).
%
% MC: The number of Monte Carlo iterations to perform in the
%       spatial-interdependence test.
%
% noiseType: A flag for the type of noise to use in the Monte Carlo test
%       'white': White gaussian noise
%       'red': Lag-1 autocorrelated noise with added white noise
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
% lagArgs: A cell containing additional parameters for a lagged field
%       correlation. See the documentation of the "lagFieldcorr" 'fixedN' 
%       and 'restrictBounds' flags
%
% d: A scalar specifying the dimension of the field on which to perform the
%       correlation.
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
%
% ----- Outputs -----
% 
% s: a structure containing the following fields
%   corrmaps: The correlation coefficients at each field point. Separate
%       lags are along the last dimension.
%
%   pmaps: The p-value for each correlation coefficient.
%
%   corrSize: The number of points used in the calculation of the
%       correlation coefficients for each lag.
%
%   tsPoints: The points in the time series used for each lag. Each column
%       records the points for one lag.
%
%   fPoints: The points for the field used for each lag. Each column
%       records the points for one lag.
%
%   p: The tested significance levels.
%    
%   nPass: The number of p-values exceeding the significance levels for
%       each lag. Rows are significance levels, and columns are lags.
%
%   N: The total number of hypothesis tests for each significance level and
%       lag. As previously, rows are significance levels, and columns are lags.
%
%   nPassG: The number of hypothesis tests globally (all lags), that pass
%       each significance level.
%
%   NG: The global (all lags) number of hypothesis tests for each
%       significance level.
%
%   areSigM: 
%
%   nNeededM
%
%   areSigS: A logical indicating whether each lagged correlation maintains
%       significance at each significance level.
%
%   nNeededS: The number of passed tests needed to maintain signficance at
%       each significance level and lag, given spatial interdependence.
%
%   iterNPassed: The number of passed tests needed to mainain significance
%       at each Monte Carlo iteration.
%
%   iterTrueConf: The true confidence level of each Monte Carlo iteration.
%
%   areSig: 
%
%   sigP: The indices of p-values that pass the significance tests
%
%   q: The tested false detection rates
%
%   fdrSigP: The indices of p-values that pass the false detection rate test
%
%   metadata: Metadata for the analysis%   
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
[dim, corrArgs, iTS, iField, lagArgs, spatialTest, MC, noiseType, ...
    convergeFlag, fdrTest, q, fdrType, noLags] = parseInputs(varargin{:});

% Reshape field to 2D
[field, dSize, dOrder] = dimNTodim2(field, dim);

% Run the analysis with no lags
if noLags
    % Run the field correlation
    [s.corrmaps, s.pmaps] = fieldcorr(ts, field, 'corrArgs', corrArgs);
    
    % Get the initial number of passed significance tests
    [s.nPass, s.N] = getNPassN(p, s.pmaps);
    
    % Account for mutliple comparisons
    [s.areSigM, s.nNeededM] = multCompSig(p, s.N, s.nPass);
    
    % Account for spatial interdependence
    if spatialTest
        [s.areSigS, s.nNeededS, s.iterNPassed, s.iterTrueConf] = ...
            MC_fieldcorr(ts, field, MC, noiseType, p, s.nPass, convergeFlag, 'corrArgs', corrArgs);
    end
    

% Run the analysis with lags
else
    % Run the lagged field correlation
    [s.corrmaps, s.pmaps, s.corrSize, s.tsPoints, s.fPoints] = ...
        lagFieldcorr(ts, iTS, field, iField, lagArgs{:}, 'corrArgs', corrArgs);
    
    % Get the global number of passed tests
    [s.nPassG, s.NG] = getNPassN(p, s.pmaps);
    % Get the number of passed significance tests for each lag
    [s.nPass, s.N] = getNPassN(p, s.pmaps, 3);
      
    % Test significance while accounting for multiple comparisons
    if any( s.N ~= s.N(1) )
        [s.areSigM, s.nNeededM] = multCompSig(p, s.N, s.nPass);
    else
        % All N are the same, can be used to improve runtime
        [s.areSigM, s.nNeededM] = multCompSig(p, s.N(1), s.nPass);
    end
    
    % Account for spatial interdependence
    if spatialTest
        % Preallocate at all lags
        nlags = size(s.tsPoints,2);
        s.areSigS = NaN(length(p), nlags);
        s.nNeededS = NaN(length(p), nlags);
        s.iterNPass = NaN(MC, length(p), nlags);
        s.iterTrueConf = NaN(MC, length(p), nlags);
        
        % Run MC_fieldcorr at each lag
        for k = 1:nlags
            [s.areSigS(:,k), s.nNeededS(:,k), s.iterNPassed(:,:,k), s.iterTrueConf(:,:,k)] = ...
                MC_fieldcorr(ts(s.tsPoints(:,k)), field(s.fPoints(:,k),:), MC, noiseType, p, s.nPass(:,k), convergeFlag, 'corrArgs', corrArgs);
        end
    end
end

% Get the sets of p-values that remain significant under all tests
if spatialTest
    s.areSig =  ( s.areSigM & s.areSigS );
else
    s.areSig = s.areSigM;
end

% Get the indices of significant p-values in significant sets
s.sigP = cell( size(s.areSig) );
for j = 1:size(s.pmaps, 3)
    for k = 1:length(p)
        if s.areSig(k,j)
            s.sigP{k,j} =  ( s.pmaps(:,:,j) <= p(k) );
            s.sigP{k,j} = dim2TodimN( s.sigP{k,j}, [1 dSize(2:end)], dOrder);
        end
    end
end

% Apply a false detection rate procedure
if fdrTest
    s.fdrSigP = cell(length(q), size(s.pmaps,3),1);
    for j = 1:size(s.pmaps,3)
        for k = 1:length(q)
            [~, s.fdrSigP{k,j}] = fdr( s.pmaps(:,:,j), q(k), fdrType);
            s.fdrSigP{k,j} = dim2TodimN( s.fdrSigP{k,j}, [1 dSize(2:end)], dOrder);
        end
    end
end
    
% Reshape to original dimensions
s.corrmaps = dim2TodimN(s.corrmaps, [1 dSize(2:end) size(s.pmaps,3)], [dOrder, max(dOrder)+1] );
s.pmaps = dim2TodimN(s.pmaps, [1 dSize(2:end) size(s.pmaps,3)], [dOrder, max(dOrder)+1] );

% Get metadata
s.metadata = [{'Spatial Test', 'MC', 'Noise Type', 'MC Convergence Test', ...
    'FDR Test', 'FDR Type', 'Lag Args', 'corr Args'};...
    {spatialTest, MC, noiseType, convergeFlag, fdrTest, fdrType, lagArgs, corrArgs}];
s.p = p;
s.q = q;

end

% ----- Helper Functions -----

function[dim, corrArgs, iTS, iField, lagArgs, spatialTest, MC, noiseType, ...
    convergeFlag, fdrTest, q, fdrType, noLags] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
dim = 1;
corrArgs = {'type', 'Pearson'};

iTS = NaN;
iField = NaN;
tsLags = [];
fieldLags = [];

spatialTest = true;
MC = NaN;
noiseType = NaN;
p = NaN;
convergeFlag = 'convergeTest';

fdrTest = false;
q = NaN;
fdrType = NaN;

noLags = false;
otherLagArg = {};

if ~isempty(inArgs)
   
    % Set flag switches
    isdim = false;
    iscorrargs = false;
    
    isits = false;
    isifield = false;

    isnoisetype = false;

    isq = false;
    isfdrtype = false;

    istslags = false;
    isfieldlags = false;
    bothLags = false;
    
    islagarg = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        % True switches
        if isnoisetype
            noiseType = arg;
            isnoisetype = false;
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
        elseif islagarg
            otherLagArg = arg;
            islagarg = false;
        elseif isdim
            dim = arg;
            isdim = false;
        elseif iscorrargs
            corrArgs = arg;
            iscorrargs = false;
               
        % Switch triggers
        elseif k==1 && ~strcmpi(arg, 'noSpatial')
            if length(inArgs)>= k+1
                MC = arg;
                isnoisetype = true;
                spatialTest = true;
            else 
                error('NoiseType is not specified');
            end
            
        elseif k==1 && strcmpi(arg, 'noSpatial')
            spatialTest = false;
            
        elseif strcmpi(arg, 'lagArgs')
            if length(inArgs)>=k+1
                islagarg = true;
            else
                error('lagArgs not specified');
            end
            
        elseif strcmpi(arg, 'fdr')
            fdrTest = true;
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
            convergeFlag = arg;
            
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

% Get the full set of lag arguments
lagArgs = {};
if isempty( [tsLags, fieldLags])
    noLags = true;
else
    
    if ~isempty(tsLags)
        lagArgs(end+1:end+2) = {'lagTS', tsLags};
    end
    if ~isempty(fieldLags)
        lagArgs(end+1:end+2) = {'lagF', fieldLags};
    end
    lagArgs = [lagArgs, otherLagArg];   
end

end

