% Copied from https://github.com/TOSSHtoolbox/TOSSH/blob/master/TOSSH_code/utility_functions/util_DataCheck.m
% Gnann, S.J., Coxon, G., Woods, R.A., Howden, N.J.K., McMillan H.K., 2021. TOSSH: A Toolbox for Streamflow Signatures in Hydrology. Environmental Modelling & Software. https://doi.org/10.1016/j.envsoft.2021.104983

function [DataCheck_Result, timestep, t] = util_DataCheck(Q, t, varargin)
%util_DataCheck Checks data for various things (e.g. unrealistic values).
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%   OPTIONAL
%   P: precipitation [mm/timestep]
%   PET: potential evapotranspiration [mm/timestep]
%
%   OUTPUT
%   DataCheck_Result: 0: Good Data, 1: Warning, 2: Error
%   timestep: (median) timestep of data
%   t: time in Matlab datetime format
%
%   EXAMPLE
%   % example data [matlab_date flow]
%   data = load('./TOSSH_development/example/example_data/Q_test.mat');
%   Q = data.Q_test(:,2);
%   t = data.Q_test(:,1);
%   data = load('./TOSSH_development/example/example_data/P_test.mat');
%   P = data.P_test(:,2);
%   dataCheck = util_DataCheck(Q, t, 'P', P);
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

% check input parameters
if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
% date time series has to be numeric or datetime and either a (n,1) or a (1,n) vector
addRequired(ip, 't', @(t) (isnumeric(t) || isdatetime(t)) && (size(t,1)==1 || size(t,2)==1))

% optional input arguments
addParameter(ip, 'P', [], @(P) isnumeric(P) && (size(P,1)==1 || size(P,2)==1)) % P has to be numeric and either a (n,1) or a (1,n) vector
addParameter(ip, 'PET', [], @(PET) isnumeric(PET) && (size(PET,1)==1 || size(PET,2)==1)) % PET has to be numeric and either a (n,1) or a (1,n) vector

parse(ip, Q, t, varargin{:})
P = ip.Results.P;
PET = ip.Results.PET;

% default setting reads as good data
DataCheck_Result = 0;

% data checks
if min(Q)<0
    DataCheck_Result = 2;
    error('Negative values in flow series.')
end

if length(Q) ~= length(t)
    DataCheck_Result = 2;
    error('Flow series and time vector have different lengths.')
end

if any(isnan(Q))
    DataCheck_Result = 1;
    % warning('Ignoring NaNs in streamflow data.')
end

if length(Q) < 30
    DataCheck_Result = 1;
    % warning('Extremely short time series.')
end

% timestep checks
if isnumeric(t)
    DataCheck_Result = 1;
    t = datetime(t,'ConvertFrom','datenum');
    % warning('Converted datenum to datetime.')
end

timesteps = diff(t);
timestep = median(timesteps);
if any(diff(timesteps)~=0)
     warning('Record is not continuous (some timesteps are missing).')
end

% optionally check P
if ~isempty(P)
    
    if any(isnan(P))
        DataCheck_Result = 1;
        % warning('Ignoring NaNs in precipitation data.')
    end
    
    if length(Q) ~= length(P)
        DataCheck_Result = 2;
        % error('Precipitation and flow series have different lengths.')
    end
    
    if min(P)<0
        DataCheck_Result = 2;
        % error('Negative values in precipitation series.')
    end
    
end

% optionally check PET
if ~isempty(PET)
    
    if any(isnan(PET))
        DataCheck_Result = 1;
        % warning('Ignoring NaNs in potential evpotranspiration data.')
    end
    
    if length(Q) ~= length(PET)
        DataCheck_Result = 2;
        % error('Potential evpotranspiration and flow series have different lengths.')
    end
    
    if min(PET)<0
        DataCheck_Result = 2;
        % error('Negative values in potential evpotranspiration series.')
    end
    
end

end