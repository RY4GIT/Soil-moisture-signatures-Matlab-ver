% Copied from https://github.com/TOSSHtoolbox/TOSSH/blob/master/TOSSH_code/utility_functions/util_RisingLimbs.m
% Gnann, S.J., Coxon, G., Woods, R.A., Howden, N.J.K., McMillan H.K., 2021. TOSSH: A Toolbox for Streamflow Signatures in Hydrology. Environmental Modelling & Software. https://doi.org/10.1016/j.envsoft.2021.104983

function [flow_section] = util_RisingLimbs(Q, t, varargin)
%util_RecessionSegments Identify all rising limbs.
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%	OPTIONAL
%   rising_limb_length: length of rising limbs (days), default = 1
%   eps: allowed decrease in flow during rising limb, default = 0
%   plot_results: whether to plot results, default = false
%   minimum_peak: minimum peak to be counted as rising limb
%
%   OUTPUT
%   flow_section: n-by-2 array where n is the number of recession segments
%   columns are the indices into the flow array of the start and end of the recession segments
%
%   EXAMPLE
%   % example data [matlab_date flow]
%   data = load('./TOSSH_development/example/example_data/Q_test.mat'); 
%   Q = data.Q_test(:,2); 
%   t = data.Q_test(:,1); t = datetime(t,'ConvertFrom','datenum');
%   flow_section = util_RisingLimbs(Q, t);
%   flow_section = util_RisingLimbs(Q, t, 'rising_limb_length', 2, 'plot_results', true);
%
%   References
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
addParameter(ip, 'rising_limb_length', 1, @isnumeric) % length of increasing flow in days to be declared a rising limb
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'minimum_peak', median(Q,'omitnan'), @isnumeric) % minimum peak to be counted as rising limb
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
rising_limb_length = ip.Results.rising_limb_length;
plot_results = ip.Results.plot_results;
eps = ip.Results.eps;
minimum_peak = ip.Results.minimum_peak;

% Identify all individual rising limbs with length > rising_limb_length days.
% how many increasing timesteps depends on length of timestep
len_increase = rising_limb_length/days(t(2)-t(1));
% hind timesteps with increasing flow
increasing_flow = Q(2:end)>(Q(1:end-1)-eps);
% start on a non-increasing point 
start_point = find(increasing_flow==0,1);
increasing_flow = increasing_flow(start_point:end);
% find start and end of increasing sections
flow_change = find(increasing_flow(1:end-1) ~= increasing_flow(2:end));
% reshape into x by 2 array (columns = start, end of decrease)
flow_change = flow_change(1:(2*floor(size(flow_change,1)./2)));
flow_change = reshape(flow_change,2,[]).';
% find sections
flow_section = flow_change((flow_change(:,2)-flow_change(:,1))>=len_increase,:);
flow_section = flow_section+start_point;
flow_section(:,1) = flow_section(:,1); % move start point n days
% remove rising limbs which have a peak lower than minimum_peak
% flow_section((Q(flow_section(:,2)) < minimum_peak),:) = [];
% remove rising limbs which have a peak smaller than minimum_peak
nonsig_peaks = find( (Q(flow_section(:,2)) - Q(flow_section(:,1))) < minimum_peak);
flow_section(nonsig_peaks,:) = [];

% if numel(flow_section)==0
%     % error('No long enough rising limbs, consider setting eps parameter > 0')
% end
% 
% % if plot_results
%     figure('Position',[100 100 1000 400]); hold on;
%     h1=plot(t,Q);
%     for i = 1:size(flow_section,1)
%         h2=plot(t(flow_section(i,1):flow_section(i,2)),Q(flow_section(i,1):flow_section(i,2)),'r-');
%     end
%     h3=plot(t,minimum_peak.*ones(size(t)),'k--');
%     title('Selected rising limbs')
%     legend([h1 h2 h3],{'Full flow series', 'Selected rising limbs', 'Minimum peak threshold'})
% %     datetick('x')
%     ylabel('Flow')
% % end

end