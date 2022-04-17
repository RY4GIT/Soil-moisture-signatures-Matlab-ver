% Copied from https://github.com/TOSSHtoolbox/TOSSH/blob/master/TOSSH_code/signature_functions/sig_RisingLimbDensity.m
% Gnann, S.J., Coxon, G., Woods, R.A., Howden, N.J.K., McMillan H.K., 2021. TOSSH: A Toolbox for Streamflow Signatures in Hydrology. Environmental Modelling & Software. https://doi.org/10.1016/j.envsoft.2021.104983

function [RLD] = sig_RLD(Q, t, varargin)
%sig_RLD Calculates rising limb density (RLD).
%   Calculates the rising limb density, the ratio between the number of
%   rising limbs and the total amount of timesteps the hydrograph is
%   rising.
%
%	INPUT
%   Q: streamflow [mm/timestep]
%   t: time [Matlab datetime]
%	OPTIONAL
%   rising_limb_length: length of increasing flow section (days) to be 
%   declared a rising limb
%   eps: allowed increase in flow during recession period
%   plot_results: whether to plot results
%   minimum_peak: added
%
%   OUTPUT
%   RLD: rising limb density [1/timestep]
%   % TODO: rising_limb_month: approx. month of rising limb
%
%   EXAMPLE
%   % example data [matlab_date flow]
%   data = load('./TOSSH_development/example/example_data/Q_test.mat'); 
%   Q = data.Q_test(:,2); 
%   t = data.Q_test(:,1); 
%   RLD = sig_RLD(Q,t);
%   RLD = sig_RLD(Q,t,'plot_results',true,'rising_limb_length',2);
%
%   References
%   Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
%   2011. Catchment classification: empirical analysis of hydrologic
%   similarity based on catchment function in the eastern USA. Hydrology
%   and Earth System Sciences, 15(9), pp.2895-2911.
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
addParameter(ip, 'rising_limb_length', 1, @isnumeric)% length of increasing flow section (days) to be declared a rising limb
addParameter(ip, 'eps', 0, @isnumeric) % allowed increase in flow during recession period
addParameter(ip, 'minimum_peak', median(Q,'omitnan'), @isnumeric) % minimum peak to be counted as rising limb
addParameter(ip, 'plot_results', false, @islogical) % whether to plot results (2 graphs)

parse(ip, Q, t, varargin{:})
rising_limb_length = ip.Results.rising_limb_length;
plot_results = ip.Results.plot_results;
eps = ip.Results.eps;
minimum_peak = ip.Results.minimum_peak;

% data checks
[~,~,t] = util_DataCheck(Q, t);

% calculate signature
[flow_section] = util_RisingLimbs(Q, t, 'eps', eps, 'minimum_peak', minimum_peak, ...
    'rising_limb_length', rising_limb_length, 'plot_results', false);
if numel(flow_section)==0
    RLD = NaN;
end
RLD = 1./mean(flow_section(:,2)-flow_section(:,1));

% TODO: return individual rising limb lengths
% get rising limb month
% date_tmp = datevec(floor(mean(flow_section,2)));
% rising_limb_month = date_tmp(:,2);

% if plot_results
%     figure
%     hold on
%      colour_mat_seasons = [...
%         0 0 1;  0 0 1;...
%         0 1 0; 0 1 0; 0 1 0;...
%         1 0 0; 1 0 0; 1 0 0;...
%         1 1 0; 1 1 0; 1 1 0; ...
%         0 0 1];
%     p1=plot(0,0,'.','Color',[0 1 0]); 
%     p2=plot(0,0,'.','Color',[1 0 0]);  
%     p3=plot(0,0,'.','Color',[1 1 0]); 
%     p4=plot(0,0,'.','Color',[0 0 1]);
%     for i = 1:size(flow_section,1)
%         RL = [flow_section(i,1):flow_section(i,2)]'; % get rising limb
%         Q_tmp = Q(RL);
%         t_tmp = t(RL);
%         date_vec = datevec(t_tmp);
%         ind = floor(median(date_vec(:,2))); % get approx. month
%         plot(1:length(Q_tmp),Q_tmp,'-','color',colour_mat_seasons(ind,:),'linewidth',1)
%     end
%     xlabel('Time since start of rising limb') % [d]
%     ylabel('Q') % [mm/d]
%     set(gca,'YScale','log')
%     legend([p1 p2 p3 p4],{'MAM','JJA','SON','DJF'},'box','off','Location','best');
% end

% old
%{
differences = diff(Q); % differences between one value and the next

isRisingLimb = false(size(differences));
isPeak = false(size(differences));

for i = min_RL+1:length(differences)-1
    if all(differences(i-min_RL+1:i)>0)
        isRisingLimb(i-min_RL+1:i) = true;
    elseif all(differences(i-min_RL:i-1)>0) && differences(i)<=0
        isRisingLimb(i) = false;
        isPeak(i) = true;
    else
    end
end

up = differences(isRisingLimb);
peaks = differences(isPeak);

% test plot to check rising limbs
%{
figure; plot(t-min(t),Q); hold on
plot(t(isRisingLimb)-min(t),Q(isRisingLimb),'.')
plot(t(isPeak)-min(t),Q(isPeak),'.')
xlabel('Days')
ylabel('Flow [mm/h]')
xlim([554 574])
%}

try % try & catch structure in case tmp_peak is undefined or tmp_up is 0
    RLD = length(peaks)/length(up); % number of peaks / by number of time steps that flow goes up
catch
    RLD = NaN;
end

%}

end
