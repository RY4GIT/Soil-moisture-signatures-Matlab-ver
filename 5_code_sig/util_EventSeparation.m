% Copied from https://github.com/TOSSHtoolbox/TOSSH/blob/master/TOSSH_code/utility_functions/util_EventSeparation.m
% Gnann, S.J., Coxon, G., Woods, R.A., Howden, N.J.K., McMillan H.K., 2021. TOSSH: A Toolbox for Streamflow Signatures in Hydrology. Environmental Modelling & Software. https://doi.org/10.1016/j.envsoft.2021.104983


function [stormarray] = util_EventSeparation(...
    dates,P,timestep,min_termination,min_duration,min_intensity_hour,...
    min_intensity_day,min_intensity_hour_during,min_intensity_day_during,...
    max_recessiondays,plot_yn)

% util_EventSeparation Take rainfall data and pick out storm periods.
% Hilary McMillan 06/11/2019
% Dept. of Geography, San Diego State University
%
% stormarray = util_EventSeparation(dates,P,timestep,min_termination,...
%    min_duration,min_intensity_hour,min_intensity_day,min_intensity_during,plot_yn)
%
% Inputs
% dates = column vector of datenums
% P = column vector of precipitation
% timestep = time step of precipitation array (hour) 1 = hourly, 24 = daily
% min_termination = Minimum termination time (time between storms) in hour
% min_duration = Minimum duration of storm in hour
% min_intensity_hour = Minimum intensity (per hour)
% min_intensity_day = Minimum intensity (per day)
% min_intensity_hour_during = Minimum timestep intensity allowed during storm
% event without contributing to termination time
% min_intensity_day_during = Minimum timestep intensity allowed during storm
% event without contributing to termination time
% max_recessiondays = maximum number of days to allow recession after rain
% ends
% plot_yn (0/1) = whether to produce graph showing storm events
%
% Outputs
% stormarray = 2-column array with start and end locations of each storm
% event as indices into the precipitation array
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.


% data checks--------------------------------------------------------
if ~ismember(timestep,[0.25,1,24])
    warning('Caution: The event separation function was designed for timesteps of 15min, 1hr or 1day.')
end

if and(timestep<1,1/timestep ~= floor(1/timestep))
    error('Timestep must divide into 1 hour')
end

%------------------------------------------------------------------------

%P is rainfall with timestep (hr) - 15 min (0.25), hr (1), day (24)

%Create moving average series to check hourly/daily intensities ---------
%Create hourly moving average series
if timestep == 1
    P_hr = P;
elseif timestep < 1
    %Calculate size of moving average window for 1 hr
    hr_window = 1/timestep;
    %Append zeros, filter, remove zeros to center filter on each timestep
    P_hr = [P; zeros(floor(hr_window/2),1)];
    P_hr = filter((1/hr_window)*ones(1,hr_window),1,P_hr);
    P_hr = P_hr(1+floor(hr_window/2):end);
end
%Create daily moving average series
if timestep == 24
    P_day = P;
elseif timestep < 24
    day_window = 24/timestep;
    P_day = [P; zeros(floor(day_window/2),1)];
    P_day = filter((1/day_window)*ones(1,day_window),1,P_day);
    P_day = P_day(1+floor(day_window/2):end);
end


% Find gaps between storm events ----------------------------------------
%Storm gaps when hourly rainfall below threshold for time
%greater than min_termination
if timestep <= 1
    %Find all timesteps with hourly rainfall below threshold
    P_lowrain = P_hr <= min_intensity_hour_during;
else
    P_lowrain = P_day <= min_intensity_day_during;
end
P_lowrain(1) = 0;
%Find beginning and end of runs of hourly rainfall below threshold
P_lowrain_change = P_lowrain(2:end)-P_lowrain(1:end-1);
begin_gap = find(P_lowrain_change == 1)+1;
end_gap = find(P_lowrain_change == -1);
%Get complete gaps only
begin_gap = begin_gap(1:length(end_gap));
%Get length of gaps in hours
length_gap = end_gap-begin_gap+1;
%Identify too short gaps
short_gaps = find(length_gap < min_termination/timestep);
for i = 1:length(short_gaps) %delete these short gaps
    P_lowrain(begin_gap(short_gaps(i)):end_gap(short_gaps(i)))=0;
    
end

%Check potential storm periods meet criteria -----------------------------
%Get potential storm periods
%All timesteps where intensity is high enough at hourly/daily timescale
potential_storms = (P_lowrain == 0);
potential_storms(1) = 0;
%Identify runs (consecutive timesteps of rainfall intensity)
potential_storms_change = potential_storms(2:end)-potential_storms(1:end-1);
%Identify beginning and end of storm periods
begin_storm = find(potential_storms_change == 1)+1;
end_storm = find(potential_storms_change == -1);
%Remove last 'beginning of storm' if it does not complete within timeseries
begin_storm = begin_storm(1:length(end_storm));

valid_storm = zeros(length(begin_storm),3);
%Cycle through potential storms and check if valid
for i = 1:length(begin_storm)
    
    %Check duration
    if end_storm(i) - begin_storm(i) + 1 >= min_duration/timestep
        valid_storm(i,1) = 1;
    end
    
    %Check hourly intensity (if timestep <= 1 hr)
    if timestep <= 1
        if max(P_hr(begin_storm(i):end_storm(i))) >= min_intensity_hour
            valid_storm(i,2) = 1;
        end
    end
    
    %Check daily intensity
    if max(P_day(begin_storm(i):end_storm(i))) >= min_intensity_day/(24/timestep)
        valid_storm(i,3) = 1;
    end
end

%Valid storm should have long enough duration and high enough intensity at
%either hourly or daily timescale
valid_overall = and(valid_storm(:,1),or(valid_storm(:,2),valid_storm(:,3)));

%Output array records beginning and end of storm periods -----------------
stormarray = [begin_storm(valid_overall),end_storm(valid_overall)];
%disp(['Number of storm events: ', num2str(size(stormarray,1))]);


%Get suitable end of storm event for signatures that use flow. An event---
%goes for either 5 days after rain, or until the next event starts, or----
%until rainfall is greater than the min_intensity_hour_during-------------

%5-day criterion or until next event
stormarray(1:end-1,3) = min(stormarray(1:end-1,2)+max_recessiondays*24/timestep,stormarray(2:end,2)-1);
stormarray(end,3) = min(stormarray(end,2)+max_recessiondays*24/timestep,length(P));

%Until rainfall is over the maximum
if timestep <= 1
    %Find all timesteps with hourly rainfall below threshold
    P_lowrain = P_hr > min_intensity_hour_during;
else
    P_lowrain = P_day > min_intensity_day_during;
end

recession_rain = zeros(size(stormarray,1),1);
for i = 1:size(stormarray,1)
    rain_index = find(P_lowrain((stormarray(i,2)+1):stormarray(i,3))==1,1,'first');
    if numel(rain_index) > 0    
        recession_rain(i) = rain_index-1;
    else
        recession_rain(i) = inf;
    end
end
stormarray(:,3) = min(stormarray(:,3),stormarray(:,2)+recession_rain);


%Plotting ----------------------------------------------------------------
%If plotting requested, show rainfall with overlaid storm events
%Warning - plotting is slow
if plot_yn
    disp(['Please wait a few seconds: plotting is slow'])
    dates_dt = datetime(dates,'ConvertFrom','datenum');
    figure('Position',[300 400 1400 400])
    h3=plot(dates_dt(:),P(:),'b-'); %Sets axes to correct length
    hold on
    for i = 1:length(stormarray)
        h1=fill([dates_dt(stormarray(i,1)),dates_dt(stormarray(i,1)),...
            dates_dt(stormarray(i,2)),dates_dt(stormarray(i,2))],...
            [0, 60, 60, 0],'g');
        h2=fill([dates_dt(stormarray(i,2)),dates_dt(stormarray(i,2)),...
            dates_dt(stormarray(i,3)),dates_dt(stormarray(i,3))],...
            [0, 60, 60, 0],'c');        
    end
    h3=plot(dates_dt(:),P(:),'b-');
    legend([h3,h1,h2],{'Rainfall (mm)','Storm Period','Recession Period'})
end