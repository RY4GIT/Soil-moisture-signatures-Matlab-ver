%% Signature code to characterize three event-scale dynamics in soil moisture timeseries
%% Originally written by Ryoko Araki (raraki8159@sdsu.edu) at San Diego State University
%% April 17, 2022

function [amplitude, risingtime, noresrate] = sig_event(stormtt, ptt, wp, fc, sm, smtt)

% INPUT
% stormtt: the result of storm separation calculated using util_EventSeparation
% ptt: precipitation [mm/timestep] in timetable format
% wp: wilting point calculated using sig_fcwp()
% fc: field capacity value calculated using sig_fcwp()
% sm: timeseries of soil moisture data (only soil moisture values)
% smtt: timeseries of soil moisture data in timetable format

% OUTPUT
% amplitude: normalized event response amplitude signature using field capacity and wilting point
% risingtime: event response time signature
% noresrate: no response rate signature

% Convert timetable -> datenum
smttnum = datenum(smtt.Properties.RowTimes);
pttnum = datenum(ptt.Properties.RowTimes);

%initialize matrix
cumrain = NaN(size(stormtt,1),1);
eventduration = NaN(size(stormtt,1),1);
initialvalue = NaN(size(stormtt,1),1);
amplitude = NaN(size(stormtt,1),1);
maxtime = NaN(size(stormtt,1),1);
risingtime = NaN(size(stormtt,1),1);

% loop for all events
for event = 1:size(stormtt,1)
    
    % obtain the storm features
    sevent = datenum(stormtt(event,1)); % start of the event
    eevent = datenum(stormtt(event,3)); % end of the event
    t_ps = find(pttnum == sevent); % start row of the event in precip time series
    t_pe = find(pttnum == eevent);% end row of the event in precip time series
    cumrain(event,1) = sum(ptt.Var1(t_ps:t_pe));
    eventduration(event,1) = hours(stormtt(event,3)-stormtt(event,1)); % because hourly data
    
    % get the indices of start and the end of event in sm time series, and
    % max time in precip time series
    t_sms = find(smttnum == sevent); % start row of the event in SM time series
    t_sme = find(smttnum == eevent); % end row of the event in SM time series
    maxtime0 = find(ptt.Var1(t_ps:t_pe) == max(ptt.Var1(t_sms:t_sme)));
    maxtime(event,1) = maxtime0(1,1);
    
    % calculate event signatures only if less than 20% of missing values
    if sum(isnan(sm(t_sms:t_sme,:))) < size(sm(t_sms:t_sme,:),1)*0.2
        initialvalue(event,1) = sm(t_sms,:);
        
        % event amplitude normalized by fc-wp
        if (fc-wp) ~= 0
            amplitude(event,1) = (max(sm(t_sms:t_sme,:)) - sm(t_sms,1)) / (fc-wp);
        end
        
        % event rising time
        risingtime0 = find(sm(t_sms:t_sme,:) == max(sm(t_sms:t_sme,:)));
        if ~isempty(risingtime0)
            risingtime(event,1) = risingtime0(1,1) - maxtime(event,1);
        end
        
        % if the SM is still increasing during the event window, expand
        % the window by 12hr until it gets the maxium value (i.e. reach to the event peak),|
        % or it gets to the end of the time series, or it get's more than three days
        addwin = 0;
        if (risingtime0(1,1)-1 == t_sme-t_sms+addwin)
            while (risingtime0(1,1)-1 == t_sme-t_sms+addwin) && (addwin < 24*3) && ((t_sme+addwin+5) < length(sm))
                addwin = addwin + 12;
                maxsm = max(sm(t_sms:t_sme + addwin,:));
                risingtime0 = find(sm(t_sms:t_sme+addwin,:) == maxsm(1,1));
            end
            
            if addwin < 24*3 || ((t_sme+addwin) < length(sm))
                risingtime(event,1) = risingtime0(1,1) - maxtime(event,1);
                amplitude(event,1) = (max(sm(t_sms:t_sme + addwin,:)) - sm(t_sms,1)) / (fc-wp);
            else
                risingtime(event,1) = NaN;
                amplitude(event,1) = NaN;
            end
            
        end
        
        %         % plot per event
        %         figure(event);
        %         subplot(2,1,1)
        %         plot(ptt(eventrange,:).Properties.RowTimes,ptt(eventrange,:).(pvarname));
        %         subplot(2,1,2)
        %         plot(smtt.Properties.RowTimes(stime:etime+addwin),sm(stime:etime+addwin));
        %         if ~isnan(amplitude(event,1))
        %         yline(sm(stime,1)+(fc-wp)*amplitude(event,1));
        %         end
        %         if ~isnan(risingtime(event,1))
        %         xline(smtt.Properties.RowTimes(stime + risingtime0(1,1)-1));
        %         end
        
    end
end

% Count event with no-responses
% no-response is defined as response time = 0 and amplitude = 0
% get the rate of no-response event out of all events
risingtime(risingtime < 0) = 0;
norestime = (risingtime == 0);
noresamp = (amplitude == 0);
noresrate = sum(norestime .* noresamp)/sum(~isnan(amplitude));
risingtime((norestime .* noresamp)==1) = NaN;
amplitude((norestime .* noresamp)==1) = NaN;

end

