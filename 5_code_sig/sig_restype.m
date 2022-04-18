%% Signature code to calculate the response type of the timeseries
%% Originally written by Inge Wiekenkamp (inge.wiekenkamp@gfz-potsdam.de), and
%% simplified for this analysis by Ryoko Araki (raraki8159@sdsu.edu) at San Diego State University
%% April 17, 2022

function   [restype, timestamp] = sig_restype(sm_n, sm_n_ts, ptt, stormtt, noise_level, increment_level, plot_yn)

% INPUT
% sm_n: an array of soil moisture timeseries (soil moisture values, n depth combined together)
% sm_n_ts: an array of soil moisture timeseries (datenum values, n depth combined together)
% ptt: precipitation timeseries in datetime format
% stormtt: precipitation [mm/timestep] in timetable format
% noise_level: not used 
% increment_level: the level at which this function reject the small reponse
% plot_yn: wheather to plot the results (boolean)

% OUTPUT
% restype: response type of an event
% timestamp: timestamp of the event (sequential, non-sequential, no response)

close all;

% Timetable -> datenum
pttnum = datenum(ptt.Properties.RowTimes);

% 1. load settings 
c = [8,88,158;
43,140,190;
78,179,211;
123,204,196;
168,221,181;
204,235,197;
    ]./255;

% 2. initialization
cumrain = nan(size(stormtt,1),1);
eventduration = nan(size(stormtt,1),1);
restime = nan(size(stormtt,1),size(sm_n,2));
amplitude = nan(size(stormtt,1),size(sm_n,2));
addwin = nan(size(stormtt,1),size(sm_n,2));
restype = strings([size(stormtt,1) 1]);
timestamp = nan(size(stormtt,1),2);
maxtime = nan(size(stormtt,1),1);

% 3. calculate response time
for event = 1:length(stormtt)
    % 3.1. crop the SM and time array
    % obtain the storm features
    sevent = datenum(stormtt(event,1)); % start of the event
    eevent = datenum(stormtt(event,3)); % end of the event
    t_ps = find(pttnum == sevent); % start row of the event in precip time series
    t_pe = find(pttnum == eevent);% end row of the event in precip time series
    cumrain(event,1) = sum(ptt.Var1(t_ps:t_pe));
    eventduration(event,1) = hours(stormtt(event,3)-stormtt(event,1)); % because hourly data
    
    % get the indices of start and the end of event in sm time series, and
    t_sms = find(sm_n_ts == sevent); % start row of the event in SM time series
    t_sme = find(sm_n_ts == eevent); % end row of the event in SM time series
    
    % find max time in precip time series
    maxtime0 = find(ptt.Var1(t_ps:t_pe) == max(ptt.Var1(t_sms:t_sme)));
    maxtime(event,1) = maxtime0(1,1);
 
    % 3.2. calculate the response time
    for k = 1:size(sm_n,2)
        sm = sm_n(:,k);
        
        % if the data is enough, execute the signature
        if sum(isnan(sm(t_sms:t_sme,:))) < size(sm(t_sms:t_sme,:),1)*0.7 
            % event rising time
            amplitude(event,k) = (max(sm(t_sms:t_sme,:)) - sm(t_sms,1));
            restime0 = find(sm(t_sms:t_sme,:) == max(sm(t_sms:t_sme,:)));
            if ~isempty(restime0)
                restime(event,k) = restime0(1,1) - maxtime(event,1);
            end
            
            addwin0 = 0;
            % if the SM is still increasing during the event window (=1d),
            if (restime0(1,1)-1 == t_sme-t_sms+addwin0)      
                % expand the window by 12hr until it gets the maxium value (i.e. reach to the event peak),
                % or additional window size reach to five days 
                % or length of the time series exceeds the soil moisture time series
                while (restime0(1,1)-1 == t_sme-t_sms+addwin0) && (addwin0 < 24*5) && ((t_sme+addwin0+5) < length(sm))
                    addwin0 = addwin0 + 12;
                    maxsm = max(sm(t_sms:t_sme + addwin0,:));
                    restime0 = find(sm(t_sms:t_sme+addwin0,:) == maxsm(1,1));
                end

                % if response time is calculated, get the response time
                if addwin0 < 24*5 || ((t_sme+addwin0) < length(sm))
                    restime(event,k) = restime0(1,1) - maxtime(event,1);
                    amplitude(event,k) = maxsm;
                else
                    restime(event,k) = NaN;
                    amplitude(event,k) = NaN;
                end
                
            end
            addwin(event,k) = addwin0;
            
            % if the event response amplitude is too small, reject the response
            if amplitude(event,k) < increment_level
                restime(event,k) = NaN;
            end
            
        % if data is not enough, do nothing
        end
        clear sm
    end

    % 4. calculate the order of sequence
    resseq = restime(event,:);
    % if there is more than two data, check the sequence
    if sum(~isnan(resseq)) >= 2
        if issorted(resseq(~isnan(resseq)),2,'ascend') 
            restype0 = "sequential";
        else
            restype0 = "nonsequential";
        end
    % if not, categorize as no-response
    else
        restype0 = "noresponse";
    end      
    restype(event,1) = restype0; 
    timestamp(event,1) = sevent;
    timestamp(event,2) = eevent;
    
    % 5. visualization
    if plot_yn
        % plot per event
        if restype0 ~= "noresponse"
            figure(event);
            set(gcf, 'Position', [100 100 300 300])
            subplot(3,1,1)
            plot(pttnum(t_ps:t_pe),ptt(t_ps:t_pe,:).Var1,'color','k');
            subplot(3,1,[2 3])
            % plot sm for each depth
            for k = 1:size(sm_n,2)
                if ~isnan(addwin(event,k))
                    plot(datetime(sm_n_ts(t_sms:t_sme+addwin(event,k)),'ConvertFrom','datenum'),sm_n(t_sms:t_sme+addwin(event,k),k), ...
                        'color',c(k,:), 'LineWidth', 2); hold on
                    ylim([0.1 0.9])
                    % plot the peak timing
                    title(restype(event,1))
                    if ~isnan(restime(event,k)) && maxtime(event,1) + restime(event,k)-1 >= 0
                        maxtimerow = t_sms + maxtime(event,1)+restime(event,k)-1;
                        xline(datetime(sm_n_ts(maxtimerow),'ConvertFrom','datenum'), 'color',c(k,:), 'LineWidth', 2); hold on;
                        plot(datetime(sm_n_ts(maxtimerow),'ConvertFrom','datenum'), sm_n(maxtimerow,k),'o', 'color',c(k,:), 'LineWidth', 2)
                    end
                end
            end
            hold off;
        end
    end
    
end

end