% To detrend using moving average

function [smdtr, smttdtr] = util_detrend(sm,smtt)

y = sm;
T = length(y);

% 
% figure
% plot(y)
% h1 = gca;
% h1.XLim = [0,T];
% title 'Soil moisture depth5cm-station42';
% hold on

% Detrend the data using moving average using 1-year window (8760 hrs)
yS = movmean(y,8760*2,'omitnan');
% Repeating the first and last observations six times to prevent data loss. 
yS(1:6) = yS(7); yS(T-5:T) = yS(T-6);

% Divide the original series by the smoothed series to detrend the data.
% Add the moving average trend estimate to the observed time series plot. 

% multiplicative decomposition
% smdtr = y./yS;

% additive decomposition
smdtr = y - yS;

% add 0.5 to avoid negative SM value (do not trust the absolute value of SM
% anymore
smdtr = smdtr + 0.5;

% return timetable, too
smttdtr = timetable(smtt.Properties.RowTimes, smdtr);

end

% h = plot(yS,'r','LineWidth',2);
% legend(h,'8760-Term Moving Average')
% hold off

% yS ... trend
% xt ... seasonality + noise

% you can apply improved estimate of the trend coponent by applying
% Henderson filter for each month
% % Create seasonal indices
% s = 8760;
% sidx = cell(s,1); % Preallocation
% for i = 1:s
%  sidx{i,1} = i:s:T;
% end
% sidx{1:2}

