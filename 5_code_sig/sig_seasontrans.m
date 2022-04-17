%% Signature code to obtain seasonal transition metrics
%% Originally written by Flora Branger, and
%% the algorithm was further developed by Ryoko Araki (raraki8159@sdsu.edu) at San Diego State University

% 'sig_seasontrans' estimates seasonal transition duration and date(start/end)
% both for dry-to-wet and wet-to-dry transition
% based on linear piecewise model and logistic models

% applicable only for soil moisture time series with bimodal distribution PDF (i.e. clear seasonality)
% this code does not include the judgement of bimodal distribution

% INPUT
% smtt: time series of soil moisture in timetable format. can only have one column of data
% e.g.    |         Time        |   Var1  |
%         |2006-01-05 00:00:00: |   0.25  |
%         |       ...           |   ...   |
% wp: wilting point (dry equilibrium value) in double format
% fc: field capacity (wet equilibrium value) in double format
% format_date: format date in datenum or day of the year

% OUTPUT
% seasonal trans sdate & edate: start and the end of the day of seasonal transition
% seasonal trans duration: number of days it took to complete the seasonal transition
% dry2wet: dry to wet seasoon transition
% wet2dry: wet to dry season transition
% l: logit model
% p: piecewise linear model

% requires other function files
% util_piecewise_constraint.m
% util_piecewise_constraint2.m

function [seasontrans_sdate_wet2dry_p, seasontrans_edate_wet2dry_p, seasontrans_sdate_dry2wet_p, seasontrans_edate_dry2wet_p, ...
    seasontrans_duration_wet2dry_p, seasontrans_duration_dry2wet_p, seasontrans_sdate_wet2dry_l, seasontrans_edate_wet2dry_l, ...
    seasontrans_sdate_dry2wet_l, seasontrans_edate_dry2wet_l, seasontrans_duration_wet2dry_l, seasontrans_duration_dry2wet_l] = sig_seasontrans(smtt, wp, fc, plot_results, format_date)

%% Preparation of time series
% aggregate into daily to remove diel signals
smtt = retime(smtt,'daily','mean');

% fit the sine curve and find the appropriate seasonal crop window
% initial models are set after manual testing
% [shift in y-dir, amplitude, shift in x-dir]
S0 = [0.2 0.1 150];   % initial parameter
lb = [0 0.01 0]; % upper boundary
ub = [1.5 10 365];  % lower boundary
SineModel = @(S,x) S(1) * ( 1 + S(2)*sin(2*pi*x/365 + S(3)));
x = datenum(smtt.Properties.RowTimes); %[1:size(smtt,1)]';
y = table2array(smtt);
A = find(~isnan(y)); % to skip date with NaN data
options = optimoptions('lsqcurvefit','Display','off'); % turn off the optimization results
Sfit = lsqcurvefit(SineModel, S0, x(A), y(A),lb,ub, options);

% plot the fitted sine curve
if plot_results
    sinemodelpred = SineModel(Sfit,sort(x(A)));
    figure;plot(x(A),y(A),'o',sort(x(A)),sinemodelpred,'r-');
end

% obtain the start and end of the sine fitting
x2 = x(A);
sine_start = fix(x2(1)/365 + Sfit(3)/2/pi);
sine_end = fix(x2(end)/365 + Sfit(3)/2/pi)-1;

%% Main execution. Loop for water years
%% initialization

record_seasontrans_date = datetime(NaT);
% store following variables
%(1) dry2wet_start_piecewise; (2) dry2wet_end_piecewise; (3) dry2wet_start_logit; (4) dry2wet_end_logit;
%(5) wet2dry_start_piecewise; (6) wet2dry_end_piecewise; (7) wet2dry_start_logit; (8) wet2dry_end_logit;
record_seasontrans_duration = days([]);
% store following variables
%(1) dry2wet_duration_piecewise; (2) dry2wet_duration_logit;
%(3) wet2dry_duration_piecewise; (4) wet2dry_duration_logit

Pfit = zeros(100,6); % store parameter sets for linear piecewise eq.
Lfit = zeros(100,4);
Pfit2 = zeros(100,6); % store parameter sets for linear piecewise eq.
Lfit2 = zeros(100,4);
i = 0;

for sinecycle =  sine_start:1:sine_end %cut out at every two peaks
    i = i+1; % display #seasonal cycle
    
    % specify the dryest & wettest points from the sine curve
    dry2wet_start = datetime( 365/2/pi*(2*sinecycle*pi -pi/2 - Sfit(3)), 'ConvertFrom','datenum');
    dry2wet_end = datetime( 365/2/pi*(2*sinecycle*pi +pi/2 -Sfit(3)), 'ConvertFrom','datenum');
    wet2dry_start = datetime( 365/2/pi*(2*sinecycle*pi +pi/2 -Sfit(3)), 'ConvertFrom','datenum');
    wet2dry_end = datetime( 365/2/pi*(2*sinecycle*pi +3*pi/2 -Sfit(3)), 'ConvertFrom','datenum');
    
    %% Dry2wet procedure
    
    % ==crop out the dry2wet time series== %
    % crop a season
    TR = timerange(dry2wet_start-days(30),dry2wet_end+days(30));
    seasonsm = smtt(TR,:);
    seasonsmvalue = table2array(seasonsm);
    seasonsmtime = seasonsm.Properties.RowTimes;
    
    % check the season has how much NaN sequence. If there is too much
    % NaN, skip the season
    if sum(isnan(seasonsmvalue))/size(seasonsmvalue,1) > 0.2 || isempty(seasonsm)
        % return NaT
        record_seasontrans_date(i,1:8) = NaT;
        record_seasontrans_duration(1,1:4) = NaN;
    else
        
        % find the actual dryest & wettest point during a season using max and min
        % dryest point
        dryestpoint = find(seasonsmvalue == min(seasonsmvalue));
        dry2wet_start2 = seasonsmtime(dryestpoint(1));
        
        % dryest period is short, do not use wettest point for cutting
        % out the TS
        %             % wettest point
        %             % should be happening later than dryest point
        %             dry2wet_end2 = seasonsmtime(seasonsmvalue == max(seasonsmvalue(dryestpoint:end)));
        %             dry2wet_end2 = dry2wet_end2(1);
        
        % specify seasonal time window
        % 1-month buffer from the dryest/wettest point
        TR = timerange(dry2wet_start2-days(15),dry2wet_end+days(30));
        seasonsm = smtt(TR,:);
        
        % if the VWC is smaller than the minimum VWC after the wettest point, remove the data
        if length(seasonsm.Var1) >= 30
            seasonsm.Var1(seasonsm.Var1(end-30:end)< min(seasonsmvalue)) = NaN;
        end
        
        % if the time series is too short, skip the time series
        if length(seasonsm.Var1) < 100
            % return NaT
            record_seasontrans_date(i,1:8) = NaT;
            record_seasontrans_duration(1,1:4) = NaN;
        else
            
            seasonsmvalue = table2array(seasonsm);
            % seasonsmtime = seasonsm.Properties.RowTimes;
            
            % ==Piecewise linear regression == %
            % define the model
            y2 = fillmissing(seasonsmvalue,'linear'); %soil moisture data
            x2 = [1:size(y2,1)]'; % #data (daily) from the cropping point
            
            P0 = [0.1       0.001      50      10     wp   fc]; %[P1 P2 P3 P4 wilting_point field_capacity]
            Plb = [0        0        0       1    wp-0.5   fc-0.5];
            Pub = [1        1        150     250    wp+0.5   fc+0.5];
            
            
            I = ones(size(y2,1),1);
            
            % Minimize the least square of piecewise model
            plusfun = @(x2) max(x2,0);
            piecewisemodel2 =  @(P) sum((y2 - (P(1) + I.*P(2).*x2 + I.*P(2).*plusfun(I.*P(3) - x2) + (-I.*P(2)).*plusfun(x2-(I.*P(3)+I.*P(4))))).^2);
            nonlcon = @util_piecewise_constraint;
            A = []; % No other constraints
            b = [];
            Aeq = [];
            beq = [];
            options = optimoptions('fmincon','Display','off'); % turn off the optimization results
            try
                Pfit(i,:) = fmincon(piecewisemodel2,P0,A,b,Aeq,beq,Plb,Pub,nonlcon,options);
            catch
                Pfit(i,:) = NaN(1,6);
            end
            
            % if the parameters reached to upper boundary or lower boundary,
            % reject it because the fitting has likely to have failed
            % case 1: the slope maxed out
            % case 2: the x shift reached to upper boundary
            % case 3: transition duration reached to lower boundary
            if Pfit(i,1) > Pub(1) - 0.01 || ...
                    (Pfit(i,3) > Pub(3)-0.01) || ...
                    (Pfit(i,4) < Plb(4)+0.01*(Pub(4)-Plb(4))) || ...
                    (Pfit(i,4) > Pub(4)-0.01*(Pub(4)-Plb(4))) || ...
                    (abs(Pfit(i,6) - Pfit(i,5)) < 1.0E-04)
                Pfit(i,:) = NaN(1,6);
            end
            
            piecewisemode1 = @(P,x2) P(1) + P(2)*x2 + P(2)*plusfun(P(3)-x2) + (-P(2))*plusfun(x2-(P(3)+P(4))) + 0*P(5) + 0*P(6);
            modelpred = piecewisemode1(Pfit(i,:),[1:250]');
            
            % get the seasonal transition date
            record_seasontrans_date(i,1) =  dry2wet_start2-days(30)+ days(Pfit(i,3)); % dry2wet_startdate
            record_seasontrans_date(i,2) =  dry2wet_start2-days(30)+ days(Pfit(i,3)+Pfit(i,4)); % dry2wet_enddate
            % get the seasonal transition duration
            record_seasontrans_duration(i,1) = days(Pfit(i,4)); % dry2wet_duration
            % record all / take an average?
            
            % plot
            %     figure;plot(x2,y2,'o',sort(x2),modelpred,'r-');
            %     xline(Pfit(i,3),'LineWidth',2); xline(Pfit(i,3)+Pfit(i,4),'LineWidth',2);
            
            % == Logit model == %
            L0 =  [wp  fc-wp     1E-1    50]; %[y-shift amplitude k_steepness x-shift]
            Llb = [0  0.01       0       0];
            Lub = [1  1.5     20       120];
            
            % Not constrained
            logit = @(L,x2) L(1) + L(2)./( 1 + exp(1).^(-L(3).*(x2-L(4))) );
            options = optimoptions('lsqcurvefit','Display','off'); % turn off the optimization results
            try
                Lfit(i,:) = lsqcurvefit(logit, L0 , x2, y2, Llb, Lub, options);
            catch
                Lfit(i,:) = NaN(1,4);
            end
            
            % if the parameters reached to upper boundary or lower boundary,
            % reject it because the fitting has likely to have failed
            % case 1: steepness got to be zero
            if (Lfit(i,3) == 0) ||  ...
                    ( Lfit2(i,4)> Lub(4)+0.01*(Lub(4)-Llb(4))) || ...
                    ( Lfit2(i,4)> Lub(4)-0.01*(Lub(4)-Llb(4)))
                Lfit(i,:) = NaN(1,4);
            end
            
            x250d = [1:250]';
            logitpred = logit(Lfit(i,:),x250d);
            %     figure;plot(x2,y2,'o',sort(x2),logitpred,'r-')
            
            % ==Get signature== %
            % get the seasonal transition date
            percent005 = find(logitpred >= wp+0.05*(fc-wp));
            if isempty(percent005) == 0
                record_seasontrans_date(i,3) = dry2wet_start2-days(30)+ days(x250d(percent005(1)));
            else
                record_seasontrans_date(i,3) = NaT;
            end
            percent095 = find(logitpred >= wp+0.95*(fc-wp));
            
            % get the seasonal transition duration
            if isempty(percent095) == 0
                record_seasontrans_date(i,4) = dry2wet_start2-days(30)+ days(x250d(percent095(1)));
                record_seasontrans_duration(i,2) = days(percent095(1) - percent005(1));
            else
                percent005 = [];
                record_seasontrans_date(i,3) = NaT;
                record_seasontrans_date(i,4) = NaT;
                record_seasontrans_duration(i,2) = NaN;
            end
            % record all / take an average?
            
            % plot of both functions
            if plot_results
                figure; plot(x2,y2,'o','DisplayName','obs');
                hold on; plot(x250d,logitpred,'r-','DisplayName','logit','LineWidth',2);
                hold on; plot(x250d,modelpred,'g-','DisplayName','piecewise linear','LineWidth',2);
                if ~isnan(Pfit(i,3))
                    hold on; xline(Pfit(i,3),'g','LineWidth',1.5); xline(Pfit(i,3)+Pfit(i,4),'g','LineWidth',1.5);
                end
                hold on;
                if isempty(percent005) == 0
                    xline(percent005(1),'r','LineWidth',1.5)
                end
                if isempty(percent095) == 0
                    xline(percent095(1),'r','LineWidth',1.5)
                end
                legend('Location','northwest'); xlabel('days from the dryest day of the year'); ylabel('VWC [m^3/m^3]');
                %             titlename = sprintf('%s-depth%cm-station%d', site, depth(k), n);
                %             title(titlename);
            end
        end
        
    end
    
    %% Wet2dry procedure
    % ==crop out the wet2dry time series== %
    % crop a season
    TR2 = timerange(wet2dry_start-days(30),wet2dry_end+days(30));
    seasonsm2 = smtt(TR2,:);
    
    % check the season has how much NaN sequence. If there is too much
    % NaN, skip the season
    if sum(isnan(seasonsm2.Var1))/size(seasonsm2.Var1,1) > 0.2 || isempty(seasonsm2)
        % return NaN
        record_seasontrans_date(i,1:8) = NaT;
        record_seasontrans_duration(1,1:4) = NaN;
    else
        
        seasonsmvalue2 = table2array(seasonsm2);
        seasonsmtime2 = seasonsm2.Properties.RowTimes;
        
        % find the actual dryest & wettest point during a season using max and min
        % wettest point
        wettestpoint = find(seasonsmvalue2 == max(seasonsmvalue2));
        wet2dry_start2 = seasonsmtime2(wettestpoint(1));
        % driest point
        % should be happening later than wettest point
        wet2dry_end2 = seasonsmtime2(seasonsmvalue2 == min(seasonsmvalue2(wettestpoint:end)));
        wet2dry_end2 = wet2dry_end2(1);
        
        % specify seasonal time window
        % Buffer from the dryest/wettest point
        TR2 = timerange(wet2dry_start2-days(30),wet2dry_end2(end)+days(30));
        seasonsm2 = smtt(TR2,:);
        
        % if the VWC is larger than the max VWC after the dryest point, remove the data
        % find(seasonsm2.Var1(end-30:end)> max(seasonsmvalue2));
        if length(seasonsm2.Var1) > 30
            seasonsm2.Var1(seasonsm2.Var1(end-30:end)> max(seasonsmvalue2)) = NaN;
        end
        
        % if the time series is too short, skip the time series
        if length(seasonsm2.Var1) < 100
            % return NaN
            record_seasontrans_date(i,1:8) = NaT;
            record_seasontrans_duration(1,1:4) = NaN;
        else
            
            
            seasonsmvalue2 = table2array(seasonsm2);
            
            % seasonsmtime2 = seasonsm2.Properties.RowTimes;
            
            % ==Piecewise linear regression == %
            % define the model
            y3 = fillmissing(seasonsmvalue2,'linear'); %soil moisture data
            x3 = [1:size(y3,1)]'; % #data (daily) from the cropping point
            
            P02 =   [0.5     -0.0001          50      10     wp   fc]; %[P1 P2 P3 P4 wilting_point field_capacity]
            Plb2 =  [0.1     -1           0       1    wp-0.5   fc-0.5];
            Pub2 =  [2.0    0            150      250    wp+0.5   fc+0.5];
            
            I2 = ones(size(y3,1),1);
            
            % Minimize the least square of piecewise model
            plusfun = @(x3) max(x3,0);
            piecewisemodel4 =  @(P2) sum((y3 - (P2(1) + I2.*P2(2).*x3 + I2.*P2(2).*plusfun(I2.*P2(3) - x3) + (-I2.*P2(2)).*plusfun(x3-(I2.*P2(3)+I2.*P2(4))))).^2);
            nonlcon2 = @util_piecewise_constraint2;
            A = []; % No other constraints
            b = [];
            Aeq = [];
            beq = [];
            options = optimoptions('fmincon','Display','off'); % turn off the optimization results
            try
                Pfit2(i,:) = fmincon(piecewisemodel4,P02,A,b,Aeq,beq,Plb2,Pub2,nonlcon2,options);
            catch
                Pfit2(i,:) = NaN(1,6);
            end
            
            % if the parameters reached to upper boundary or lower boundary,
            % reject it because the fitting has likely to have failed\
            % case 1: the slope got flat
            % case 2: the x shift reached to upper boundary
            % case 3: transition duration reached to lower boundary
            % case 4: ... to upper boundary
            if Pfit2(i,1) < Plb2(1) + 0.01 ...
                    || (Pfit2(i,3) > Pub2(3)-0.01) ...
                    || (Pfit2(i,4) < Plb2(4)+0.01*(Pub2(4)-Plb2(4))) ...
                    || (Pfit2(i,4) > Pub2(4)-0.01*(Pub2(4)-Plb2(4))) ...
                    || (abs(Pfit2(i,6) - Pfit2(i,5)) < 1.0E-04)
                Pfit2(i,:) = NaN(1,6);
            end
            
            % get the seasonal transition date
            record_seasontrans_date(i,5) =  wet2dry_start2-days(30)+ days(Pfit2(i,3)); % dry2wet_startdate
            record_seasontrans_date(i,6) =  wet2dry_start2-days(30)+ days(Pfit2(i,3)+Pfit2(i,4)); % dry2wet_enddate
            % get the seasonal transition duration
            record_seasontrans_duration(i,3) = days(Pfit2(i,4)); % dry2wet_duration
            % record all / take an average?
            
            % plot
            piecewisemodel3 = @(P2,x3) P2(1) + P2(2)*x3 + P2(2)*plusfun(P2(3)-x3) + (-P2(2))*plusfun(x3-(P2(3)+P2(4)));
            modelpred2 = piecewisemodel3(Pfit2(i,:),[1:250]');
            %     figure;plot(x3,y3,'o',sort(x3),modelpred2,'r-');
            %     xline(Pfit2(i,3),'LineWidth',2); xline(Pfit2(i,3)+Pfit2(i,4),'LineWidth',2);
            
            % == Logit model == %
            L02 =  [wp  fc-wp     -0.1  50]; %[y-shift amplitude k_steepness x-shift wilting_point field_capacity]
            Llb2 = [0  0.01    -20      0];
            Lub2 = [1.5  1.5     0       200];
            
            % Optimization, nonconstrained (constrained in parameter sets)
            logit3 = @(L2,x3) L2(1) + L2(2)./( 1 + exp(1).^(-L2(3).*(x3-L2(4))) );
            options = optimoptions('lsqcurvefit','Display','off'); % turn off the optimization results
            try
                Lfit2(i,:) = lsqcurvefit(logit3, L02, x3, y3, Llb2, Lub2, options);
            catch
                Lfit2(i,:) = NaN(1,4);
            end
            
            % if the parameters reached to upper boundary or lower boundary,
            % reject it because the fitting has likely to have failed
            % case 1: steepness got to be almost zero
            % case 2: x-shift maxed out
            % case 3: x-shift reached min-out
            if (Lfit2(i,3) > -1.0E-05) || ...
                    ( Lfit2(i,4)> Lub2(4)+0.01*(Lub2(4)-Llb2(4))) || ...
                    ( Lfit2(i,4)> Lub2(4)-0.01*(Lub2(4)-Llb2(4))) || ...
                    ( Lfit2(i,4)< Llb2(4)+1.0E-05)
                Lfit2(i,:) = NaN(1,4);
            end
            
            x250d = [1:250]';
            logitpred2 = logit3(Lfit2(i,:),x250d);
            % logitpred2 = logit3(Lfit2(i,:),sort(x3));
            %     figure;plot(x3,y3,'o',sort(x3),logitpred2,'r-')
            
            % ==Get signature== %
            % get the seasonal transition date
            % 5% percentile
            percent0052 = find(logitpred2 <= wp+0.95*(fc-wp));
            if isempty(percent0052) == 0
                record_seasontrans_date(i,7) = wet2dry_start2-days(30)+ days(x250d(percent0052(1)));
            else
                record_seasontrans_date(i,7) = NaT;
            end
            percent0952 = find(logitpred2 <= wp+0.05*(fc-wp));
            
            % get the seasonal transition duration
            if isempty(percent0952) == 0
                record_seasontrans_date(i,8) = wet2dry_start2-days(30)+ days(x250d(percent0952(1)));
                record_seasontrans_duration(i,4) = days(percent0952(1) - percent0052(1));
            else
                percent0052 = [];
                record_seasontrans_date(i,7) = NaT;
                record_seasontrans_date(i,8) = NaT;
                record_seasontrans_duration(i,4) = NaN;
            end
            
            
            % plot the results from both models
            if plot_results
                figure;plot(x3,y3,'o','DisplayName','obs');
                hold on; plot(x250d,logitpred2,'r-','DisplayName','logit','LineWidth',2);
                hold on; plot(x250d,modelpred2,'g-','DisplayName','piecewise linear','LineWidth',2);
                if ~isnan(Pfit2(i,3))
                    hold on; xline(Pfit2(i,3),'g','LineWidth',1.5); xline(Pfit2(i,3)+Pfit2(i,4),'g','LineWidth',1.5);
                end
                hold on;
                if isempty(percent0052) == 0
                    xline(percent0052(1),'r','LineWidth',1.5); hold on;
                end
                if isempty(percent0952) == 0
                    xline(percent0952(1),'r','LineWidth',1.5)
                end
                legend('Location','northeast'); xlabel('days from the wettest day of the year'); ylabel('VWC [m^3/m^3]');
                %             titlename = sprintf('%s-depth%cm-station%d', site, depth(k), n);
                %             title(titlename);
            end
            
        end
        
    end
end

%% return signatures
% 0days --> NaN
record_seasontrans_duration(record_seasontrans_duration == 0) = days(NaN);

if size(record_seasontrans_date, 2) ~= 1 && size(record_seasontrans_date, 2) ~= 0
    % return the array
    if format_date == "dayofyear"
        seasontrans_sdate_dry2wet_p = day(record_seasontrans_date(:,1),'dayofyear');
        seasontrans_edate_dry2wet_p = day(record_seasontrans_date(:,2),'dayofyear');
        
        seasontrans_sdate_dry2wet_l = day(record_seasontrans_date(:,3),'dayofyear');
        seasontrans_edate_dry2wet_l = day(record_seasontrans_date(:,4),'dayofyear');
        
        seasontrans_sdate_wet2dry_p = day(record_seasontrans_date(:,5),'dayofyear');
        seasontrans_edate_wet2dry_p = day(record_seasontrans_date(:,6),'dayofyear');
        
        seasontrans_sdate_wet2dry_l = day(record_seasontrans_date(:,7),'dayofyear');
        seasontrans_edate_wet2dry_l = day(record_seasontrans_date(:,8),'dayofyear');
    elseif format_date == "date"
        seasontrans_sdate_dry2wet_p = record_seasontrans_date(:,1);
        seasontrans_edate_dry2wet_p = record_seasontrans_date(:,2);
        
        seasontrans_sdate_dry2wet_l = record_seasontrans_date(:,3);
        seasontrans_edate_dry2wet_l = record_seasontrans_date(:,4);
        
        seasontrans_sdate_wet2dry_p = record_seasontrans_date(:,5);
        seasontrans_edate_wet2dry_p = record_seasontrans_date(:,6);
        
        seasontrans_sdate_wet2dry_l = record_seasontrans_date(:,7);
        seasontrans_edate_wet2dry_l = record_seasontrans_date(:,8);
    end
else
    seasontrans_sdate_dry2wet_p = NaN;
    seasontrans_edate_dry2wet_p = NaN;
    seasontrans_sdate_dry2wet_l = NaN;
    seasontrans_edate_dry2wet_l = NaN;
    seasontrans_sdate_wet2dry_p = NaN;
    seasontrans_edate_wet2dry_p = NaN;
    seasontrans_sdate_wet2dry_l = NaN;
    seasontrans_edate_wet2dry_l = NaN;
end

if size(record_seasontrans_duration, 2) ~= 1 && size(record_seasontrans_duration, 2) ~= 0
    seasontrans_duration_dry2wet_p = days(record_seasontrans_duration(:,1));
    seasontrans_duration_dry2wet_l = days(record_seasontrans_duration(:,2));
    
    seasontrans_duration_wet2dry_p = days(record_seasontrans_duration(:,3));
    seasontrans_duration_wet2dry_l = days(record_seasontrans_duration(:,4));
else
    seasontrans_duration_dry2wet_p = NaN;
    seasontrans_duration_dry2wet_l = NaN;
    seasontrans_duration_wet2dry_p = NaN;
    seasontrans_duration_wet2dry_l = NaN;
end

% save to check the parameters
writematrix(Pfit(~isnan(Pfit)), './sig_seasontrans_para/season_Pfit.txt', 'WriteMode','append');
writematrix(Pfit2(~isnan(Pfit2)), './sig_seasontrans_para/season_Pfit2.txt','WriteMode','append');
writematrix(Lfit(~isnan(Lfit)), './sig_seasontrans_para/season_Lfit.txt','WriteMode','append');
writematrix(Lfit2(~isnan(Lfit2)), './sig_seasontrans_para/season_Lfit2.txt','WriteMode','append');

end

