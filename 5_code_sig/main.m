%% A main execusion code to run the all soil moisture signatures
%% Originally written by Ryoko Araki (raraki8159@sdsu.edu) at San Diego State University
%% April 17, 2022
% =========== BEGINNING OF THE CODE ============
clear all;
slCharacterEncoding('UTF-8');

%% Select the functions to run
run_disttype = true;
run_fcwp = true;
run_seasontrans = true;
run_eventres = true;
run_RLD = true;
run_restype = true;

% if you want to clear the previous results and save new results
save_results = true;

%% Preparation
% Set your file path
cd("your_file_path");
in_path = "..\4_data_after_qc";
out_path = "..\6_out_signatures";

% Site information
site = ["HB"]; % Hamburg data is uploaded as an example
% Observation tyle (F: Field observation)
obs = "F";

%% Main execution
for i = 1:length(site)
    for j = 1
        
        % create/read new file for recording results
        % timeseries-based signatures
        fn1 = fullfile(out_path, sprintf('disttype_%s_%s.txt',site(i,:), obs(j,:)));
        fn2 = fullfile(out_path, sprintf('fc_%s_%s.txt',site(i,:), obs(j,:)));
        fn3 = fullfile(out_path, sprintf('wp_%s_%s.txt',site(i,:), obs(j,:)));
        
        % seasonal signatures
        fn5 = fullfile(out_path, sprintf('seasontrans_sdate_wet2dry_p_%s_%s.txt',site(i,:), obs(j,:)));
        fn5_2 = fullfile(out_path, sprintf('seasontrans_sdate_wet2dry_l_%s_%s.txt',site(i,:), obs(j,:)));
        fn6 = fullfile(out_path, sprintf('seasontrans_edate_wet2dry_p_%s_%s.txt',site(i,:), obs(j,:)));
        fn6_2 = fullfile(out_path, sprintf('seasontrans_edate_wet2dry_l_%s_%s.txt',site(i,:), obs(j,:)));
        fn7 = fullfile(out_path, sprintf('seasontrans_sdate_dry2wet_p_%s_%s.txt',site(i,:), obs(j,:)));
        fn7_2 = fullfile(out_path, sprintf('seasontrans_sdate_dry2wet_l_%s_%s.txt',site(i,:), obs(j,:)));
        fn8 = fullfile(out_path, sprintf('seasontrans_edate_dry2wet_p_%s_%s.txt',site(i,:), obs(j,:)));
        fn8_2 = fullfile(out_path, sprintf('seasontrans_edate_dry2wet_l_%s_%s.txt',site(i,:), obs(j,:)));
        fn9 = fullfile(out_path, sprintf('seasontrans_duration_wet2dry_p_%s_%s.txt',site(i,:), obs(j,:)));
        fn9_2 = fullfile(out_path, sprintf('seasontrans_duration_wet2dry_l_%s_%s.txt',site(i,:), obs(j,:)));
        fn10 = fullfile(out_path, sprintf('seasontrans_duration_dry2wet_p_%s_%s.txt',site(i,:), obs(j,:)));
        fn10_2 = fullfile(out_path, sprintf('seasontrans_duration_dry2wet_l_%s_%s.txt',site(i,:), obs(j,:)));
        
        % event signatures
        fn11 = fullfile(out_path, sprintf('amplitude_%s_%s.txt',site(i,:), obs(j,:)));
        fn12 = fullfile(out_path, sprintf('risingtime_%s_%s.txt',site(i,:), obs(j,:)));
        fn13 = fullfile(out_path, sprintf('noresrate_%s_%s.txt',site(i,:), obs(j,:)));
        fn14 = fullfile(out_path, sprintf('RLD_%s_%s.txt',site(i,:), obs(j,:)));
        fn15 = fullfile(out_path, sprintf('restype_%s_%s.txt',site(i,:), obs(j,:)));
        
        % delete existing files
        if save_results
            if run_disttype
                delete(fn1);
            end
            if run_fcwp
                delete(fn2); delete(fn3);
            end
            if run_seasontrans
                delete(fn5);delete(fn5_2); delete(fn6); delete(fn6_2); ...
                    delete(fn7); delete(fn7_2); delete(fn8);delete(fn8_2); ...
                    delete(fn9); delete(fn9_2); delete(fn10); delete(fn10_2);
            end
            if run_eventres
                delete(fn11);delete(fn12);delete(fn13);
            end
            if run_RLD
                delete(fn14);
            end
            if run_restype
                delete(fn15);
            end
        end
        
        % Read site information (depth of sensors, number of sm stations, and
        % number of precip stations)
        [depth, ~, nstation, npstation] = io_siteinfo(site(i));
        
        % Read stationflag (something like quality flag for a station)
        filename2 = fullfile(in_path, sprintf('\\%s_%s_cleaned_csv\\stationflag.csv', site(i,:), obs(j,:)));
        fid = fopen(filename2, 'r');
        stationflag0 = textscan(fid,repmat('%f', 1, length(depth)),'HeaderLines',0,'Delimiter',',');
        fclose(fid);
        stationflag = cell2mat(stationflag0);
        
        % Loop for station
        for n = 1:nstation
            statement = sprintf('Currently processing the station #%d of network %s', n, site(i,:));
            disp(statement)
            
            % Storm event separation
            % if there is only one station in the network, read precip data only once
            if npstation == 1 && n ~=1
            % else, read precipitation data corresponding to the soil moisture station
            else
                % read precipitation time table
                pfilename = fullfile(in_path, sprintf('\\%s_%s_cleaned_csv\\ppt_s%02d.csv', site(i,:), obs(j,:), n));
                fid = fopen(pfilename, 'r');
                ptt0 = textscan(fid,'%{dd-MMM-yyyy HH:mm:ss}D %f','HeaderLines',1,'DateLocale','en_US','Delimiter',',');
                fclose(fid);
                ptt = timetable(ptt0{1},ptt0{2});
                clear ptt0
                
                % run the event separation
                stormarray = util_EventSeparation(...
                    datenum(ptt.Properties.RowTimes), ptt.Var1, ...
                    1, 12, 5, 0.2, 1.0, 0.2, 1.0, 5, 0);
                % The output is stormarray = 2-column array with start and end locations of each storm
                % event as indices into the precipitation array
                % 3rd column ... 5-day criterion or until next event
                stormtt = [...
                    ptt.Properties.RowTimes(stormarray(:,1)), ...
                    ptt.Properties.RowTimes(stormarray(:,2)), ...
                    ptt.Properties.RowTimes(stormarray(:,3))];
                
                % save the event separation results
                fnstorm = fullfile(out_path, sprintf('storm_%s_%s', site(i,:), obs(j,:)));
                delete(fnstorm);
                fid = fopen(fnstorm,'a');
                for p = 1:size(stormarray,1)
                    fprintf(fid, '%f %f \n', ...
                        datenum(ptt.Properties.RowTimes(stormarray(p,1))), ...
                        datenum(ptt.Properties.RowTimes(stormarray(p,2))) );
                end
                fclose(fid);
            end
            
            % Loop for sensor depth
            for k = 1:length(depth)
                depthfield = sprintf('depth%dcm', depth(k));
                
                % Read SM data
                filename = fullfile(in_path, sprintf('\\%s_%s_cleaned_csv\\sm_d%02d_s%02d.csv', site(i,:), obs(j,:), depth(k), n));
                
                if exist(filename, 'file') == 2
                    
                    fid = fopen(filename, 'r');
                    smtt0 = textscan(fid,'%q %f','HeaderLines',1,'DateLocale','en_US','Delimiter',',');
                    smtt = timetable(datetime(smtt0{1}),smtt0{2});
                    fclose(fid);
                    sm = smtt0{2};
                    clear smtt0
                    
                    %% Calculate signatures
                    %% if the data is problematic, put NaN in results
                    if stationflag(n,k) == 101 || stationflag(n,k) == 107 || stationflag(n,k) == 106 || stationflag(n,k) == 108
                        
                        disttype = 'NaN';
                        fc = NaN; wp = NaN;
                        seasontrans_sdate_wet2dry_p = NaN;seasontrans_edate_wet2dry_p = NaN;
                        seasontrans_sdate_dry2wet_p = NaN;seasontrans_edate_dry2wet_p = NaN;
                        seasontrans_duration_wet2dry_p = NaN;seasontrans_duration_dry2wet_p = NaN;
                        seasontrans_sdate_wet2dry_l = NaN;seasontrans_edate_wet2dry_l = NaN;
                        seasontrans_sdate_dry2wet_l = NaN;seasontrans_edate_dry2wet_l = NaN;
                        seasontrans_duration_wet2dry_l = NaN;seasontrans_duration_dry2wet_l = NaN;
                        amplitude = NaN; risingtime = NaN; noresrate = NaN;
                        RLD = NaN;
                        restype = 'NaN';
                        
                    %% if the data has no issues, execute the signatures
                    elseif stationflag(n,k) == 100 || stationflag(n,k) == 102 || stationflag(n,k) == 103 ...
                            || stationflag(n,k) == 104 || stationflag(n,k) == 105
                        
                        % Detrend the time series
                        [smdtr, smttdtr] = util_detrend(sm,smtt);
                        
                        %  [SMSignature] Distribution type
                        if run_disttype
                            disttype = sig_pdf(smdtr);
                        end
                        
                        %  [SMSignature] Field capacity and wilting points
                        if run_fcwp || run_seasontrans
                            [fc, wp] = sig_fcwp(sm, smtt, false);
                        end
                        
                        %  [SMSignature] Seasonal transition duration and dates
                        if run_seasontrans
                            [seasontrans_sdate_wet2dry_p, seasontrans_edate_wet2dry_p, ...
                                seasontrans_sdate_dry2wet_p, seasontrans_edate_dry2wet_p, ...
                                seasontrans_duration_wet2dry_p, seasontrans_duration_dry2wet_p, ...
                                seasontrans_sdate_wet2dry_l, seasontrans_edate_wet2dry_l, ...
                                seasontrans_sdate_dry2wet_l, seasontrans_edate_dry2wet_l, ...
                                seasontrans_duration_wet2dry_l, seasontrans_duration_dry2wet_l] ...
                                = sig_seasontrans(smtt, wp, fc, false, "dayofyear");
                        end
                        
                        %  [SMSignature] Event signatures
                        if run_eventres
                            [amplitude, risingtime, noresrate] ...
                                = sig_event(stormtt, ptt, wp, fc, sm, smtt);
                        end
                        
                        %  [SMSignature] Rising Limb Density
                        if run_RLD
                            RLD = ...
                                sig_RLD(sm,smtt.Properties.RowTimes,'plot_results',false,...
                                'eps',0.0001,'rising_limb_length',0.04167, 'minimum_peak',0.01);
                            % rising_limb_length ~= 1h, minimum peak = 1percent of VWC, eps is arbitrary
                            
                        end
                        
                        
                    end
                    
                    % save the results
                    if save_results
                        if run_disttype
                            fid = fopen(fn1,'a');
                            fprintf(fid, '%d %d %s \n', depth(k), n, disttype);
                            fclose(fid);
                        end
                        
                        if run_fcwp
                            fid = fopen(fn2,'a');
                            fprintf(fid, '%d %d %f \n', depth(k), n, fc);
                            fclose(fid);
                            
                            fid = fopen(fn3,'a');
                            fprintf(fid, '%d %d %f \n', depth(k), n, wp);
                            fclose(fid);
                        end
                        
                        if run_seasontrans
                            fid = fopen(fn5,'a');
                            for p = 1:size(seasontrans_sdate_wet2dry_p,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_sdate_wet2dry_p(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn5_2,'a');
                            for p = 1:size(seasontrans_sdate_wet2dry_l,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_sdate_wet2dry_l(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn6,'a');
                            for p = 1:size(seasontrans_edate_wet2dry_p,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_edate_wet2dry_p(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn6_2,'a');
                            for p = 1:size(seasontrans_edate_wet2dry_l,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_edate_wet2dry_l(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn7,'a');
                            for p = 1:size(seasontrans_sdate_dry2wet_p,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_sdate_dry2wet_p(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn7_2,'a');
                            for p = 1:size(seasontrans_sdate_dry2wet_l,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_sdate_dry2wet_l(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn8,'a');
                            for p = 1:size(seasontrans_edate_dry2wet_p,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_edate_dry2wet_p(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn8_2,'a');
                            for p = 1:size(seasontrans_edate_dry2wet_l,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, seasontrans_edate_dry2wet_l(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn9,'a');
                            for p = 1:size(seasontrans_duration_wet2dry_p,1)
                                fprintf(fid, '%d %d %d \n', depth(k), n, seasontrans_duration_wet2dry_p(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn9_2,'a');
                            for p = 1:size(seasontrans_duration_wet2dry_l,1)
                                fprintf(fid, '%d %d %d \n', depth(k), n, seasontrans_duration_wet2dry_l(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn10,'a');
                            for p = 1:size(seasontrans_duration_dry2wet_p,1)
                                fprintf(fid, '%d %d %d \n', depth(k), n, seasontrans_duration_dry2wet_p(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn10_2,'a');
                            for p = 1:size(seasontrans_duration_dry2wet_l,1)
                                fprintf(fid, '%d %d %d \n', depth(k), n, seasontrans_duration_dry2wet_l(p));
                            end
                            fclose(fid);
                        end
                        
                        if run_eventres
                            fid = fopen(fn11,'a');
                            for p = 1:size(amplitude,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, amplitude(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn12,'a');
                            for p = 1:size(risingtime,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, risingtime(p));
                            end
                            fclose(fid);
                            
                            fid = fopen(fn13,'a');
                            for p = 1:size(noresrate,1)
                                fprintf(fid, '%d %d %f \n', depth(k), n, noresrate(p));
                            end
                            fclose(fid);
                        end
                        
                        if run_RLD
                            fid = fopen(fn14,'a');
                            fprintf(fid, '%d %d %f \n', depth(k), n, RLD);
                            fclose(fid);
                        end
                        
                    end
                    
                    % Clear variables for house keeping
                    if run_disttype
                        clear disttype; 
                    end
                    if run_fcwp
                        clear fc wp; 
                    end
                    if run_seasontrans
                        clear seasontrans_sdate_wet2dry_p seasontrans_edate_wet2dry_p seasontrans_sdate_dry2wet_p seasontrans_edate_dry2wet_p seasontrans_duration_wet2dry_p seasontrans_duration_dry2wet_p ...
                            seasontrans_sdate_wet2dry_l seasontrans_edate_wet2dry_l seasontrans_sdate_dry2wet_l seasontrans_edate_dry2wet_l seasontrans_duration_wet2dry_l seasontrans_duration_dry2wet_l; 
                    end
                    if run_eventres
                        clear amplitude risingtime noresrate; 
                    end
                    if run_RLD
                        clear RLD; 
                    end
                    
                else
                    % if the soil moisture data file does not exist, do nothing
                    smtt = timetable(ptt.Properties.RowTimes, NaN(length(ptt.Properties.RowTimes),1));
                end
                
                % create a big table of SM data for ResponseType analysis
                if k == 1
                    smtt_n = smtt;
                else
                    smtt_n = synchronize(smtt_n,smtt);
                end
                
                clear smtt sm smttdtr smdtr
            end 
            
            sm_n0 = timetable2table(smtt_n);
            sm_n = sm_n0{:,2:end};
            sm_n_ts = datenum(smtt_n.Properties.RowTimes);
            
            % run signature that requires multiple depth data
            %  [SMSignature] Response type
            if run_restype
                [restype, timestamp] = sig_restype(sm_n, sm_n_ts, ptt, stormtt, 0.004, 0.02, false);
            end
            
            % Save results 
            if save_results
                if run_restype
                    fid = fopen(fn15,'a');
                    restype(restype == "") = "NaN";
                    for p = 1:size(restype,1)
                        fprintf(fid, '%d %d %s \n', NaN, n, restype(p));
                    end
                    fclose(fid);
                end
            end
            
            % if there is only one precipitation data, use the data for the rest of the loop
            if npstation == 1
                continue
            else
                % otherwise, clear precipitation data
                clear ppt stormarray stormtt
            end
            
            
        end 
        
        clear depth nstation npstation stationflag
        
    end
end

load handel
sound(y,Fs)

%% =========== END OF THE CODE ============
