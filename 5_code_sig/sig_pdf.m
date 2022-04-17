%% Signature code to calculate the modality of PDF of soil moisture timeseries
%% Originally written by Flora Branger, and
%% the algorithm was further developed by Ryoko Araki (raraki8159@sdsu.edu) at San Diego State University
%% April 17, 2022

function [disttype] = sig_pdf(smdtr)

% INPUT
% smdtr: detrended soil moisture timeseries (only soil moisture values) using
% util_detrend

% OUTPUT
% disttype: number of peaks in the soil moisture PDF

% Get the PDF 
[f,xi, bw] = ksdensity(smdtr);

% Remove unsignificant peak (<10% of the largest peak)
[pks,lks] = findpeaks(ksdensity(smdtr, 'BandWidth', bw*2.5),'MinPeakHeight',0.2*max(f));

% Get PDF again just for plotting purpose
[f2,xi2] = ksdensity(smdtr ,'BandWidth', bw*2.5);

% Judge the modality
if size(pks,2) == 1
    disttype = "unimodal";
elseif size(pks,2) == 2
    disttype = "bimodal";
elseif size(pks,2) >= 3
    disttype = "multimodal";
else
    disttype = "NaN";
end

%     figure;
%     x0=0;y0=0;width=250;height=200;set(gcf,'position',[x0,y0,width,height]);
%     plot(xi2, f2, xi2(lks), pks, 'o');
%     hold on
%     plot(xi, f)
%     title(modality);
%     xlabel('VWC(%)');ylabel('Probability');
%
%     figure;
%     x0=0;y0=0;width=250;height=200;set(gcf,'position',[x0,y0,width,height]);
%     plot(smdtr);
%     xlabel('Year');ylabel('VWC(%)');

end