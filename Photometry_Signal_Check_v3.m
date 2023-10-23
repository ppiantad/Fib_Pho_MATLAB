
% use this script to quickly check if mice have signal. should see
% relatively small fluctuations in 405 channel, and larger ones in 465

close all; clear all; clc
%% Make some nice colors
tic
green = [0.4660 0.6740 0.1880];
red = [0.6350 0.0780 0.1840];
%% Import Data
root = upper('k:');

BLOCKPATH = 'H:\MATLAB\TDTbin2mat\Photometry\D1-eOP-5\D1-eOP-5-210216-114052'; 

% BLOCKPATH = strcat(root,'K:\MATLAB\TDTbin2mat\Photometry\RRD340\RRD340-230328-144152'); 

data = TDTbin2mat(BLOCKPATH);


%work behavioral data into photo data, generate descriptive statistics
[BehavData,ABETfile,Descriptives, block_end, largeRewSide, smallRewSide]=ABET2TableFn_Chamber_A_v6('D1-eOP-5 02162021.csv',[]);

%%
Sensor = input('Please enter name of sensor used (C or D): \n','s'); %input which sensor the mouse was recorded from 
Sensor = upper(Sensor);
if Sensor =='A'
    Sensor = 'A';
elseif Sensor =='B'
    Sensor = 'B';
elseif Sensor == 'C'
    Sensor = 'C';
elseif Sensor == 'D'
    Sensor = 'D';
% elseif Sensor == 'A_2'
%     Sensor = 'A_2';
% elseif Sensor == 'A_3'
%     Sensor = 'A_3';
% elseif Sensor == 'RCAMP'
%     Sensor = 'RCAMP';
elseif Sensor == 'RG'
    Sensor = 'RG';
end

%find time of first TTL.  This value will be passed to Behavioral Data so
%that the timestamps of the ABET data will match the timestamps of Synapse

% if Sensor == 'C';
% timeStart=data.epocs.TTL_.onset(1);
% elseif Sensor == 'D';
%         timeStart=data.epocs.TL2_.onset(1);
% end

streamNames=fieldnames(data.streams);

if Sensor == 'A';
    STREAM_STORE1 = 'x405S'
    STREAM_STORE2 = 'x470S'
elseif Sensor =='B';
    STREAM_STORE1 = 'x405G'
    STREAM_STORE2 = 'x470G'
elseif Sensor == 'D';
    STREAM_STORE1 = streamNames{2};%'d5sP' 
    STREAM_STORE2 = streamNames{1}; %'d0sP'
elseif Sensor == 'C';
    STREAM_STORE1 = streamNames{3}; % 'x45sP'
    STREAM_STORE2 = streamNames{4}; % 'x40sP'
% elseif Sensor == 'A_2'
%     STREAM_STORE1 = 'x405A'
%     STREAM_STORE2 = 'x470A'
elseif Sensor == 'G';
    STREAM_STORE1 = 'x05A'
    STREAM_STORE2 = 'x70A'
% elseif Sensor == 'RCAMP';
%     STREAM_STORE1 = streamNames{3}; % 'x405C'
%     STREAM_STORE2 = streamNames{4}; % 'x565G'
elseif Sensor == 'RG';
    STREAM_STORE1 = streamNames{2}; % 'x405C' %405
    STREAM_STORE2 = streamNames{1}; % 'x565G' %470
    STREAM_STORE3 = streamNames{3}; % %565
end


if Sensor == 'A';
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'B';
    timeStart=data.epocs.TTL_.onset(1);  
elseif Sensor == 'C';
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'D' && contains(BLOCKPATH, 'D1-eOP-4') && contains(BLOCKPATH, '210215');
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'D' && contains(BLOCKPATH, 'D1-eOP-5') && contains(BLOCKPATH, '210210') || contains(BLOCKPATH, '210216');
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'D';
    timeStart=data.epocs.TL2_.onset(1);
elseif Sensor == 'RG';
%     timeStart=data.epocs.TTL_.onset(1);
    timeStart=data.epocs.TL2_.onset(1);
end



%align timestamps of behavioral data to timestamps of photometry by adding
%the time elapsed from when Synapse recording began to when the ABET
%program was started back to each relevant time column (choiceTime,
%collectionTime, stTime)
% timeShift=timeStart*ones(numel(BehavData.choiceTime(:)),1);
% BehavData.choiceTime(:)=BehavData.choiceTime(:)+timeShift;
% BehavData.collectionTime(:)=BehavData.collectionTime(:)+timeShift;
% BehavData.stTime(:)=BehavData.stTime(:)+timeShift;
% block_end = block_end +timeShift;

%% Make Time and DC Vectors


time = (1:length(data.streams.(STREAM_STORE2).data))/data.streams.(STREAM_STORE2).fs;


idx_start = time >= timeStart;
trim_470 = data.streams.(STREAM_STORE2).data(idx_start);
data.streams.(STREAM_STORE2).data = [];
data.streams.(STREAM_STORE2).data = trim_470;
trim_405 = data.streams.(STREAM_STORE1).data(idx_start);
data.streams.(STREAM_STORE1).data = [];
data.streams.(STREAM_STORE1).data = trim_405;

time2 = time(idx_start)-timeStart;
% find last valid timestamp (sometimes mice omit final bunch of trials
% leaving no valid collectionTime
if BehavData.collectionTime(end,:) ~= 0
    idx_end = time2 <= (BehavData.collectionTime(end,:)+10);
elseif BehavData.collectionTime(end,:) == 0
    idx_end = time2 <= (BehavData.choiceTime(end,:)+10);
end


trim_470 = data.streams.(STREAM_STORE2).data(idx_end);
data.streams.(STREAM_STORE2).data = [];
data.streams.(STREAM_STORE2).data = trim_470;
trim_405 = data.streams.(STREAM_STORE1).data(idx_end);
data.streams.(STREAM_STORE1).data = [];
data.streams.(STREAM_STORE1).data = trim_405;
% test3 = test1(idx_end);
time3 = time2(idx_end);





% test1 = data.streams.(STREAM_STORE2).data(idx_start & idx_end);
% time1 = time(idx_start1 & idx_end);
% time2 = time1-data.epocs.TTL_.onset; 
% time3 = (1:length(test1))/data.streams.(STREAM_STORE2).fs;
% 
% indx_test = time >= data.epocs.TTL_.onset & <= BehavData.collectionTime(end,:);
% time2 = time(idx_start)-data.epocs.TTL_.onset;  %time2 = time(idx_start)-data.epocs.TTL_.onset;
% block1_time(1,1:sum(idx_start)) = data.streams.(STREAM_STORE2).data(idx_start);
%%
GCamp_DC = data.streams.(STREAM_STORE2).data - mean(data.streams.(STREAM_STORE2).data);
UV_DC = data.streams.(STREAM_STORE1).data - mean(data.streams.(STREAM_STORE1).data); 

if Sensor == 'RG'
    RCaMP_DC = data.streams.(STREAM_STORE3).data - mean(data.streams.(STREAM_STORE3).data); 
elseif Sensor ~= 'RG'
    RCaMP_DC = UV_DC;
end

minLength_array = floor(min([length(GCamp_DC), length(UV_DC), length(RCaMP_DC)])/10);
GCamp_DC = GCamp_DC(1:min([length(GCamp_DC), length(UV_DC), length(RCaMP_DC)]));
UV_DC = UV_DC(1:min([length(GCamp_DC), length(UV_DC), length(RCaMP_DC)]));
RCaMP_DC = RCaMP_DC(1:min([length(GCamp_DC), length(UV_DC), length(RCaMP_DC)]));
time = time(1:min([length(GCamp_DC), length(UV_DC), length(RCaMP_DC)]));


% downsample 10x and average 405 signal
allSignals = data.streams.(STREAM_STORE1).data;
N = 10; %10
F405 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength1 = size(F405,2);

F405 = F405(1:minLength_array);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405);
stdSignal1 = std(double(F405))/sqrt(size(F405,1));
dcSignal1 = mean(meanSignal1);


% downsample 10x and average 465 signal
allSignals = data.streams.(STREAM_STORE2).data;
F465 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength2 = size(F465,2);

F465 = F465(1:minLength_array);

% Create mean signal, standard error of signal, and DC offset of 465 signal
meanSignal2 = mean(F465);
stdSignal2 = std(double(F465))/sqrt(size(F465,1));
dcSignal2 = mean(meanSignal2);

% downsample 10x and average 565 signal



% downsample time 10x for the dF analysis
timedownsample = zeros(size(time(:,1:N:end-N+1)));
for ii = 1:size(time,1)
    timedownsample(ii,:) = arrayfun(@(i) mean(time(ii,i:i+N-1)),1:N:length(time)-N+1);
end
minLength1 = size(timedownsample,2);



% fit 405 to 465 to detrend signal (Y_dF_all)
bls = polyfit(F405(1:end), F465(1:end), 1); %polyfit(F465(1:end), F405(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F465 - Y_fit_all;
Y_dF_all_delta_F_over_F = Y_dF_all./Y_fit_all;

%testing a moving average with a window size of 100 (roughly 1 s at
%downsampled resolution). Definitely smooths things a bit. 
movmean_Y_dF_all = movmean(Y_dF_all(1,:),100);



if Sensor == 'RG';
    allSignals = data.streams.(STREAM_STORE3).data;
    F565 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F565(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength3 = size(F565,2);
    
    F565 = F565(1:minLength_array);


% Create mean signal, standard error of signal, and DC offset of 465 signal
meanSignal3 = mean(F565);
stdSignal3 = std(double(F565))/sqrt(size(F565,1));
dcSignal3 = mean(meanSignal3);

%fit 405 to 565 to detrend signal (Y_dF_RCAMP)
bls_rcamp = polyfit(F405(1:end), F565(1:end), 1);
Y_fit_all_rcamp = bls_rcamp(1) .* F405 + bls_rcamp(2);
Y_dF_all_rcamp = F565 - Y_fit_all_rcamp;

end




%% Plot Data
figure; 
plot(time,GCamp_DC,'Color',green,'LineWidth',2); hold on; 
plot(time, UV_DC,'b','LineWidth',2); hold on;
plot(time, RCaMP_DC,'r','LineWidth',2);

title('DC Plots')
xlabel('Seconds')
legend('GCamMP','UV','RCaMP');


% figure; 
% plot(time, data.streams.(STREAM_STORE2).data,'Color',green,'LineWidth',2);
% xlabel('Seconds')
% ylim([0 1000])
% title('GCaMP');
% 
% figure; 
% plot(time, data.streams.(STREAM_STORE1).data,'b','LineWidth',2);
% title('UV');
% xlabel('Seconds')
% ylim([0 1000])

% plot de-trended dF signal
figure;
plot(timedownsample, Y_dF_all,'Color',green,'LineWidth',2);
xlabel('Seconds')
title('Detrended GCaMP');
%% Myles's normalized dFF of GCaMP signal
dFF = (data.streams.(STREAM_STORE2).data - median(data.streams.(STREAM_STORE1).data))/median(data.streams.(STREAM_STORE1).data);
figure; plot(time, dFF,'Color',red,'LineWidth',2);
xlabel('Seconds')
title('Normalized dFF');

%% some stuff that i was doing to show distribution of UV signal
% UV_mean = mean(data.streams.(STREAM_STORE1).data); 
% UV_median = median(data.streams.(STREAM_STORE1).data); 
% UV_std = std(data.streams.(STREAM_STORE1).data); 
% [mu,s,muci,sci] = normfit(data.streams.(STREAM_STORE1).data);
% dist = (1/sqrt(2*pi*sci(1)))*exp(-((data.streams.(STREAM_STORE1).data - mu).^2)/(2*sci(1)));
% figure; plot(data.streams.(STREAM_STORE1).data,dist);
% norm = normpdf(data.streams.(STREAM_STORE1).data, mu, s);
% figure; plot(data.streams.(STREAM_STORE1).data, norm)

%% compile sensor info in table (modulation freq, current info, etc
% test = data.scalars.FiPi.data';
% 
% for ii = 1:8
%     evens(ii) = (ii*1)+ii;
% end
% 
% for ii = 1:numel(evens)
%     sensor_tbl = array2table(test(:,evens),...
%     'VariableNames',{'DAC1','modulation_freq1','current1','DC_offset1','DAC2'...
%     'modulation_freq2','current2','DC_offset2'});
% end
% 
toc
