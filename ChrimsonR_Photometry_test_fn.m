function [ts1, zall] = ChrimsonR_Photometry_test_fn(animalName, blockpath, whichStreams, Pu1_epoc)


data = TDTbin2mat(blockpath);

%%
Pulse_epoc=fieldnames(data.epocs);
if Pu1_epoc==5
   epoc1=5; 
elseif Pu1_epoc==3
   epoc1=3;
end

streamNames=fieldnames(data.streams);
if whichStreams==12
    stream1=2; stream2=1;
elseif whichStreams==34
    stream1=3; stream2=4;
end


%%
%create epoc for trials if using SYNAPSE TO CONTROL LASER (40 pulses in
%sypanse for 20 Hz, collects onset for every laser by skipping all the
%interleved pulses


% %create epoc for trials
% data.epocs.Session.name = 'Session';
%     
% data.epocs.Session.onset = data.epocs.(Pulse_epoc{epoc1}).onset(1:40:end);
% data.epocs.Session.offset = data.epocs.(Pulse_epoc{epoc1}).offset(1:40:end);
% 
% data.epocs.Session.data = ones(1,size(data.epocs.Session.onset,1));
    


%create epoc for trials
data.epocs.Session.name = 'Session';

if size(data.epocs.(Pulse_epoc{epoc1}).onset,1) > 20
    
    data.epocs.Session.onset = data.epocs.(Pulse_epoc{epoc1}).onset(1:40:end);
    data.epocs.Session.offset = data.epocs.(Pulse_epoc{epoc1}).offset(1:40:end);

    data.epocs.Session.data = ones(1,size(data.epocs.Session.onset,1));
    
elseif size(data.epocs.(Pulse_epoc{epoc1}).onset,1) <= 20
    data.epocs.Session.onset = data.epocs.(Pulse_epoc{epoc1}).onset(1:end);
    data.epocs.Session.offset = data.epocs.(Pulse_epoc{epoc1}).offset(1:end);
    data.epocs.Session.data = ones(1,size(data.epocs.Session.onset,1));
end

%create epoc for trials if using PULSEPAL VIA ABET TO CONTROL LASER

% data.epocs.Session.name = 'Session';
% data.epocs.Session.onset =data.epocs.TL2_.onset(2:end);
% data.epocs.Session.offset = data.epocs.TL2_.offset(2:end);
% 
% data.epocs.Session.data = ones(1,size(data.epocs.Session.onset,1));




%%
%FOR SENSOR G (temporary, add code to select which STREAMSTORE later)
%use TDT filter to get just the trials you want
EPOC = 'Session'; %set which epoc to look at
STREAM_STORE1 = streamNames{stream1}; %put here what your 405 channel is It's ax05A for others! 'x05A'
STREAM_STORE2 = streamNames{stream2}; %put here what your 465 channel is ax70A 'x70A'
STREAM_STORE3 = 'Fi1r'; %aFi1r
TRANGE = [-5 15]; %  window si ze [start time relative to epoc onset, window duration]
BASELINE_PER = [-5 0]; % baseline period within our window
ARTIFACT = Inf; % optionally set an artifact rejection level
data = TDTfilter(data,EPOC,'TIME',TRANGE); %perform TDTfilter 

% %for Sensor D
% %use TDT filter to get just the trials you want
% EPOC = 'Session'; %set which epoc to look at
% STREAM_STORE1 = 'd5sP'; %put here what your 405 channel is It's ax05A for others!
% STREAM_STORE2 = 'd0sP'; %put here what your 465 channel is ax70A
% % STREAM_STORE3 = 'Fi1r'; %aFi1r
% TRANGE = [-5 35]; %  window si ze [start time relative to epoc onset, window duration]
% BASELINE_PER = [-5 0]; % baseline period within our window
% ARTIFACT = Inf; % optionally set an artifact rejection level
% data = TDTfilter(data,EPOC,'TIME',TRANGE); %perform TDTfilter 



%%
% Optionally remove artifacts. If any waveform is above ARTIFACT level, or
% below -ARTIFACT level, remove it from the data set.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
good2 = ~art1 & ~art2;
data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);

numArtifacts = sum(~good) + sum(~good2);

%%
% Applying a time filter to a uniformly sampled signal means that the
% length of each segment could vary by one sample.  Let's find the minimum
% length so we can trim the excess off before calculating the mean.
minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
% minLength3 = min(cellfun('prodofsize', data.streams.(STREAM_STORE3).filtered));
data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);
% data.streams.(STREAM_STORE3).filtered = cellfun(@(x) x(1:minLength3), data.streams.(STREAM_STORE3).filtered, 'UniformOutput',false);


allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');

% downsample 10x and average 405 signal
N = 10;
F405 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405,1);
stdSignal1 = std(double(F405))/sqrt(size(F405,1));
dcSignal1 = mean(meanSignal1);

% downsample 10x and average 465 signal
allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
F465 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength2 = size(F465,2);

% Create mean signal, standard error of signal, and DC offset of 465 signal
meanSignal2 = mean(F465,1);
stdSignal2 = std(double(F465))/sqrt(size(F465,1));
dcSignal2 = mean(meanSignal2);

% % downsample 10x and average 465 signal
% allSignals = cell2mat(data.streams.(STREAM_STORE3).filtered');
% Fraw = zeros(size(allSignals(:,1:N:end-N+1)));
% for ii = 1:size(allSignals,1)
%     Fraw(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
% end
% minLength3 = size(Fraw,2);





% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
% ts3=TRANGE(1)+(1:minLength3) / data.streams.(STREAM_STORE3).fs*N;

% Subtract DC offset to get signals on top of one another
meanSignal1 = meanSignal1 - dcSignal1;
meanSignal2 = meanSignal2 - dcSignal2;



%% 


% Fitting 405 channel onto 465 channel to detrend signal bleaching
% Scale and fit data
% Algorithm sourced from Tom Davidson's Github:
% https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m

bls = polyfit( F405(1:end),F465(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F465 - Y_fit_all;

zall = zeros(size(Y_dF_all));
tmp = 0;
for i = 1:size(Y_dF_all,1)
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean (-6 to -1 seconds)
    zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    for j = 1:size(Y_dF_all,2) % Z score per bin
        tmp = tmp + 1;
        zall(i,tmp)=(Y_dF_all(i,j) - zb)/zsd;
    end
    tmp=0;
end

% Standard error of the z-score
zerror = std(zall)/size(zall,1);


%BASELINE FOR US
% 
% zBL4US = zeros(size(Y_dF_all));
% tmp = 0;
% for i = 1:size(Y_dF_all,1)
%     ind = ts2(1,:) < 28 & ts2(1,:) > 23;
%     zb2 = mean(Y_dF_all(i,ind)); % baseline period mean (-6 to -1 seconds)
%     zsd2 = std(Y_dF_all(i,ind)); % baseline period stdev
%     for j = 1:size(Y_dF_all,2) % Z score per bin
%         tmp = tmp + 1;
%         zBL4US(i,tmp)=(Y_dF_all(i,j) - zb2)/zsd2;
%     end
%     tmp=0;
% end
% 



end
