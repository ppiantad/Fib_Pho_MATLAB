%% Created and modified by Vaibhav Konanur
%   Roitman Lab
%   1007 W Harrison St.
%   Room 1066
%   Chicago, IL 60607
%   Phone: 312-996-8658
%   email: vkonan2@uic.edu
%   Edits by Patrick Piantadosi (Holmes lab @ NIAAA)

%%  LOAD DATA
tic
clear all
reset(gca)
reset(gcf)
close all
clc
root = upper('e:');

BLOCKPATH =  strcat(root,'\MATLAB\TDTbin2mat\Photometry\143 & D1-eOP-1\143_D1-eOP-1-201023-102841'); 



v1 = strcat(BLOCKPATH,'\D1-eOP-1_RM_CHOICE_5mW2020-10-23T10_28_43.avi'); %set path for video file

video_name_split = strsplit(v1, ["\", " ", "-", "_"]);

data = TDTbin2mat(BLOCKPATH);



clip = 0; % seconds to clip off the front end of the data; FOR USB CAM vids value should be 0; for whatever reason PointGrey seem to line up better @ 2-3

%%


%%

Sensor = input('Please enter name of sensor used (A,B, C, or D): \n','s'); %input which sensor the mouse was recorded from 
Sensor = upper(Sensor);

%find time of first TTL.  This value will be passed to Behavioral Data so
%that the timestamps of the ABET data will match the timestamps of Synapse.
%This does not work perfect for mice < #RRD47, because Sensor D was in Box
%1 (with TTL_), and Sensor C was in box 2 (with TL2_)
if Sensor == 'A';
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'B';
    timeStart=data.epocs.TTL_.onset(1);  
elseif Sensor == 'C';
    timeStart=data.epocs.TTL_.onset(1);
elseif Sensor == 'D';
    timeStart=data.epocs.TL2_.onset(1);
end

%work behavioral data into photo data, generate descriptive statistics
[BehavData,ABETfile,Descriptives]=ABET2TableFn_Chamber_Av3('D1-eOP-1 10232020.csv');



%align timestamps of behavioral data to timestamps of photometry by adding
%the time elapsed from when Synapse recording began to when the ABET
%program was started back to each relevant time column (choiceTime,
%collectionTime, stTime)
timeShift=timeStart*ones(numel(BehavData.choiceTime(:)),1);
BehavData.choiceTime(:)=BehavData.choiceTime(:)+timeShift;
BehavData.collectionTime(:)=BehavData.collectionTime(:)+timeShift;
BehavData.stTime(:)=BehavData.stTime(:)+timeShift;

%%
%pick trial types

[BehavData, ~, filtervars] =TrialFilter(BehavData,'ALL',1);

%Convert filtervars to str so that we can use in the video file name
out=cellfun(@num2str,filtervars,'un',0);



%%
%create epoc for trials aligned to choice
data.epocs.Choice.name = 'Choice';
data.epocs.Choice.onset = BehavData.choiceTime; %
data.epocs.Choice.offset = BehavData.choiceTime+0.1*ones(numel(BehavData.choiceTime),1);
data.epocs.Choice.data = [1:numel(BehavData.choiceTime)]';


%create epoc for trials aligned to collection time
data.epocs.Collect.name = 'Collect';
data.epocs.Collect.onset = BehavData.collectionTime; %
data.epocs.Collect.offset = BehavData.collectionTime+0.1*ones(numel(BehavData.collectionTime),1);
data.epocs.Collect.data = [1:numel(BehavData.collectionTime)]';


%create epoc for trials aligned to Trial Start time
data.epocs.start.name = 'start';
data.epocs.start.onset = BehavData.stTime; %
data.epocs.start.offset = BehavData.stTime+0.1*ones(numel(BehavData.stTime),1);
data.epocs.start.data = [1:numel(BehavData.stTime)];

%%
%use TDT filter to select just the trials you want
EPOC = 'Choice'; %set which epoc to look at
StreamFields = fieldnames(data.streams); 
% check for known names that we have used for 405 / 470 channels here,
% selects the one that is used for this particular session. if name
% stream names change, must be added to this loop manually
if Sensor == 'A';
    STREAM_STORE1 = 'x405S'
    STREAM_STORE2 = 'x470S'
elseif Sensor =='B';
    STREAM_STORE1 = 'x405G'
    STREAM_STORE2 = 'x470G'
elseif Sensor == 'D';
    STREAM_STORE1 = 'd5sP'
    STREAM_STORE2 = 'd0sP'
elseif Sensor == 'C';
    STREAM_STORE1 = 'x45sP'
    STREAM_STORE2 = 'x40sP'
end
% STREAM_STORE1 = input('Please enter name of 405 stream: \n','s'); %check data.streams 
% STREAM_STORE2 = input('Please enter name of 470 stream: \n','s'); %check data.streams
% if ismember('x405G',StreamFields);
%     STREAM_STORE1 = 'x405G';
% elseif ismember('x405S',StreamFields);
%     STREAM_STORE1 = 'x405S';
% elseif ismember('x45sP',StreamFields);
%     STREAM_STORE1 = 'x45sP';
% end;
% 
% if ismember('x470G',StreamFields);
%     STREAM_STORE2 = 'x470G';
% elseif ismember('x470S',StreamFields);
%     STREAM_STORE2 = 'x470S';
% elseif ismember('x40sP',StreamFields);
%     STREAM_STORE2 = 'x40sP';
% end;
% if ~any(strcmp(StreamFields,'x405G'));
%      STREAM_STORE1 = 'x405G';
% elseif ~any(strcmp(StreamFields, 'x405S'));
%     STREAM_STORE1 = 'x405S';
% elseif ~any(strcmp(StreamFields, 'x45sP'));
%     STREAM_STORE1 = 'x45sP';
% end
% STREAM_STORE1 = isfields(data.streams('x405S')); %UV stream? can we not add to code to check if x405G / x405S / x45sP exist, and if so select that as the STREAM_STORE1 name?
% STREAM_STORE2 = isfields(data.streams('x470G','x470S','x45sP')); %GCamp stream?
%TRANGE sets the window size, data will be filtered around EPOC (whatever event), size of data will be window duration (2nd # in TRANGE) x sampling
%rate (1017) approximately
TRANGE = [-33 43]; % window size [start time relative to epoc onset, window duration] -10 20 default -4 12 -10 18 -31 41
BASELINE_PER = [-3 -1]; % baseline period within our window -6 -1 default; THE FIRST VALUE CAN'T BE < THE FIRST TRANGE VALUE -8 0 is ITI
ARTIFACT = Inf; % optionally set an artifact rejection level
data = TDTfilter(data,EPOC,'TIME',TRANGE); %perform TDTfilter 
%%
% if EPOC == 'Trial';
%     BL_time = BehavData.stTime(:)-BehavData.choiceTime(:);
%     time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
% elseif EPOC == 'Collect';
%     BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
% end

time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);

if strcmpi('start',EPOC)
    time2EPOC = BehavData.stTime(:) - BehavData.choiceTime(:);

elseif strcmpi('Collect', EPOC)
    BL_time = BehavData.collectionTime(:)-BehavData.choiceTime(:);
    time2EPOC = BehavData.choiceTime(:)-BehavData.collectionTime(:); %time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
elseif strcmpi('Choice', EPOC)
    BL_time = BehavData.stTime(:)-BehavData.choiceTime(:);
%     time2Collect = BehavData.collectionTime(:)-BehavData.choiceTime(:);
    time2EPOC = BehavData.stTime(:)-BehavData.choiceTime(:);
end


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
data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);

allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');

% downsample 10x and average 405 signal
N = 10;
F405 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405);
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
meanSignal2 = mean(F465);
stdSignal2 = std(double(F465))/sqrt(size(F465,1));
dcSignal2 = mean(meanSignal2);

% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
% ts3 = (0:1/fsg1:((size(F465,2)-1)/fsg1))';


%% SET VARIABLES

trial_select = 44; %select which trial to view

timestep = 300; % Smaller number = higher resolution data and video

%Gather trial number for video file name
TrialNum = BehavData.Trial(trial_select);

fsg1=floor(data.streams.(STREAM_STORE2).fs); % sampling rate
data.sig = double(data.streams.(STREAM_STORE2).filtered{1, trial_select}'); % 465nm signal vector
data.sig = data.sig - mean(data.sig);
% data.sig = double(data.sig((fsg1):(end-(fsg1)))); % 465nm signal vector
data.baq = double(data.streams.(STREAM_STORE1).filtered{1, trial_select}'); % 405nm signal vector
data.baq = data.baq - mean(data.baq);
% data.baq = double(data.baq((fsg1):(end-(fsg1))));% 405nm signal vector
% data.sig = double(data.streams.(STREAM_STORE2).data((clip*fsg1):(end-(clip*fsg1)))'); % 465nm signal vector
% data.baq = double(data.streams.(STREAM_STORE1).data((clip*fsg1):(end-(clip*fsg1)))'); % 405nm signal vector

ws = (((BehavData.choiceTime(trial_select)))*30)+(TRANGE(1)*30);
we = ws + TRANGE(2)*30;


% ws = 1; % desired window startts1(1)
% we = ts1(end); % desired window end
% data.entry = ceil(((data.epocs.entry.onset).*fsg1)-(clip*fsg1)); % head entries vector
% data.trainOn = ceil(((data.epocs.Lver.onset).*fsg1)-(clip*fsg1)); % event vector
% pletno = 1:length(data.trainOn)-1;

info = aviinfo(v1) % display info for avi

obj1 = VideoReader(v1); % read in avi as object

OGdatalength = length(data.streams.(STREAM_STORE2).data); 
OGvidlength = info.NumFrames;
% OGvidlength = obj1.FrameRate(1); %read frames for .mp4s
moddatalength = length(data.sig);
nspf = OGdatalength/OGvidlength;

tg1 = (0:1/fsg1:((numel(data.sig)-1)/fsg1))'; % time vector


%%
% Filter 465nm data (if you have your own way to normalize data, use your
% own method

% [b,a] = butter(3, 0.00003 , 'high'); % filter out low freq oscillations
% [c,d] = butter(3, 0.001,   'low'); % filter out high freq noise, 0.01 for less filtering, 0.001 for more etc.
% 
%     data.sig = filtfilt(b,a,data.sig);
%     data.sig = filtfilt(c,d,data.sig);
%     data.baq = filtfilt(b,a,data.baq);
%     data.baq = filtfilt(c,d,data.baq);




%% SET WINDOW BOUNDS TO VIEW 

wt = we-ws; % range of window


viewstart = 1; % (ceil(ws*fsg1)); %start 
viewend = size(data.sig,1); %(ceil(we*fsg1)); %end %
wlength = viewend - viewstart; %window length
% data window


wdata_470 = data.sig; %zeros(length(data.sig(viewstart(1):viewend(1))), length(wlength));
wdata_405 = data.baq;


% for a = 1:length(wlength)
%     wdata_470(:,a) = data.sig(viewstart(a):viewend(a)); % snip out windows of 465 data
%     wdata_405(:,a) = data.baq(viewstart(a):viewend(a)); % snip out windows of 405 data
%     tw_test(:,a) = tg1(viewstart(a):viewend(a));
% end
frame_sample = size(read(obj1,1)); % video frame size

tw = (0:(frame_sample(2))/(size(wdata_470,1)-1):frame_sample(2))'; % time vector of window


ratio = viewend(1) / wt;
f_idx = 1;
wd_idx = 1;
SP = ((wt-we)/wt)*frame_sample(2);
tracestep = -200;
entrystep = tracestep;
imstep = 0;
tracefactor = 1; %scale the signal in frame

%THESE CALCULATIONS ASSUME TRANGE -33 43; MUST BE CHANGED IF NECESSARY
cho_time1 = abs(ceil(TRANGE(1)*fsg1));
st_time1 = ceil(abs(TRANGE(1)+ (BehavData.choiceTime(trial_select) - BehavData.stTime(trial_select)))*fsg1); %calc diff from stTime, then scale based on TRANGE (puts on same scale as tw array)
con_time1 = ceil(abs(TRANGE(1)- (BehavData.collectionTime(trial_select) - BehavData.choiceTime(trial_select)))*fsg1);

%%
% wdata_dFF = (wdata_470 - median(wdata_405))/median(wdata_405);

%%
% line 116 changed from viewend to moddatalength
%line 130 changed from viewstart(a):timestep:viewend(a); to
%viewstart(a):timestep:moddatalength(a)
for a = 1:length(wlength) %viewend
    figure(1)
    hold on;
    curve = animatedline('LineWidth',1,'Color','w');
    curve2 = animatedline('LineWidth',1,'Color','c');
    stepup = 0;
%     st_txt = text(tw(st_time1), entrystep-75,'St','color','b','VerticalAlignment','top','HorizontalAlignment','center');
%     cho_txt = text(tw(cho_time1),entrystep-75,'Cho','color','r','VerticalAlignment','top','HorizontalAlignment','center');
%     con_txt = text(tw(con_time1),entrystep-75,'Con','color','g','VerticalAlignment','top','HorizontalAlignment','center');
    
    line([tw(st_time1) tw(st_time1)],[entrystep-125 imstep-25],'Color',[0 1 0],'LineWidth',2); %0 0 1
    line([tw(cho_time1) tw(cho_time1)],[entrystep-125 imstep-25],'Color',[1 0 0],'LineWidth',2); %tw(33562) corresponds to 33s, which is the "0" point of our -33 43 window
%     line([tw(37634) tw(37634)],[entrystep-50 imstep-25],'Color',[0 1 0],'LineWidth',1);
%     line([tw(con_time1) tw(con_time1)],[entrystep-50 imstep-25],'Color',[0 1 0],'LineWidth',1);
    
%     line([tw(37634) tw(37634)],[entrystep-50 imstep-25],'Color',[0 1 0],'LineWidth',1);
%     line([tw(33562) tw(33562)],[entrystep-50 imstep-25],'Color',[1 0 0],'LineWidth',1); %tw(33562) corresponds to 33s, which is the "0" point of our -33 43 window
    line([tw(end-fsg1) tw(end)],[entrystep-50 entrystep-50],'Color','w','LineWidth',3);
    lejtxt = text((tw(end)-((tw(end)-tw(end-fsg1))/2)),entrystep-50,'1 sec','color','w','VerticalAlignment','top','HorizontalAlignment','center');
    for idx = viewstart(a):timestep:viewend(a);
        figure(1)
        hold on
        set(gca,'XLim',[tw(1) tw(end)],'YLim',[tracestep*2 frame_sample(1)],'FontSmoothing', 'off');
        set(gcf,'Position',[680   0   560   624],'color','k'); %680   354   560   624
        frame = ceil(ws + (idx / ratio)); %take the WINDOW START time (first frame in the video to play), then increment by an idx that is proportional to the size of the photometry signal / the window size
        frame1 = read(obj1,frame);
        J = imtranslate(flipud(frame1), [0 -imstep], 'FillValues',255,'OutputView','full');
        image(J)
        axis off
        idx = idx-viewstart(a)+1;
        addpoints(curve,tw(idx),(wdata_470(idx,a)*tracefactor)+entrystep); %addpoints(curve,tw(idx),(wdata_470(idx,a)*tracefactor)+tracestep);
        addpoints(curve2,tw(idx),(wdata_405(idx,a)*tracefactor)+entrystep); %addpoints(curve2,tw(idx),(wdata_405(idx,a)*tracefactor)+tracestep);
        drawnow
        figure(1)
        F(f_idx) = getframe(gcf);
        stepup = stepup + timestep;
        f_idx = f_idx + 1;
    end
    delete(lejtxt)
%     delete(enttxt1)
%     delete(enttxt2)
%     delete(cuetxt1)
%     delete(cuetxt2)
    plot(tw,(wdata_470(:,1:wd_idx)*tracefactor)+entrystep, 'Color','w'); %plot(tw,(wdata_470(:,1:wd_idx)*tracefactor)+tracestep, 'Color','w');
    plot(tw,(wdata_405(:,1:wd_idx)*tracefactor)+entrystep, 'Color','c'); %plot(tw,(wdata_405(:,1:wd_idx)*tracefactor)+tracestep, 'Color','c');
    clearpoints(curve)
    clearpoints(curve2)
    wd_idx = wd_idx + 1;
end
%%
%uncomment if you want to write video
%if traces are clipping into video, the video will report a cdata mismatch
%and not be able to write the vid; edit "tracefactor" to re-scale

video = VideoWriter([(video_name_split{13}),'_','RDT','_',(video_name_split{11}),'_', (out{1,:}),'_','Trial',num2str(TrialNum(1))],'MPEG-4');%datestr(now,'HHMMSS''.avi')%['pletdrop',num2str(pletno(1)),'_', filename(end-9:end),'.avi']
video.FrameRate = ceil(fsg1/timestep);
video.Quality = 100;
open(video)
writeVideo(video,F)
close(video)


