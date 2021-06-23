%% Created and modified by Vaibhav Konanur
%   Roitman Lab
%   1007 W Harrison St.
%   Room 1066
%   Chicago, IL 60607
%   Phone: 312-996-8658
%   email: vkonan2@uic.edu
%% LOAD DATA
tic
clear all
reset(gca)
reset(gcf)
close all
clc

clip = 1; % seconds to clip off the front end of the data
filename = '<filepath>'; % insert file path here
pathname = filename(1:(end-11)) % replace this with the path for the folder containing the block
data=TDT2mat(filename); % load tdt file
timestep = 50; % Smaller number = higher resolution data and video
pletno = [2 4 7]; % these are the event numbers of the trails that want to be viewed
% time of window
ws = -3; % window around the event (seconds before)
we = 10; % seconds after event
%% SET VARIABLES

fsg1=floor(data.streams.x65n.fs); % sampling rate
data.sig = double(data.streams.x65n.data((clip*fsg1):(end-(clip*fsg1))))'; % 465nm signal vector
data.baq = double(data.streams.x05n.data((clip*fsg1):(end-(clip*fsg1))))'; % 405nm signal vector
data.entry = ceil(((data.epocs.etry.onset).*fsg1)-(clip*fsg1)); % head entries vector
data.trainOn = ceil(((data.epocs.Lver.onset).*fsg1)-(clip*fsg1)); % event vector
% pletno = 1:length(data.trainOn)-1;
v1 = ([pathname,filename(end-10:end),'\test_experiment-',filename(end-17:end-12),'_',filename(end-10:end),'_Cam1.avi']) % set path of video file

info = aviinfo(v1) % display info for avi
obj1 = VideoReader(v1); % read in avi as object

OGdatalength = length(data.streams.x65n.data); 
OGvidlength = info.NumFrames;
moddatalength = length(data.sig);
nspf = OGdatalength/OGvidlength;

tg1 = (0:1/fsg1:((numel(data.sig)-1)/fsg1))'; % time vector


%%
% Filter 465nm data (if you have your own way to normalize data, use your
% own method

[b,a] = butter(3, 0.00003 , 'high');
[c,d] = butter(3, 0.01,   'low');

    data.sigfilt = filtfilt(b,a,data.sig);
    data.sigfilt = filtfilt(c,d,data.sigfilt);
    data.baqfilt = filtfilt(b,a,data.baq);
    data.baqfilt = filtfilt(c,d,data.baqfilt);




%% SET WINDOW BOUNDS TO VIEW 

wt = we-ws; % range of window

% pellet window

viewstart = data.trainOn(pletno)+(ceil(ws*fsg1));
viewend = data.trainOn(pletno)+(ceil(we*fsg1));

% data of window
wdata = zeros(length(data.sigfilt(viewstart(1):viewend(1))), length(pletno));
for a = 1:length(pletno)
    wdata(:,a) = data.sigfilt(viewstart(a):viewend(a)); % snip out windows of data
end
frame_sample = size(read(obj1,1)); % video frame size
tw = (0:(frame_sample(2))/(size(wdata,1)-1):frame_sample(2))'; % time vector of window

% entry
entry_log = zeros(length(data.sigfilt),1);
entry_log(data.entry) = 1;

eap = zeros(length(data.sigfilt(viewstart(1):viewend(1))), length(pletno));
for a = 1:length(pletno)
    eap(:,a) = entry_log(viewstart(a):viewend(a));
end
eap_shft = (min(data.sigfilt)-1);
eap = (eap.*(2))+eap_shft;
eap(eap==+eap_shft) = NaN;
if isempty(eap(~isnan(eap)))==1
else
    eap1=eap(~isnan(eap));
    eap1=eap1(1);
    asdf=find(~isnan(eap)==1);
    eap2=eap;
    for abc = 1:numel(asdf)-1
        eap2(asdf(abc):asdf(abc)+111)=eap1;
    end
    eap=eap2;
end


f_idx = 1;
wd_idx = 1;
SP = ((wt-we)/wt)*frame_sample(2);
tracestep = -100;
entrystep = tracestep;
imstep = 0;
tracefactor = 2;


for a = 1:length(pletno)
    figure(1)
    hold on;
    curve = animatedline('LineWidth',2,'Color','w');
    curve2 = animatedline('LineWidth',3,'Color','c');
    stepup = 0;
    enttxt1 = text(0,entrystep-75,'Head','color','c','VerticalAlignment','bottom');
    enttxt2 = text(0,entrystep-75,'Entry','color','c','VerticalAlignment','top');
    cuetxt1 = text(SP,entrystep-75,'Light','color','r','VerticalAlignment','bottom','HorizontalAlignment','center');
    cuetxt2 = text(SP,entrystep-75,'Cue','color','r','VerticalAlignment','top','HorizontalAlignment','center');
    line([SP SP],[entrystep-50 imstep-25],'Color',[1 0 0],'LineWidth',2);
    line([tw(end-fsg1) tw(end)],[entrystep-50 entrystep-50],'Color','w','LineWidth',3);
    lejtxt = text((tw(end)-((tw(end)-tw(end-fsg1))/2)),entrystep-50,'1 sec','color','w','VerticalAlignment','top','HorizontalAlignment','center');
    for idx = viewstart(a):timestep:viewend(a);
        figure(1)
        hold on
        set(gca,'XLim',[tw(1) tw(end)],'YLim',[tracestep*2 frame_sample(1)],'FontSmoothing', 'off');
        set(gcf,'Position',[680   354   560   624],'color','k');
        frame = ceil((idx+(clip*fsg1)-1)/nspf);
        frame1 = read(obj1,frame);
        J = imtranslate(flipud(frame1), [0 -imstep], 'FillValues',255,'OutputView','full');
        image(J)
        axis off
        idx = idx-viewstart(a)+1;
        addpoints(curve,tw(idx),(wdata(idx,a)*tracefactor)+tracestep);
        addpoints(curve2,tw(idx),(eap(idx,a)+entrystep-10));
        drawnow
        figure(1)
        F(f_idx) = getframe(gcf);
        stepup = stepup + timestep;
        f_idx = f_idx + 1;
    end
    delete(lejtxt)
    delete(enttxt1)
    delete(enttxt2)
    delete(cuetxt1)
    delete(cuetxt2)
    plot(tw,(wdata(:,1:wd_idx)*tracefactor)+tracestep, 'Color',[0.5 0.5 0.5]);
    clearpoints(curve)
    clearpoints(curve2)
    wd_idx = wd_idx + 1;
end
%%
video = VideoWriter(['pletdrop',num2str(pletno(1)),'_', filename(end-9:end),'_','.avi'],'MPEG-4');%datestr(now,'HHMMSS''.avi')%['pletdrop',num2str(pletno(1)),'_', filename(end-9:end),'.avi']
video.FrameRate = ceil(fsg1/timestep);
video.Quality = 100;
open(video)
writeVideo(video,F)
close(video)


