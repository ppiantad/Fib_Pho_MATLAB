%% Run PhotometryAnalysis_new_Chamber_A_BL2nsp_v4_motion_integrated.m first
%%

v1 = strcat(BLOCKPATH,'\68_RDT2019-10-31T14_38_10.avi'); %set path for video file

video_name_split = strsplit(v1, ["\", " ", "-", "_"]);

clip = 0; % seconds to clip off the front end of the data; FOR USB CAM vids value should be 0; for whatever reason PointGrey seem to line up better @ 2-3

%%

trial_select = 44; %select which trial to view

timestep = 1; % Smaller number = higher resolution data and video

%Gather trial number for video file name
TrialNum = BehavData.Trial(trial_select);

ws = (((BehavData.choiceTime(trial_select)))*30)+(TRANGE(1)*30);
we = ws + TRANGE(2)*30;

info = aviinfo(v1); % display info for avi

obj1 = VideoReader(v1); % read in avi as object

%THESE CALCULATIONS ASSUME TRANGE -33 43; MUST BE CHANGED IF NECESSARY
cho_time1 = abs(ceil(TRANGE(1)*30));
st_time1 = ceil(abs(TRANGE(1)+ (BehavData.choiceTime(trial_select) - BehavData.stTime(trial_select)))*30); %calc diff from stTime, then scale based on TRANGE (puts on same scale as tw array)
con_time1 = ceil(abs(TRANGE(1)- (BehavData.collectionTime(trial_select) - BehavData.choiceTime(trial_select)))*30);


%% SET WINDOW BOUNDS TO VIEW 

wt = we-ws; % range of window

viewstart = 1; % (ceil(ws*fsg1)); %start 
viewend = size(zall(trial_select,:),2); %(ceil(we*fsg1)); %end %
wlength = viewend - viewstart; %window length

wdata_zall = zall(trial_select,:)';
wdata_zmotion = zall_motion(trial_select,:)';

%%
frame_sample = size(read(obj1,1)); % video frame size

tw = (0:(frame_sample(2))/(size(wdata_zall,1)-1):frame_sample(2))'; % time vector of window

ratio = viewend(1) / wt;
f_idx = 1;
wd_idx = 1;
SP = ((wt-we)/wt)*frame_sample(2);
tracestep = -200;
entrystep = tracestep;
imstep = 0;
tracefactor = 50; %scale the signal in frame

%%
for a = 1:length(wlength) %viewend
    figure(1)
    hold on;
    curve = animatedline('LineWidth',1,'Color','w');
    curve2 = animatedline('LineWidth',1,'Color','c');
    stepup = 0;
%     st_txt = text(tw(st_time1), entrystep-75,'St','color','b','VerticalAlignment','top','HorizontalAlignment','center');
%     cho_txt = text(tw(cho_time1),entrystep-75,'Cho','color','r','VerticalAlignment','top','HorizontalAlignment','center');
%     con_txt = text(tw(con_time1),entrystep-75,'Con','color','g','VerticalAlignment','top','HorizontalAlignment','center');
    
%     line([tw(st_time1) tw(st_time1)],[entrystep-125 imstep-25],'Color',[0 1 0],'LineWidth',2); %0 0 1
%     line([tw(cho_time1) tw(cho_time1)],[entrystep-125 imstep-25],'Color',[1 0 0],'LineWidth',2); %tw(33562) corresponds to 33s, which is the "0" point of our -33 43 window
% 
    line([tw(end-30) tw(end)],[entrystep-50 entrystep-50],'Color','w','LineWidth',3);
    lejtxt = text((tw(end)-((tw(end)-tw(end-1))/2)),entrystep-50,'1 sec','color','w','VerticalAlignment','top','HorizontalAlignment','center');
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
        addpoints(curve,tw(idx),(wdata_zall(idx,a)*tracefactor)+entrystep); %addpoints(curve,tw(idx),(wdata_470(idx,a)*tracefactor)+tracestep);
        addpoints(curve2,tw(idx),(wdata_zmotion(idx,a)*tracefactor)+entrystep); %addpoints(curve2,tw(idx),(wdata_405(idx,a)*tracefactor)+tracestep);
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
    plot(tw,(wdata_zall(:,1:wd_idx)*tracefactor)+entrystep, 'Color','w'); %plot(tw,(wdata_470(:,1:wd_idx)*tracefactor)+tracestep, 'Color','w');
    plot(tw,(wdata_zmotion(:,1:wd_idx)*tracefactor)+entrystep, 'Color','c'); %plot(tw,(wdata_405(:,1:wd_idx)*tracefactor)+tracestep, 'Color','c');
    clearpoints(curve)
    clearpoints(curve2)
    wd_idx = wd_idx + 1;
end


%%
%uncomment if you want to write video
%if traces are clipping into video, the video will report a cdata mismatch
%and not be able to write the vid; edit "tracefactor" to re-scale

video = VideoWriter([(video_name_split{13}),'_','RDT','_',(video_name_split{11}),'_', (out{1,:}),'_','Trial',num2str(TrialNum(1))],'MPEG-4');%datestr(now,'HHMMSS''.avi')%['pletdrop',num2str(pletno(1)),'_', filename(end-9:end),'.avi']
video.FrameRate = ceil(30/timestep);
video.Quality = 100;
open(video)
writeVideo(video,F)
close(video)

