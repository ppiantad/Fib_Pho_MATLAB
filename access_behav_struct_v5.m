%%
% to access BehavData from inside behav_struct for trial filtering


iter = 0;
% These .mat files contain 4 files processed using AllAnimalFPanalysis
% with the Analysis_BL2initiation_fn_beta.m function. The files are as
% follows:

    %behav_struct
    %corrTable
    %ts1_struct
    %zall_struct

% As of 2022 - the .mat files should contain the following: 

    %behav_struct
    %corrTable
    %ts1_struct
    %zall_struct
    %zall_motion_struct
    %Y_dF_all_struct

% Uncomment the particular file below to load the .mat file (will create a list of inputs later in a function to tidy this up)  
% Can uncomment and load diff. sessions to graph data from separate groups,
% etc (e.g. vmPFC RM Early REW 1.2 vs vmPFC RM Late REW 1.2)

% load('BLA-PL_RM_EARLY_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_05062023.mat')

% load('BLA-PL_RM_LATE_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_05062023.mat')

% load('BLA-PL_RDT_SESSIONS_ALL_EPOC_start_TRANGE_-8_38_Z-from_window_05062023.mat')

% load('BLA-PL_RDT_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_05062023.mat')

% load('BLA-PL_RDT_SESSIONS_ALL_EPOC_Collect_TRANGE_-8_38_Z-from_window_05062023.mat')

% load('vmPFC-NAcSh_RM_EARLY_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_12112022_includes_aborts_and_motion.mat')

% load('vmPFC-NAcSh_RM_LATE_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_12112022_includes_aborts_and_motion.mat')

% load('vmPFC-NAcSh_RDT_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_12112022_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RM_EARLY_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_08312023_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RM_Early_SESSIONS_ALL_EPOC_Collect_TRANGE_-4_12_Z-from_window_05082023_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RM_LATE_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_08312023_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RM_LATE_SESSIONS_ALL_EPOC_Collect_TRANGE_-4_12_Z-from_window_05082023_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RDT_SESSIONS_ALL_EPOC_start_TRANGE_-8_38_Z-from_window_04272023.mat')

load('BLA-NAcSh_RDT_SESSIONS_ALL_EPOC_Choice_10232023.mat')

% load('BLA-NAcSh_RDT_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_-10_to_-5_02022023_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RDT_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_-10_to_-5_adjusted_to_ITI_02022023_includes_aborts_and_motion.mat')

% load('BLA-NAcSh_RDT_SESSIONS_ALL_EPOC_collect_TRANGE_-33_43_Z-from_window_07062022_includes_aborts_and_motion.mat')

% load('vHPC-NAcSh_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_Z-from_window_07242022_includes_motion_and_Y_df_struct.mat')

% load('vHPC-NAcSh_RDT_SESSIONS_ALL_EPOC_collect_TRANGE_-33_43_Z-from_window_07242022_includes_motion_and_Y_df_struct.mat')

% load('GFP_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5_11192021.mat')

% load('GFP_RM_LATE_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_09302022_includes_aborts_and_motion.mat')

% load('All_GFP_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_Z-from_window_07062022_includes_motion_and_Y_df_struct.mat')

% load('D1_RewMag_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D2_RewMag_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D1_RM_EARLY_1_or_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('D2_RM_EARLY_1_or_2_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5.mat')

% load('D1_RM_LATE_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_09302022_includes_aborts_and_motion.mat')

% load('D2_RM_LATE_SESSIONS_ALL_EPOC_Choice_TRANGE_-33_43_Z-from_window_09302022_includes_aborts_and_motion.mat')

% load('D1_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5_01212022.mat')

% load('D1_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_Z-from_window_04252022_includes_motion_and_Y_df_struct.mat')

% load('D2_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_Z-from_window_04252022_includes_motion_and_Y_df_struct.mat')

% load('D1_RDT_SESSIONS_ALL_EPOC_collect_TRANGE_-33_43_Z-from_window_07122022_includes_aborts_and_motion.mat')

% load('D2_RDT_SESSIONS_ALL_EPOC_collect_TRANGE_-33_43_Z-from_window_07122022_includes_aborts_and_motion.mat')

% load('D2_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5_01212022.mat')

% load('D1_RDT_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D2_RDT_ALL_TRANGE_-33_43_EPOC_Choice_BL_-3_0_BEFORE_INITIATION')

% load('D1-iOP_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5_10162021.mat')

% load('D1-iOP_RDT_SESSIONS_ALL_EPOC_start_TRANGE_-33_43_Z-from_window_07052022_includes_motion_and_Y_df_struct.mat')

% load('D1-iOP_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_Z-from_window_07022022_includes_motion_and_Y_df_struct.mat')

% load('D1-eOP_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_BL_-10_-5_10162021.mat')

% load('D1-eOP_RDT_SESSIONS_ALL_EPOC_start_TRANGE_-33_43_Z-from_window_07052022_includes_motion_and_Y_df_struct.mat')

% load('D1-eOP_RDT_SESSIONS_ALL_EPOC_CHOICE_TRANGE_-33_43_Z-from_window_07022022_includes_motion_and_Y_df_struct.mat')

% load('D2-eOP_RDT_SESSIONS_ALL_EPOC_Choice_10222023.mat')


%because depending on when the data were analyzed, the size of the
%ABET2Table generated behavioral data structure differs (contains slightly
%diff columns), check the size of the table and adjust the titles
%accordingly. This is hard coded and should be verified if anything odd
%happens, or any unusual analyses need to be conducted
if size(behav_struct(1).data,2) == 20
    titles = {'Trial','Block','ForceFree','bigSmall','RewSelection', 'TrialPossible', 'stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','smallRew','bigRew', 'type_binary'};
elseif size(behav_struct(1).data,2) == 19
    titles = {'Trial','Block','ForceFree','bigSmall','RewSelection', 'TrialPossible', 'stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','smallRew','bigRew'};
elseif size(behav_struct(1).data,2) == 17
    titles = {'Trial','Block','ForceFree','bigSmall','stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','smallRew','bigRew'};
elseif size(behav_struct(1).data,2) == 21
    titles = {'Trial','Block','ForceFree','bigSmall','RewSelection', 'TrialPossible', 'stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','lose_stay', 'smallRew','bigRew', 'type_binary'};
end


% Sort table by Risk % (> 80 = Risky, < 80 = not risky)
%calculate the % Large Rew Choices as a function of all free choice trials (66).
%This is more accurate than the original RiskPercent column, because it
%factors in that mice who did not complete the session were erroneously
%listed as having a high risk %, because it was originally calculated as
%Large Rew Choices / Num of Completed Free Choice Trials

for ii = 1:size(corrTable)
    corrTable.RiskPercent_91Trials(ii) = ((corrTable.TotalWins(ii) + corrTable.TotalLosses(ii))/66)*100;
end



for ii = 1:size(corrTable)
    if corrTable.RiskPercent_91Trials(ii) >= 70
        corrTable.RiskSeeking(ii) = 1;
%         behav_struct_risky(ii) = behav_struct(ii);
%         zall_struct_risky(ii) = zall_struct.data(ii);
    elseif corrTable.RiskPercent_91Trials(ii) < 70
        corrTable.RiskSeeking(ii) = 2;
%         behav_struct_NOTrisky(ii) = behav_struct(ii);
%         zall_struct_NOTrisky(ii) = zall_struct.data(ii);
    end  
end



%% Filter by SESSION #

include_session = [2];

for cc = 1:size(corrTable);
    if isempty(include_session);
        disp(['Using all sessions']);
    elseif include_session == 1
        idx = corrTable.session == 1;
        corrTable= corrTable(idx,:);
        behav_struct = behav_struct(idx);
        zall_struct.data = zall_struct.data(idx);
        zall_motion_struct.data = zall_motion_struct.data(idx);
        disp(['Using session 1 only']);
        
    elseif include_session == 2
        idx = corrTable.session == 2;
        corrTable= corrTable(idx,:);
        behav_struct = behav_struct(idx);
        zall_struct.data = zall_struct.data(idx);
        zall_motion_struct.data = zall_motion_struct.data(idx);
        disp(['Using session 2 only']);
    end
end
 
%% Filter by IMPLANT SIDE & Large Rew CONTRA or IPSI

% use side_match = 'ipsi' or side_match = 'contra' or side_match = [];

side_match = 'contra';


for cc = 1:size(corrTable);
    if isempty(side_match)
        disp(['Using all sessions']);
    elseif strcmp(side_match,'ipsi')
        idx = corrTable.implant_side == corrTable.large_rew_side;
        corrTable= corrTable(idx,:);
        behav_struct = behav_struct(idx);
        zall_struct.data = zall_struct.data(idx);
        zall_motion_struct.data = zall_motion_struct.data(idx);
        disp(['Using only IPSI fiber + large reward screen'])
    elseif strcmp(side_match,'contra')
        idx = corrTable.implant_side ~= corrTable.large_rew_side;
        corrTable= corrTable(idx,:);
        behav_struct = behav_struct(idx);
        zall_struct.data = zall_struct.data(idx);
        zall_motion_struct.data = zall_motion_struct.data(idx);
        disp(['Using only CONTRA fiber + large reward screen'])
    end
end




%% limit analyses to just mice with SLEAP motion data

only_motion_sessions = [2];



for cc = 1:size(corrTable);
    if isempty(only_motion_sessions);
        disp(['Using all sessions']);
    elseif only_motion_sessions == 1
        emptyCells = cellfun(@isempty,corrTable.SLEAP_files);
        corrTable= corrTable(~emptyCells,:);
        behav_struct = behav_struct(~emptyCells);
        zall_struct.data = zall_struct.data(~emptyCells);
        zall_motion_struct.data = zall_motion_struct.data(~emptyCells);
        disp(['Using only sessions with valid motion data']);
        
    end
end

%%
clear zall_filtered; clear trials; clear ii; clear jj; clear ZallMean; clear collect_mice; clear ZallMean_motion; clear zall_motion_filtered; clear behav_struct_filtered; clear xcf_struct; clear lags_struct;
            

collect_times = [];
choice_times = [];
struct_AUC = [];
struct_choice_lat = [];
struct_collect_lat = [];
mean_choice_lat = [];
mean_collect_lat = [];
num_trials = [];
% collect_mice  = [];


ts1 = cell2mat(ts1_struct.data);
N = 33.8822618125484; %this is the current downsample factor thanks to needing to downsample the motion and calcium data similarly.

%change this range to set the crosscorrelation window
xcor_range = [-10 10];
%set the number of lags equal to the range of the xcor window
numLags = abs(xcor_range(1)-xcor_range(2));
xcor_ind = ts1 > xcor_range(1) & ts1 < xcor_range(2);
xcor_lag_bins = sum(xcor_ind)/N;


ZallMean = NaN(size(behav_struct,2), size(ts1, 2));
ZallMean_motion = NaN(size(behav_struct,2), size(ts1, 2));


row = 1; 
for ii = 1:size(behav_struct,2)
%     arg1 = {'REW' 1.2}
%     arg2 = {'BLOCK' 3}
    BehavData = cell2mat(behav_struct(ii).data);
    BehavData = array2table(BehavData);
    
%     if exist('behav_struct_risky') == 1
%         BehavDataRisky = cell2mat(behav_struct_risky(ii).data);
%         BehavDataRisky = array2table(BehavDataRisky);
%         BehavDataNOTRisky = cell2mat(behav_struct_NOTrisky(ii).data);
%         BehavDataNOTRisky = array2table(BehavDataNOTRisky);
%     end
    
    %New version of ABET2Tbl add a few columns, RewSelection (BehavData5) and TrialPossible (BehavData6)
    %Delete these columns to proceed with analysis - can be added back in later
%     if size(BehavData,2) > 17
%         BehavData.BehavData5 = [];
%         BehavData.BehavData6 = [];
%     end
    BehavData.Properties.VariableNames = titles;
    
%     if exist('behav_struct_risky') == 1
%         BehavDataRisky.Properties.VariableNames = titles;
%         BehavDataRisky.Properties.VariableNames = titles;
%         BehavDataNOTRisky.Properties.VariableNames = titles;
%         BehavDataNOTRisky.Properties.VariableNames = titles;
%     end
    

    %Trim size of BehavData to match zall_struct (sometimes BehavData has a
    %few extra fake trials)
%     keep = size(zall_struct.data{1, (ii)}, 1);
%     BehavData = BehavData(1:keep, :);
%     for kk = 1:size(BehavData, 1)
%         if BehavData.stTime(kk) == BehavData.choiceTime(kk);
%             BehavData(kk,:)=[];
%         end
%     end
    % Changed to BehavData.stTime(kk) - BehavData.choiceTime(kk) == 0; from
    % BehavData.stTime(kk) & BehavData.choiceTime(kk) == 0; to account for
    % incidences where ABET adds an extra row at the end of the session.
    % Did this 1/20/2022
    for kk = 1:size(BehavData, 1)
        if BehavData.stTime(kk) && BehavData.choiceTime(kk) == 0;
            BehavData(kk,:)=[];
        end
    end
    [BehavData, trials, varargin] = TrialFilter(BehavData, 'ALL', 1);
    
    behav_struct_filtered(ii).data = BehavData;
    
     
    if isempty(BehavData)
        disp(['Animal has none of the specified trials'])
        ZallMean(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2));
        ZallMean_motion(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2));
   

        zerror(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2));
        zerror_motion(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2));

%         xcf_mean(ii,:) = NaN(1, size)

        choice_lat(ii) = NaN;
        collect_lat(ii) = NaN;
        struct_choice_lat(ii).data = NaN;
        struct_collect_lat(ii).data = NaN;

    continue
    else
    for jj = 1:size(trials,1)
        zall_filtered(ii).data(jj,:)= zall_struct.data{1, (ii)}(trials{jj},:);
    
        
        
        if ~isempty(zall_motion_struct.data{1, (ii)})
            %changed on 7/14/2022
             %to make zall_filtered and zall_motion_filtered equal, take
             %only the length of zall_struct
            zall_motion_filtered(ii).data(jj,:)=zall_motion_struct.data{1, (ii)}(trials{jj},1:size(zall_struct.data{1,(ii)},2));
            
            % get crosscorrelations for the window around the 0 point for
            % all trials for all mice (will end up as NaNs if mice do not
            % have video/movement data

            
            % Victoria recommends doing a mean subtracted calculation, but
            % why? 
            [xcf, lags] = xcorr(zall_motion_filtered(ii).data(jj, xcor_ind),zall_filtered(ii).data(jj, xcor_ind), 'coeff'); % numLags=150;
%             [xcf, lags] = xcorr((zall_motion_filtered(ii).data(jj, xcor_ind) - mean(zall_motion_filtered(ii).data(jj, xcor_ind))) , (zall_filtered(ii).data(jj, xcor_ind) - mean(zall_filtered(ii).data(jj, xcor_ind))));
%             
            xcf_struct(ii).data(jj,:) = xcf;
            %each "lag" corresponds to a time bin of the xcor_ind window,
            %with each bin being equal to 1 sample from ts1 (e.g., 0.0333)
            lags_struct(ii).data(jj,:) = lags;
            
            
            %because some mice don't have videos (motion data), need to fill
        %these zall_motion_filtered cells with NaNs
        elseif isempty(zall_motion_struct.data{1, (ii)})
            zall_motion_filtered(ii).data(jj,1:size(ts1,2)) = NaN;

            
        end
        jj = jj+1;
    end
    jj = 1;

    


    if size(zall_filtered(ii).data, 1) == 1
%         disp(['Animal only has one trial, ignoring']) % comment out if you want to omit data if animal only has 1 trial of a given type
%         ZallMean(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2)); % comment out if you want to omit data if animal only has 1 trial of a given type
%         zerror(ii,:) = NaN(1, size(zall_struct.data{1, (ii)},2)); % comment out if you want to omit data if animal only has 1 trial of a given type
        ZallMean(ii,:) = zall_filtered(ii).data(jj,:);
        
        zerror(ii,:) = zeros(1, size(zall_struct.data{1, (ii)},2));
        
        ZallMean_motion(ii,:) = zall_motion_filtered(ii).data(jj,:);
        zerror_motion(ii,:) = zeros(1, size(zall_struct.data{1, (ii)},2));
        


        num_trials(ii) = size(behav_struct_filtered(ii).data,1);
        collect_mice(ii,:) = [corrTable.animalNames(ii), corrTable.behavFiles(ii)];
        
        %need to add these otherwise it misses out on animals with 1 trial
        choice_lat(ii) = (behav_struct_filtered([ii]).data.choiceTime - behav_struct_filtered([ii]).data.stTime)';
        collect_lat(ii) = (behav_struct_filtered([ii]).data.collectionTime - behav_struct_filtered([ii]).data.choiceTime)';
        struct_choice_lat(ii).data = (behav_struct_filtered([ii]).data.choiceTime - behav_struct_filtered([ii]).data.stTime)';
        mean_choice_lat(ii) = struct_choice_lat(ii).data(jj,:);
        struct_collect_lat(ii).data = (behav_struct_filtered([ii]).data.collectionTime - behav_struct_filtered([ii]).data.choiceTime)';
        mean_collect_lat(ii) = struct_collect_lat(ii).data(jj,:);


    elseif size(zall_filtered(ii).data,1) ~= 1
        ZallMean(ii,:) = mean(zall_filtered(ii).data);
        ZallMean_motion(ii,:) = mean(zall_motion_filtered(ii).data);
        zerror(ii,:) = std(zall_filtered(ii).data)/size(zall_filtered(ii).data,1);
        zerror_motion(ii,:) = std(zall_motion_filtered(ii).data)/size(zall_motion_filtered(ii).data,1);



        choice_lat(ii) = mean(behav_struct_filtered([ii]).data.choiceTime - behav_struct_filtered([ii]).data.stTime)';
        collect_lat(ii) = mean(behav_struct_filtered([ii]).data.collectionTime - behav_struct_filtered([ii]).data.choiceTime)';
        struct_choice_lat(ii).data = (behav_struct_filtered([ii]).data.choiceTime - behav_struct_filtered([ii]).data.stTime)';
        mean_choice_lat(ii) = mean(struct_choice_lat(ii).data);
        struct_collect_lat(ii).data = (behav_struct_filtered([ii]).data.collectionTime - behav_struct_filtered([ii]).data.choiceTime)';
        mean_collect_lat(ii) = mean(struct_collect_lat(ii).data);
        num_trials(ii) = size(behav_struct_filtered(ii).data,1); 
        collect_mice(ii,:) = [corrTable.animalNames(ii), corrTable.behavFiles(ii)];
    end
%     clear zall_filtered;
    
%     collect_times(ii,:) = (behav_struct_filtered([ii]).data.collectionTime) - (behav_struct_filtered([ii]).data.stTime);
%     choice_times(ii,:) = (behav_struct_filtered([ii]).data.choiceTime) - (behav_struct_filtered([ii]).data.stTime);
     
    ii = ii+1;
    jj = 1;

    end
end    

[~,cols]=size(ZallMean);

iter = iter+1;

for hh=1:cols
    zAvg(hh)=nanmean(ZallMean(:,hh));
    zAvg_motion(hh)=nanmean(ZallMean_motion(:,hh));
%     values=ZallMean(:,hh);
%     values_ZallMean_motion = ZallMean_motion(:,hh);
%     stdDev(hh)=nanstd(values);
%     stdDev_motion(hh)=nanstd(values_ZallMean_motion);
end

varargin_array(iter,:) = varargin;
ZallMean_no_nan = ZallMean((all((~isnan(ZallMean)),2)),:);
ZallMean_for_perm_test{iter} = ZallMean_no_nan;
% ZallMean_for_perm_test{iter} = ZallMean;
% SEM_real(iter,:)= nanstd(ZallMean,1)/(sqrt(size(ZallMean_no_nan, 1)));
SEM_real = nanstd(ZallMean,1)/(sqrt(size(ZallMean_no_nan, 1)));

ZallMean_motion_no_nan = ZallMean_motion((all((~isnan(ZallMean_motion)),2)),:);
ZallMean_motion_for_perm_test{iter} = ZallMean_motion_no_nan;
% ZallMean_motion_for_perm_test{iter} = ZallMean_motion;

% SEM__motion_real(iter,:)= nanstd(ZallMean_motion,1)/(sqrt(size(ZallMean_motion_no_nan, 1)));
SEM_motion_real= nanstd(ZallMean_motion,1)/(sqrt(size(ZallMean_motion_no_nan, 1)));


num_trials_sum{iter} = sum(num_trials);
num_sessions_sum{iter} = numel(ZallMean_no_nan(:,1));
disp(collect_mice);
disp(num_trials_sum);
% 
% SEM=stdDev./(sqrt(size(ZallMean, 1)));
% lo=zAvg-SEM;
% hi=zAvg+SEM;

% SEM_motion=stdDev./(sqrt(size(ZallMean_motion, 1)));
% lo_motion=zAvg_motion-SEM_motion;
% hi_motion=zAvg_motion+SEM_motion;


ind = ts1(1,:) > -10 & ts1(1,:) < 10;
ts1forprism = ts1(:,ind);
zAvgforprism = nanmean(ZallMean(:,ind));
zAvgforprism_motion = nanmean(ZallMean_motion(:,ind));
semforprism_motion = SEM_motion_real(:,ind);
sem_real_for_prism = SEM_real(:,ind);





%change index values for different AUC calculations as necessary
AUC=[]; % cue, shock
for qq = 1:size(ZallMean,1);
    AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -5)); % -0 -2 %proxy for pre-choice
    AUC(qq,2)=trapz(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 3)); % 0 2 %proxy for post-choice
    AUC(qq,3)=trapz(ZallMean(qq,ts1(1,:) > 4 & ts1(1,:) < 8)); % 0 2 %proxy for collectionTime
    qq=qq+1;
end
AUC_mean(:,iter) = nanmean(AUC,1);
% AUC_sem(:,iter) = nanstd(AUC,1)/sqrt(numel(AUC(:,1)));
%changed 8/30/2022 to reflect that numel includes NaN values, which it
%shouldnt
AUC_sem(:,iter) = nanstd(AUC,1)/sqrt(sum(~isnan(AUC(:,1))));
AUC_SD(:,iter) = nanstd(AUC,1);




%change index values for different AUC calculations as necessary
AUC_motion=[]; % cue, shock
for qq = 1:size(ZallMean_motion_no_nan,1);
    AUC_motion(qq,1)=trapz(ZallMean_motion_no_nan(qq,ts1(1,:) < 0 & ts1(1,:) > -4)); % -0 -2 %proxy for pre-choice
    AUC_motion(qq,2)=trapz(ZallMean_motion_no_nan(qq,ts1(1,:) > 0 & ts1(1,:) < 4)); % 0 2 %proxy for post-choice
    AUC_motion(qq,3)=trapz(ZallMean_motion_no_nan(qq,ts1(1,:) > 4 & ts1(1,:) < 8)); % 0 2 %proxy for collectionTime
    qq=qq+1;
end

AUC_motion_iters(iter) = {AUC_motion};

AUC_iters(iter) = {AUC};



AUC_mean_motion(:,iter) = nanmean(AUC_motion,1);
% AUC_sem_motion(:,iter) = nanstd(AUC_motion)/sqrt(numel(AUC_motion(:,1)));
%changed 8/30/2022 to reflect that numel includes NaN values, which it
%shouldnt
%need to have this if statement because for unfrequent events, sometimes
%there can be only 1 trial, which mangles the SEM/SD calculation. if there
%is only one trial - that trial has 0 mean or SD
if size(AUC_motion, 1) <=1
    AUC_sem_motion(:,iter) = zeros(1, size(AUC_motion,2));
    AUC_SD_motion(:,iter) = zeros(1, size(AUC_motion,2));
elseif size(AUC_motion, 1) > 1
    AUC_sem_motion(:,iter) = nanstd(AUC_motion,1)/sqrt(sum(~isnan(AUC_motion(:,1))));
    AUC_SD_motion(:,iter) = nanstd(AUC_motion,1);
end

% AUC_pasted = [];
% choice_lat_pasted = [];
% for ll = 1:size(zall_filtered, 2)
%     for nn = 1:size(zall_filtered(ll).data, 1);
%         struct_AUC(ll).data(nn,1) = trapz(zall_filtered(ll).data(nn,ts1(1,:) < 0 & ts1(1,:) > -5));
%         struct_AUC(ll).data(nn,2) = trapz(zall_filtered(ll).data(nn,ts1(1,:) > 0 & ts1(1,:) < 4));
%         AUC_pasted = [AUC_pasted; (struct_AUC(ll).data(nn,1))];
%         choice_lat_pasted = [choice_lat_pasted; (struct_choice_lat(ll).data(nn))];
% %         AUC_pasted_2 = [AUC_pasted_2; {struct_AUC(ll).data(nn,1)}];
%     end
% end



arg_string = string(varargin);

arg_string_combine = join(arg_string, '=');

% Get labels for each comparison based on the events that the data are filtered on


for ii = 1:size(arg_string, 2)
    for qq = 1:size(arg_string, 1)
        if arg_string(qq,ii) == 'BLOCK';
            arg_string_block(qq,:) = join([arg_string(qq,ii) arg_string(qq, ii+1)], '=');
        end
    end
end

cc = 1;

for ii = 1:size(varargin_array, 2)-1
    if isstring(varargin_array(qq,ii))
        arg_string_other = join([arg_string(:,ii) arg_string(:, ii+1)], '=')
        cc = cc+1;
    elseif ~isstring(varargin_array(qq,ii))
        arg_string_other = join([arg_string(:,ii) arg_string(:, ii+1)], '=')
        cc = cc+1;
    end
end

num_trials_sum_string = string(num_trials_sum)';
num_sessions_sum_strong = string(num_sessions_sum)';
cmb_strings = append(arg_string_other, ' ','num events=', num_trials_sum_string, ' ', 'num sessions=', num_sessions_sum_strong);





% hold on;
% figure(1)
% hold on



% Plot vertical line at epoch onset, time = 0

%for loop for 5 tone retrieval- retrieval is 50 tones though??

% a=0;
% for ii=1:50
%     p1=patch([a a+30 a+30 a], [-5 -5 20 20], [.8 1 1], 'EdgeColor','none');
%     a=a+35;
% end 

% vline = plot([0 0],[min(lo) max(hi)],'--','LineWidth',2)

% y = ylim; % current y-axis limits
% vline = plot([0 0],[min(lo) max(hi)]);
% vline.Color = 'k';
% vline.LineStyle = '--';
% vline.LineWidth = 2;

% xxx=[ts1, ts1(end:-1:1)];
% yyy=[lo, hi(end:-1:1)];
% hp= fill(xxx,yyy,[ .8 .8 .8]); hold on;
% set(hp,'FaceColor', 'k','EdgeColor','none','FaceAlpha',0.10,'EdgeAlpha',0.10);hold on;
% p3=plot(ts1, zAvg, 'k', 'LineWidth', 2,'DisplayName','z-score'); hold on;
% xlim([-8 8]); hold on;



% xlim([-30 5]); hold on;
% xlim([-4 8]); hold on;
% ylim([-4 5]);

%from scientific color maps: Crameri, F. (2018), Scientific colour-maps, Zenodo, doi:10.5281/zenodo.1243862
load('batlow.mat');

% change color of lines used in plots based on # of iterations of analysis
if iter == 1
    color_iter = iter;
elseif iter == 2
    color_iter = iter+98;
elseif iter == 3
    color_iter = iter+197;
elseif iter > 3
    color_iter = iter+20;
end

hold on;
figure(1); 
% shadedErrorBar(ts1, zAvg_motion(1,:), SEM_motion_real(1,:),'lineProps', 'r');
plot(ts1, zAvg_motion(1,:),'color', batlow(color_iter,:),'DisplayName',arg_string_combine);
errorplot3(zAvg_motion(1,:)-SEM_motion_real(1,:),zAvg_motion(1,:)+SEM_motion_real(1,:),[-33 10],batlow(color_iter,:),.15)
% shadedErrorBar(ts1, zAvg_motion(1,:), SEM_motion_real(1,:),'lineProps', {'color', batlow(color_iter,:)});
ylabel('z-scored velocity', 'FontSize', 12);
xlabel('Time from choice (s)');
xlim([-8 8]); 
% legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')
legend('-DynamicLegend')


hold off;

hold on;
figure(2); 
% shadedErrorBar(ts1, zAvg(1,:), SEM_real(1,:), 'lineProps', {'color', [255, 83, 26]/255});
shadedErrorBar(ts1, zAvg(1,:), SEM_real(1,:), 'lineProps', {'color', batlow(color_iter,:)});
ylabel('z-scored dF/F', 'FontSize', 12);
xlabel('Time from choice (s)');
xlim([-8 8]);
legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')



hold off;

hold on;
figure(3);
bar([iter], AUC_mean(1,iter));
hold on;
er = errorbar([iter],AUC_mean(1,iter),AUC_sem(1,iter),AUC_sem(1,iter));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

hold off;

hold on;
figure(4);
bar([iter], AUC_mean(2,iter));
hold on;
er = errorbar([iter],AUC_mean(2,iter),AUC_sem(2,iter),AUC_sem(2,iter));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off





% xxx=[ts1, ts1(end:-1:1)];
% yyy_motion=[lo_motion, hi_motion(end:-1:1)];
% hp= fill(xxx,yyy_motion,[ .8 .8 .8]); hold on;
% set(hp,'FaceColor', 'k','EdgeColor','none','FaceAlpha',0.10,'EdgeAlpha',0.10);hold on;
% p3=plot(ts1, zAvg_motion, 'k', 'LineWidth', 2,'DisplayName','z-score'); hold on;
% xlim([-8 8]); hold on;


% Make a legend
% legend([p3 hp],{'Block 1','SEM'}, 'AutoUpdate', 'off');

% set(gca,'Layer','top')



% loop through xcf struct, created earlier by doing the crosscorr between motion & photometry signal for each
% trial, for each mouse. will contain NaNs because some mice have no
% video/motion data
for i = 1:length(xcf_struct)
    if length(xcf_struct(i).data) > 0
        xcf_means(i,:) = nanmean(xcf_struct(i).data);
        xcf_error(i,:) = nanstd(xcf_struct(i).data)/size(xcf_struct(i).data,1);
    elseif length(xcf_struct(i).data) == 0
        xcf_means(i,:) = NaN(1,length(lags));
        xcf_error(i,:) = NaN(1,length(lags));
    end
end

% get NaNmean of crosscorrelation array. 
% iter refers to how many times the code is run - for example, if you want
% to do blocks 1, 2, 3, run code 3 times

xcf_means_no_nan = xcf_means((all((~isnan(xcf_means)),2)),:);

xcf_array_mean(iter,:) = nanmean(xcf_means_no_nan, 1);
% xcf_array_avg_SEM(iter,:)= nanmean(xcf_error,1);

xcf_array_avg_SEM(iter,:)= nanstd(xcf_means,1)/(sqrt(size(xcf_means_no_nan, 1)));


%%
%Heatmap each animal
figure;
% start=[];fin=[];
load('batlowW.mat');

for an = 1:size(ZallMean_no_nan,1)
%     start=30*(an-1)+1;fin=start+29;
    subplot(ceil(size(ZallMean_no_nan,1)/2),2,an);
    imagesc(ts1, 1, zall_filtered(an).data(:,:));
    colormap(batlowW); % c1 = colorbar;
    title(sprintf(corrTable.animalNames{an,1}),'FontSize',7);
    xlim([-10 10]);
    
end
sgtitle('Trial Based Heat Map (Z-Scored Ca)') 

%heatmap of motion for each mouse
figure;
for an = 1:size(ZallMean_no_nan,1)
%     start=30*(an-1)+1;fin=start+29;
    subplot(ceil(size(ZallMean_no_nan,1)/2),2,an);
    imagesc(ts1, 1, zall_motion_filtered(an).data(:,:));
    colormap(batlowW); % c1 = colorbar;
    title(sprintf(corrTable.animalNames{an,1}),'FontSize',7);
    xlim([-10 10]);
    
end
sgtitle('Trial Based Heat Map (Z-Scored Motion)') 

%line graph of individual mean+SEM for each mouse
figure;
for an = 1:size(ZallMean,1)
%     start=30*(an-1)+1;fin=start+29;
    subplot(ceil(size(ZallMean,1)/2),2,an);
    shadedErrorBar(ts1, ZallMean(an,:), zerror(an,:),'lineProps', {'color', batlow(color_iter,:)});
%     ylabel('z-scored Ca', 'FontSize', 12);
%     xlabel('Time from choice (s)');
    colormap(batlowW); % c1 = colorbar;
    title(sprintf(corrTable.animalNames{an,1}),'FontSize',7);
    xlim([-10 10]);
    
end
sgtitle('Avg Z-Scored Ca') 

%line graph of individual mean+SEM for each mouse
figure;
for an = 1:size(ZallMean,1)
%     start=30*(an-1)+1;fin=start+29;
    subplot(ceil(size(ZallMean,1)/2),2,an);
    shadedErrorBar(ts1, ZallMean_motion(an,:), zerror_motion(an,:),'lineProps', {'color', batlow(color_iter,:)});
%     ylabel('z-scored Ca', 'FontSize', 12);
%     xlabel('Time from choice (s)');
    colormap(batlowW); % c1 = colorbar;
    title(sprintf(corrTable.animalNames{an,1}),'FontSize',7);
    xlim([-10 10]);
    
end
sgtitle('Avg Z-Scored Motion') 


%heatmap Y axis each individual mouse Calcium
figure;
imagesc(ts1, 1, ZallMean_no_nan);
sgtitle('Mean Ca Each Mouse') 

%heatmap Y axis each individual mouse Motion
figure;
imagesc(ts1, 1, ZallMean_motion_no_nan);
sgtitle('Mean Motion Each Mouse') 

%%
figure; 
shadedErrorBar(lags, xcf_array_mean(1,:), xcf_array_avg_SEM(1,:));
hold on;
shadedErrorBar(lags, xcf_array_mean(2,:), xcf_array_avg_SEM(2,:),'lineProps', '-b');
hold on;
shadedErrorBar(lags, xcf_array_mean(3,:), xcf_array_avg_SEM(3,:),'lineProps', '-r');
legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')

%%
% would be useful to get AUC from the entire pre-choice period for each
% trial, so if a trial takes 15.5 s, take entire period AUC
%%
% would be useful to get AUC from the entire pre-choice period for each
% trial, so if a trial takes 15.5 s, take entire period AUC
AUC_pasted = [];
choice_lat_pasted = [];
collect_lat_pasted = [];
zall_pasted = [];
corrcoef_struct = []
for ll = 1:size(zall_filtered, 2)
    for nn = 1:size(zall_filtered(ll).data, 1)
        struct_AUC(ll).data(nn,1) = trapz(zall_filtered(ll).data(nn,ts1(1,:) < 0 & ts1(1,:) > -5));
        struct_AUC(ll).data(nn,2) = trapz(zall_filtered(ll).data(nn,ts1(1,:) > 0 & ts1(1,:) < 4));
        %this should be AUC for the entire pre-choice period
        struct_AUC(ll).data(nn,3) = trapz(zall_filtered(ll).data(nn,ts1(1,:) < 0 & ts1(1,:) > behav_struct_filtered(ll).data.stTime(nn) - behav_struct_filtered(ll).data.choiceTime(nn)));
        AUC_pasted = [AUC_pasted; (struct_AUC(ll).data(nn,1)),(struct_AUC(ll).data(nn,2)),(struct_AUC(ll).data(nn,3))];
        choice_lat_pasted = [choice_lat_pasted; (struct_choice_lat(ll).data(nn))];
        collect_lat_pasted = [collect_lat_pasted; (struct_collect_lat(ll).data(nn))];
        zall_pasted = [zall_pasted; (zall_filtered(ll).data(nn,:))];
        %get correlations for trial from trial start till the end of the
        %chopped up trace
        corrcoef_mouse = corrcoef(zall_filtered(ll).data(nn, ts1(1,:) > behav_struct_filtered(ll).data.stTime(nn) - behav_struct_filtered(ll).data.choiceTime(nn)), zall_motion_filtered(ll).data(nn, ts1(1,:) > behav_struct_filtered(ll).data.stTime(nn) - behav_struct_filtered(ll).data.choiceTime(nn)));
        corrcoef_struct(ll).data(nn,1) = corrcoef_mouse(1,2); 
        
%         AUC_pasted_2 = [AUC_pasted_2; {struct_AUC(ll).data(nn,1)}];
    end
end

for pp = 1:size(corrcoef_struct, 2)
    corrcoef_mean_mouse(pp) = nanmean(corrcoef_struct(pp).data);
end



%%
c = distinguishable_colors(numel(corrTable.animalNames));

figure;
set(gca, 'ColorOrder', c, 'NextPlot', 'replacechildren');

plot(ts1, ZallMean, 'LineWidth',2);

ID_legend = corrTable.animalNames;
legend(ID_legend);


%%
%%sum aborts for comparison in Excel or GraphPad
for ii = 1:size(behav_struct,2)
    aa_sum(ii,1) = sum(behav_struct(ii).data{1, 20} == 1, 'omitnan');
    aa_sum(ii,2) = sum(behav_struct(ii).data{1, 20} == 2, 'omitnan');
end

names = corrTable.animalNames;
sessions = corrTable.session;
AA_table = table(names, sessions, aa_sum);

AA_table_sorted = sortrows(AA_table, 'sessions');

%%  Calculate the Z-scored data for RISKY or NON-RISKY mice (defined above)
%   Works best when only 1 session is included (see code to select sessions
%   above
clear ZallMeanRisky zerrorRisky ZallMeanNOTRisky zerrorNOTRisky MeanZallRisky meanzerrorRisky MeanZallNOTRisky MeanzerrorNOTRisky

qq = 1; 
zz = 1;
for ii = 1:size(corrTable,1)
    if corrTable.RiskSeeking(ii) == 1
        ZallMeanRisky(qq,:) = ZallMean(ii,:);
        zerrorRisky(qq,:) = zerror(ii,:);

        qq = qq+1;
    elseif corrTable.RiskSeeking(ii) == 2
        ZallMeanNOTRisky(zz,:) = ZallMean(ii,:);
        zerrorNOTRisky(zz,:) = zerror(ii,:);

        zz=zz+1;
    end

end
% 
        MeanZallRisky = nanmean(ZallMeanRisky,1);
        MeanzerrorRisky = nanmean(zerrorRisky,1);
        MeanZallNOTRisky = nanmean(ZallMeanNOTRisky,1);
        MeanzerrorNOTRisky = nanmean(zerrorNOTRisky,1);
        
        
figure;
% vline = plot([0 0],[min(lo) max(hi)]);
% vline.Color = 'k';
% vline.LineStyle = '--';
% vline.LineWidth = 2;

% xxx=[ts1, ts1(end:-1:1)];
% yyy=[lo, hi(end:-1:1)];
% hp= fill(xxx,yyy,[ .8 .8 .8]); hold on;
set(hp,'FaceColor', [ .8 .8 .8],'EdgeColor','none','FaceAlpha',0.8,'EdgeAlpha',0.8);hold on;
p3=plot(ts1, MeanZallNOTRisky, 'k', 'LineWidth', 2,'DisplayName','z-score'); hold on;
plot(ts1, MeanZallRisky, 'b', 'LineWidth', 2,'DisplayName','z-score');


% Make a legend
% legend([p3 hp],{'Block 1','SEM'}, 'AutoUpdate', 'off');

set(gca,'Layer','top')
%%

% Plot heat map for a given variable (change _filtered(#) etc to correspond
% to variable / mouse of interest
% subplot(4,1,2)

animalNum_var = 18;

[numTrials,~]=size(behav_struct_filtered(animalNum_var).data.collectionTime(:));
Tris=[1:numTrials]';

time2Collect = behav_struct_filtered(animalNum_var).data.collectionTime(:)-behav_struct_filtered(animalNum_var).data.choiceTime(:);


figure;
IM=imagesc(ts1, 1, zall_filtered(animalNum_var).data); hold on;

load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 


scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
%scatter(time2choose,Tris,'Marker','>','MarkerFaceColor','k')
% title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 12);
hold off;


figure;
IM=imagesc(ts1, 1, zall_motion_filtered(animalNum_var).data); hold on;

load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 


scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
%scatter(time2choose,Tris,'Marker','>','MarkerFaceColor','k')
% title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 12);
hold off;


figure; shadedErrorBar(ts1, ZallMean_motion(animalNum_var,:), zerror_motion(animalNum_var,:),'lineProps', 'k'); hold on; shadedErrorBar(ts1, ZallMean(animalNum_var, :), zerror(animalNum_var,:), 'lineProps', {'color', [255, 83, 26]/255});

%% Calculate AUC for different time windows, store in AUC array
% AUC=[]; % cue, shock
% for qq = 1:size(ZallMean,1);
%     AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
%     AUC(qq,2)=trapz(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 5));
%     qq=qq+1;
% end
% AUC_mean(:,iter) = nanmean(AUC);
% AUC_std(:,iter) = nanstd(AUC)/sqrt(numel(AUC(:,1)));
%%
means=[]; % cue, shock
for qq = 1:size(ZallMean,1);
    means(qq,1)=mean(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    means(qq,2)=mean(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 5));
    qq=qq+1;
end
means_mean = mean(means);
means_std = std(means)/sqrt(numel(means(:,1)));

%%
%AUC correlation w/ choice latency
figure;
%create a scatterplot w/ variables of interest
scatter(choice_lat,AUC(:,1));

%create idx to omit NaN values for idxs to be correlated
idx_choice_lat = isnan(choice_lat);
idx_AUC = isnan(AUC);

% get line of best fit while omitting NaN values using idx
coefficients = polyfit(choice_lat(~idx_choice_lat), AUC(~idx_AUC(:,1),1)', 1); %added a AUC(:,1)' here to make the arrays the same orientation?
xFit = linspace(0, max(choice_lat), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Choice Latency % vs Pre-Choice -5 to 0 AUC')
xlabel('Choice Latency')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;
choice_lat_AUC_corr = corrcoef(choice_lat, AUC(:,1),'rows','complete'); %correlate choice latency w/ pre-choice AUC



%AUC for entire pre-choice period correlated w/ choice latency
figure;
scatter(choice_lat_pasted,AUC_pasted(:,3));
% get line of best fit while omitting NaN values using idx
coefficients = polyfit(choice_lat_pasted, AUC_pasted(:,3), 1); %added a AUC(:,1)' here to make the arrays the same orientation?
xFit = linspace(0, max(choice_lat_pasted), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Choice Latency % vs Pre-Choice -5 to 0 AUC')
xlabel('Choice Latency')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;
choice_lat_AUC_corr_pasted = corrcoef(choice_lat_pasted, AUC_pasted(:,3),'rows','complete'); %correlate choice latency w/ pre-choice AUC


%% Separate slow and fast choice trials

%create index for choice latencies that are longer than 4s
slow_choice_lat_ind = choice_lat_pasted > 5;


%plot slow vs fast
figure;
plot(ts1, mean(zall_pasted(slow_choice_lat_ind,:)));
hold on; plot(ts1, mean(zall_pasted(~slow_choice_lat_ind,:)));

%%
%if a given animal has 0 trials, will need to delete rows to deal with
%mismatch
figure;

scatter(corrTable.RiskPercent_91Trials,AUC(:,1));
coefficients = polyfit(corrTable.RiskPercent_91Trials, AUC(:,1), 1);
xFit = linspace(0, max(corrTable.RiskPercent_91Trials), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
% fontsize(gca,18,"pixels")
title('Scatter: Risk % vs Pre-Choice')
xlabel('Proportion Risky Choices')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
% grid on;

pre_choice_riskiness_mdl = fitlm(corrTable.RiskPercent_91Trials, AUC(:,1));

figure;
% fontsize(gca,100,"pixels")
scatter(corrTable.RiskPercent_91Trials,AUC(:,2));
coefficients = polyfit(corrTable.RiskPercent_91Trials, AUC(:,2), 1);
xFit = linspace(0, max(corrTable.RiskPercent_91Trials), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
% fontsize(gca,18,"pixels")
title('Scatter: Risk % vs Post-Choice')
xlabel('Proportion Risky Choices')
ylabel('Normalized Post-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
% grid on;

post_choice_riskiness_mdl = fitlm(corrTable.RiskPercent_91Trials,AUC(:,2));

figure;
% fontsize(gca,18,"pixels")
scatter(corrTable.LoseShiftPercent,AUC(:,1));
coefficients = polyfit(corrTable.LoseShiftPercent, AUC(:,1), 1);
xFit = linspace(0, max(corrTable.LoseShiftPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
% fontsize(gca,18,"pixels")
title('Scatter: LoseShift vs Pre-Choice')
xlabel('Proportion of Lose Shifts')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
% grid on;

%%
figure;
hold on
C = [.8 .2 .8;
     .2 .75 .2];
CE = [.5 .1 .5];
superbar(AUC_mean(1,:), 'BarFaceColor', C, 'BarWidth', .5, 'E', AUC_std(1,:), 'ErrorbarColor', CE);
xlim([0.5 3.5]);

% title('superbar');

%% IF USING ALL SESSIONS + Need to double check that any mice with 1 session only are removed

LoseShift_D(:,1) = corrTable.LoseShiftPercent(1:2:end);
LoseShift_D(:,2) = corrTable.LoseShiftPercent(2:2:end);
WinStay_D(:,1) = corrTable.WinStayPercent(1:2:end);
WinStay_D(:,2) = corrTable.WinStayPercent(2:2:end);
RiskyPercent_D(:,1) = corrTable.RiskPercent(1:2:end);
RiskyPercent_D(:,2) = corrTable.RiskPercent(2:2:end);
meanRiskyPercent_D = mean(RiskyPercent_D);
SEMRiskyPercent_D = std(RiskyPercent_D)/sqrt(numel(RiskyPercent_D(:,1)));

%%
figure;
C = [.8 .2 .8;
     .2 .75 .2];
CE = [.5 .1 .5];
superbar(meanRiskyPercent_D, 'BarFaceColor', C, 'BarWidth', .5, 'E', SEMRiskyPercent_D, 'ErrorbarColor', CE);
xlim([0.5 2.5]);
% title('superbar');

%%
num_animals = size(zAvgforprism,1);

rc = num_animals/5;

sp_rows = floor(rc);
sp_cols = ceil(rc);

for i = 1:num_animals
    axi = subplot(sp_rows,sp_cols,i);
    plot(ts1forprism,zAvgforprism(i,:));
    title(ID_legend{i});
end

%%

for qq = 1:size(zall_filtered, 1)
    for pp = 1:size(zall_filtered(qq).data, 1)
        BL_shifted(pp,qq)=[behav_struct_filtered(qq).data.stTime(pp) - behav_struct_filtered(qq).data.choiceTime(pp)]; 
        peaks(pp, qq) = findchangepts(zall_filtered(qq).data(pp,ts1(1,:) > BL_shifted(pp) & ts1(1,:) < 0));    
        peaks_norm(pp, qq) = peaks(pp, qq) / abs(BL_shifted(pp,qq));
    end
    clear pp
end
%     ind = ts1(1,:) < BL_shifted(2) & ts1(1,:) > BL_shifted(1);
%     findchangepts(zall_filtered(qq).data(:,ts1(1,:) > 0 & ts1(1,:) < 2))
%%
for qq = 1:size(ZallMean,1)
    premean(qq,1)=mean(ZallMean(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    test_change(qq,1)=findchangepts(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -2));
    qq=qq+1;
end
increasing_means = ZallMean([1:6,13:14,19:20],:);
decreasing_means = ZallMean([7:12,15:18,21:26],:);

%% start to think about how to do cross correlations w/signal & velocity around choice time
% https://www.mathworks.com/help/econ/crosscorr.html#btzkul8-3

% ts1_lags = xcor_range(1):0.0333:(xcor_range(2)+5);


sum_xcor_window = abs(xcor_range(1))+abs(xcor_range(2));

ts1_lags = -(sum_xcor_window):0.033333333:(sum_xcor_window);
ts1_lags_correct = ts1_lags(2:end-1);



% % loop through xcf struct, created earlier by doing the crosscorr between motion & photometry signal for each
% % trial, for each mouse. will contain NaNs because some mice have no
% % video/motion data
% for i = 1:length(xcf_struct)
%     if length(xcf_struct(i).data) > 0
%         xcf_means(i,:) = nanmean(xcf_struct(i).data);
%         xcf_error(i,:) = nanstd(xcf_struct(i).data)/size(xcf_struct(i).data,1);
%     elseif length(xcf_struct(i).data) == 0
%         xcf_means(i,:) = NaN(1,length(lags));
%         xcf_error(i,:) = NaN(1,length(lags));
%     end
% end
% 
% % get NaNmean of crosscorrelation array. 
% % iter refers to how many times the code is run - for example, if you want
% % to do blocks 1, 2, 3, run code 3 times
% 
% xcf_means_no_nan = xcf_means((all((~isnan(xcf_means)),2)),:);
% 
% xcf_array_mean(iter,:) = nanmean(xcf_means_no_nan, 1);
% % xcf_array_avg_SEM(iter,:)= nanmean(xcf_error,1);
% 
% xcf_array_avg_SEM(iter,:)= nanstd(xcf_means,1)/(sqrt(size(xcf_means_no_nan, 1)));
% 




figure; 
shadedErrorBar(ts1_lags_correct, xcf_array_mean(1,:), xcf_array_avg_SEM(1,:));
hold on;
shadedErrorBar(ts1_lags_correct, xcf_array_mean(2,:), xcf_array_avg_SEM(2,:),'lineProps', '-b');
hold on;
shadedErrorBar(ts1_lags_correct, xcf_array_mean(3,:), xcf_array_avg_SEM(3,:),'lineProps', '-r');
legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')



