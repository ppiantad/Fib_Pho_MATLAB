%%
% to access BehavData from inside behav_struct for trial filtering


iter = 0;
% These .mat files contain 4 files processed using AllAnimalFPanalysis
% with the Analysis_BL2initiation_fn_beta.m function. The files are as
% follows:


% Uncomment the particular file below to load the .mat file (will create a list of inputs later in a function to tidy this up)  
% Can uncomment and load diff. sessions to graph data from separate groups,

load 'BLA-NAcSh_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_01032023_no_motion.mat'

% load 'BLA-PL_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_07112023_no_motion.mat'

% load 'vmPFC-NAcSh_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_01032023_no_motion.mat'

% load 'vHPC-NAcSh_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_01032023_no_motion.mat'

% load 'GFP_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_01032023_no_motion.mat'

% load 'D2_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_01032023_no_motion.mat'

% load 'D1_Shock_Test_EPOC_Trial_TRANGE_-5_10_zscored_01032023_no_motion.mat'

if size(behav_struct(1).data,2) == 4
    titles = {'Trial','shock','shockIntensity','choiceTime'};
elseif size(behav_struct(1).data,2) == 19
    titles = {'Trial','Block','ForceFree','bigSmall','RewSelection', 'TrialPossible', 'stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','smallRew','bigRew'};
elseif size(behav_struct(1).data,2) == 17
    titles = {'Trial','Block','ForceFree','bigSmall','stTime','choiceTime'...
    'collectionTime','shock','omission','omissionALL','WL','WSLScode','win_stay','lose_shift','lose_omit','smallRew','bigRew'};
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

only_motion_sessions = [1];



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


struct_AUC = [];
num_trials = [];


ts1 = cell2mat(ts1_struct.data);
N = 33.8822618125484; %this is the current downsample factor thanks to needing to downsample the motion and calcium data similarly.

%change this range to set the crosscorrelation window
xcor_range = [-5 5];
%set the number of lags equal to the range of the xcor window
numLags = abs(xcor_range(1)-xcor_range(2));
xcor_ind = ts1 > xcor_range(1) & ts1 < xcor_range(2);
xcor_lag_bins = sum(xcor_ind)/N;

for ii = 1:size(behav_struct,2)

    BehavData = cell2mat(behav_struct(ii).data);
    BehavData = array2table(BehavData);

    
    BehavData.Properties.VariableNames = titles;

    [BehavData, trials, varargin] = TrialFilter(BehavData, 'SHK',1);
    
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
        

    elseif size(zall_filtered(ii).data,1) ~= 1
        ZallMean(ii,:) = mean(zall_filtered(ii).data);
        ZallMean_motion(ii,:) = mean(zall_motion_filtered(ii).data);
        zerror(ii,:) = std(zall_filtered(ii).data)/size(zall_filtered(ii).data,1);
        zerror_motion(ii,:) = std(zall_motion_filtered(ii).data)/size(zall_motion_filtered(ii).data,1);

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
% SEM_real(iter,:)= nanstd(ZallMean,1)/(sqrt(size(ZallMean_no_nan, 1)));
SEM_real = nanstd(ZallMean,1)/(sqrt(size(ZallMean_no_nan, 1)));

ZallMean_motion_no_nan = ZallMean_motion((all((~isnan(ZallMean_motion)),2)),:);
ZallMean_motion_for_perm_test{iter} = ZallMean_motion;

% SEM__motion_real(iter,:)= nanstd(ZallMean_motion,1)/(sqrt(size(ZallMean_motion_no_nan, 1)));
SEM_motion_real= nanstd(ZallMean_motion,1)/(sqrt(size(ZallMean_motion_no_nan, 1)));


num_trials_sum = sum(num_trials);
disp(collect_mice);
disp(num_trials_sum);
% 
% SEM=stdDev./(sqrt(size(ZallMean, 1)));
% lo=zAvg-SEM;
% hi=zAvg+SEM;

% SEM_motion=stdDev./(sqrt(size(ZallMean_motion, 1)));
% lo_motion=zAvg_motion-SEM_motion;
% hi_motion=zAvg_motion+SEM_motion;


ind = ts1(1,:) > -5 & ts1(1,:) < 5;
ts1forprism = ts1(:,ind);
zAvgforprism = nanmean(ZallMean(:,ind));
zAvgforprism_motion = nanmean(ZallMean_motion(:,ind));
semforprism_motion = SEM_motion_real(:,ind);
sem_real_for_prism = SEM_real(:,ind);





%change index values for different AUC calculations as necessary
AUC=[]; % cue, shock
for qq = 1:size(ZallMean,1);
    AUC(qq,1)=trapz(ZallMean(qq,ts1(1,:) < 0 & ts1(1,:) > -4)); % -0 -2 %proxy for pre-choice
    AUC(qq,2)=trapz(ZallMean(qq,ts1(1,:) > 0 & ts1(1,:) < 4)); % 0 2 %proxy for post-choice
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
AUC_mean_motion(:,iter) = nanmean(AUC_motion,1);
% AUC_sem_motion(:,iter) = nanstd(AUC_motion)/sqrt(numel(AUC_motion(:,1)));
%changed 8/30/2022 to reflect that numel includes NaN values, which it
%shouldnt

if AUC_motion ~= []
    AUC_sem_motion(:,iter) = nanstd(AUC_motion,1)/sqrt(sum(~isnan(AUC_motion(:,1))));
    AUC_SD_motion(:,iter) = nanstd(AUC_motion,1);
elseif AUC_motion == []
    disp('No motion data exists')
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
shadedErrorBar(ts1, zAvg_motion(1,:), SEM_motion_real(1,:),'lineProps', {'color', batlow(color_iter,:)});
ylabel('z-scored velocity', 'FontSize', 12);
xlabel('Time from choice (s)');
xlim([-8 8]); 
legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')

hold off;

hold on;
figure(2); 
% shadedErrorBar(ts1, zAvg(1,:), SEM_real(1,:), 'lineProps', {'color', [255, 83, 26]/255});
shadedErrorBar(ts1, zAvg(1,:), SEM_real(1,:), 'lineProps', {'color', batlow(color_iter,:)});
ylabel('z-scored dF/F', 'FontSize', 12);
xlabel('Time from choice (s)');
xlim([-8 8]);
% legend('block 1 (0%)','block 2 (50%)', 'block 3 (75%)')
legend('BLA','vmPFC', 'vHPC', 'GFP')


hold off;





% loop through xcf struct, created earlier by doing the crosscorr between motion & photometry signal for each
% trial, for each mouse. will contain NaNs because some mice have no
% video/motion data

if exist('xcf_struct','var') == 1
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
elseif exist('xcf_struct','var') == 0
    disp('No crosscorr structure exists')
end
%%


for ii = 1:size(zall_filtered, 2)
    first_five{ii} = zall_filtered(ii).data(2:6,:);
        mean_first_five(ii,:) = mean(first_five{ii});
        sem_first_five(ii,:) = nanstd(first_five{ii})/sqrt(size(first_five{ii}, 1));
        
    second_five{ii} = zall_filtered(ii).data(7:11,:);
        mean_second_five(ii,:) = mean(second_five{ii});
        sem_second_five(ii,:) = nanstd(second_five{ii})/sqrt(size(second_five{ii}, 1));
    third_five{ii} = zall_filtered(ii).data(12:16,:);
        mean_third_five(ii,:) = mean(third_five{ii});
        sem_third_five(ii,:) = nanstd(third_five{ii})/sqrt(size(third_five{ii}, 1));
    fourth_five{ii} = zall_filtered(ii).data(17:21,:);
        mean_fourth_five(ii,:) = mean(fourth_five{ii});
        sem_fourth_five(ii,:) = nanstd(fourth_five{ii})/sqrt(size(fourth_five{ii}, 1));
    fifth_five{ii} = zall_filtered(ii).data(22:26,:);
        mean_fifth_five(ii,:) = mean(fifth_five{ii});
        sem_fifth_five(ii,:) = nanstd(fifth_five{ii})/sqrt(size(fifth_five{ii}, 1));

end

mean_first_five_Avg = mean(mean_first_five);
mean_first_five_SEM = mean(sem_first_five);
mean_second_five_Avg = mean(mean_second_five);
mean_second_five_SEM = mean(sem_second_five);
mean_third_five_Avg = mean(mean_third_five);
mean_third_five_SEM = mean(sem_third_five);
mean_fourth_five_Avg = mean(mean_fourth_five);
mean_fourth_five_SEM = mean(sem_fourth_five);
mean_fifth_five_Avg = mean(mean_fifth_five);
mean_fifth_five_SEM = mean(sem_fifth_five);


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


figure; 
shadedErrorBar(ts1, mean_first_five_Avg(1,:), mean_first_five_SEM(1,:));
hold on;
shadedErrorBar(ts1, mean_second_five_Avg(1,:), mean_second_five_SEM(1,:),'lineProps', '-b');
hold on;
shadedErrorBar(ts1, mean_third_five_Avg(1,:), mean_third_five_SEM(1,:),'lineProps', '-g');
hold on; 
shadedErrorBar(ts1, mean_fourth_five_Avg(1,:), mean_fourth_five_SEM(1,:),'lineProps', '-y');
hold on; 
shadedErrorBar(ts1, mean_fifth_five_Avg(1,:), mean_fifth_five_SEM(1,:),'lineProps', '-r');
legend('0.02 - 0.1 mA','0.12 - 0.20 mA', '0.22 - 0.30 mA', '0.32 - 0.40 mA', '0.42 - 0.50 mA')



