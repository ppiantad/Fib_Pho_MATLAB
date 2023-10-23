for ii = 1:size(ZallMean)
    pre_choice_corr = corrcoef(ZallMean(ii,(ts1(1,:) < 0 & ts1(1,:) > -5)),ZallMean_motion(ii,(ts1(1,:) < 0 & ts1(1,:) > -5)));
    pre_choice_corr_entire_window = corrcoef(ZallMean(ii,:),ZallMean_motion(ii,:));
    pre_choice_corr_list(ii) = pre_choice_corr(2,1);
    pre_choice_corr_entire_window_list(ii)=pre_choice_corr_entire_window(2,1);
end




pre_choice_corr_list_mean = nanmean(pre_choice_corr_list);
pre_choice_corr_list_sem = nansem(pre_choice_corr_list);

pre_choice_corr_entire_window_list_mean = nanmean(pre_choice_corr_entire_window_list);
pre_choice_corr_entire_window_list_sem = nansem(pre_choice_corr_entire_window_list);

idx = ~isnan(pre_choice_corr_list);

pre_choice_corr_list_2 = pre_choice_corr_list(idx);

pre_choice_corr_entire_window_list_2 = pre_choice_corr_entire_window_list(idx);

corrTable_motion = corrTable.animalNames(idx);

% nanstd(pre_choice_corr_list)/sqrt(numel(pre_choice_corr_list(idx)))



%%
animalNum_var = 5;

[numTrials,~]=size(behav_struct_filtered(animalNum_var).data.collectionTime(:));
Tris=[1:numTrials]';

time2Collect = behav_struct_filtered(animalNum_var).data.collectionTime(:)-behav_struct_filtered(animalNum_var).data.choiceTime(:);
%%
figure; 
fontsize(gca,18,"pixels")
% shadedErrorBar(ts1, zAvg(1,:), SEM_real(1,:), 'lineProps', {'color', [255, 83, 26]/255});
shadedErrorBar(ts1, ZallMean(animalNum_var,:), zerror(animalNum_var,:), 'lineProps', 'r');
ylabel('z-scored dF/F', 'FontSize', 20);
xlabel('Time from choice (s)');
xlim([-8 8]);
ylim([-2 4]);


hold on;
shadedErrorBar(ts1, ZallMean_motion(animalNum_var,:), zerror_motion(animalNum_var,:), 'lineProps', 'k');
ylabel('z-scored dF/F', 'FontSize', 20);
xlabel('Time from choice (s)');
xlim([-8 8]);
legend('GCaMP','Velocity')

hold off
%%

% time2Collect = BehavData.collectionTime(:)-BehavData.stTime(:);
% time2Choice = BehavData.choiceTime(:)-BehavData.stTime(:);
% 
% large_choice = BehavData.bigSmall == 1.2;
% small_choice = BehavData.bigSmall == 0.3;
% shk_outcome = BehavData.shock == 1;
% 
% [numTrials,~]=size(BehavData.collectionTime(:));
% Tris=[1:numTrials]';


figure;

IM=imagesc(ts1, 1, zall_filtered(animalNum_var).data); hold on;
fontsize(gca,30,"pixels")
load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 


scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
%scatter(time2choose,Tris,'Marker','>','MarkerFaceColor','k')
% title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 34);
xlim([-10 10]); %xlim([-10 10]);xlim([-8 30])
colorbar;
hold off;


figure;

IM=imagesc(ts1, 1, zall_motion_filtered(animalNum_var).data); hold on;
fontsize(gca,30,"pixels")
load('batlowW.mat'); %using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
colormap(batlowW); % c1 = colorbar; 


scatter(time2Collect,Tris,'Marker','p','MarkerFaceColor','w')
plot(zeros(numTrials,1),Tris)
%scatter(time2choose,Tris,'Marker','>','MarkerFaceColor','k')
% title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 34);
xlim([-10 10]); %xlim([-10 10]);xlim([-8 30])
colorbar;
hold off;


