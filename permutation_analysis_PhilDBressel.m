%% ERT example
sig = .05;
consec_thresh = 10; % 1017.3Hz sample rate / 3Hz filter %340 PRD used 340 because his data WERE NOT DOWNSAMPLED (e.g., they were 1018 samples per sec) our data are 30 samples per sec or so.
% this means if his threshold is 1017/340 = ~3, ours should be 30/x = 3,
% which is 10

% Graphing parameters
ylims = [-1 2];
xlims = [-8 8];
sig_plot_level = linspace(4,3.2,7);

ind_2 = ts1(1,:) < xlims(2) & ts1(1,:) > xlims(1);

ZallMean_large_trunc = ZallMean_large(:, ind_2);
ZallMean_small_trunc = ZallMean_small(:, ind_2);

%%
[n_Cp,ev_win] = size(ZallMean_large_trunc);
[n_Cm,~] = size(ZallMean_small_trunc);
timeline = linspace(-8,8,ev_win);

mean_Cp = nanmean(ZallMean_large_trunc,1);
sem_Cp = sem(ZallMean_large_trunc);

mean_Cm = nanmean(ZallMean_small_trunc,1);
sem_Cm = sem(ZallMean_small_trunc);

%%
perm_p = permTest_array(ZallMean_large_trunc,ZallMean_small_trunc,1000);% permTest_array(ERT_test.Cp_off1,ERT_test.Cm_off3,1000);
diff_bCI = boot_diffCI(ZallMean_large_trunc,ZallMean_small_trunc,1000,sig);
[adjLCI,adjUCI] = CIadjust(diff_bCI(1,:),diff_bCI(2,:),[],n_Cm,2);
diff_bCIexp = [adjLCI;adjUCI];

%%
%Permutation test
perm_p_sig = NaN(1,ev_win);
sig_idx = find(perm_p < sig);
consec = consec_idx(sig_idx,consec_thresh);
perm_p_sig(sig_idx(consec)) = sig_plot_level(7);

%% Plot
figure; hold on
plot(timeline,mean_Cm,'Color',col_rep(3))
errorplot3(mean_Cm-sem_Cm,mean_Cm+sem_Cm,[-8 8],col_rep(3),.15)

plot(timeline,mean_Cp,'Color',col_rep(2))
errorplot3(mean_Cp-sem_Cp,mean_Cp+sem_Cp,[-8 8],col_rep(2),.15)

%Plot permutation test sig
plot(timeline,perm_p_sig,'Color',col_rep(1),'Marker','.')
text(xlims(1),sig_plot_level(7),'\bf Perm','Color',col_rep(1));

plot([-0.5 -0.5],ylim,'k:')
plot(xlim,[0 0],'k--')

xlim(xlims);