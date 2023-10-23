IND_no_nan = ~isnan(ZallMean_motion_B3(:,1))

pre_choice_ind = ts1(1,:) > -5 & ts1(1,:) < 0;

ZallMean_B1_no_nan = ZallMean_B1(IND_no_nan, :);
ZallMean_B2_no_nan = ZallMean_B2(IND_no_nan, :);
ZallMean_B3_no_nan = ZallMean_B3(IND_no_nan, :);

ZallMean_motion_B1_no_nan = ZallMean_motion_B1(IND_no_nan, :);
ZallMean_motion_B2_no_nan = ZallMean_motion_B2(IND_no_nan, :);
ZallMean_motion_B3_no_nan = ZallMean_motion_B3(IND_no_nan, :);
ZallMean_concat = vertcat(ZallMean_B1_no_nan(:, pre_choice_ind), ZallMean_B2_no_nan(:, pre_choice_ind),ZallMean_B3_no_nan(:, pre_choice_ind));
ZallMean_motion_concat = vertcat(ZallMean_motion_B1_no_nan(:, pre_choice_ind), ZallMean_motion_B2_no_nan(:, pre_choice_ind), ZallMean_motion_B3_no_nan(:, pre_choice_ind));


ZallMean_concat_mean = mean(ZallMean_concat,2);
ZallMean_motion_concat_mean = mean(ZallMean_motion_concat,2);





blocks(13:24,1) = 2;
blocks(1:12,1) = 1;
blocks(25:36, 1) = 3;

tbl_data = table(ZallMean_concat_mean,ZallMean_motion_concat_mean,blocks);

ZallMean_B1_no_nan_mean = mean(ZallMean_B1_no_nan);
ZallMean_B2_no_nan_mean = mean(ZallMean_B2_no_nan);
ZallMean_B3_no_nan_mean = mean(ZallMean_B3_no_nan);
figure; plot(ts1, ZallMean_B1_no_nan_mean); hold on; plot(ts1, ZallMean_B2_no_nan_mean);hold on; plot(ts1, ZallMean_B3_no_nan_mean);

ZallMean_motion_B1_no_nan_mean = mean(ZallMean_motion_B1_no_nan);
ZallMean_motion_B2_no_nan_mean = mean(ZallMean_motion_B2_no_nan);
ZallMean_motion_B3_no_nan_mean = mean(ZallMean_motion_B3_no_nan);
figure; plot(ts1, ZallMean_motion_B1_no_nan_mean); hold on; plot(ts1, ZallMean_motion_B2_no_nan_mean);hold on; plot(ts1, ZallMean_motion_B3_no_nan_mean);


X = [ZallMean_concat_mean ZallMean_motion_concat_mean];

blocks = ordinal(blocks)
% X = [ZallMean_concat blocks];

[B,dev,stats] = mnrfit(X, blocks,'model','ordinal');
% [B,dev,stats] = mnrfit(X, ZallMean_motion_concat);

s = scatter(tbl_data,'ZallMean_concat_mean','ZallMean_motion_concat_mean', 'filled', 'ColorVariable', 'blocks');


%% attempts since talking to Ruairi. run PhotometryAnalysis_new_Chamber_A_BL2nsp_v4_motion_integrated.m
pre_choice_ind = ts1(1,:) > -5 & ts1(1,:) < 0;
B1_ind = BehavData.Block == 1;
B2_ind = BehavData.Block == 2;
B3_ind = BehavData.Block == 3;

B1_mean = mean(zall(B1_ind, pre_choice_ind))';
B2_mean = mean(zall(B2_ind, pre_choice_ind))';
B3_mean = mean(zall(B3_ind, pre_choice_ind))';

B1_dF_mean = mean(Y_dF_all(B1_ind, pre_choice_ind))';
B2_dF_mean = mean(Y_dF_all(B2_ind, pre_choice_ind))';
B3_dF_mean = mean(Y_dF_all(B3_ind, pre_choice_ind))';

B1_allSignals_motion = mean(allSignals_motion(B1_ind, pre_choice_ind))';
B2_allSignals_motion = mean(allSignals_motion(B2_ind, pre_choice_ind))';
B3_allSignals_motion = mean(allSignals_motion(B3_ind, pre_choice_ind))';

B1_motion_mean = mean(zall_motion(B1_ind, pre_choice_ind))';
B2_motion_mean = mean(zall_motion(B2_ind, pre_choice_ind))';
B3_motion_mean = mean(zall_motion(B3_ind, pre_choice_ind))';


blocks(1:150,1) = 1;
blocks(151:300,1) = 2;
blocks(301:450, 1) = 3;

blocks = categorical(blocks);

large_rew_Ca_data = vertcat(B1_mean,B2_mean,B3_mean);

large_rew_Ca_data_unnormalized = vertcat(B1_dF_mean, B2_dF_mean, B3_dF_mean);
zb_Ca = mean(large_rew_Ca_data_unnormalized(:,1)); % entire window mean
zsd_Ca = std(large_rew_Ca_data_unnormalized(:,1)); % entire window stdev
zall_Ca_data=(large_rew_Ca_data_unnormalized - zb_Ca)/zsd_Ca;

large_rew_motion_data = vertcat(B1_motion_mean,B2_motion_mean,B3_motion_mean);

large_rew_motion_data_unnormalized = vertcat(B1_allSignals_motion, B2_allSignals_motion, B3_allSignals_motion);
zb_motion = mean(large_rew_motion_data_unnormalized(:,1)); % entire window mean
zsd_motion = std(large_rew_motion_data_unnormalized(:,1)); % entire window stdev
zall_motion_data=(large_rew_motion_data_unnormalized - zb_motion)/zsd_motion;

X_tbl = table(large_rew_motion_data, blocks, large_rew_Ca_data);

X_tbl_2 = table(zall_motion_data, blocks, zall_Ca_data);

Y = large_rew_Ca_data;

mdl = fitlm(X_tbl, 'interactions')

mdl_2 = fitlm(X_tbl_2, 'linear')

%%
writetable(XY_tbl_2, 'XY_tbl_concat_normalized.csv');