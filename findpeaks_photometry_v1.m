function [block1_pks_per_min_filtered, block2_pks_per_min_filtered, block3_pks_per_min_filtered] = findpeaks_photometry_v1(timedownsample, Y_dF_all, block_end, animalName);
%%
%use photometry_signal_check_v3.m to get whole session timeseries data

% [BehavData,ABETfile,Descriptives, block_end]=ABET2TableFn_Chamber_Av4('RRD23_03192019.csv',[]);

%assuming data were downsampled by a factor of 10
% Fs = data.streams.FiPr.fs/10;
if length(block_end) < 2
    disp(['Animal' animalName 'did not complete all three blocks, skipping'])
    block1_pks_per_min_filtered = nan;
    block2_pks_per_min_filtered = nan;
    block3_pks_per_min_filtered = nan;

elseif length(block_end) >=2

    % findpeaks(Y_dF_all,timedownsample,'MinPeakDistance',6)
    idx_1 = timedownsample < block_end(1);
    idx_2 = timedownsample > block_end(1) & timedownsample < block_end(2);
    idx_3 = timedownsample > block_end(2);

    block1_time(1,1:sum(idx_1)) = timedownsample(idx_1);
    block2_time(1,1:sum(idx_2)) = timedownsample(idx_2);
    block3_time(1,1:sum(idx_3)) = timedownsample(idx_3);

    Block1_Y_dF_all(1,1:sum(idx_1)) = Y_dF_all(idx_1);
    Block2_Y_dF_all(1,1:sum(idx_2)) = Y_dF_all(idx_2);
    Block3_Y_dF_all(1,1:sum(idx_3)) = Y_dF_all(idx_3);

    Block1_Y_dF_all_filtered = sgolayfilt(Block1_Y_dF_all, 1, 71);
    Block2_Y_dF_all_filtered = sgolayfilt(Block2_Y_dF_all, 1, 71);
    Block3_Y_dF_all_filtered = sgolayfilt(Block3_Y_dF_all, 1, 71);

    % figure;
    % plot(block1_time, Block1_Y_dF_all);
    %
    % figure;
    % plot(block2_time, Block2_Y_dF_all);
    %
    % figure;
    % plot(block3_time, Block3_Y_dF_all);

    [pks_filtered,locs_filtered, proms_filtered] = findpeaks(Block1_Y_dF_all_filtered,block1_time,'MinPeakDistance',3);
    [pks2_filtered,locs2_filtered,proms2_filtered] = findpeaks(Block2_Y_dF_all_filtered,block2_time,'MinPeakDistance',3);
    [pks3_filtered,locs3_filtered,proms3_filtered] = findpeaks(Block3_Y_dF_all_filtered,block3_time,'MinPeakDistance',3);





    [pks,locs, proms] = findpeaks(Block1_Y_dF_all,block1_time,'MinPeakDistance',3);
    [pks2,locs2,proms2] = findpeaks(Block2_Y_dF_all,block2_time,'MinPeakDistance',3);
    [pks3,locs3,proms3] = findpeaks(Block3_Y_dF_all,block3_time,'MinPeakDistance',3);

    
    figure;
    plot(block1_time, Block1_Y_dF_all);
    hold on
    plot(locs, pks, 'pg', 'MarkerFaceColor', 'g')
    hold off
    
    figure;
    plot(block2_time, Block2_Y_dF_all);
    hold on
    plot(locs2, pks2, 'pg', 'MarkerFaceColor', 'g')
    hold off
    
    figure;
    plot(block3_time, Block3_Y_dF_all);
    hold on
    plot(locs3, pks3, 'pg', 'MarkerFaceColor', 'g')
    hold off
    
    
    figure;
    plot(block1_time, Block1_Y_dF_all_filtered);
    xlim([block1_time(1) block1_time(end)]);
    ylim([-50 100]);
    hold on
    plot(locs_filtered, pks_filtered, 'pg', 'MarkerFaceColor', 'g')
    hold off
    
    figure;
    plot(block2_time, Block2_Y_dF_all_filtered);
    xlim([block2_time(1) block2_time(end)]);
    ylim([-50 100]);
    hold on
    plot(locs2_filtered, pks2_filtered, 'pg', 'MarkerFaceColor', 'g')
    hold off
    
    figure;
    plot(block3_time, Block3_Y_dF_all_filtered);
    xlim([block3_time(1) block3_time(end)]);
    ylim([-50 100]);
    hold on
    plot(locs3_filtered, pks3_filtered, 'pg', 'MarkerFaceColor', 'g')
    hold off



    block1_len = block1_time(end) - block1_time(1);
    block1_min = block1_len/60;

    block2_len = block2_time(end) - block2_time(1);
    block2_min = block2_len/60;


    block3_len = block3_time(end) - block3_time(1);
    block3_min = block3_len/60;

    block1_pks_per_min = numel(pks)/block1_len;
    block2_pks_per_min = numel(pks2)/block2_len;
    block3_pks_per_min = numel(pks3)/block3_len;


    block1_pks_per_min_filtered = numel(pks_filtered)/block1_len;
    block2_pks_per_min_filtered = numel(pks2_filtered)/block2_len;
    block3_pks_per_min_filtered = numel(pks3_filtered)/block3_len;

    pks_per_min_filtered = [block1_pks_per_min_filtered block2_pks_per_min_filtered block3_pks_per_min_filtered];
end
end

