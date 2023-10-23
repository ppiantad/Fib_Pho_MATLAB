% Use this code after first running access_behav_struct_v5 (or latest
% iteration) using TrialFilter ALL, 1



% Initialize cell arrays to store latencies for different trial types, blocks, and choice types
latencies_12_reward_free = cell(1, 3); % 3 blocks
latencies_12_reward_forced = cell(1, 3); % 3 blocks
latencies_03_reward_free = cell(1, 3); % 3 blocks
latencies_03_reward_forced = cell(1, 3); % 3 blocks
trial_blocks = [1, 2, 3];

% Initialize variables to store maximum frequency and latency limits
max_frequency = 0;
latency_limits = [0, 30];

% Set bin width for histogram
bin_width = 0.5;

% Loop through your struct
for i = 1:numel(behav_struct_filtered)
    data = behav_struct_filtered(i).data;
    
    % Extract the relevant data for choiceTime, stTime, and ForceFree
    choiceTime = data.choiceTime;
    stTime = data.stTime;
    forceFree = data.ForceFree;
    
    % Calculate the latencies for each trial
    latencies = choiceTime - stTime;
    
    % Determine the trial type (1.2 or 0.3 reward)
    trial_type = data.bigSmall;
    
    % Loop through trial blocks (1, 2, 3)
    for block = trial_blocks
        % Filter the latencies for the current trial block
        latencies_block = latencies(data.Block == block);
        
        % Find the indices for 1.2 and 0.3 reward trials within this block
        indices_12_reward_free = find(trial_type == 1.2 & data.Block == block & forceFree == 0);
        indices_12_reward_forced = find(trial_type == 1.2 & data.Block == block & forceFree == 1);
        indices_03_reward_free = find(trial_type == 0.3 & data.Block == block & forceFree == 0);
        indices_03_reward_forced = find(trial_type == 0.3 & data.Block == block & forceFree == 1);
        
        % Separate latencies based on ForceFree (Forced or Free Choice)
        latencies_12_reward_free{block} = [latencies_12_reward_free{block}; latencies(indices_12_reward_free)];
        latencies_12_reward_forced{block} = [latencies_12_reward_forced{block}; latencies(indices_12_reward_forced)];
        latencies_03_reward_free{block} = [latencies_03_reward_free{block}; latencies(indices_03_reward_free)];
        latencies_03_reward_forced{block} = [latencies_03_reward_forced{block}; latencies(indices_03_reward_forced)];

        % Update the maximum frequency
        max_frequency = max([max_frequency, ...
            max(histcounts(latencies_block, 'BinWidth', bin_width))]);
    end
    clear data indices_12_reward_free indices_12_reward_forced indices_03_reward_free indices_03_reward_forced choiceTime stTime latencies trial_type latencies_block
end

% Create histograms for each combination of block and reward type
for block = trial_blocks
    % For 1.2 reward trials - Free Choice
    figure;
    subplot(2, 1, 1);
    histogram(latencies_12_reward_free{block}, 'BinWidth', bin_width);
    title(['Latency Histogram for 1.2 Reward Trials - Free Choice - Block ' num2str(block)]);
    xlabel('Latency');
    ylabel('Frequency');
    ylim([0, max_frequency]);
    xlim(latency_limits);
    
    % Calculate and plot the mean and median values as vertical dashed lines
    mean_latencies_free = mean(latencies_12_reward_free{block});
    median_latencies_free = median(latencies_12_reward_free{block});
    hold on;
    plot([mean_latencies_free, mean_latencies_free], [0, max_frequency], 'r--', 'LineWidth', 1.5);
    plot([median_latencies_free, median_latencies_free], [0, max_frequency], 'k--', 'LineWidth', 1.5);
    legend({'Latencies', ['Mean (' num2str(mean_latencies_free) ')'], ['Median (' num2str(median_latencies_free) ')']});
    
    % For 1.2 reward trials - Forced Choice
    subplot(2, 1, 2);
    histogram(latencies_12_reward_forced{block}, 'BinWidth', bin_width);
    title(['Latency Histogram for 1.2 Reward Trials - Forced Choice - Block ' num2str(block)]);
    xlabel('Latency');
    ylabel('Frequency');
    ylim([0, max_frequency]);
    xlim(latency_limits);
    
    % Calculate and plot the mean and median values as vertical dashed lines
    mean_latencies_forced = mean(latencies_12_reward_forced{block});
    median_latencies_forced = median(latencies_12_reward_forced{block});
    hold on;
    plot([mean_latencies_forced, mean_latencies_forced], [0, max_frequency], 'r--', 'LineWidth', 1.5);
    plot([median_latencies_forced, median_latencies_forced], [0, max_frequency], 'k--', 'LineWidth', 1.5);
    legend({'Latencies', ['Mean (' num2str(mean_latencies_forced) ')'], ['Median (' num2str(median_latencies_forced) ')']});
end

for block = trial_blocks
    % For 0.3 reward trials - Free Choice
    figure;
    subplot(2, 1, 1);
    histogram(latencies_03_reward_free{block}, 'BinWidth', bin_width);
    title(['Latency Histogram for 0.3 Reward Trials - Free Choice - Block ' num2str(block)]);
    xlabel('Latency');
    ylabel('Frequency');
    ylim([0, max_frequency]);
    xlim(latency_limits);
    
    % Calculate and plot the mean and median values as vertical dashed lines
    mean_latencies_free = mean(latencies_03_reward_free{block});
    median_latencies_free = median(latencies_03_reward_free{block});
    hold on;
    plot([mean_latencies_free, mean_latencies_free], [0, max_frequency], 'r--', 'LineWidth', 1.5);
    plot([median_latencies_free, median_latencies_free], [0, max_frequency], 'k--', 'LineWidth', 1.5);
    legend({'Latencies', ['Mean (' num2str(mean_latencies_free) ')'], ['Median (' num2str(median_latencies_free) ')']});

    % For 0.3 reward trials - Forced Choice
    subplot(2, 1, 2);
    histogram(latencies_03_reward_forced{block}, 'BinWidth', bin_width);
    title(['Latency Histogram for 0.3 Reward Trials - Forced Choice - Block ' num2str(block)]);
    xlabel('Latency');
    ylabel('Frequency');
    ylim([0, max_frequency]);
    xlim(latency_limits);
    
    % Calculate and plot the mean and median values as vertical dashed lines
    mean_latencies_forced = mean(latencies_03_reward_forced{block});
    median_latencies_forced = median(latencies_03_reward_forced{block});
    hold on;
    plot([mean_latencies_forced, mean_latencies_forced], [0, max_frequency], 'r--', 'LineWidth', 1.5);
    plot([median_latencies_forced, median_latencies_forced], [0, max_frequency], 'k--', 'LineWidth', 1.5);
    legend({'Latencies', ['Mean (' num2str(mean_latencies_forced) ')'], ['Median (' num2str(median_latencies_forced) ')']});
end




%%

% % Create histograms for each combination of block and reward type
% for block = trial_blocks
%     % For 1.2 reward trials
%     create_histogram(latencies_12_reward_free{block}, latencies_12_reward_forced{block}, '1.2 Reward Trials', block);
% 
%     % For 0.3 reward trials
%     create_histogram(latencies_03_reward_free{block}, latencies_03_reward_forced{block}, '0.3 Reward Trials', block);
% end

% % Define the create_histogram function
% function create_histogram(latencies_free, latencies_forced, reward_type, block, bin_width, max_frequency, latency_limits)
%     figure;
%     
%     subplot(2, 1, 1); % for Free Choice
%     histogram(latencies_free, 'BinWidth', bin_width);
%     title(['Latency Histogram for ' reward_type ' - Free Choice - Block ' num2str(block)]);
%     xlabel('Latency');
%     ylabel('Frequency');
%     ylim([0, max_frequency]);
%     xlim(latency_limits);
% 
%     subplot(2, 1, 2); % for Forced Choice
%     histogram(latencies_forced, 'BinWidth', bin_width);
%     title(['Latency Histogram for ' reward_type ' - Forced Choice - Block ' num2str(block)]);
%     xlabel('Latency');
%     ylabel('Frequency');
%     ylim([0, max_frequency]);
%     xlim(latency_limits);
% end