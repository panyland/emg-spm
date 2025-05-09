%% EMG envelope and statistical parametric mapping

% Parameters
n_subjects = 13;
n_points = 1000;
window_size = 100;
condition_labels = {'Norm', 'Exo'};
colors = {'b', 'r'};

% Store per-subject averages
valid_subjects = setdiff(1:13, [2,6]); % Skip subjects if necessary
n_valid = numel(valid_subjects);
all_subject_avg = cell(1, 2); % norm, exo
all_subject_avg{1} = zeros(n_valid, n_points);
all_subject_avg{2} = zeros(n_valid, n_points);

subj_idx = 1;
for subj = valid_subjects
    
    for cond = 1:2
        suffix = {'norm', 'exo'};
        filename = sprintf('data/%d_%s.txt', subj, suffix{cond});
        
        data = readtable(filename, DecimalSeparator=",");
        time = data{:, 1};
        markers = data{:, "Markers"}; % Column name for markers
        emg = data{:,"LUMBARESLT_uV"}; % Column name for particular electrode
     
        % Find reps
        rep_starts = find(markers == 1);
        rep_ends = zeros(size(rep_starts));
        for j = 1:length(rep_starts)
            next4 = find(markers(rep_starts(j):end) == 4, 1, 'first');
            if ~isempty(next4)
                rep_ends(j) = rep_starts(j) + next4 - 1;
            end
        end
        valid = rep_ends > rep_starts;
        rep_starts = rep_starts(valid);
        rep_ends = rep_ends(valid);

        % Interpolate reps
        emg_interp = zeros(length(rep_starts), n_points);
        for j = 1:length(rep_starts)
            t_seg = time(rep_starts(j):rep_ends(j));
            emg_seg = emg(rep_starts(j):rep_ends(j));
            t_norm = linspace(t_seg(1), t_seg(end), n_points);
            emg_interp(j, :) = interp1(t_seg, emg_seg, t_norm);
        end

        % Average reps first, then smooth
        emg_avg = mean(emg_interp, 1);
        emg_avg_smooth = movmean(emg_avg, window_size);

        all_subject_avg{cond}(subj_idx, :) = emg_avg_smooth;
    end
    subj_idx = subj_idx + 1;
end

% Group averages
group_avg_norm = mean(all_subject_avg{1}, 1);
group_avg_exo  = mean(all_subject_avg{2}, 1);

% Compute std across subjects
std_norm = std(all_subject_avg{1}, 0, 1);
std_exo  = std(all_subject_avg{2}, 0, 1);

x = linspace(0, 100, n_points);

% Add spm1d path
addpath('C:\Users\Paavo\Documents\MATLAB\0todd0000-spm1dmatlab-b80f01b');
addpath('C:\Users\Paavo\Documents\MATLAB\0todd0000-spm1dmatlab-b80f01b\spm8\');

% Run paired t-test across time
spm = spm1d.stats.ttest_paired(all_subject_avg{1}, all_subject_avg{2});
spmi = spm.inference(0.05, 'two_tailed', true, 'interp', true);
figure;
spmi.plot()

% Plot results with shaded significance
figure;
hold on;

% Std fill
fill([x, fliplr(x)], ...
     [group_avg_norm + std_norm, fliplr(group_avg_norm - std_norm)], ...
     colors{1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');
fill([x, fliplr(x)], ...
     [group_avg_exo + std_exo, fliplr(group_avg_exo - std_exo)], ...
     colors{2}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

% Plot means
plot(x, group_avg_norm, 'Color', colors{1}, 'LineWidth', 2, 'DisplayName', 'Norm');
plot(x, group_avg_exo,  'Color', colors{2}, 'LineWidth', 2, 'DisplayName', 'Exo');

% Highlight significant regions
if spmi.h0reject
    clusters = spmi.clusters;
    for i = 1:length(clusters)
        cluster = clusters{i};  
        endpoints = cluster.endpoints;

        % Find the closest indices in x to the cluster's start and end points
        [~, start_idx] = min(abs(x - endpoints(1)/10));
        [~, end_idx]   = min(abs(x - endpoints(2)/10));

        sig_x = x(start_idx:end_idx);

        y_min = min([group_avg_norm, group_avg_exo], [], 'all') - 0.05;
        y_max = y_min + 0.05;
        
        if i == 1
            display_name = 'Significant Region';
        else
            display_name = '';
        end
        bar_y_low = 0.01;
        bar_y_high = 0.03;
        % Plot horizontal significance bar
        fill([sig_x, fliplr(sig_x)], ...
             [bar_y_low * ones(size(sig_x)), bar_y_high * ones(size(sig_x))], ...
             'k', 'EdgeColor', 'none', 'DisplayName', 'SPM paired t-test');
    
        % Add "p < 0.05" text label near the middle of the bar
        text_x = mean(sig_x);
        text_y = bar_y_high + 0.02;  % Slightly above the bar
        text(text_x, text_y, 'p < 0.05', ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 9, ...
             'FontWeight', 'bold', ...
             'Color', 'k');

    end
else
    disp('No significant regions at this threshold.');
end

xlabel('% Movement Cycle');
ylabel('RMS EMG (uV)');
title('Grand average Left Lumbar EMG (Stretcher lift)');
legend show;


