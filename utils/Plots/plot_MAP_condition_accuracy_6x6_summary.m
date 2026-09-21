function plot_MAP_condition_accuracy_6x6_summary(area, savepaths, varargin)
% Plot 6x6 decoder accuracy and prediction-error heatmaps across sessions.
%
% Outputs:
%   - ME decoder accuracy/error heatmaps
%   - SE-add decoder accuracy/error heatmaps
%   - Split ME and SE congruent/incongruent bar plots
%   - Single-neuron mean FR scatter: congruent vs incongruent ME trials

    p = inputParser;
    addParameter(p, 'NBins', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(p, 'PlotBestSessionToo', true, @(x)islogical(x)||ismember(x,[0 1]));
    parse(p, varargin{:});

    nbins = p.Results.NBins;
    plot_best = logical(p.Results.PlotBestSessionToo);
    bin_deg = 360 / nbins;

    chance_acc = 1 / (nbins^2);
    chance_err = expected_joint_chance_error_deg(nbins, bin_deg);

    counts_true_ME = zeros(nbins, nbins);
    counts_corr_ME = zeros(nbins, nbins);
    counts_true_SE = zeros(nbins, nbins);
    counts_corr_SE = zeros(nbins, nbins);

    counts_err_ME = zeros(nbins, nbins);
    sum_err_ME    = zeros(nbins, nbins);
    counts_err_SE = zeros(nbins, nbins);
    sum_err_SE    = zeros(nbins, nbins);

    used = strings(0,1);
    skipped = strings(0,1);

    best_score = -Inf;
    best_file = "";
    best_ME_map = [];
    best_SE_map = [];
    best_ME_err_map = [];
    best_SE_err_map = [];

    sess_cong_ME = [];
    sess_incong_ME = [];
    sess_cong_SE = [];
    sess_incong_SE = [];

    sess_err_cong_ME = [];
    sess_err_incong_ME = [];
    sess_err_cong_SE = [];
    sess_err_incong_SE = [];

    % ---------- NEW: single-neuron activity accumulators ----------
    all_neuron_FR_cong = [];
    all_neuron_FR_incong = [];
    all_neuron_session = [];

    for s = 1:numel(savepaths)
        sp = savepaths{s};
        f = fullfile(sp, 'models', sprintf('%s_MAP_Smoothed.mat', area));

        if ~exist(f, 'file')
            warning('Missing file: %s', f);
            skipped(end+1,1) = string(f); %#ok<AGROW>
            continue;
        end

        S = load(f);
        if ~isfield(S, 'Results')
            warning('No Results struct in: %s', f);
            skipped(end+1,1) = string(f); %#ok<AGROW>
            continue;
        end

        R = S.Results;

        if ~all(isfield(R, {'true_flat','idx_hat_ME','idx_hat_SE'}))
            warning(['Skipping %s because it does not contain joint prediction fields ', ...
                     'true_flat, idx_hat_ME, idx_hat_SE.'], f);
            skipped(end+1,1) = string(f); %#ok<AGROW>
            continue;
        end

        true_flat = R.true_flat(:);
        idx_ME    = R.idx_hat_ME(:);
        idx_SE    = R.idx_hat_SE(:);

        ok = isfinite(true_flat) & isfinite(idx_ME) & isfinite(idx_SE) & ...
             true_flat >= 1 & idx_ME >= 1 & idx_SE >= 1 & ...
             true_flat <= nbins^2 & idx_ME <= nbins^2 & idx_SE <= nbins^2;

        true_flat = true_flat(ok);
        idx_ME    = idx_ME(ok);
        idx_SE    = idx_SE(ok);

        if isempty(true_flat)
            warning('No valid decoded trials in: %s', f);
            skipped(end+1,1) = string(f); %#ok<AGROW>
            continue;
        end

        [Htrue, Etrue] = ind2sub([nbins nbins], true_flat);
        [Hhat_ME, Ehat_ME] = ind2sub([nbins nbins], idx_ME);
        [Hhat_SE, Ehat_SE] = ind2sub([nbins nbins], idx_SE);

        correct_ME = (Hhat_ME == Htrue) & (Ehat_ME == Etrue);
        correct_SE = (Hhat_SE == Htrue) & (Ehat_SE == Etrue);

        true_counts_sess = accumarray([Htrue, Etrue], 1, [nbins nbins], @sum, 0);
        corr_ME_sess     = accumarray([Htrue, Etrue], double(correct_ME), [nbins nbins], @sum, 0);
        corr_SE_sess     = accumarray([Htrue, Etrue], double(correct_SE), [nbins nbins], @sum, 0);

        counts_true_ME = counts_true_ME + true_counts_sess;
        counts_corr_ME = counts_corr_ME + corr_ME_sess;
        counts_true_SE = counts_true_SE + true_counts_sess;
        counts_corr_SE = counts_corr_SE + corr_SE_sess;

        acc_ME_sess = safe_divide(corr_ME_sess, true_counts_sess);
        acc_SE_sess = safe_divide(corr_SE_sess, true_counts_sess);

        err_ME = joint_circular_bin_error_deg(Htrue, Etrue, Hhat_ME, Ehat_ME, nbins, bin_deg);
        err_SE = joint_circular_bin_error_deg(Htrue, Etrue, Hhat_SE, Ehat_SE, nbins, bin_deg);

        sum_err_ME_sess = accumarray([Htrue, Etrue], err_ME, [nbins nbins], @sum, 0);
        sum_err_SE_sess = accumarray([Htrue, Etrue], err_SE, [nbins nbins], @sum, 0);

        sum_err_ME = sum_err_ME + sum_err_ME_sess;
        sum_err_SE = sum_err_SE + sum_err_SE_sess;
        counts_err_ME = counts_err_ME + true_counts_sess;
        counts_err_SE = counts_err_SE + true_counts_sess;

        err_ME_sess_map = safe_divide(sum_err_ME_sess, true_counts_sess);
        err_SE_sess_map = safe_divide(sum_err_SE_sess, true_counts_sess);

        diag_mask_sess = eye(nbins) == 1;
        off_mask_sess  = ~diag_mask_sess;

        sess_cong_ME(end+1,1)   = mean(acc_ME_sess(diag_mask_sess), 'omitnan'); %#ok<AGROW>
        sess_incong_ME(end+1,1) = mean(acc_ME_sess(off_mask_sess),  'omitnan'); %#ok<AGROW>
        sess_cong_SE(end+1,1)   = mean(acc_SE_sess(diag_mask_sess), 'omitnan'); %#ok<AGROW>
        sess_incong_SE(end+1,1) = mean(acc_SE_sess(off_mask_sess),  'omitnan'); %#ok<AGROW>

        sess_err_cong_ME(end+1,1)   = mean(err_ME_sess_map(diag_mask_sess), 'omitnan'); %#ok<AGROW>
        sess_err_incong_ME(end+1,1) = mean(err_ME_sess_map(off_mask_sess),  'omitnan'); %#ok<AGROW>
        sess_err_cong_SE(end+1,1)   = mean(err_SE_sess_map(diag_mask_sess), 'omitnan'); %#ok<AGROW>
        sess_err_incong_SE(end+1,1) = mean(err_SE_sess_map(off_mask_sess),  'omitnan'); %#ok<AGROW>

        score_sess = mean([acc_ME_sess(:); acc_SE_sess(:)], 'omitnan');
        if isfinite(score_sess) && score_sess > best_score
            best_score = score_sess;
            best_file = string(f);
            best_ME_map = acc_ME_sess;
            best_SE_map = acc_SE_sess;
            best_ME_err_map = err_ME_sess_map;
            best_SE_err_map = err_SE_sess_map;
        end

        used(end+1,1) = string(f); %#ok<AGROW>

        % ---------- NEW: single-neuron congruent vs incongruent mean FR ----------
        data_file = processed_data_file(sp, area);
        if exist(data_file, 'file')
            D = load(data_file);
            if isfield(D, 'data') && all(isfield(D.data, {'FR','TP1','TP2','TPi1','TPi2'}))
                data = D.data;

                idxME_activity = find(data.TPi2 ~= 7 & data.TPi1 ~= 7);
                FR_me_activity = data.FR(idxME_activity, :);
                TP1_me_activity = data.TP1(idxME_activity, :);
                TP2_me_activity = data.TP2(idxME_activity, :);

                thetaH_activity = atan2(TP1_me_activity(:,2), TP1_me_activity(:,1));
                thetaE_activity = atan2(TP2_me_activity(:,2), TP2_me_activity(:,1));

                theta_edges = linspace(-pi, pi, nbins+1) + 0.1;
                Hbin_activity = discretize(thetaH_activity, theta_edges);
                Ebin_activity = discretize(thetaE_activity, theta_edges);

                valid_activity = isfinite(Hbin_activity) & isfinite(Ebin_activity);
                cong_activity = valid_activity & (Hbin_activity == Ebin_activity);
                incong_activity = valid_activity & (Hbin_activity ~= Ebin_activity);

                if any(cong_activity) && any(incong_activity)
                    fr_cong = mean(FR_me_activity(cong_activity, :), 1, 'omitnan')';
                    fr_incong = mean(FR_me_activity(incong_activity, :), 1, 'omitnan')';

                    keep_neurons = isfinite(fr_cong) & isfinite(fr_incong);
                    all_neuron_FR_cong = [all_neuron_FR_cong; fr_cong(keep_neurons)]; %#ok<AGROW>
                    all_neuron_FR_incong = [all_neuron_FR_incong; fr_incong(keep_neurons)]; %#ok<AGROW>
                    all_neuron_session = [all_neuron_session; s * ones(sum(keep_neurons),1)]; %#ok<AGROW>
                else
                    warning('No congruent or incongruent ME trials found for activity plot in: %s', data_file);
                end
            else
                warning('Processed data file missing required fields for activity plot: %s', data_file);
            end
        else
            warning('Missing processed data file for activity plot: %s', data_file);
        end
    end

    if isempty(used)
        error(['No usable session files found. Exact 6x6 condition accuracy/error requires ', ...
               'Results.true_flat, Results.idx_hat_ME, and Results.idx_hat_SE.']);
    end

    acc_ME = safe_divide(counts_corr_ME, counts_true_ME);
    acc_SE = safe_divide(counts_corr_SE, counts_true_SE);
    err_ME_map = safe_divide(sum_err_ME, counts_err_ME);
    err_SE_map = safe_divide(sum_err_SE, counts_err_SE);

    diag_mask = eye(nbins) == 1;
    off_mask  = ~diag_mask;

    root_for_output = savepaths{1};
    plotdir = fullfile(root_for_output, '..', 'plots', 'summary', 'decoding');
    modeldir = fullfile(root_for_output, 'models');

    if ~exist(plotdir, 'dir'), mkdir(plotdir); end
    if ~exist(modeldir, 'dir'), mkdir(modeldir); end

    Summary = struct();
    Summary.area = area;
    Summary.usedFiles = used;
    Summary.skippedFiles = skipped;
    Summary.nbins = nbins;
    Summary.bin_deg = bin_deg;
    Summary.chance_acc = chance_acc;
    Summary.chance_err_deg = chance_err;

    Summary.counts_true = counts_true_ME;
    Summary.counts_correct_ME = counts_corr_ME;
    Summary.counts_correct_SE = counts_corr_SE;

    Summary.acc_ME_6x6 = acc_ME;
    Summary.acc_SE_6x6 = acc_SE;
    Summary.err_ME_6x6_deg = err_ME_map;
    Summary.err_SE_6x6_deg = err_SE_map;

    Summary.congruent_ME = mean(acc_ME(diag_mask), 'omitnan');
    Summary.incongruent_ME = mean(acc_ME(off_mask), 'omitnan');
    Summary.congruent_SE = mean(acc_SE(diag_mask), 'omitnan');
    Summary.incongruent_SE = mean(acc_SE(off_mask), 'omitnan');

    Summary.err_congruent_ME_deg = mean(err_ME_map(diag_mask), 'omitnan');
    Summary.err_incongruent_ME_deg = mean(err_ME_map(off_mask), 'omitnan');
    Summary.err_congruent_SE_deg = mean(err_SE_map(diag_mask), 'omitnan');
    Summary.err_incongruent_SE_deg = mean(err_SE_map(off_mask), 'omitnan');

    Summary.sessionwise.congruent_ME = sess_cong_ME;
    Summary.sessionwise.incongruent_ME = sess_incong_ME;
    Summary.sessionwise.congruent_SE = sess_cong_SE;
    Summary.sessionwise.incongruent_SE = sess_incong_SE;

    Summary.sessionwise.err_congruent_ME_deg = sess_err_cong_ME;
    Summary.sessionwise.err_incongruent_ME_deg = sess_err_incong_ME;
    Summary.sessionwise.err_congruent_SE_deg = sess_err_cong_SE;
    Summary.sessionwise.err_incongruent_SE_deg = sess_err_incong_SE;

    Summary.activity.congruent_FR = all_neuron_FR_cong;
    Summary.activity.incongruent_FR = all_neuron_FR_incong;
    Summary.activity.session_index = all_neuron_session;
    Summary.activity.delta_cong_minus_incong = all_neuron_FR_cong - all_neuron_FR_incong;

    Summary.bestSession.file = best_file;
    Summary.bestSession.score = best_score;
    Summary.bestSession.acc_ME_6x6 = best_ME_map;
    Summary.bestSession.acc_SE_6x6 = best_SE_map;
    Summary.bestSession.err_ME_6x6_deg = best_ME_err_map;
    Summary.bestSession.err_SE_6x6_deg = best_SE_err_map;

    % ---------- accuracy heatmaps ----------
    cmap_acc = make_light_to_dark_green_colormap(256);

    f = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6.3 11.2]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot_condition_heatmap(acc_ME, [0 1], cmap_acc, nbins, chance_acc, ...
        sprintf('%s ME decoder: condition accuracy', area), 'Accuracy', 'accuracy');

    nexttile;
    plot_condition_heatmap(acc_SE, [0 1], cmap_acc, nbins, chance_acc, ...
        sprintf('%s SE-add decoder: condition accuracy', area), 'Accuracy', 'accuracy');

    outpdf_acc = fullfile(plotdir, sprintf('%s_MAP_ConditionAccuracy_6x6_AcrossSessions.pdf', area));
    print(f, outpdf_acc, '-dpdf', '-painters');
    close(f);

    % ---------- error heatmaps ----------
    max_err = max([err_ME_map(:); err_SE_map(:); chance_err], [], 'omitnan');
    if ~isfinite(max_err) || max_err <= 0
        max_err = 180;
    end
    max_err = max(180, ceil(max_err/10)*10);

    f_err = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6.3 11.2]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot_condition_heatmap(err_ME_map, [0 max_err], parula(256), nbins, chance_err, ...
        sprintf('%s ME decoder: prediction error', area), 'Error (deg)', 'error');

    nexttile;
    plot_condition_heatmap(err_SE_map, [0 max_err], parula(256), nbins, chance_err, ...
        sprintf('%s SE-add decoder: prediction error', area), 'Error (deg)', 'error');

    outpdf_err = fullfile(plotdir, sprintf('%s_MAP_ConditionError_6x6_AcrossSessions.pdf', area));
    print(f_err, outpdf_err, '-dpdf', '-painters');
    close(f_err);

    % ---------- optional best-session heatmaps ----------
    if plot_best && ~isempty(best_ME_map)
        f2 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6.3 11.2]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile;
        plot_condition_heatmap(best_ME_map, [0 1], cmap_acc, nbins, chance_acc, ...
            sprintf('%s BEST ME decoder accuracy', area), 'Accuracy', 'accuracy');

        nexttile;
        plot_condition_heatmap(best_SE_map, [0 1], cmap_acc, nbins, chance_acc, ...
            sprintf('%s BEST SE-add decoder accuracy', area), 'Accuracy', 'accuracy');

        outpdf_best = fullfile(plotdir, sprintf('%s_MAP_ConditionAccuracy_6x6_BestSession.pdf', area));
        print(f2, outpdf_best, '-dpdf', '-painters');
        close(f2);

        f2e = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6.3 11.2]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile;
        plot_condition_heatmap(best_ME_err_map, [0 max_err], parula(256), nbins, chance_err, ...
            sprintf('%s BEST ME decoder error', area), 'Error (deg)', 'error');

        nexttile;
        plot_condition_heatmap(best_SE_err_map, [0 max_err], parula(256), nbins, chance_err, ...
            sprintf('%s BEST SE-add decoder error', area), 'Error (deg)', 'error');

        outpdf_best_err = fullfile(plotdir, sprintf('%s_MAP_ConditionError_6x6_BestSession.pdf', area));
        print(f2e, outpdf_best_err, '-dpdf', '-painters');
        close(f2e);
    end

    % ---------- bar summaries ----------
    acc_me_means = [mean(sess_cong_ME,'omitnan'), mean(sess_incong_ME,'omitnan')];
    acc_me_sems  = [sem_local(sess_cong_ME), sem_local(sess_incong_ME)];
    acc_se_means = [mean(sess_cong_SE,'omitnan'), mean(sess_incong_SE,'omitnan')];
    acc_se_sems  = [sem_local(sess_cong_SE), sem_local(sess_incong_SE)];

    err_me_means = [mean(sess_err_cong_ME,'omitnan'), mean(sess_err_incong_ME,'omitnan')];
    err_me_sems  = [sem_local(sess_err_cong_ME), sem_local(sess_err_incong_ME)];
    err_se_means = [mean(sess_err_cong_SE,'omitnan'), mean(sess_err_incong_SE,'omitnan')];
    err_se_sems  = [sem_local(sess_err_cong_SE), sem_local(sess_err_incong_SE)];

    Summary.barplot_accuracy.ME.means = acc_me_means;
    Summary.barplot_accuracy.ME.sems  = acc_me_sems;
    Summary.barplot_accuracy.SE.means = acc_se_means;
    Summary.barplot_accuracy.SE.sems  = acc_se_sems;
    Summary.barplot_error_deg.ME.means = err_me_means;
    Summary.barplot_error_deg.ME.sems  = err_me_sems;
    Summary.barplot_error_deg.SE.means = err_se_means;
    Summary.barplot_error_deg.SE.sems  = err_se_sems;
    Summary.barplot_labels = {'Congruent','Incongruent'};

    % ---------- split accuracy bars ----------
    f3 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot_cong_incong_bar(acc_me_means, acc_me_sems, chance_acc, ...
        sprintf('%s ME accuracy', area), 'Joint decoding accuracy', [0 1], 'higher');

    nexttile;
    plot_cong_incong_bar(acc_se_means, acc_se_sems, chance_acc, ...
        sprintf('%s SE-add accuracy', area), 'Joint decoding accuracy', [0 1], 'higher');

    outpdf_bar_acc = fullfile(plotdir, sprintf('%s_MAP_Congruent_Incongruent_Accuracy_Bar_Split_ME_SE.pdf', area));
    print(f3, outpdf_bar_acc, '-dpdf', '-painters');
    close(f3);

    % ---------- split error bars ----------
    ytop_err = max([err_me_means + err_me_sems, err_se_means + err_se_sems, chance_err], [], 'omitnan');
    if ~isfinite(ytop_err), ytop_err = 180; end
    ytop_err = max(180, ytop_err * 1.15);

    f4 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot_cong_incong_bar(err_me_means, err_me_sems, chance_err, ...
        sprintf('%s ME prediction error', area), 'Mean prediction error (deg)', [0 ytop_err], 'lower');

    nexttile;
    plot_cong_incong_bar(err_se_means, err_se_sems, chance_err, ...
        sprintf('%s SE-add prediction error', area), 'Mean prediction error (deg)', [0 ytop_err], 'lower');

    outpdf_bar_err = fullfile(plotdir, sprintf('%s_MAP_Congruent_Incongruent_Error_Bar_Split_ME_SE.pdf', area));
    print(f4, outpdf_bar_err, '-dpdf', '-painters');
    close(f4);

    % ---------- NEW: single-neuron congruent vs incongruent activity scatter ----------
    outpdf_activity = "";
    if ~isempty(all_neuron_FR_cong)
        f5 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 5.2 5.2]);
        hold on; box on; grid on;

        scatter(all_neuron_FR_cong, all_neuron_FR_incong, 22, ...
            'filled', ...
            'MarkerFaceColor', [0.35 0.35 0.35], ...
            'MarkerFaceAlpha', 0.45, ...
            'MarkerEdgeColor', 'none');

        mn = min([all_neuron_FR_cong(:); all_neuron_FR_incong(:)], [], 'omitnan');
        mx = max([all_neuron_FR_cong(:); all_neuron_FR_incong(:)], [], 'omitnan');
        if ~isfinite(mn), mn = 0; end
        if ~isfinite(mx), mx = 1; end
        pad = 0.05 * max(mx - mn, eps);
        lims = [max(0, mn - pad), mx + pad];

        plot(lims, lims, 'k--', 'LineWidth', 1.3);
        xlim(lims); ylim(lims);
        axis square;

        xlabel('Mean FR on congruent ME trials');
        ylabel('Mean FR on incongruent ME trials');
        title(sprintf('%s single-neuron activity: congruent vs incongruent', area), ...
            'Interpreter','none');

        delta = all_neuron_FR_cong - all_neuron_FR_incong;
        nAbove = sum(all_neuron_FR_incong > all_neuron_FR_cong);
        nBelow = sum(all_neuron_FR_incong < all_neuron_FR_cong);
        nTotal = numel(delta);

        text(0.05, 0.95, sprintf('n = %d neurons', nTotal), ...
            'Units','normalized', 'FontSize', 10, 'FontWeight','bold', ...
            'BackgroundColor','w', 'Margin', 2);
        text(0.05, 0.88, sprintf('mean(cong - incong) = %.2f', mean(delta,'omitnan')), ...
            'Units','normalized', 'FontSize', 10, ...
            'BackgroundColor','w', 'Margin', 2);
        text(0.05, 0.81, sprintf('above unity = %d, below = %d', nAbove, nBelow), ...
            'Units','normalized', 'FontSize', 10, ...
            'BackgroundColor','w', 'Margin', 2);

        outpdf_activity = fullfile(plotdir, sprintf('%s_ME_NeuronMeanFR_Congruent_vs_Incongruent.pdf', area));
        print(f5, outpdf_activity, '-dpdf', '-painters');
        close(f5);
    else
        warning('No single-neuron congruent/incongruent activity values were available for plotting.');
    end

    % ---------- save ----------
    outmat = fullfile(modeldir, sprintf('%s_MAP_ConditionAccuracy_6x6_Summary.mat', area));
    save(outmat, 'Summary', '-v7.3');

    fprintf('\nSaved condition accuracy/error summary:\n  %s\n', outmat);
    fprintf('Saved accuracy heatmap plot:\n  %s\n', outpdf_acc);
    fprintf('Saved error heatmap plot:\n  %s\n', outpdf_err);
    fprintf('Saved split accuracy bar plot:\n  %s\n', outpdf_bar_acc);
    fprintf('Saved split error bar plot:\n  %s\n', outpdf_bar_err);
    if strlength(outpdf_activity) > 0
        fprintf('Saved single-neuron activity plot:\n  %s\n', outpdf_activity);
    end

    fprintf('\n[%s] Chance levels:\n', area);
    fprintf('  Joint accuracy chance: %.2f%%\n', 100*chance_acc);
    fprintf('  Mean circular prediction error chance: %.1f deg\n', chance_err);

    fprintf('\n[%s] Across-session pooled joint decoding accuracy:\n', area);
    fprintf('  ME decoder congruent:    %.1f%%\n', 100*Summary.congruent_ME);
    fprintf('  ME decoder incongruent:  %.1f%%\n', 100*Summary.incongruent_ME);
    fprintf('  SE decoder congruent:    %.1f%%\n', 100*Summary.congruent_SE);
    fprintf('  SE decoder incongruent:  %.1f%%\n', 100*Summary.incongruent_SE);

    fprintf('\n[%s] Across-session pooled mean prediction error:\n', area);
    fprintf('  ME decoder congruent:    %.1f deg\n', Summary.err_congruent_ME_deg);
    fprintf('  ME decoder incongruent:  %.1f deg\n', Summary.err_incongruent_ME_deg);
    fprintf('  SE decoder congruent:    %.1f deg\n', Summary.err_congruent_SE_deg);
    fprintf('  SE decoder incongruent:  %.1f deg\n', Summary.err_incongruent_SE_deg);

    if ~isempty(all_neuron_FR_cong)
        fprintf('\n[%s] Single-neuron congruent/incongruent ME activity:\n', area);
        fprintf('  n neurons: %d\n', numel(all_neuron_FR_cong));
        fprintf('  mean congruent FR:   %.2f\n', mean(all_neuron_FR_cong, 'omitnan'));
        fprintf('  mean incongruent FR: %.2f\n', mean(all_neuron_FR_incong, 'omitnan'));
        fprintf('  mean cong - incong:  %.2f\n', mean(all_neuron_FR_cong - all_neuron_FR_incong, 'omitnan'));
    end

    if ~isempty(skipped)
        fprintf('\nSkipped %d files because they lacked joint predictions or were invalid.\n', numel(skipped));
    end
end

% ============================================================
% helpers
% ============================================================

function A = safe_divide(num, den)
    A = nan(size(num));
    ok = den > 0 & isfinite(den);
    A(ok) = num(ok) ./ den(ok);
end

function err_deg = joint_circular_bin_error_deg(Htrue, Etrue, Hhat, Ehat, nbins, bin_deg)
    dH = abs(Hhat(:) - Htrue(:));
    dE = abs(Ehat(:) - Etrue(:));

    dH = min(dH, nbins - dH);
    dE = min(dE, nbins - dE);

    errH_deg = dH * bin_deg;
    errE_deg = dE * bin_deg;

    err_deg = 0.5 * (errH_deg + errE_deg);
end

function e = expected_joint_chance_error_deg(nbins, bin_deg)
    true_bin = 1;
    pred_bins = (1:nbins)';
    d = abs(pred_bins - true_bin);
    d = min(d, nbins - d);
    e = mean(d * bin_deg, 'omitnan');
end

function s = sem_local(x)
    x = x(:);
    x = x(isfinite(x));
    if numel(x) <= 1
        s = NaN;
    else
        s = std(x, 0) ./ sqrt(numel(x));
    end
end

function plot_condition_heatmap(M, climVals, cmap, nbins, chance_val, title_text, cb_label, metric_type)
    ax = gca;
    imagesc(ax, M, climVals);
    axis(ax, 'xy');
    axis(ax, 'square');
    colormap(ax, cmap);

    title(ax, title_text, 'Interpreter','none');
    xlabel(ax, 'True eye target');
    ylabel(ax, 'True hand target');
    set(ax, 'XTick', 1:nbins, 'YTick', 1:nbins, 'FontSize', 15);
    hold(ax, 'on');

    cb = colorbar(ax);
    ylabel(cb, cb_label);
    add_colorbar_chance_tick(cb, climVals, chance_val, metric_type);

    cong_color = [0 0 0];
    incong_color = [0.55 0.55 0.55];
    anti_color = [0.65 0 0.85];

    for h = 1:nbins
        for e = 1:nbins
            if h == e
                rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
                    'EdgeColor', cong_color, 'LineWidth', 1.45, 'LineStyle', '-');
            else
                rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
                    'EdgeColor', incong_color, 'LineWidth', 0.55, 'LineStyle', ':');
            end
        end
    end

    if mod(nbins,2) == 0
        half = nbins / 2;
        for h = 1:nbins
            e = mod(h - 1 + half, nbins) + 1;
            rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
                'EdgeColor', anti_color, 'LineWidth', 1.65, 'LineStyle', '--');
        end
    end

    x_left = 0.72;
    x_right = 1.18;
    y_top = nbins - 0.18;
    dy = 0.38;

    plot(ax, [x_left x_right], [y_top y_top], '-', ...
        'Color', cong_color, 'LineWidth', 1.45);
    text(ax, x_right + 0.08, y_top, 'congruent', ...
        'FontSize', 9, 'FontWeight','bold', ...
        'VerticalAlignment','middle', ...
        'BackgroundColor','w', 'Margin', 1);

    plot(ax, [x_left x_right], [y_top-dy y_top-dy], ':', ...
        'Color', incong_color, 'LineWidth', 1.0);
    text(ax, x_right + 0.08, y_top-dy, 'incongruent', ...
        'FontSize', 9, ...
        'VerticalAlignment','middle', ...
        'BackgroundColor','w', 'Margin', 1);

    if mod(nbins,2) == 0
        plot(ax, [x_left x_right], [y_top-2*dy y_top-2*dy], '--', ...
            'Color', anti_color, 'LineWidth', 1.65);
        text(ax, x_right + 0.08, y_top-2*dy, 'anti-parallel', ...
            'FontSize', 9, ...
            'VerticalAlignment','middle', ...
            'BackgroundColor','w', 'Margin', 1);
    end
end

function add_colorbar_chance_tick(cb, climVals, chance_val, metric_type)
    if ~isfinite(chance_val), return; end

    lo = climVals(1);
    hi = climVals(2);
    if chance_val < lo || chance_val > hi || hi <= lo
        return;
    end

    old_ticks = cb.Ticks(:)';
    ticks = unique([old_ticks, chance_val]);
    ticks = ticks(ticks >= lo & ticks <= hi);
    ticks = sort(ticks);

    labels = cell(size(ticks));
    for i = 1:numel(ticks)
        if abs(ticks(i) - chance_val) < 1e-10
            if strcmpi(metric_type, 'accuracy')
                labels{i} = sprintf('chance %.1f%%', 100*chance_val);
            else
                labels{i} = sprintf('chance %.0f°', chance_val);
            end
        else
            if strcmpi(metric_type, 'accuracy')
                labels{i} = sprintf('%.2g', ticks(i));
            else
                labels{i} = sprintf('%.0f', ticks(i));
            end
        end
    end

    cb.Ticks = ticks;
    cb.TickLabels = labels;
end

function plot_cong_incong_bar(means, sems, chance_val, title_text, ylabel_text, ylimVals, better_direction)
    hold on; box on; grid on;

    bar(1:2, means, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'k');
    errorbar(1:2, means, sems, 'k.', 'LineWidth', 1.5, 'CapSize', 12);

    yline(chance_val, 'k--', 'LineWidth', 1.3);

    if strcmpi(better_direction, 'higher')
        chance_label = sprintf('chance = %.1f%%', 100*chance_val);
    else
        chance_label = sprintf('chance = %.0f°', chance_val);
    end

    text(1.5, chance_val, [' ' chance_label], ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'FontSize', 10, ...
        'FontWeight','bold', ...
        'BackgroundColor','w');

    set(gca, 'XTick', 1:2, 'XTickLabel', {'Congruent','Incongruent'}, 'FontSize', 13);
    ylabel(ylabel_text);
    ylim(ylimVals);
    title(title_text, 'Interpreter','none');

    if strcmpi(better_direction, 'higher')
        text(0.02, 0.95, 'higher is better', 'Units','normalized', ...
            'FontSize', 9, 'FontWeight','bold', 'BackgroundColor','w', 'Margin', 2);
    else
        text(0.02, 0.95, 'lower is better', 'Units','normalized', ...
            'FontSize', 9, 'FontWeight','bold', 'BackgroundColor','w', 'Margin', 2);
    end
end

function cmap = make_light_to_dark_green_colormap(n)
    if nargin < 1, n = 256; end
    t = linspace(0,1,n)';

    c_light = [0.97, 0.99, 0.97];
    c_dark  = [0.18, 0.60, 0.22];

    cmap = (1-t).*c_light + t.*c_dark;
    cmap = max(0, min(1, cmap));
end
