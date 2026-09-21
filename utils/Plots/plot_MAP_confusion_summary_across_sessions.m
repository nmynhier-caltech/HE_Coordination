function plot_MAP_confusion_summary_across_sessions(area, savepaths, varargin)
% Plot MAP decoding summary across sessions.
%
% CHANGE FROM PREVIOUS VERSION:
%   The confusion heatmaps are NO LONGER averaged across sessions.
%   Instead, they are taken from the SINGLE RIGHTMOST session.
%
% Rightmost session definition:
%   The session with the largest full-population neuron count in
%   Results.acc_vs_neurons.neu_levels. Combined decoding accuracy is used
%   only to break a tie in neuron count.
%
% Robust to sessions that:
%   - save raw predictions (true_flat, idx_hat_ME/SE)
%   - only save 6x6 matrices (C_*_row or C_*_counts)
%
% Across-session plots still retained:
%   1) Final accuracy vs total neurons per session (scatter + fit, 2x extrap)
%   2) All-points neurons scatter (ME-only) + fit
%   3) All-points neurons scatter (SE-only) + fit
%   4) All-points trials scatter (ME-only) + fit
%   5) All-points trials scatter (SE-only) + fit
%
% Inputs
%   area      : string, e.g., 'M1' or 'SPL'
%   savepaths : cellstr of session root folders
%
% Name-Value
%   'NBins' : default 6 — target number of angular bins per effector
%
% Saves
%   plots/summary/decoding/<area>_MAP_Confusions_BestSession_HE.pdf
%   plots/summary/decoding/<area>_FinalAcc_vs_N_Scatter.pdf
%   plots/summary/decoding/<area>_AllPoints_Accuracy_vs_Neurons_ME.pdf
%   plots/summary/decoding/<area>_AllPoints_Accuracy_vs_Neurons_SE.pdf
%   plots/summary/decoding/<area>_AllPoints_Accuracy_vs_Trials_ME.pdf
%   plots/summary/decoding/<area>_AllPoints_Accuracy_vs_Trials_SE.pdf
%   models/<area>_MAP_Confusion_Summary.mat

    % ---------- args ----------
    p = inputParser;
    addParameter(p, 'NBins', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    parse(p, varargin{:});
    nbins_target = p.Results.NBins;

    % ---------- accumulators for non-heatmap summary/plots ----------
    % Neurons — all group points (k neurons, accuracy)
    allpts_neur_ME_H = [];  allpts_neur_ME_E = [];
    allpts_neur_SE_H = [];  allpts_neur_SE_E = [];

    % Neurons — final acc vs session N
    sess_N = [];  sess_ME_H = []; sess_ME_E = []; sess_SE_H = []; sess_SE_E = [];

    % Trials — raw pairs (count, accuracy), separate per model+effector
    raw_trials_ME_H = [];  raw_trials_ME_E = [];
    raw_trials_SE_H = [];  raw_trials_SE_E = [];

    used = strings(0,1);

    % ---------- RIGHTMOST SESSION tracking ----------
    bestN = -Inf;
    bestScore = -Inf;
    bestFile = "";
    bestSessionIdx = NaN;

    best_C_ME_H_row = [];
    best_C_SE_H_row = [];
    best_C_ME_E_row = [];
    best_C_SE_E_row = [];

    best_acc_ME_H = NaN;
    best_acc_SE_H = NaN;
    best_acc_ME_E = NaN;
    best_acc_SE_E = NaN;

    % ---------- per-session loop ----------
    for d = 1:numel(savepaths)
        sp = savepaths{d};

        f = fullfile(sp, 'models', sprintf('%s_MAP_Smoothed.mat', area));
        if ~exist(f, 'file')
            warning('Missing file: %s', f);
            continue;
        end

        S = load(f);
        if ~isfield(S,'Results')
            warning('No Results in %s', f);
            continue;
        end
        R = S.Results;

        % --- determine session nbins ---
        nbins_sess = nbins_target;
        if isfield(R,'theta_centers') && ~isempty(R.theta_centers)
            nbins_sess = numel(R.theta_centers);
        elseif all(isfield(R,{'true_flat','idx_hat_ME','idx_hat_SE'}))
            max_state = max([R.true_flat(:); R.idx_hat_ME(:); R.idx_hat_SE(:)], [], 'omitnan');
            if isfinite(max_state) && max_state > 0
                k = round(sqrt(double(max_state)));
                if k >= 2
                    nbins_sess = k;
                end
            end
        elseif isfield(R,'C_ME_H_row')
            nbins_sess = size(R.C_ME_H_row,1);
        end

        % ============================================================
        % Build SESSION-SPECIFIC confusion matrices (for best-session selection)
        % ============================================================
        session_has_confusion = false;
        C_ME_H_row_sess = [];
        C_SE_H_row_sess = [];
        C_ME_E_row_sess = [];
        C_SE_E_row_sess = [];

        % --- preferred: reconstruct from raw predictions ---
        if all(isfield(R, {'true_flat','idx_hat_ME','idx_hat_SE'}))
            tf  = R.true_flat(:);
            im  = R.idx_hat_ME(:);
            is  = R.idx_hat_SE(:);

            [h_true_s, e_true_s] = ind2sub([nbins_sess nbins_sess], tf);
            [h_pred_ME_s, e_pred_ME_s] = ind2sub([nbins_sess nbins_sess], im);
            [h_pred_SE_s, e_pred_SE_s] = ind2sub([nbins_sess nbins_sess], is);

            m_ok = isfinite(h_true_s) & isfinite(e_true_s) & ...
                   isfinite(h_pred_ME_s) & isfinite(e_pred_ME_s) & ...
                   isfinite(h_pred_SE_s) & isfinite(e_pred_SE_s);

            h_true_s    = h_true_s(m_ok);
            e_true_s    = e_true_s(m_ok);
            h_pred_ME_s = h_pred_ME_s(m_ok);
            e_pred_ME_s = e_pred_ME_s(m_ok);
            h_pred_SE_s = h_pred_SE_s(m_ok);
            e_pred_SE_s = e_pred_SE_s(m_ok);

            if ~isempty(h_true_s)
                if nbins_sess ~= nbins_target
                    h_true    = remap_bins(h_true_s,    nbins_sess, nbins_target);
                    e_true    = remap_bins(e_true_s,    nbins_sess, nbins_target);
                    h_pred_ME = remap_bins(h_pred_ME_s, nbins_sess, nbins_target);
                    e_pred_ME = remap_bins(e_pred_ME_s, nbins_sess, nbins_target);
                    h_pred_SE = remap_bins(h_pred_SE_s, nbins_sess, nbins_target);
                    e_pred_SE = remap_bins(e_pred_SE_s, nbins_sess, nbins_target);
                else
                    h_true    = h_true_s;    e_true    = e_true_s;
                    h_pred_ME = h_pred_ME_s; e_pred_ME = e_pred_ME_s;
                    h_pred_SE = h_pred_SE_s; e_pred_SE = e_pred_SE_s;
                end

                C_ME_H_counts_sess = accumarray([h_true, h_pred_ME], 1, [nbins_target nbins_target], @sum, 0);
                C_SE_H_counts_sess = accumarray([h_true, h_pred_SE], 1, [nbins_target nbins_target], @sum, 0);
                C_ME_E_counts_sess = accumarray([e_true, e_pred_ME], 1, [nbins_target nbins_target], @sum, 0);
                C_SE_E_counts_sess = accumarray([e_true, e_pred_SE], 1, [nbins_target nbins_target], @sum, 0);

                C_ME_H_row_sess = rownorm_local(C_ME_H_counts_sess);
                C_SE_H_row_sess = rownorm_local(C_SE_H_counts_sess);
                C_ME_E_row_sess = rownorm_local(C_ME_E_counts_sess);
                C_SE_E_row_sess = rownorm_local(C_SE_E_counts_sess);

                session_has_confusion = true;
            end

        % --- next best: raw count matrices ---
        elseif isfield(R,'C_ME_H_counts') && isfield(R,'C_SE_H_counts') && ...
               isfield(R,'C_ME_E_counts') && isfield(R,'C_SE_E_counts')

            C_ME_H_counts_sess = resize_confusion(R.C_ME_H_counts, nbins_target, false);
            C_SE_H_counts_sess = resize_confusion(R.C_SE_H_counts, nbins_target, false);
            C_ME_E_counts_sess = resize_confusion(R.C_ME_E_counts, nbins_target, false);
            C_SE_E_counts_sess = resize_confusion(R.C_SE_E_counts, nbins_target, false);

            C_ME_H_row_sess = rownorm_local(C_ME_H_counts_sess);
            C_SE_H_row_sess = rownorm_local(C_SE_H_counts_sess);
            C_ME_E_row_sess = rownorm_local(C_ME_E_counts_sess);
            C_SE_E_row_sess = rownorm_local(C_SE_E_counts_sess);

            session_has_confusion = true;

        % --- fallback: row-normalized matrices directly ---
        elseif isfield(R,'C_ME_H_row') && isfield(R,'C_SE_H_row') && ...
               isfield(R,'C_ME_E_row') && isfield(R,'C_SE_E_row')

            C_ME_H_row_sess = resize_confusion(R.C_ME_H_row, nbins_target, true);
            C_SE_H_row_sess = resize_confusion(R.C_SE_H_row, nbins_target, true);
            C_ME_E_row_sess = resize_confusion(R.C_ME_E_row, nbins_target, true);
            C_SE_E_row_sess = resize_confusion(R.C_SE_E_row, nbins_target, true);

            session_has_confusion = true;
        end

        % Session-specific confusion accuracies
        acc_ME_H_sess = NaN; acc_SE_H_sess = NaN; acc_ME_E_sess = NaN; acc_SE_E_sess = NaN;
        if session_has_confusion
            acc_ME_H_sess = mean(diag(C_ME_H_row_sess), 'omitnan');
            acc_SE_H_sess = mean(diag(C_SE_H_row_sess), 'omitnan');
            acc_ME_E_sess = mean(diag(C_ME_E_row_sess), 'omitnan');
            acc_SE_E_sess = mean(diag(C_SE_E_row_sess), 'omitnan');
        end

        % ============================================================
        % Session neuron count for choosing the RIGHTMOST SESSION
        % ============================================================
        sessionN = -Inf;
        sessionScore = -Inf;

        if isfield(R,'acc_vs_neurons')
            B = R.acc_vs_neurons;
            needB = {'neu_levels','ME_H','ME_E','SE_H','SE_E'};
            if all(isfield(B, needB)) && ~isempty(B.neu_levels)
                sessionN = B.neu_levels(end);
                vals = [B.ME_H(end), B.SE_H(end), B.ME_E(end), B.SE_E(end)];
                if any(isfinite(vals))
                    sessionScore = mean(vals, 'omitnan');
                end
            end
        end

        % fallback to confusion-derived score
        if ~isfinite(sessionScore) && session_has_confusion
            sessionScore = mean([acc_ME_H_sess, acc_SE_H_sess, acc_ME_E_sess, acc_SE_E_sess], 'omitnan');
        end

        % Update rightmost session; use accuracy only to break an N tie.
        is_rightmost = sessionN > bestN || ...
            (sessionN == bestN && sessionScore > bestScore);
        if isfinite(sessionN) && is_rightmost && session_has_confusion
            bestN = sessionN;
            bestScore = sessionScore;
            bestFile = string(f);
            bestSessionIdx = d;

            best_C_ME_H_row = C_ME_H_row_sess;
            best_C_SE_H_row = C_SE_H_row_sess;
            best_C_ME_E_row = C_ME_E_row_sess;
            best_C_SE_E_row = C_SE_E_row_sess;

            best_acc_ME_H = acc_ME_H_sess;
            best_acc_SE_H = acc_SE_H_sess;
            best_acc_ME_E = acc_ME_E_sess;
            best_acc_SE_E = acc_SE_E_sess;
        end

        % ============================================================
        % Collect plotting data across sessions (unchanged)
        % ============================================================
        if isfield(R,'acc_vs_neurons')
            B = R.acc_vs_neurons;
            needB = {'neu_levels','ME_H','ME_E','SE_H','SE_E'};
            if all(isfield(B, needB))
                xN = B.neu_levels(:);

                allpts_neur_ME_H = [allpts_neur_ME_H; [xN, B.ME_H(:)]]; %#ok<AGROW>
                allpts_neur_ME_E = [allpts_neur_ME_E; [xN, B.ME_E(:)]]; %#ok<AGROW>
                allpts_neur_SE_H = [allpts_neur_SE_H; [xN, B.SE_H(:)]]; %#ok<AGROW>
                allpts_neur_SE_E = [allpts_neur_SE_E; [xN, B.SE_E(:)]]; %#ok<AGROW>

                if ~isempty(xN)
                    sess_N(end+1,1)    = xN(end); %#ok<AGROW>
                    sess_ME_H(end+1,1) = B.ME_H(end); %#ok<AGROW>
                    sess_ME_E(end+1,1) = B.ME_E(end); %#ok<AGROW>
                    sess_SE_H(end+1,1) = B.SE_H(end); %#ok<AGROW>
                    sess_SE_E(end+1,1) = B.SE_E(end); %#ok<AGROW>
                end
            end
        end

        if isfield(R,'acc_vs_trials')
            A = R.acc_vs_trials;
            needA = {'x_me_trials','x_se_trials','ME_H','ME_E','SE_H','SE_E'};
            if all(isfield(A, needA))
                xm  = A.x_me_trials(:);
                xs  = A.x_se_trials(:);
                ymh = A.ME_H(:);
                yme = A.ME_E(:);
                ysh = A.SE_H(:);
                yse = A.SE_E(:);

                if ~isempty(xm)
                    raw_trials_ME_H = [raw_trials_ME_H; [xm, ymh]]; %#ok<AGROW>
                    raw_trials_ME_E = [raw_trials_ME_E; [xm, yme]]; %#ok<AGROW>
                end
                if ~isempty(xs)
                    raw_trials_SE_H = [raw_trials_SE_H; [xs, ysh]]; %#ok<AGROW>
                    raw_trials_SE_E = [raw_trials_SE_E; [xs, yse]]; %#ok<AGROW>
                end
            end
        end

        used(end+1,1) = string(f);
    end

    % ---------- output dir ----------
    plotdir = fullfile(savepaths{1}, '..', 'plots', 'summary', 'decoding');
    if ~exist(plotdir, 'dir')
        mkdir(plotdir);
    end

    % ---------- require best session ----------
    if isempty(best_C_ME_H_row) || isempty(best_C_SE_H_row) || isempty(best_C_ME_E_row) || isempty(best_C_SE_E_row)
        error('No valid session-level confusion matrices found. Check your session files and NBins setting.');
    end

    % ---------- green colormap for heatmaps (light -> dark) ----------
    cmap_green = make_light_to_dark_green_colormap(256);

    % ---------- Confusion plots: RIGHTMOST SESSION ONLY ----------
    fprintf('\n[%s] Rightmost session chosen for heatmaps (N=%d):\n', area, round(bestN));
    fprintf('  File : %s\n', bestFile);
    mean_acc = [mean(sess_ME_H, 'omitnan'), mean(sess_SE_H, 'omitnan'), ...
                mean(sess_ME_E, 'omitnan'), mean(sess_SE_E, 'omitnan')];
    accuracy_table = array2table(100 * [best_acc_ME_H, best_acc_SE_H, best_acc_ME_E, best_acc_SE_E, bestScore; ...
                                       mean_acc, mean(mean_acc, 'omitnan')], ...
        'RowNames', {'Rightmost session', 'Across-session mean'}, ...
        'VariableNames', {'ME_Hand_pct', 'SE_Hand_pct', 'ME_Eye_pct', 'SE_Eye_pct', 'Overall_pct'});
    disp(accuracy_table);

    f = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 8]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(best_C_ME_H_row,[0 1]); axis xy; axis square;
    colormap(gca, cmap_green); colorbar;
    title(sprintf('%s — RIGHTMOST session ME: Hand (acc=%.1f%%)', area, 100*best_acc_ME_H));
    set(gca,'FontSize',20);
    xlabel('Pred H bin'); ylabel('True H bin');

    nexttile;
    imagesc(best_C_SE_H_row,[0 1]); axis xy; axis square;
    colormap(gca, cmap_green); colorbar;
    title(sprintf('%s — RIGHTMOST session SE-add: Hand (acc=%.1f%%)', area, 100*best_acc_SE_H));
    set(gca,'FontSize',20);
    xlabel('Pred H bin'); ylabel('True H bin');

    nexttile;
    imagesc(best_C_ME_E_row,[0 1]); axis xy; axis square;
    colormap(gca, cmap_green); colorbar;
    title(sprintf('%s — RIGHTMOST session ME: Eye (acc=%.1f%%)', area, 100*best_acc_ME_E));
    set(gca,'FontSize',20);
    xlabel('Pred E bin'); ylabel('True E bin');

    nexttile;
    imagesc(best_C_SE_E_row,[0 1]); axis xy; axis square;
    colormap(gca, cmap_green); colorbar;
    title(sprintf('%s — RIGHTMOST session SE-add: Eye (acc=%.1f%%)', area, 100*best_acc_SE_E));
    set(gca,'FontSize',20);
    xlabel('Pred E bin'); ylabel('True E bin');

    print(f, fullfile(plotdir, sprintf('%s_MAP_Confusions_BestSession_HE.pdf', area)), '-dpdf','-painters');
    close(f);

    % ---------- Metric 1: mean neurons per session ----------
    mean_session_N = mean(sess_N, 'omitnan');

    % ---------- Metric 2: RIGHTMOST points (max N) per model+effector ----------
    [ME_H_rightN, ME_H_rightAcc] = rightmost_point(allpts_neur_ME_H);
    [ME_E_rightN, ME_E_rightAcc] = rightmost_point(allpts_neur_ME_E);
    [SE_H_rightN, SE_H_rightAcc] = rightmost_point(allpts_neur_SE_H);
    [SE_E_rightN, SE_E_rightAcc] = rightmost_point(allpts_neur_SE_E);

    % ---------- Plot 1: final accuracy vs total neurons (per session) ----------
    if ~isempty(sess_N)
        xext = linspace(0, max(sess_N)*2, 300);
        [xf_meh, yf_meh] = fit_saturating(sess_N, sess_ME_H, xext);
        [xf_seh, yf_seh] = fit_saturating(sess_N, sess_SE_H, xext);
        [xf_mee, yf_mee] = fit_saturating(sess_N, sess_ME_E, xext);
        [xf_see, yf_see] = fit_saturating(sess_N, sess_SE_E, xext);

        sfig = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 10]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        scatter(sess_N, sess_ME_H, 40, 'r', 'filled');
        scatter(sess_N, sess_SE_H, 40, 'b', 'filled');
        plot(xf_meh, yf_meh, 'r-','LineWidth',1.4);
        plot(xf_seh, yf_seh, 'b-','LineWidth',1.4);
        if isfinite(mean_session_N), xline(mean_session_N, 'k--','LineWidth',1.2); end
        axis square;
        xlabel('# Neurons (per session)');
        ylabel('Final Accuracy (H)');
        ylim([0 1]);
        title(sprintf('%s: Final accuracy vs N (per session)', area));
        legend({'ME','SE','ME fit','SE fit','mean N'},'Location','southeast');

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        scatter(sess_N, sess_ME_E, 40, 'r', 'filled');
        scatter(sess_N, sess_SE_E, 40, 'b', 'filled');
        plot(xf_mee, yf_mee, 'r-','LineWidth',1.4);
        plot(xf_see, yf_see, 'b-','LineWidth',1.4);
        if isfinite(mean_session_N), xline(mean_session_N, 'k--','LineWidth',1.2); end
        axis square;
        xlabel('# Neurons (per session)');
        ylabel('Final Accuracy (E)');
        ylim([0 1]);

        print(sfig, fullfile(plotdir, sprintf('%s_FinalAcc_vs_N_Scatter.pdf', area)), '-dpdf','-painters');
        close(sfig);
    end

    % ---------- Plot 2: ALL-POINTS neurons scatter (ME-only) ----------
    if ~isempty(allpts_neur_ME_H)
        x_me_h = allpts_neur_ME_H(:,1); y_me_h = allpts_neur_ME_H(:,2);
        x_me_e = allpts_neur_ME_E(:,1); y_me_e = allpts_neur_ME_E(:,2);
        xextH = linspace(0, 2*max([x_me_h; 1]), 300);
        xextE = linspace(0, 2*max([x_me_e; 1]), 300);
        [xf_meH, yf_meH] = fit_saturating(x_me_h, y_me_h, xextH);
        [xf_meE, yf_meE] = fit_saturating(x_me_e, y_me_e, xextE);

        fN_ME = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 10]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        scatter(x_me_h, y_me_h, 16, 'r', 'filled', 'MarkerFaceAlpha',0.3);
        plot(xf_meH, yf_meH, 'r-','LineWidth',1.3);
        hline(1/6, 'r-');
        if isfinite(ME_H_rightAcc) && isfinite(ME_H_rightN)
            plot(ME_H_rightN, ME_H_rightAcc, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
            xline(ME_H_rightN, 'k--', 'LineWidth',1.2);
            yline(ME_H_rightAcc,'k--', 'LineWidth',1.2);
            text(ME_H_rightN, ME_H_rightAcc, sprintf('  N=%d, %.1f%%', ME_H_rightN, 100*ME_H_rightAcc), ...
                'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','left');
        end
        axis square;
        xlabel('# Neurons used');
        ylabel('Accuracy (H)');
        ylim([0 1]);
        title(sprintf('%s: Accuracy vs Neurons (ME only, all points)', area));

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        scatter(x_me_e, y_me_e, 16, 'r', 'filled', 'MarkerFaceAlpha',0.3);
        plot(xf_meE, yf_meE, 'r-','LineWidth',1.3);
        hline(1/6, 'r-');
        if isfinite(ME_E_rightAcc) && isfinite(ME_E_rightN)
            plot(ME_E_rightN, ME_E_rightAcc, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
            xline(ME_E_rightN, 'k--', 'LineWidth',1.2);
            yline(ME_E_rightAcc,'k--', 'LineWidth',1.2);
            text(ME_E_rightN, ME_E_rightAcc, sprintf('  N=%d, %.1f%%', ME_E_rightN, 100*ME_E_rightAcc), ...
                'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','left');
        end
        axis square;
        xlabel('# Neurons used');
        ylabel('Accuracy (E)');
        ylim([0 1]);

        print(fN_ME, fullfile(plotdir, sprintf('%s_AllPoints_Accuracy_vs_Neurons_ME.pdf', area)), '-dpdf','-painters');
        close(fN_ME);
    end

    % ---------- Plot 2b: ALL-POINTS neurons scatter (SE-only) ----------
    if ~isempty(allpts_neur_SE_H)
        x_se_h = allpts_neur_SE_H(:,1); y_se_h = allpts_neur_SE_H(:,2);
        x_se_e = allpts_neur_SE_E(:,1); y_se_e = allpts_neur_SE_E(:,2);
        xextH = linspace(0, 2*max([x_se_h; 1]), 300);
        xextE = linspace(0, 2*max([x_se_e; 1]), 300);
        [xf_seH, yf_seH] = fit_saturating(x_se_h, y_se_h, xextH);
        [xf_seE, yf_seE] = fit_saturating(x_se_e, y_se_e, xextE);

        fN_SE = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 10]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        scatter(x_se_h, y_se_h, 16, 'b', 'filled', 'MarkerFaceAlpha',0.3);
        plot(xf_seH, yf_seH, 'b-','LineWidth',1.3);
        hline(1/6, 'b-');
        if isfinite(SE_H_rightAcc) && isfinite(SE_H_rightN)
            plot(SE_H_rightN, SE_H_rightAcc, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
            xline(SE_H_rightN, 'k--', 'LineWidth',1.2);
            yline(SE_H_rightAcc,'k--', 'LineWidth',1.2);
            text(SE_H_rightN, SE_H_rightAcc, sprintf('  N=%d, %.1f%%', SE_H_rightN, 100*SE_H_rightAcc), ...
                'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','left');
        end
        axis square;
        xlabel('# Neurons used');
        ylabel('Accuracy (H)');
        ylim([0 1]);
        title(sprintf('%s: Accuracy vs Neurons (SE only, all points)', area));

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        scatter(x_se_e, y_se_e, 16, 'b', 'filled', 'MarkerFaceAlpha',0.3);
        plot(xf_seE, yf_seE, 'b-','LineWidth',1.3);
        hline(1/6, 'b-');
        if isfinite(SE_E_rightAcc) && isfinite(SE_E_rightN)
            plot(SE_E_rightN, SE_E_rightAcc, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
            xline(SE_E_rightN, 'k--', 'LineWidth',1.2);
            yline(SE_E_rightAcc,'k--', 'LineWidth',1.2);
            text(SE_E_rightN, SE_E_rightAcc, sprintf('  N=%d, %.1f%%', SE_E_rightN, 100*SE_E_rightAcc), ...
                'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','left');
        end
        axis square;
        xlabel('# Neurons used');
        ylabel('Accuracy (E)');
        ylim([0 1]);

        print(fN_SE, fullfile(plotdir, sprintf('%s_AllPoints_Accuracy_vs_Neurons_SE.pdf', area)), '-dpdf','-painters');
        close(fN_SE);
    end

    % ---------- Plot 3: ALL-POINTS trials scatter (ME-only) ----------
    if ~isempty(raw_trials_ME_H) || ~isempty(raw_trials_ME_E)
        fT_ME = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 10]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        if ~isempty(raw_trials_ME_H)
            x_me = raw_trials_ME_H(:,1); y_me = raw_trials_ME_H(:,2);
            scatter(x_me, y_me, 16, 'r', 'filled', 'MarkerFaceAlpha',0.3);
            xmaxH = max([x_me; 1]);
            xextH = linspace(0, 2*xmaxH, 300);
            [xf_me, yf_me] = fit_saturating(x_me, y_me, xextH);
            plot(xf_me, yf_me, 'r-','LineWidth',1.3);
        end
        axis square;
        xlabel('# Training trials');
        ylabel('Accuracy (H)');
        ylim([0 1]);
        title(sprintf('%s: Accuracy vs Training Trials (ME only, all points)', area));

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        if ~isempty(raw_trials_ME_E)
            x_me = raw_trials_ME_E(:,1); y_me = raw_trials_ME_E(:,2);
            scatter(x_me, y_me, 16, 'r', 'filled', 'MarkerFaceAlpha',0.3);
            xmaxE = max([x_me; 1]);
            xextE = linspace(0, 2*xmaxE, 300);
            [xf_me, yf_me] = fit_saturating(x_me, y_me, xextE);
            plot(xf_me, yf_me, 'r-','LineWidth',1.3);
        end
        axis square;
        xlabel('# Training trials');
        ylabel('Accuracy (E)');
        ylim([0 1]);

        print(fT_ME, fullfile(plotdir, sprintf('%s_AllPoints_Accuracy_vs_Trials_ME.pdf', area)), '-dpdf','-painters');
        close(fT_ME);
    end

    % ---------- Plot 3b: ALL-POINTS trials scatter (SE-only) ----------
    if ~isempty(raw_trials_SE_H) || ~isempty(raw_trials_SE_E)
        fT_SE = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 10]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        if ~isempty(raw_trials_SE_H)
            x_se = raw_trials_SE_H(:,1); y_se = raw_trials_SE_H(:,2);
            scatter(x_se, y_se, 16, 'b', 'filled', 'MarkerFaceAlpha',0.3);
            xmaxH = max([x_se; 1]);
            xextH = linspace(0, 2*xmaxH, 300);
            [xf_se, yf_se] = fit_saturating(x_se, y_se, xextH);
            plot(xf_se, yf_se, 'b-','LineWidth',1.3);
        end
        axis square;
        xlabel('# Training trials');
        ylabel('Accuracy (H)');
        ylim([0 1]);
        title(sprintf('%s: Accuracy vs Training Trials (SE only, all points)', area));

        nexttile; hold on; grid on; box on; set(gca,'FontSize',20);
        if ~isempty(raw_trials_SE_E)
            x_se = raw_trials_SE_E(:,1); y_se = raw_trials_SE_E(:,2);
            scatter(x_se, y_se, 16, 'b', 'filled', 'MarkerFaceAlpha',0.3);
            xmaxE = max([x_se; 1]);
            xextE = linspace(0, 2*xmaxE, 300);
            [xf_se, yf_se] = fit_saturating(x_se, y_se, xextE);
            plot(xf_se, yf_se, 'b-','LineWidth',1.3);
        end
        axis square;
        xlabel('# Training trials');
        ylabel('Accuracy (E)');
        ylim([0 1]);

        print(fT_SE, fullfile(plotdir, sprintf('%s_AllPoints_Accuracy_vs_Trials_SE.pdf', area)), '-dpdf','-painters');
        close(fT_SE);
    end

    % ---------- Save numeric summary ----------
    Summary = struct();
    Summary.area = area;
    Summary.usedFiles = used;

    % BEST-session confusion matrices
    Summary.bestSession = struct();
    Summary.bestSession.file = bestFile;
    Summary.bestSession.sessionIndex = bestSessionIdx;
    Summary.bestSession.N = bestN;
    Summary.bestSession.score = bestScore;

    Summary.C_ME_H_row = best_C_ME_H_row;
    Summary.C_SE_H_row = best_C_SE_H_row;
    Summary.C_ME_E_row = best_C_ME_E_row;
    Summary.C_SE_E_row = best_C_SE_E_row;

    Summary.acc_ME_H = best_acc_ME_H;
    Summary.acc_SE_H = best_acc_SE_H;
    Summary.acc_ME_E = best_acc_ME_E;
    Summary.acc_SE_E = best_acc_SE_E;

    Summary.mean_session_N = mean_session_N;
    Summary.rightmostPoint = struct( ...
        'ME_H', struct('N', ME_H_rightN, 'acc', ME_H_rightAcc), ...
        'ME_E', struct('N', ME_E_rightN, 'acc', ME_E_rightAcc), ...
        'SE_H', struct('N', SE_H_rightN, 'acc', SE_H_rightAcc), ...
        'SE_E', struct('N', SE_E_rightN, 'acc', SE_E_rightAcc) );

    if ~isempty(sess_N)
        Summary.FinalAccVsN = table(sess_N, sess_ME_H, sess_SE_H, sess_ME_E, sess_SE_E, ...
            'VariableNames', {'N','ME_H','SE_H','ME_E','SE_E'});
    end

    if ~isempty(allpts_neur_ME_H)
        Summary.AllPointsNeurons = struct( ...
            'ME_H', allpts_neur_ME_H, 'ME_E', allpts_neur_ME_E, ...
            'SE_H', allpts_neur_SE_H, 'SE_E', allpts_neur_SE_E );
    end

    if ~isempty(raw_trials_ME_H) || ~isempty(raw_trials_SE_H) || ~isempty(raw_trials_ME_E) || ~isempty(raw_trials_SE_E)
        Summary.AllPointsTrials = struct( ...
            'ME_H', raw_trials_ME_H, 'ME_E', raw_trials_ME_E, ...
            'SE_H', raw_trials_SE_H, 'SE_E', raw_trials_SE_E );
    end

    outmat = fullfile(savepaths{1}, 'models', sprintf('%s_MAP_Confusion_Summary.mat', area));
    save(outmat, 'Summary', '-v7.3');
    fprintf('\nSaved: %s\n', outmat);
end

% ===== helpers =====
function idx2 = remap_bins(idx1, nb1, nb2)
    idx2 = 1 + floor((double(idx1)-1) * (double(nb2)/double(nb1)));
    idx2 = max(1, min(nb2, idx2));
end

function M2 = resize_confusion(M, nbins_target, isRowNorm)
    if nargin < 3, isRowNorm = false; end
    [r,c] = size(M);
    if r == nbins_target && c == nbins_target
        M2 = M;
        return;
    end

    map_r = remap_bins((1:r)', r, nbins_target);
    map_c = remap_bins((1:c)', c, nbins_target);

    M2 = zeros(nbins_target);
    for i = 1:r
        for j = 1:c
            M2(map_r(i), map_c(j)) = M2(map_r(i), map_c(j)) + M(i,j);
        end
    end

    if isRowNorm
        M2 = M2 ./ max(sum(M2,2), 1);
    end
end

function Mrow = rownorm_local(M)
    Mrow = M ./ max(sum(M,2), 1);
end

function [xfit, yfit] = fit_saturating(x, y, xfit)
    x = x(:);
    y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    if nargin < 3 || isempty(xfit)
        if numel(x) < 2
            xfit = x;
            yfit = y;
            return;
        end
        xfit = linspace(min(x), max(x), 200);
    end

    if numel(x) < 2
        yfit = nan(size(xfit));
        return;
    elseif numel(x) == 2
        yfit = interp1(x, y, xfit, 'linear', 'extrap');
        yfit = max(0, min(1, yfit));
        return;
    end

    y0 = max(0, min(1, y(1)));
    L0 = max(y);
    k0 = 1 / max(x);

    obj = @(p) sum(((p(1) - (p(1) - p(2))*exp(-p(3)*x)) - y).^2, 'omitnan');
    p = fminsearch(obj, [L0, y0, k0], optimset('Display','off'));

    yfit = p(1) - (p(1) - p(2))*exp(-p(3)*xfit);
    yfit = max(0, min(1, yfit));
end

function hline(y, ls)
    if nargin < 2, ls = 'k-'; end
    xl = xlim;
    plot(xl, [y y], ls, 'LineWidth', 1.0);
end

function [Nstar, accstar] = rightmost_point(X)
% X is [N, acc] rows. Returns N at maximum N (empirical), and the accuracy at that N.
% If multiple points share max N, take the maximum accuracy among them.
    Nstar = NaN;
    accstar = NaN;

    if isempty(X) || size(X,2) < 2
        return;
    end

    N = X(:,1);
    acc = X(:,2);
    m = isfinite(N) & isfinite(acc);
    N = N(m);
    acc = acc(m);

    if isempty(N)
        return;
    end

    Nmax = max(N);
    ii = (N == Nmax);
    acc_at = acc(ii);

    if isempty(acc_at)
        return;
    end

    Nstar = round(Nmax);
    accstar = max(acc_at);
end

function cmap = make_light_to_dark_green_colormap(n)
% Light-to-dark green colormap, matching the style of the attached figure.
    if nargin < 1, n = 256; end
    t = linspace(0,1,n)';

    c_light = [0.97, 0.99, 0.97];
    c_dark  = [0.18, 0.60, 0.22];

    cmap = (1-t).*c_light + t.*c_dark;
    cmap = max(0, min(1, cmap));
end
