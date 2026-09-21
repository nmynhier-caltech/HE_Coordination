%% added color coding of heatmaps
function plot_FR_amplitude_summary_across_sessions(area, savepaths)
% Plot FR amplitude summary across sessions using ONLY the median z-scores
% across bootstraps saved by the compute script (z_med_multi, z_med_single).
% No fallbacks are used.
%
% Requires fields in each session file:
%   - z_med_multi.H, z_med_multi.E
%   - z_med_single.H, z_med_single.E
%   - BOOT.MAXDIR_multi.H, BOOT.MAXDIR_multi.E
%   - BOOT.MAXDIR_single.H, BOOT.MAXDIR_single.E
%   - class_labels, class_labels_single
%
% Saves:
%   *_Amp_Selectivity_Across_Sessions_MULTI_vs_SINGLE.pdf
%   *_Amp_Tuning_Classification_Bar_MULTI_vs_SINGLE.pdf
%   *_Single_vs_Multi_Amplitude_HandEye_byClass_4col.pdf
%   *_Single_vs_Multi_TuningDir_HandEye_byClass_4col.pdf
%   *_Class_Agreement_and_Suppression.pdf
%   *_PreferredDirection_Distribution_SignificantOnly.pdf

    % -----------------------
    % Accumulators
    % -----------------------
    all_amp_H = [];
    all_amp_E = [];
    all_amp_single_H = [];
    all_amp_single_E = [];
    all_dir_multi_H  = [];
    all_dir_multi_E  = [];
    all_dir_single_H = [];
    all_dir_single_E = [];
    all_class_labels = [];
    all_class_labels_single = [];

    % NEW: accumulators for FR-by-dir draws (for amplitude mean/variance & CV)
    multi_FR_H_cells = {};   % each cell: [N x nbins x B]
    multi_FR_E_cells = {};
    single_FR_H_cells = {};
    single_FR_E_cells = {};
    missing_multi_sessions  = {};
    missing_single_sessions = {};
    % NEW: keep labels per session for class-conditioned summaries
    class_multi_cells  = {};
    class_single_cells = {};

    % -----------------------
    % Load per-session files (strict requirements)
    % -----------------------
    for d = 1:length(savepaths)
        savepath = savepaths{d};
        dat = load(fullfile(savepath, 'models', sprintf('%s_FR_amplitude_selectivity.mat', area)));

        % --- Hard requirements (no fallbacks)
        assert(isfield(dat, 'z_med_multi')  && isfield(dat.z_med_multi,  'H') && isfield(dat.z_med_multi,  'E'), ...
               'Missing z_med_multi.(H|E) in %s', savepath);
        assert(isfield(dat, 'z_med_single') && isfield(dat.z_med_single, 'H') && isfield(dat.z_med_single, 'E'), ...
               'Missing z_med_single.(H|E) in %s', savepath);
        assert(isfield(dat, 'BOOT') && isfield(dat.BOOT, 'MAXDIR_multi') && ...
               isfield(dat.BOOT.MAXDIR_multi, 'H') && isfield(dat.BOOT.MAXDIR_multi, 'E'), ...
               'Missing BOOT.MAXDIR_multi.(H|E) in %s', savepath);
        assert(isfield(dat, 'BOOT') && isfield(dat.BOOT, 'MAXDIR_single') && ...
               isfield(dat.BOOT.MAXDIR_single, 'H') && isfield(dat.BOOT.MAXDIR_single, 'E'), ...
               'Missing BOOT.MAXDIR_single.(H|E) in %s', savepath);
        assert(isfield(dat, 'class_labels'), 'Missing class_labels in %s', savepath);
        assert(isfield(dat, 'class_labels_single'), 'Missing class_labels_single in %s', savepath);

        % --- Pull median z across bootstraps
        aH  = dat.z_med_multi.H(:);
        aE  = dat.z_med_multi.E(:);
        aHs = dat.z_med_single.H(:);
        aEs = dat.z_med_single.E(:);

        % --- Circular means of PDs across bootstraps
        dmH = circ_mean_dim(dat.BOOT.MAXDIR_multi.H, 2);
        dmE = circ_mean_dim(dat.BOOT.MAXDIR_multi.E, 2);
        dsH = circ_mean_dim(dat.BOOT.MAXDIR_single.H, 2);
        dsE = circ_mean_dim(dat.BOOT.MAXDIR_single.E, 2);

        % --- Accumulate (for plotting)
        all_amp_H        = [all_amp_H;        aH(:)];
        all_amp_E        = [all_amp_E;        aE(:)];
        all_amp_single_H = [all_amp_single_H; aHs(:)];
        all_amp_single_E = [all_amp_single_E; aEs(:)];
        all_dir_multi_H  = [all_dir_multi_H;  dmH(:)];
        all_dir_multi_E  = [all_dir_multi_E;  dmE(:)];
        all_dir_single_H = [all_dir_single_H; dsH(:)];
        all_dir_single_E = [all_dir_single_E; dsE(:)];
        all_class_labels        = [all_class_labels;        dat.class_labels];
        all_class_labels_single = [all_class_labels_single; dat.class_labels_single];

        % --- Keep per-session class labels for later grouping
        class_multi_cells{end+1}  = dat.class_labels(:);         %#ok<AGROW>
        class_single_cells{end+1} = dat.class_labels_single(:);  %#ok<AGROW>

        % --- NEW: collect bootstrapped tuning curves (FR-by-direction draws)
        % ME (multi) expected at BOOT.FR_by_dir.(H|E) : [N x nbins x B]
        if isfield(dat, 'BOOT') && isfield(dat.BOOT, 'FR_by_dir') && ...
           isfield(dat.BOOT.FR_by_dir, 'H') && isfield(dat.BOOT.FR_by_dir, 'E')
            multi_FR_H_cells{end+1} = dat.BOOT.FR_by_dir.H; %#ok<AGROW>
            multi_FR_E_cells{end+1} = dat.BOOT.FR_by_dir.E; %#ok<AGROW>
        else
            missing_multi_sessions{end+1} = savepath; %#ok<AGROW>
        end

        % SE (single) expected at BOOT.FR_by_dir_single.(H|E) : [N x nbins x B]
        if isfield(dat, 'BOOT') && isfield(dat.BOOT, 'FR_by_dir_single') && ...
           isfield(dat.BOOT.FR_by_dir_single, 'H') && isfield(dat.BOOT.FR_by_dir_single, 'E')
            single_FR_H_cells{end+1} = dat.BOOT.FR_by_dir_single.H; %#ok<AGROW>
            single_FR_E_cells{end+1} = dat.BOOT.FR_by_dir_single.E; %#ok<AGROW>
        else
            missing_single_sessions{end+1} = savepath; %#ok<AGROW>
        end
    end

    %% Setup
    categories = ["H", "E", "Both", "None"];
    category_colors = containers.Map( ...
        {'H','E','Both','None'}, ...
        {[0.8, 0.2, 0.2], [0.2, 0.2, 0.8], [0.2, 0.6, 0.2], [0.5, 0.5, 0.5]} ...
    );

    plotdir = fullfile(savepaths{1}, '..', 'plots', 'summary', ...
        'selectivity_generalizability');
    if ~exist(plotdir, 'dir'), mkdir(plotdir); end

    class_labels        = all_class_labels;         % multi-based
    class_labels_single = all_class_labels_single;  % single-based

    amp_label_multi  = 'Amplitude (z) – Multi [median across bootstraps]';
    amp_label_single = 'Amplitude (z) – Single [median across bootstraps]';

    %% Plot 1: Amp(H) vs Amp(E): Multi (left) AND Single (right)
    h1 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);

    % ---- Left: Multi-effector
    subplot(1,2,1);
    set(gca, 'FontSize', 20); hold on;
    for i = 1:3  % Only 'H', 'E', 'Both' colored by MULTI class
        cat = char(categories(i));
        idx = strcmp(class_labels, cat);
        if any(idx)
            scatter(all_amp_H(idx), all_amp_E(idx), 30, ...
                'filled', 'MarkerFaceColor', category_colors(cat));
        end
    end
    xlabel(['Hand ' amp_label_multi]);
    ylabel(['Eye '  amp_label_multi]);
    xlim([0,15]); ylim([0,15]); axis square;
    title('Selectivity: Amp(H) vs Amp(E) – Multi');
    grid on; legend(categories(1:3), 'Location', 'northwest');

    % ---- Right: Single-effector (colored by MULTI class to keep same scheme)
    subplot(1,2,2);
    set(gca, 'FontSize', 20); hold on;
    for i = 1:3
        cat = char(categories(i));
        idx = strcmp(class_labels, cat);
        if any(idx)
            scatter(all_amp_single_H(idx), all_amp_single_E(idx), 30, ...
                'filled', 'MarkerFaceColor', category_colors(cat));
        end
    end
    xlabel(['Hand ' amp_label_single]);
    ylabel(['Eye '  amp_label_single]);
    xlim([0,15]); ylim([0,15]); axis square;
    title('Selectivity: Amp(H) vs Amp(E) – Single');
    grid on;

    saveas(h1, fullfile(plotdir, sprintf('%s_Amp_Selectivity_Across_Sessions_MULTI_vs_SINGLE.pdf', area)));

    %% Plot 2: Bar chart of classification counts – Multi (left) AND Single (right)
    h2 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);

    bar_colors = [
        0.8 0.2 0.2;  % H
        0.2 0.2 0.8;  % E
        0.2 0.6 0.2;  % Both
        0.5 0.5 0.5   % None
    ];

    % ---- Left: Multi counts
    subplot(1,2,1);
    set(gca, 'FontSize', 20);
    counts_multi = [sum(strcmp(class_labels, 'H')), ...
                    sum(strcmp(class_labels, 'E')), ...
                    sum(strcmp(class_labels, 'Both')), ...
                    sum(strcmp(class_labels, 'None'))];
    hold on;
    for i = 1:length(counts_multi)
        bar(i, counts_multi(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'k');
    end
    xticks(1:4); xticklabels(categories);
    ylabel('Neuron Count');
    title('Tuning Classification – Multi-effector');
    grid on; box on;

    ymax = max(counts_multi);
    ylim([0, max(1, ymax)*1.15]);
    dy = max(1, 0.03*ymax);
    for i = 1:length(counts_multi)
        text(i, counts_multi(i) + dy, sprintf('%d', counts_multi(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
    end

    % ---- Right: Single counts
    subplot(1,2,2);
    set(gca, 'FontSize', 20);
    counts_single = [sum(strcmp(class_labels_single, 'H')), ...
                     sum(strcmp(class_labels_single, 'E')), ...
                     sum(strcmp(class_labels_single, 'Both')), ...
                     sum(strcmp(class_labels_single, 'None'))];
    hold on;
    for i = 1:length(counts_single)
        bar(i, counts_single(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'k');
    end
    xticks(1:4); xticklabels(categories);
    ylabel('Neuron Count');
    title('Tuning Classification – Single-effector');
    grid on; box on;

    ymax2 = max(counts_single);
    ylim([0, max(1, ymax2)*1.15]);
    dy2 = max(1, 0.03*ymax2);
    for i = 1:length(counts_single)
        text(i, counts_single(i) + dy2, sprintf('%d', counts_single(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
    end

    saveas(h2, fullfile(plotdir, sprintf('%s_Amp_Tuning_Classification_Bar_MULTI_vs_SINGLE.pdf', area)));

    %% Plot 3: Single vs Multi Amplitudes split by class (2x4 grid)
    h3 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);

    col_classes = ["H","E","Both","Both∩Both"];
    class_masks = cell(1,4);
    class_masks{1} = strcmp(class_labels, 'H');
    class_masks{2} = strcmp(class_labels, 'E');
    class_masks{3} = strcmp(class_labels, 'Both');
    class_masks{4} = strcmp(class_labels, 'Both') & strcmp(class_labels_single, 'Both');

    col_colors = {category_colors('H'), category_colors('E'), category_colors('Both'), [0.1 0.6 0.1]};

    for c = 1:4
        mask = class_masks{c};
        clr  = col_colors{c};
        label_text = char(col_classes(c));

        % ----- Row 1: Hand -----
        subplot(2,4,c); hold on; grid on; set(gca,'FontSize',16);
        x = all_amp_single_H(mask);
        y = all_amp_H(mask);
        scatter(x, y, 25, 'filled', 'MarkerFaceColor', clr);
        xlim([-5,15]); ylim([-5,15]); axis square;
        xl = xlim; yl = ylim; mm = [min([xl yl]) max([xl yl])];
        plot(mm, mm, 'k-', 'LineWidth', 1.5, 'Color', [0 0 0 0.25]);

        plot_cov_ellipse(x, y, clr, 0.15);
        v = ~isnan(x) & ~isnan(y);
        if any(v)
            r = corr(x(v), y(v));
            text(0.05, 0.93, sprintf('r = %.2f', r), 'Units','normalized','FontWeight','bold');
        end
        xlabel('Single-Eff Hand Amp (z)'); if c==1, ylabel('Multi-Eff Hand Amp (z)'); end
        title(sprintf('Hand – %s', label_text));

        % ----- Row 2: Eye -----
        subplot(2,4,4+c); hold on; grid on; set(gca,'FontSize',16);
        x = all_amp_single_E(mask);
        y = all_amp_E(mask);
        scatter(x, y, 25, 'filled', 'MarkerFaceColor', clr);
        xlim([-5,15]); ylim([-5,15]); axis square;
        xl = xlim; yl = ylim; mm = [min([xl yl]) max([xl yl])];
        plot(mm, mm, 'k-', 'LineWidth', 1.5, 'Color', [0 0 0 0.25]);

        plot_cov_ellipse(x, y, clr, 0.15);
        v = ~isnan(x) & ~isnan(y);
        if any(v)
            r = corr(x(v), y(v));
            text(0.05, 0.93, sprintf('r = %.2f', r), 'Units','normalized','FontWeight','bold');
        end
        xlabel('Single-Eff Eye Amp (z)'); if c==1, ylabel('Multi-Eff Eye Amp (z)'); end
        title(sprintf('Eye – %s', label_text));
    end

    saveas(h3, fullfile(plotdir, sprintf('%s_Single_vs_Multi_Amplitude_HandEye_byClass_4col.pdf', area)));

    %% Plot 4: Tuning Direction Heatmaps split by class (2x4 grid)
    h4 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);

    nbins = 6;
    edges0 = linspace(-pi, pi, nbins+1);
    bin_centers = edges0(1:end-1) + diff(edges0)/2;

    base_colors = containers.Map( ...
        {'H','E','Both','Both∩Both','None'}, ...
        {category_colors('H'), category_colors('E'), category_colors('Both'), category_colors('Both'), category_colors('None')} ...
    );

    for c = 1:4
        mask = class_masks{c};
        label_text = char(col_classes(c));

        Hs = wrapToPi(all_dir_single_H(mask));
        Hm = wrapToPi(all_dir_multi_H(mask));
        Es = wrapToPi(all_dir_single_E(mask));
        Em = wrapToPi(all_dir_multi_E(mask));

        Hs_bin = circ_bin_nearest(Hs, bin_centers);
        Hm_bin = circ_bin_nearest(Hm, bin_centers);
        Es_bin = circ_bin_nearest(Es, bin_centers);
        Em_bin = circ_bin_nearest(Em, bin_centers);

        h_bins = accumarray([Hs_bin, Hm_bin], 1, [nbins nbins], @sum, 0);
        e_bins = accumarray([Es_bin, Em_bin], 1, [nbins nbins], @sum, 0);

        vmax = max([h_bins(:); e_bins(:); 0]);
        if ~isfinite(vmax) || vmax <= 0, vmax = 1; end

        base = base_colors(label_text);
        cmap = mono_cmap(base, 256);

        axH = subplot(2,4,c);
        imagesc(axH, bin_centers, bin_centers, h_bins);
        axis(axH, 'xy'); axis(axH, 'square');
        colormap(axH, cmap);
        colorbar(axH);
        set(axH, 'FontSize', 16);
        xlabel(axH, 'Single-Eff Pref Dir (rad)');
        if c==1, ylabel(axH, 'Multi-Eff Pref Dir (rad)'); end
        title(axH, sprintf('Hand Dir – %s', label_text));
        xticks(axH, [-pi+.01, 0, pi-.01]); yticks(axH, [-pi+.01, 0, pi-.01]);
        xticklabels(axH, {'-\pi','0','\pi'}); yticklabels(axH, {'-\pi','0','\pi'});
        caxis(axH, [0 vmax]);

        rhoH = circ_corrcc(Hs, Hm);
        if isfinite(rhoH), tH = sprintf('\\rho_{circ}=%.2f', rhoH);
        else,              tH = '\rho_{circ}=n/a'; end
        text(axH, 0.5, 0.02, tH, 'Units','normalized', 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontWeight','bold', 'BackgroundColor','w', ...
             'Margin', 2, 'Interpreter','tex');

        axE = subplot(2,4,4+c);
        imagesc(axE, bin_centers, bin_centers, e_bins);
        axis(axE, 'xy'); axis(axE, 'square');
        colormap(axE, cmap);
        colorbar(axE);
        set(axE, 'FontSize', 16);
        xlabel(axE, 'Single-Eff Pref Dir (rad)');
        if c==1, ylabel(axE, 'Multi-Eff Pref Dir (rad)'); end
        title(axE, sprintf('Eye Dir – %s', label_text));
        xticks(axE, [-pi+.01, 0, pi-.01]); yticks(axE, [-pi+.01, 0, pi-.01]);
        xticklabels(axE, {'-\pi','0','\pi'}); yticklabels(axE, {'-\pi','0','\pi'});
        caxis(axE, [0 vmax]);

        rhoE = circ_corrcc(Es, Em);
        if isfinite(rhoE), tE = sprintf('\\rho_{circ}=%.2f', rhoE);
        else,              tE = '\rho_{circ}=n/a'; end
        text(axE, 0.5, 0.02, tE, 'Units','normalized', 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontWeight','bold', 'BackgroundColor','w', ...
             'Margin', 2, 'Interpreter','tex');
    end

    saveas(h4, fullfile(plotdir, sprintf('%s_Single_vs_Multi_TuningDir_HandEye_byClass_4col.pdf', area)));

    %% Plot 5: Agreement and suppression summary
    h5 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 4]);

    multi_core = ["H","E","Both"];
    same_pct   = zeros(1,3);
    diff_pct   = zeros(1,3);

    for i = 1:3
        m = multi_core(i);
        idx = strcmp(class_labels, m);
        total = sum(idx);
        if total == 0
            same_pct(i) = NaN; diff_pct(i) = NaN;
            continue;
        end
        same = sum(strcmp(class_labels_single(idx), m));
        same_pct(i) = 100 * same / total;
        diff_pct(i) = 100 * (total - same) / total;
    end

    axes('FontSize',16); hold on; grid on; box on;

    X = 1:3;
    Y = [same_pct(:), diff_pct(:)];
    h = bar(X, Y, 'stacked');

    set(h(1), 'FaceColor', [0.3 0.7 0.3], 'EdgeColor','k');
    set(h(2), 'FaceColor', [0.7 0.3 0.3], 'EdgeColor','k');

    xticks(X); xticklabels(multi_core);
    ylabel('% of neurons');
    title('Single vs Multi Class Agreement (by Multi Class)');
    ylim([0 100]);
    legend(h, {'Same','Different'}, 'Location','northeast');

    saveas(h5, fullfile(plotdir, sprintf('%s_Class_Agreement_and_Suppression.pdf', area)));

    %% Plot 6: Preferred direction distribution, significant neurons only
    % Hand plots use neurons classified as H or Both.
    % Eye plots use neurons classified as E or Both.
    % This avoids assigning interpretive weight to peak directions from neurons
    % with no significant directional tuning.

    sig_multi_H  = strcmp(class_labels, 'H') | strcmp(class_labels, 'Both');
    sig_multi_E  = strcmp(class_labels, 'E') | strcmp(class_labels, 'Both');
    sig_single_H = strcmp(class_labels_single, 'H') | strcmp(class_labels_single, 'Both');
    sig_single_E = strcmp(class_labels_single, 'E') | strcmp(class_labels_single, 'Both');

    pd_multi_H  = all_dir_multi_H(sig_multi_H);
    pd_multi_E  = all_dir_multi_E(sig_multi_E);
    pd_single_H = all_dir_single_H(sig_single_H);
    pd_single_E = all_dir_single_E(sig_single_E);

    nbins_pd = 6;
    pd_edges = linspace(-pi, pi, nbins_pd+1);
    pd_centers = pd_edges(1:end-1) + diff(pd_edges)/2;
    pd_labels = {'-\pi', '-2\pi/3', '-\pi/3', '0', '\pi/3', '2\pi/3'};

    h6 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 6]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot_prefdir_hist_panel(pd_multi_H, pd_edges, pd_centers, pd_labels, ...
        'Multi hand PDs', [0.8 0.2 0.2]);

    nexttile;
    plot_prefdir_hist_panel(pd_multi_E, pd_edges, pd_centers, pd_labels, ...
        'Multi eye PDs', [0.2 0.2 0.8]);

    nexttile;
    plot_prefdir_hist_panel(pd_single_H, pd_edges, pd_centers, pd_labels, ...
        'Single hand PDs', [0.8 0.2 0.2]);

    nexttile;
    plot_prefdir_hist_panel(pd_single_E, pd_edges, pd_centers, pd_labels, ...
        'Single eye PDs', [0.2 0.2 0.8]);

    sgtitle(sprintf('%s preferred direction distribution, significant neurons only', area), ...
        'Interpreter','none', 'FontSize', 16, 'FontWeight','bold');

    saveas(h6, fullfile(plotdir, sprintf('%s_PreferredDirection_Distribution_SignificantOnly.pdf', area)));

    fprintf('\n--- Preferred direction counts, significant neurons only (%s) ---\n', area);
    print_prefdir_counts('Multi hand',  pd_multi_H,  pd_edges);
    print_prefdir_counts('Multi eye',   pd_multi_E,  pd_edges);
    print_prefdir_counts('Single hand', pd_single_H, pd_edges);
    print_prefdir_counts('Single eye',  pd_single_E, pd_edges);

    %% NEW: CV summaries (per-curve) as classwise tables, incl. 'All' and 'None'
    fprintf('--- Median per-curve CV by class and effector (ME) ---\n');
    if ~isempty(missing_multi_sessions)
        fprintf('Cannot compute ME (multi) CV: missing BOOT.FR_by_dir.(H|E) in:\n');
        for i = 1:numel(missing_multi_sessions)
            fprintf('  - %s\n', missing_multi_sessions{i});
        end
    else
        ME = compute_classwise_CVs(multi_FR_H_cells, multi_FR_E_cells, class_multi_cells);
        print_cv_table(ME);
    end

    fprintf('--- Median per-curve CV by class and effector (SE) ---\n');
    if ~isempty(missing_single_sessions)
        fprintf('Cannot compute SE (single) CV: missing BOOT.FR_by_dir_single.(H|E) in:\n');
        for i = 1:numel(missing_single_sessions)
            fprintf('  - %s\n', missing_single_sessions{i});
        end
    else
        SE = compute_classwise_CVs(single_FR_H_cells, single_FR_E_cells, class_single_cells);
        print_cv_table(SE);
    end
end

% =========================
% Helper: covariance ellipse
% =========================
function plot_cov_ellipse(x, y, color, faceAlpha)
    v = ~isnan(x) & ~isnan(y);
    X = [x(v), y(v)];
    if size(X,1) < 2, return; end
    mu = mean(X,1);
    S  = cov(X,1);
    if any(~isfinite(S(:))) || rank(S) < 2, return; end
    k = sqrt(chi2inv(0.95, 2));
    [V,D] = eig(S);
    t = linspace(0, 2*pi, 200);
    e = (V * (k*sqrt(D)) * [cos(t); sin(t)])';
    ex = e(:,1) + mu(1);
    ey = e(:,2) + mu(2);
    patch(ex, ey, color, 'FaceAlpha', faceAlpha, 'EdgeColor', color, 'LineWidth', 1.2);
end

% =========================
% Helper: circular mean along a dimension
% =========================
function mu = circ_mean_dim(Theta, dim)
    if dim == 2
        S = mean(sin(Theta), 2, 'omitnan');
        C = mean(cos(Theta), 2, 'omitnan');
        mu = atan2(S, C);
    elseif dim == 1
        S = mean(sin(Theta), 1, 'omitnan');
        C = mean(cos(Theta), 1, 'omitnan');
        mu = atan2(S, C);
    else
        error('Only dim=1 or 2 supported.');
    end
end

% =========================
% Helper: circular nearest-center binning
% =========================
function idx = circ_bin_nearest(theta, centers)
    theta = wrapToPi(theta);
    centers = wrapToPi(centers);
    [T, C] = ndgrid(theta(:), centers(:)');
    D = atan2(sin(T - C), cos(T - C));
    [~, idx] = min(abs(D), [], 2);
    idx = reshape(idx, size(theta));
end

function rho = circ_corrcc(alpha, beta)
    alpha = wrapToPi(alpha(:));
    beta  = wrapToPi(beta(:));
    v = isfinite(alpha) & isfinite(beta);
    alpha = alpha(v);  beta = beta(v);
    n = numel(alpha);
    if n < 3
        rho = NaN; return;
    end
    a0 = atan2(mean(sin(alpha)), mean(cos(alpha)));
    b0 = atan2(mean(sin(beta)),  mean(cos(beta)));
    sa = sin(alpha - a0);
    sb = sin(beta  - b0);
    num = sum(sa .* sb);
    den = sqrt(sum(sa.^2) * sum(sb.^2));
    if den <= 0
        rho = NaN;
    else
        rho = num / den;
    end
end

function y = wrapToPi(x)
    y = mod(x + pi, 2*pi) - pi;
end

% =========================
% Preferred direction histogram helpers
% =========================
function plot_prefdir_hist_panel(theta, edges, centers, labels, title_text, color)
    theta = theta(:);
    theta = theta(isfinite(theta));
    theta = wrapToPi(theta);

    counts = histcounts(theta, edges);
    n = sum(counts);

    bar(1:numel(counts), counts, 'FaceColor', color, 'EdgeColor', 'k');
    hold on; box on; grid on;

    if n > 0
        expected = n / numel(counts);
        yline(expected, 'k--', 'LineWidth', 1.2);
    end

    xticks(1:numel(counts));
    xticklabels(labels);
    xlabel('Preferred direction');
    ylabel('Neuron count');
    title(sprintf('%s (n=%d)', title_text, n), 'Interpreter','none');

    ymax = max([counts(:); 1]);
    ylim([0 ymax * 1.25]);

    for i = 1:numel(counts)
        text(i, counts(i) + 0.03*ymax, sprintf('%d', counts(i)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', 10, 'FontWeight','bold');
    end
end

function print_prefdir_counts(name, theta, edges)
    theta = theta(:);
    theta = theta(isfinite(theta));
    theta = wrapToPi(theta);
    counts = histcounts(theta, edges);
    fprintf('%-12s: n=%d | counts = [%s]\n', name, sum(counts), num2str(counts));
end

% =========================
% NEW helper: amplitude stats from bootstrapped FR-by-dir draws (pooled H+E)
% =========================
function [mu_all, var_all, ncurves_all] = amp_stats_from_cells(H_cells, E_cells)
    mu_all  = [];
    var_all = [];

    for i = 1:numel(H_cells)
        FR = H_cells{i};
        if isempty(FR), continue; end
        A = squeeze(max(FR, [], 2) - min(FR, [], 2));
        if isvector(A), A = reshape(A, [], size(FR,3)); end
        mu_all  = [mu_all;  mean(A, 2, 'omitnan')]; %#ok<AGROW>
        var_all = [var_all; var(A, 0, 2, 'omitnan')]; %#ok<AGROW>
    end

    for i = 1:numel(E_cells)
        FR = E_cells{i};
        if isempty(FR), continue; end
        A = squeeze(max(FR, [], 2) - min(FR, [], 2));
        if isvector(A), A = reshape(A, [], size(FR,3)); end
        mu_all  = [mu_all;  mean(A, 2, 'omitnan')]; %#ok<AGROW>
        var_all = [var_all; var(A, 0, 2, 'omitnan')]; %#ok<AGROW>
    end

    ncurves_all = numel(mu_all);
end

% =========================
% NEW helper: amplitude stats from a single-effector cellset (H-only or E-only)
% =========================
function [mu_all, var_all, ncurves_all] = amp_stats_from_cellset(C_cells)
    mu_all  = [];
    var_all = [];
    for i = 1:numel(C_cells)
        FR = C_cells{i};
        if isempty(FR), continue; end
        A = squeeze(max(FR, [], 2) - min(FR, [], 2));
        if isvector(A), A = reshape(A, [], size(FR,3)); end
        mu_all  = [mu_all;  mean(A, 2, 'omitnan')]; %#ok<AGROW>
        var_all = [var_all; var(A, 0, 2, 'omitnan')]; %#ok<AGROW>
    end
    ncurves_all = numel(mu_all);
end

% =========================
% NEW helpers: classwise per-neuron CVs for H, E, and combined HE
% =========================
function OUT = compute_classwise_CVs(H_cells, E_cells, class_cells)
    muH_all = []; varH_all = [];
    muE_all = []; varE_all = [];
    labels_all = strings(0,1);

    nsess = max([numel(H_cells), numel(E_cells), numel(class_cells)]);
    for i = 1:nsess
        FRH = H_cells{i};
        FRE = E_cells{i};
        labs = class_cells{i}(:);

        if isempty(FRH) || isempty(FRE) || isempty(labs), continue; end

        AH = squeeze(max(FRH, [], 2) - min(FRH, [], 2));
        AE = squeeze(max(FRE, [], 2) - min(FRE, [], 2));
        if isvector(AH), AH = reshape(AH, [], size(FRH,3)); end
        if isvector(AE), AE = reshape(AE, [], size(FRE,3)); end

        muH = mean(AH, 2, 'omitnan');  vH = var(AH, 0, 2, 'omitnan');
        muE = mean(AE, 2, 'omitnan');  vE = var(AE, 0, 2, 'omitnan');

        muH_all = [muH_all; muH]; %#ok<AGROW>
        varH_all= [varH_all; vH];  %#ok<AGROW>
        muE_all = [muE_all; muE]; %#ok<AGROW>
        varE_all= [varE_all; vE];  %#ok<AGROW>
        labels_all = [labels_all; labs]; %#ok<AGROW>
    end

    CV_H  = sqrt(varH_all) ./ muH_all;
    CV_E  = sqrt(varE_all) ./ muE_all;

    goodH  = isfinite(CV_H);
    goodE  = isfinite(CV_E);

    OUT.classes = {'All','H','E','Both','None'};
    OUT.N = struct(); OUT.CV_H = struct(); OUT.CV_E = struct();

    for g = 1:numel(OUT.classes)
        c = OUT.classes{g};
        if strcmp(c,'All')
            idx = true(size(labels_all));
        else
            idx = strcmp(labels_all, c);
        end

        h_vals  = CV_H(idx & goodH);
        e_vals  = CV_E(idx & goodE);

        OUT.N.(c)     = [numel(h_vals), numel(e_vals)];
        OUT.CV_H.(c)  = median(h_vals,  'omitnan');
        OUT.CV_E.(c)  = median(e_vals,  'omitnan');
    end
end

function print_cv_table(OUT)
    fprintf('%-6s | %-12s | %-12s | %s\n', 'Class', 'CV_H', 'CV_E', 'N(H/E)');
    fprintf('%s\n', repmat('-',1,66));
    for g = 1:numel(OUT.classes)
        c = OUT.classes{g};
        cvh  = OUT.CV_H.(c);
        cve  = OUT.CV_E.(c);
        Ns   = OUT.N.(c);
        fprintf('%-6s | %12.6f | %12.6f | %d/%d\n', c, cvh, cve, Ns(1), Ns(2));
    end
end

% =========================
% NEW helper: monochrome colormap (white -> darker target color)
% =========================
function cmap = mono_cmap(baseColor, nLevels)
    if nargin < 2, nLevels = 256; end
    baseColor = min(max(baseColor(:).', 0), 1);
    t = linspace(0, 1, nLevels)';
    target = 0.85 * baseColor;
    cmap = (1 - t) .* ones(nLevels,3) + t .* target;
end
