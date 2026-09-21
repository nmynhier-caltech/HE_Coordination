function plot_inseparable_neuron_summary(area, savepaths, varargin)
% plot_inseparable_neuron_summary
%
% Robust plotting + inseparable/separable phenotype summary.
%
% First output file now contains 22 panels by default:
%   Top row:
%     1) Overall/raw gradient-axis distribution
%     2) Residual gradient-axis distribution
%   Then 5 rows x 4 columns:
%     1) Overall F heatmap/quiver
%     2) Overall F gradient-axis polar plot
%     3) Residual R = F - Fadd heatmap/quiver
%     4) Residual gradient-axis polar plot
%
% The 5 examples are:
%   - 5 random globally inseparable neurons
%
% Change the random examples with the RandomExampleSeed option near the top below.
%
% Other outputs:
%   2) Raw peak location heatmaps
%   3) Residual peak location heatmaps
%   4) Hand PD vs eye PD heatmaps
%   5) Angular difference distribution
%   6) Peak class fractions
%   7) Basic scalar metric comparison
%   8) Manuscript-ready amplitude/relational summaries are printed and saved

    p = inputParser;
    addParameter(p, 'NBins', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(p, 'Alpha', 0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
    addParameter(p, 'SmoothSurface', true, @(x)islogical(x)||ismember(x,[0 1]));
    addParameter(p, 'ThetaEdgeOffset', 0.1, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'NPerm', 10000, @(x)isnumeric(x)&&isscalar(x)&&x>=100);

    % ============================================================
    % Random example-neuron selection
    % ============================================================
    % Change RandomExampleSeed to regenerate a different fixed set of
    % random example neurons. NumRandomExamples controls the number of
    % example rows in the first PDF.
    addParameter(p, 'RandomExampleSeed', 678, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    addParameter(p, 'NumRandomExamples', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

    parse(p, varargin{:});

    nbins = p.Results.NBins;
    alpha = p.Results.Alpha;
    do_smooth = logical(p.Results.SmoothSurface);
    theta_edge_offset = p.Results.ThetaEdgeOffset;
    nperm = p.Results.NPerm;
    random_example_seed = round(p.Results.RandomExampleSeed);
    num_random_examples = round(p.Results.NumRandomExamples);

    % The example-neuron choice was already seeded below. Use a separate,
    % fixed stream for the permutation tests so their reported p-values are
    % reproducible as well.
    rng(random_example_seed + 1, 'twister');

    bin_deg = 360 / nbins;

    all_pHE = [];
    all_session_for_p = [];
    all_neuron_for_p = [];

    for s = 1:numel(savepaths)
        savepath = savepaths{s};
        hess_file = fullfile(savepath, 'models', sprintf('%s_hessian_results.mat', area));

        if ~exist(hess_file, 'file')
            warning('Missing hessian file: %s', hess_file);
            continue;
        end

        H = load(hess_file);
        if ~isfield(H, 'results')
            warning('No results struct in: %s', hess_file);
            continue;
        end

        r = H.results(:);

        if ~isfield(r, 'p_HE')
            warning('results does not contain p_HE in: %s', hess_file);
            continue;
        end

        pHE = arrayfun(@(x) x.p_HE, r);
        pHE = pHE(:);

        valid = isfinite(pHE) & pHE >= 0 & pHE <= 1;

        all_pHE = [all_pHE; pHE(valid)]; %#ok<AGROW>
        all_session_for_p = [all_session_for_p; s * ones(sum(valid),1)]; %#ok<AGROW>
        all_neuron_for_p = [all_neuron_for_p; find(valid(:))]; %#ok<AGROW>
    end

    if isempty(all_pHE)
        error('No valid p_HE values found for area %s.', area);
    end

    qHE = bh_fdr(all_pHE);
    is_inseparable_global = qHE < alpha;

    fprintf('\n[%s] Global HE cross-derivative FDR summary:\n', area);
    fprintf('  Valid neurons: %d\n', numel(all_pHE));
    fprintf('  Inseparable neurons, q_HE < %.3f: %d / %d = %.1f%%\n', ...
        alpha, sum(is_inseparable_global), numel(all_pHE), ...
        100 * mean(is_inseparable_global));

    class_lookup = containers.Map('KeyType','char','ValueType','logical');
    q_lookup = containers.Map('KeyType','char','ValueType','double');
    p_lookup = containers.Map('KeyType','char','ValueType','double');

    for k = 1:numel(all_pHE)
        key = sprintf('%d_%d', all_session_for_p(k), all_neuron_for_p(k));
        class_lookup(key) = is_inseparable_global(k);
        q_lookup(key) = qHE(k);
        p_lookup(key) = all_pHE(k);
    end

    insep = init_group_struct();
    sep   = init_group_struct();

    skipped_sessions = {};
    used_sessions = {};

    for s = 1:numel(savepaths)
        savepath = savepaths{s};

        data_file = processed_data_file(savepath, area);
        if ~exist(data_file, 'file')
            warning('Missing processed data file: %s', data_file);
            skipped_sessions{end+1} = savepath; %#ok<AGROW>
            continue;
        end

        D = load(data_file);
        if ~isfield(D, 'data')
            warning('No data struct in: %s', data_file);
            skipped_sessions{end+1} = savepath; %#ok<AGROW>
            continue;
        end

        data = D.data;
        required_fields = {'FR','TP1','TP2','TPi1','TPi2'};
        if ~all(isfield(data, required_fields))
            warning('Processed data missing required fields in: %s', data_file);
            skipped_sessions{end+1} = savepath; %#ok<AGROW>
            continue;
        end

        nNeurons = size(data.FR, 2);

        idxME = find(data.TPi2 ~= 7 & data.TPi1 ~= 7);
        FR_me = data.FR(idxME, :);
        TP1_me = data.TP1(idxME, :);
        TP2_me = data.TP2(idxME, :);

        thetaH = atan2(TP1_me(:,2), TP1_me(:,1));
        thetaE = atan2(TP2_me(:,2), TP2_me(:,1));

        theta_edges = linspace(-pi, pi, nbins+1) + theta_edge_offset;
        Hb = discretize(thetaH, theta_edges);
        Eb = discretize(thetaE, theta_edges);

        valid = isfinite(Hb) & isfinite(Eb);
        Hb = Hb(valid);
        Eb = Eb(valid);
        FR_me = FR_me(valid, :);

        if isempty(Hb)
            warning('No valid ME trials after binning in: %s', savepath);
            skipped_sessions{end+1} = savepath; %#ok<AGROW>
            continue;
        end

        n_insep_sess = 0;
        n_sep_sess = 0;

        for n = 1:nNeurons
            key = sprintf('%d_%d', s, n);

            if ~isKey(class_lookup, key)
                continue;
            end

            is_insep = class_lookup(key);

            F = nan(nbins, nbins);
            for h = 1:nbins
                for e = 1:nbins
                    mask = (Hb == h & Eb == e);
                    if any(mask)
                        F(h,e) = mean(FR_me(mask,n), 'omitnan');
                    end
                end
            end

            if any(~isfinite(F(:)))
                gm = mean(F(:), 'omitnan');
                F(~isfinite(F)) = gm;
            end

            if do_smooth
                F = smooth_circular_surface(F);
                F = smooth_circular_surface(F);
            end

            metrics = compute_neuron_metrics(F, nbins, bin_deg);
            metrics.session = s;
            metrics.neuron = n;
            metrics.pHE = p_lookup(key);
            metrics.qHE = q_lookup(key);

            if is_insep
                insep = append_metrics(insep, metrics);
                n_insep_sess = n_insep_sess + 1;
            else
                sep = append_metrics(sep, metrics);
                n_sep_sess = n_sep_sess + 1;
            end
        end

        fprintf('[%s session %d] globally inseparable=%d, separable=%d\n', ...
            area, s, n_insep_sess, n_sep_sess);

        used_sessions{end+1} = savepath; %#ok<AGROW>
    end

    root_for_output = savepaths{1};
    plotdir = fullfile(root_for_output, '..', 'plots', 'summary', ...
        'inseparable_neurons');
    modeldir = fullfile(root_for_output, 'models');

    if ~exist(plotdir, 'dir'), mkdir(plotdir); end
    if ~exist(modeldir, 'dir'), mkdir(modeldir); end

    if isempty(insep.raw_peak_H)
        Summary = struct();
        Summary.area = area;
        Summary.nbins = nbins;
        Summary.alpha = alpha;
        Summary.global_valid_neurons = numel(all_pHE);
        Summary.global_pHE = all_pHE;
        Summary.global_qHE = qHE;
        Summary.global_session_index = all_session_for_p;
        Summary.global_neuron_index = all_neuron_for_p;
        Summary.global_inseparable_mask = is_inseparable_global;
        Summary.n_inseparable_total = 0;

        outmat = fullfile(modeldir, sprintf('%s_Inseparable_ExtendedSummary.mat', area));
        save(outmat, 'Summary', '-v7.3');

        fprintf('\n[%s] No globally inseparable neurons to plot. Saved summary:\n  %s\n', area, outmat);
        return;
    end

    Stats = struct();
    Stats.inseparable = compute_group_stats(insep, nbins);
    Stats.separable   = compute_group_stats(sep, nbins);
    Stats.basic_metric_comparisons = compute_basic_metric_comparisons(insep, sep, nperm);

    % ============================================================
    % Plot 1: 50-panel gradient figure
    % ============================================================
    fig1_height = max(22.0, 5.0 + 4.0 * num_random_examples);
    fig1 = make_safe_figure(16.0, fig1_height);
    plot_gradient_axis_distribution_and_6_examples_full(insep, nbins, area, ...
        random_example_seed, num_random_examples);

    outpdf_axis = fullfile(plotdir, sprintf('%s_Inseparable_GradientAxis.pdf', area));
    save_safe_pdf(fig1, outpdf_axis, 16.0, fig1_height);
    close(fig1);

    % ============================================================
    % Plot 2: raw peak heatmaps
    % ============================================================
    fig2 = make_safe_figure(7.8, 13.8);
    tiledlayout(2,1,'TileSpacing','loose','Padding','loose');

    nexttile;
    plot_count_heatmap(gca, Stats.inseparable.raw_peak_counts, nbins, ...
        sprintf('%s inseparable: raw peak location', area), ...
        'Eye target bin at peak', 'Hand target bin at peak');

    nexttile;
    plot_count_heatmap(gca, Stats.separable.raw_peak_counts, nbins, ...
        sprintf('%s separable: raw peak location', area), ...
        'Eye target bin at peak', 'Hand target bin at peak');

    outpdf_raw = fullfile(plotdir, sprintf('%s_RawPeakLocation_Inseparable_vs_Separable.pdf', area));
    save_safe_pdf(fig2, outpdf_raw, 7.8, 13.8);
    close(fig2);

    % ============================================================
    % Plot 3: residual peak heatmaps
    % ============================================================
    fig3 = make_safe_figure(7.8, 13.8);
    tiledlayout(2,1,'TileSpacing','loose','Padding','loose');

    nexttile;
    plot_count_heatmap(gca, Stats.inseparable.resid_peak_counts, nbins, ...
        sprintf('%s inseparable: residual peak |F-F_{add}|', area), ...
        'Eye target bin at residual peak', 'Hand target bin at residual peak');

    nexttile;
    plot_count_heatmap(gca, Stats.separable.resid_peak_counts, nbins, ...
        sprintf('%s separable: residual peak |F-F_{add}|', area), ...
        'Eye target bin at residual peak', 'Hand target bin at residual peak');

    outpdf_resid = fullfile(plotdir, sprintf('%s_ResidualPeakLocation_Inseparable_vs_Separable.pdf', area));
    save_safe_pdf(fig3, outpdf_resid, 7.8, 13.8);
    close(fig3);

    % ============================================================
    % Plot 4: hand PD vs eye PD
    % ============================================================
    fig4 = make_safe_figure(7.8, 13.8);
    tiledlayout(2,1,'TileSpacing','loose','Padding','loose');

    nexttile;
    plot_count_heatmap(gca, Stats.inseparable.PD_pair_counts, nbins, ...
        sprintf('%s inseparable: hand PD vs eye PD', area), ...
        'Eye preferred direction bin', 'Hand preferred direction bin');

    nexttile;
    plot_count_heatmap(gca, Stats.separable.PD_pair_counts, nbins, ...
        sprintf('%s separable: hand PD vs eye PD', area), ...
        'Eye preferred direction bin', 'Hand preferred direction bin');

    outpdf_pd = fullfile(plotdir, sprintf('%s_PDPair_Inseparable_vs_Separable.pdf', area));
    save_safe_pdf(fig4, outpdf_pd, 7.8, 13.8);
    close(fig4);

    % ============================================================
    % Plot 5: angular difference distribution
    % ============================================================
    fig5 = make_safe_figure(8.0, 5.8);
    plot_angular_difference_distribution(gca, insep.raw_peak_diff_deg, sep.raw_peak_diff_deg, nbins);
    title(sprintf('%s raw peak hand-eye angular difference', area), 'Interpreter','none');

    outpdf_diff = fullfile(plotdir, sprintf('%s_RawPeak_AngularDifference_Inseparable_vs_Separable.pdf', area));
    save_safe_pdf(fig5, outpdf_diff, 8.0, 5.8);
    close(fig5);

    % ============================================================
    % Plot 6: summary fractions
    % ============================================================
    fig6 = make_safe_figure(7.8, 15.2);
    plot_summary_fraction_bars(Stats, nbins, area);

    outpdf_frac = fullfile(plotdir, sprintf('%s_PeakClassFractions_Inseparable_vs_Separable.pdf', area));
    save_safe_pdf(fig6, outpdf_frac, 7.8, 15.2);
    close(fig6);

    % ============================================================
    % Plot 7: basic scalar metrics
    % ============================================================
    fig7 = make_safe_figure(8.2, 18.5);
    plot_basic_metric_comparison_bars(insep, sep, Stats.basic_metric_comparisons, area);

    outpdf_basic = fullfile(plotdir, sprintf('%s_SeparableVsInseparable_BasicMetrics.pdf', area));
    save_safe_pdf(fig7, outpdf_basic, 8.2, 18.5);
    close(fig7);

    % ============================================================
    % Save summary
    % ============================================================
    Summary = struct();
    Summary.area = area;
    Summary.nbins = nbins;
    Summary.bin_deg = bin_deg;
    Summary.alpha = alpha;
    Summary.nperm = nperm;
    Summary.random_example_seed = random_example_seed;
    Summary.num_random_examples = num_random_examples;

    Summary.global_valid_neurons = numel(all_pHE);
    Summary.global_pHE = all_pHE;
    Summary.global_qHE = qHE;
    Summary.global_session_index = all_session_for_p;
    Summary.global_neuron_index = all_neuron_for_p;
    Summary.global_inseparable_mask = is_inseparable_global;

    Summary.used_sessions = used_sessions;
    Summary.skipped_sessions = skipped_sessions;

    Summary.inseparable = insep;
    Summary.separable = sep;
    Summary.stats = Stats;

    Summary.files.gradient_axis = outpdf_axis;
    Summary.files.raw_peak_heatmaps = outpdf_raw;
    Summary.files.residual_peak_heatmaps = outpdf_resid;
    Summary.files.PD_pair_heatmaps = outpdf_pd;
    Summary.files.raw_peak_angular_difference = outpdf_diff;
    Summary.files.summary_fraction_bars = outpdf_frac;
    Summary.files.basic_metric_comparison = outpdf_basic;

    outmat = fullfile(modeldir, sprintf('%s_Inseparable_ExtendedSummary.mat', area));
    save(outmat, 'Summary', '-v7.3');

    fprintf('\nSaved inseparable-neuron summary:\n  %s\n', outmat);
    fprintf('Saved plots:\n');
    fprintf('  %s\n', outpdf_axis);
    fprintf('  %s\n', outpdf_raw);
    fprintf('  %s\n', outpdf_resid);
    fprintf('  %s\n', outpdf_pd);
    fprintf('  %s\n', outpdf_diff);
    fprintf('  %s\n', outpdf_frac);
    fprintf('  %s\n', outpdf_basic);

    print_stats_to_console(area, Stats, nbins);
end

% ============================================================
% Safe figure/export helpers
% ============================================================
function fig = make_safe_figure(width_in, height_in)
    dpi = 120;
    fig = figure('Visible','Off', ...
        'Units','pixels', ...
        'Position',[100 100 round(width_in*dpi) round(height_in*dpi)], ...
        'Color','w', ...
        'InvertHardcopy','off', ...
        'PaperUnits','inches', ...
        'PaperSize',[width_in height_in], ...
        'PaperPosition',[0.35 0.35 width_in-0.70 height_in-0.70], ...
        'PaperPositionMode','manual');
end

function save_safe_pdf(fig, outpdf, width_in, height_in)
    drawnow;

    set(fig, 'PaperUnits','inches');
    set(fig, 'PaperSize',[width_in height_in]);
    set(fig, 'PaperPosition',[0.35 0.35 width_in-0.70 height_in-0.70]);
    set(fig, 'PaperPositionMode','manual');
    set(fig, 'InvertHardcopy','off');

    try
        exportgraphics(fig, outpdf, ...
            'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        print(fig, outpdf, '-dpdf', '-painters');
    end
end

% ============================================================
% Data containers
% ============================================================
function G = init_group_struct()
    G.raw_peak_H = [];
    G.raw_peak_E = [];
    G.resid_peak_H = [];
    G.resid_peak_E = [];
    G.PD_H = [];
    G.PD_E = [];

    G.raw_peak_diff_bins = [];
    G.raw_peak_diff_deg = [];
    G.resid_peak_diff_bins = [];
    G.resid_peak_diff_deg = [];
    G.PD_diff_bins = [];
    G.PD_diff_deg = [];

    G.axis_angle = [];
    G.axis_strength = [];
    G.resid_axis_angle = [];
    G.resid_axis_strength = [];

    G.mean_FR = [];
    G.raw_amp = [];
    G.add_amp = [];
    G.resid_amp = [];
    G.resid_max_abs = [];
    G.resid_var_frac = [];
    G.add_to_resid_amp_ratio = [];
    G.resid_to_add_amp_ratio = [];

    % Residual relational-code metrics.  The residual is R = F - Fadd.
    % Congruent axis: hand bin == eye bin.
    % Anti-congruent axis: hand bin is 180 deg from eye bin, only defined for even nbins.
    G.resid_peak_is_congruent = [];
    G.resid_peak_is_antiparallel = [];
    G.resid_peak_is_relational = [];

    G.resid_abs_mean_congruent = [];
    G.resid_abs_mean_antiparallel = [];
    G.resid_abs_mean_relational = [];
    G.resid_abs_mean_offaxis = [];
    G.resid_abs_relational_to_offaxis_ratio = [];
    G.resid_abs_congruent_to_offaxis_ratio = [];
    G.resid_abs_antiparallel_to_offaxis_ratio = [];
    G.resid_abs_relational_minus_offaxis = [];

    G.resid_signed_mean_congruent = [];
    G.resid_signed_mean_antiparallel = [];
    G.resid_signed_mean_relational = [];
    G.resid_signed_mean_offaxis = [];

    G.resid_energy_congruent = [];
    G.resid_energy_antiparallel = [];
    G.resid_energy_relational = [];
    G.resid_energy_offaxis = [];
    G.resid_energy_relational_frac = [];
    G.hand_mod_amp = [];
    G.eye_mod_amp = [];
    G.hand_vector_strength = [];
    G.eye_vector_strength = [];
    G.surface_sparsity = [];

    G.F_surface = {};
    G.Fadd_surface = {};
    G.R_surface = {};

    G.session = [];
    G.neuron = [];
    G.pHE = [];
    G.qHE = [];
end

function G = append_metrics(G, M)
    G.raw_peak_H(end+1,1) = M.raw_peak_H;
    G.raw_peak_E(end+1,1) = M.raw_peak_E;
    G.resid_peak_H(end+1,1) = M.resid_peak_H;
    G.resid_peak_E(end+1,1) = M.resid_peak_E;
    G.PD_H(end+1,1) = M.PD_H;
    G.PD_E(end+1,1) = M.PD_E;

    G.raw_peak_diff_bins(end+1,1) = M.raw_peak_diff_bins;
    G.raw_peak_diff_deg(end+1,1) = M.raw_peak_diff_deg;
    G.resid_peak_diff_bins(end+1,1) = M.resid_peak_diff_bins;
    G.resid_peak_diff_deg(end+1,1) = M.resid_peak_diff_deg;
    G.PD_diff_bins(end+1,1) = M.PD_diff_bins;
    G.PD_diff_deg(end+1,1) = M.PD_diff_deg;

    G.axis_angle(end+1,1) = M.axis_angle;
    G.axis_strength(end+1,1) = M.axis_strength;
    G.resid_axis_angle(end+1,1) = M.resid_axis_angle;
    G.resid_axis_strength(end+1,1) = M.resid_axis_strength;

    G.mean_FR(end+1,1) = M.mean_FR;
    G.raw_amp(end+1,1) = M.raw_amp;
    G.add_amp(end+1,1) = M.add_amp;
    G.resid_amp(end+1,1) = M.resid_amp;
    G.resid_max_abs(end+1,1) = M.resid_max_abs;
    G.resid_var_frac(end+1,1) = M.resid_var_frac;
    G.add_to_resid_amp_ratio(end+1,1) = M.add_to_resid_amp_ratio;
    G.resid_to_add_amp_ratio(end+1,1) = M.resid_to_add_amp_ratio;

    G.resid_peak_is_congruent(end+1,1) = M.resid_peak_is_congruent;
    G.resid_peak_is_antiparallel(end+1,1) = M.resid_peak_is_antiparallel;
    G.resid_peak_is_relational(end+1,1) = M.resid_peak_is_relational;

    G.resid_abs_mean_congruent(end+1,1) = M.resid_abs_mean_congruent;
    G.resid_abs_mean_antiparallel(end+1,1) = M.resid_abs_mean_antiparallel;
    G.resid_abs_mean_relational(end+1,1) = M.resid_abs_mean_relational;
    G.resid_abs_mean_offaxis(end+1,1) = M.resid_abs_mean_offaxis;
    G.resid_abs_relational_to_offaxis_ratio(end+1,1) = M.resid_abs_relational_to_offaxis_ratio;
    G.resid_abs_congruent_to_offaxis_ratio(end+1,1) = M.resid_abs_congruent_to_offaxis_ratio;
    G.resid_abs_antiparallel_to_offaxis_ratio(end+1,1) = M.resid_abs_antiparallel_to_offaxis_ratio;
    G.resid_abs_relational_minus_offaxis(end+1,1) = M.resid_abs_relational_minus_offaxis;

    G.resid_signed_mean_congruent(end+1,1) = M.resid_signed_mean_congruent;
    G.resid_signed_mean_antiparallel(end+1,1) = M.resid_signed_mean_antiparallel;
    G.resid_signed_mean_relational(end+1,1) = M.resid_signed_mean_relational;
    G.resid_signed_mean_offaxis(end+1,1) = M.resid_signed_mean_offaxis;

    G.resid_energy_congruent(end+1,1) = M.resid_energy_congruent;
    G.resid_energy_antiparallel(end+1,1) = M.resid_energy_antiparallel;
    G.resid_energy_relational(end+1,1) = M.resid_energy_relational;
    G.resid_energy_offaxis(end+1,1) = M.resid_energy_offaxis;
    G.resid_energy_relational_frac(end+1,1) = M.resid_energy_relational_frac;
    G.hand_mod_amp(end+1,1) = M.hand_mod_amp;
    G.eye_mod_amp(end+1,1) = M.eye_mod_amp;
    G.hand_vector_strength(end+1,1) = M.hand_vector_strength;
    G.eye_vector_strength(end+1,1) = M.eye_vector_strength;
    G.surface_sparsity(end+1,1) = M.surface_sparsity;

    G.F_surface{end+1,1} = M.F_surface;
    G.Fadd_surface{end+1,1} = M.Fadd_surface;
    G.R_surface{end+1,1} = M.R_surface;

    G.session(end+1,1) = M.session;
    G.neuron(end+1,1) = M.neuron;
    G.pHE(end+1,1) = M.pHE;
    G.qHE(end+1,1) = M.qHE;
end

% ============================================================
% Per-neuron metrics
% ============================================================
function M = compute_neuron_metrics(F, nbins, bin_deg)
    M = struct();

    [~, imax] = max(F(:));
    [M.raw_peak_H, M.raw_peak_E] = ind2sub([nbins nbins], imax);

    hand_curve = mean(F, 2, 'omitnan');
    eye_curve  = mean(F, 1, 'omitnan')';

    [~, M.PD_H] = max(hand_curve);
    [~, M.PD_E] = max(eye_curve);

    Fadd = additive_reconstruction_from_surface(F);
    R = F - Fadd;

    [~, ires] = max(abs(R(:)));
    [M.resid_peak_H, M.resid_peak_E] = ind2sub([nbins nbins], ires);

    [M.raw_peak_diff_bins, M.raw_peak_diff_deg] = circular_bin_difference(M.raw_peak_H, M.raw_peak_E, nbins, bin_deg);
    [M.resid_peak_diff_bins, M.resid_peak_diff_deg] = circular_bin_difference(M.resid_peak_H, M.resid_peak_E, nbins, bin_deg);
    [M.PD_diff_bins, M.PD_diff_deg] = circular_bin_difference(M.PD_H, M.PD_E, nbins, bin_deg);

    [M.axis_angle, M.axis_strength] = gradient_resultant_axis(F);
    [M.resid_axis_angle, M.resid_axis_strength] = gradient_resultant_axis(R);

    M.mean_FR = mean(F(:), 'omitnan');

    M.raw_amp = max(F(:), [], 'omitnan') - min(F(:), [], 'omitnan');
    M.add_amp = max(Fadd(:), [], 'omitnan') - min(Fadd(:), [], 'omitnan');
    M.resid_amp = max(R(:), [], 'omitnan') - min(R(:), [], 'omitnan');
    M.resid_max_abs = max(abs(R(:)), [], 'omitnan');

    F0 = F - mean(F(:), 'omitnan');
    denom = sum(F0(:).^2, 'omitnan');
    numer = sum(R(:).^2, 'omitnan');
    if denom > 0
        M.resid_var_frac = numer / denom;
    else
        M.resid_var_frac = NaN;
    end

    if M.resid_amp > 0
        M.add_to_resid_amp_ratio = M.add_amp / M.resid_amp;
    else
        M.add_to_resid_amp_ratio = NaN;
    end

    if M.add_amp > 0
        M.resid_to_add_amp_ratio = M.resid_amp / M.add_amp;
    else
        M.resid_to_add_amp_ratio = NaN;
    end

    rel = compute_residual_relational_metrics(R, nbins, M.resid_peak_H, M.resid_peak_E);
    rel_fields = fieldnames(rel);
    for ii = 1:numel(rel_fields)
        M.(rel_fields{ii}) = rel.(rel_fields{ii});
    end

    M.hand_mod_amp = max(hand_curve, [], 'omitnan') - min(hand_curve, [], 'omitnan');
    M.eye_mod_amp  = max(eye_curve, [], 'omitnan')  - min(eye_curve, [], 'omitnan');

    theta_centers = linspace(-pi, pi, nbins+1);
    theta_centers = theta_centers(1:end-1) + pi/nbins;

    M.hand_vector_strength = circular_vector_strength(hand_curve, theta_centers(:));
    M.eye_vector_strength  = circular_vector_strength(eye_curve,  theta_centers(:));

    M.surface_sparsity = lifetime_sparsity(F(:));

    M.F_surface = F;
    M.Fadd_surface = Fadd;
    M.R_surface = R;
end

function Fadd = additive_reconstruction_from_surface(F)
    grand = mean(F(:), 'omitnan');
    h = mean(F, 2, 'omitnan') - grand;
    e = mean(F, 1, 'omitnan') - grand;
    Fadd = grand + h + e;
end

function v = circular_vector_strength(curve, theta)
    curve = curve(:);
    theta = theta(:);

    valid = isfinite(curve) & isfinite(theta);
    curve = curve(valid);
    theta = theta(valid);

    if isempty(curve)
        v = NaN;
        return;
    end

    w = curve - min(curve);
    if sum(w) <= 0
        v = 0;
        return;
    end

    z = sum(w .* exp(1i * theta)) / sum(w);
    v = abs(z);
end

function s = lifetime_sparsity(x)
    x = x(:);
    x = x(isfinite(x));

    if isempty(x)
        s = NaN;
        return;
    end

    x = x - min(x);
    if all(x == 0)
        s = 0;
        return;
    end

    n = numel(x);
    s = (1 - (mean(x)^2 / mean(x.^2))) / (1 - 1/n);
end

% ============================================================
% Group stats
% ============================================================
function S = compute_group_stats(G, nbins)
    S = struct();

    S.n = numel(G.raw_peak_H);

    S.raw_peak_counts = accumarray([G.raw_peak_H, G.raw_peak_E], 1, [nbins nbins], @sum, 0);
    S.resid_peak_counts = accumarray([G.resid_peak_H, G.resid_peak_E], 1, [nbins nbins], @sum, 0);
    S.PD_pair_counts = accumarray([G.PD_H, G.PD_E], 1, [nbins nbins], @sum, 0);

    S.raw_peak = class_fraction_stats(G.raw_peak_H, G.raw_peak_E, nbins);
    S.resid_peak = class_fraction_stats(G.resid_peak_H, G.resid_peak_E, nbins);
    S.PD_pair = class_fraction_stats(G.PD_H, G.PD_E, nbins);

    S.raw_peak_diff_deg = G.raw_peak_diff_deg;
    S.resid_peak_diff_deg = G.resid_peak_diff_deg;
    S.PD_diff_deg = G.PD_diff_deg;

    S.basic_metrics = summarize_basic_metrics(G);
    S.manuscript_numbers = compute_manuscript_numbers(G, nbins);
end

function B = summarize_basic_metrics(G)
    metric_names = basic_metric_names();

    B = struct();
    for i = 1:numel(metric_names)
        f = metric_names{i};
        x = G.(f);
        x = x(:);
        x = x(isfinite(x));

        B.(f).n = numel(x);
        B.(f).mean = mean(x, 'omitnan');
        B.(f).median = median(x, 'omitnan');
        B.(f).std = std(x, 'omitnan');
        B.(f).sem = B.(f).std / sqrt(max(1, B.(f).n));
        B.(f).iqr = iqr_local(x);
    end
end


function Rm = compute_residual_relational_metrics(R, nbins, resid_peak_H, resid_peak_E)
    diag_mask = false(nbins, nbins);
    anti_mask = false(nbins, nbins);

    for h = 1:nbins
        diag_mask(h,h) = true;

        if mod(nbins,2) == 0
            half = nbins / 2;
            eanti = mod(h - 1 + half, nbins) + 1;
            anti_mask(h,eanti) = true;
        end
    end

    relational_mask = diag_mask | anti_mask;
    offaxis_mask = ~relational_mask;

    absR = abs(R);
    energyR = R.^2;

    Rm = struct();

    Rm.resid_peak_is_congruent = resid_peak_H == resid_peak_E;

    if mod(nbins,2) == 0
        half = nbins / 2;
        anti_E = mod(resid_peak_H - 1 + half, nbins) + 1;
        Rm.resid_peak_is_antiparallel = resid_peak_E == anti_E;
    else
        Rm.resid_peak_is_antiparallel = false;
    end

    Rm.resid_peak_is_relational = Rm.resid_peak_is_congruent || Rm.resid_peak_is_antiparallel;

    Rm.resid_abs_mean_congruent = mean(absR(diag_mask), 'omitnan');

    if any(anti_mask(:))
        Rm.resid_abs_mean_antiparallel = mean(absR(anti_mask), 'omitnan');
    else
        Rm.resid_abs_mean_antiparallel = NaN;
    end

    Rm.resid_abs_mean_relational = mean(absR(relational_mask), 'omitnan');
    Rm.resid_abs_mean_offaxis = mean(absR(offaxis_mask), 'omitnan');

    if Rm.resid_abs_mean_offaxis > 0
        Rm.resid_abs_relational_to_offaxis_ratio = Rm.resid_abs_mean_relational / Rm.resid_abs_mean_offaxis;
        Rm.resid_abs_congruent_to_offaxis_ratio = Rm.resid_abs_mean_congruent / Rm.resid_abs_mean_offaxis;
        Rm.resid_abs_antiparallel_to_offaxis_ratio = Rm.resid_abs_mean_antiparallel / Rm.resid_abs_mean_offaxis;
    else
        Rm.resid_abs_relational_to_offaxis_ratio = NaN;
        Rm.resid_abs_congruent_to_offaxis_ratio = NaN;
        Rm.resid_abs_antiparallel_to_offaxis_ratio = NaN;
    end

    Rm.resid_abs_relational_minus_offaxis = Rm.resid_abs_mean_relational - Rm.resid_abs_mean_offaxis;

    Rm.resid_signed_mean_congruent = mean(R(diag_mask), 'omitnan');

    if any(anti_mask(:))
        Rm.resid_signed_mean_antiparallel = mean(R(anti_mask), 'omitnan');
    else
        Rm.resid_signed_mean_antiparallel = NaN;
    end

    Rm.resid_signed_mean_relational = mean(R(relational_mask), 'omitnan');
    Rm.resid_signed_mean_offaxis = mean(R(offaxis_mask), 'omitnan');

    Rm.resid_energy_congruent = sum(energyR(diag_mask), 'omitnan');

    if any(anti_mask(:))
        Rm.resid_energy_antiparallel = sum(energyR(anti_mask), 'omitnan');
    else
        Rm.resid_energy_antiparallel = NaN;
    end

    Rm.resid_energy_relational = sum(energyR(relational_mask), 'omitnan');
    Rm.resid_energy_offaxis = sum(energyR(offaxis_mask), 'omitnan');

    total_energy = Rm.resid_energy_relational + Rm.resid_energy_offaxis;
    if total_energy > 0
        Rm.resid_energy_relational_frac = Rm.resid_energy_relational / total_energy;
    else
        Rm.resid_energy_relational_frac = NaN;
    end
end

function M = compute_manuscript_numbers(G, nbins)
    M = struct();
    M.n = numel(G.raw_peak_H);

    M.add_amp = summarize_vector_for_text(G.add_amp);
    M.resid_amp = summarize_vector_for_text(G.resid_amp);
    M.raw_amp = summarize_vector_for_text(G.raw_amp);
    M.resid_var_frac = summarize_vector_for_text(G.resid_var_frac);

    M.add_to_resid_amp_ratio = summarize_vector_for_text(G.add_to_resid_amp_ratio);
    M.resid_to_add_amp_ratio = summarize_vector_for_text(G.resid_to_add_amp_ratio);

    M.resid_abs_mean_relational = summarize_vector_for_text(G.resid_abs_mean_relational);
    M.resid_abs_mean_offaxis = summarize_vector_for_text(G.resid_abs_mean_offaxis);
    M.resid_abs_relational_to_offaxis_ratio = summarize_vector_for_text(G.resid_abs_relational_to_offaxis_ratio);
    M.resid_abs_relational_minus_offaxis = summarize_vector_for_text(G.resid_abs_relational_minus_offaxis);
    M.resid_energy_relational_frac = summarize_vector_for_text(G.resid_energy_relational_frac);

    is_rel = logical(G.resid_peak_is_relational(:));
    is_cong = logical(G.resid_peak_is_congruent(:));
    is_anti = logical(G.resid_peak_is_antiparallel(:));

    M.resid_peak_relational = fraction_summary_for_text(sum(is_rel), M.n);
    M.resid_peak_congruent = fraction_summary_for_text(sum(is_cong), M.n);
    M.resid_peak_antiparallel = fraction_summary_for_text(sum(is_anti), M.n);

    if mod(nbins,2) == 0
        M.chance_relational_peak_frac = 2 / nbins;
        M.chance_congruent_peak_frac = 1 / nbins;
        M.chance_antiparallel_peak_frac = 1 / nbins;
    else
        M.chance_relational_peak_frac = 1 / nbins;
        M.chance_congruent_peak_frac = 1 / nbins;
        M.chance_antiparallel_peak_frac = NaN;
    end

    M.resid_peak_relational.p_binomial_right_tail = ...
        binomial_right_tail_pvalue(M.resid_peak_relational.k, M.n, M.chance_relational_peak_frac);

    M.resid_peak_congruent.p_binomial_right_tail = ...
        binomial_right_tail_pvalue(M.resid_peak_congruent.k, M.n, M.chance_congruent_peak_frac);

    if isfinite(M.chance_antiparallel_peak_frac)
        M.resid_peak_antiparallel.p_binomial_right_tail = ...
            binomial_right_tail_pvalue(M.resid_peak_antiparallel.k, M.n, M.chance_antiparallel_peak_frac);
    else
        M.resid_peak_antiparallel.p_binomial_right_tail = NaN;
    end
end

function S = summarize_vector_for_text(x)
    x = x(:);
    x = x(isfinite(x));

    S = struct();
    S.n = numel(x);

    if isempty(x)
        S.mean = NaN;
        S.median = NaN;
        S.std = NaN;
        S.sem = NaN;
        S.iqr = NaN;
        S.p25 = NaN;
        S.p75 = NaN;
        return;
    end

    S.mean = mean(x, 'omitnan');
    S.median = median(x, 'omitnan');
    S.std = std(x, 'omitnan');
    S.sem = S.std / sqrt(max(1, S.n));
    S.iqr = iqr_local(x);
    S.p25 = prctile(x, 25);
    S.p75 = prctile(x, 75);
end

function F = fraction_summary_for_text(k, n)
    F = struct();
    F.k = k;
    F.n = n;

    if n <= 0
        F.frac = NaN;
        F.percent = NaN;
        F.sem = NaN;
        F.sem_percent = NaN;
        return;
    end

    F.frac = k / n;
    F.percent = 100 * F.frac;
    F.sem = sqrt(F.frac * (1 - F.frac) / n);
    F.sem_percent = 100 * F.sem;
end

function pval = binomial_right_tail_pvalue(k, n, p0)
    if n <= 0 || ~isfinite(p0) || p0 < 0 || p0 > 1
        pval = NaN;
        return;
    end

    k = max(0, min(n, round(k)));

    % Exact right-tail binomial probability P[X >= k], computed without
    % Statistics Toolbox dependencies.
    probs = zeros(n-k+1, 1);
    idx = 1;
    for x = k:n
        probs(idx) = nchoosek_stable(n, x) * (p0^x) * ((1-p0)^(n-x));
        idx = idx + 1;
    end

    pval = min(1, sum(probs, 'omitnan'));
end

function c = nchoosek_stable(n, k)
    if k < 0 || k > n
        c = 0;
        return;
    end

    k = min(k, n-k);

    if k == 0
        c = 1;
        return;
    end

    c = exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1));
end


function C = class_fraction_stats(H, E, nbins)
    H = H(:);
    E = E(:);
    n = numel(H);

    is_cong = H == E;
    is_anti = false(size(H));

    if mod(nbins,2) == 0
        half = nbins / 2;
        anti_E = mod(H - 1 + half, nbins) + 1;
        is_anti = E == anti_E;
    end

    is_incongruent = ~(is_cong | is_anti);

    C = struct();
    C.n = n;

    C.n_congruent = sum(is_cong);
    C.n_antiparallel = sum(is_anti);
    C.n_incongruent = sum(is_incongruent);
    C.n_other = C.n_incongruent;

    [C.frac_congruent, C.sem_congruent] = binomial_fraction_sem(C.n_congruent, n);
    [C.frac_antiparallel, C.sem_antiparallel] = binomial_fraction_sem(C.n_antiparallel, n);
    [C.frac_incongruent, C.sem_incongruent] = binomial_fraction_sem(C.n_incongruent, n);

    C.frac_other = C.frac_incongruent;
    C.sem_other = C.sem_incongruent;

    C.chance_congruent = 1 / nbins;

    if mod(nbins,2) == 0
        C.chance_antiparallel = 1 / nbins;
    else
        C.chance_antiparallel = NaN;
    end

    C.chance_incongruent = 1 - C.chance_congruent - C.chance_antiparallel;
    C.chance_other = C.chance_incongruent;
end

function [p, se] = binomial_fraction_sem(k, n)
    if n <= 0
        p = NaN;
        se = NaN;
        return;
    end

    p = k / n;
    se = sqrt(p * (1 - p) / n);
end

function [d_bins, d_deg] = circular_bin_difference(hbin, ebin, nbins, bin_deg)
    raw = abs(hbin - ebin);
    d_bins = min(raw, nbins - raw);
    d_deg = d_bins * bin_deg;
end

function metric_names = basic_metric_names()
    metric_names = { ...
        'mean_FR', ...
        'raw_amp', ...
        'add_amp', ...
        'resid_amp', ...
        'resid_max_abs', ...
        'resid_var_frac', ...
        'add_to_resid_amp_ratio', ...
        'resid_to_add_amp_ratio', ...
        'resid_abs_mean_relational', ...
        'resid_abs_mean_offaxis', ...
        'resid_abs_relational_to_offaxis_ratio', ...
        'resid_abs_relational_minus_offaxis', ...
        'resid_energy_relational_frac', ...
        'hand_mod_amp', ...
        'eye_mod_amp', ...
        'hand_vector_strength', ...
        'eye_vector_strength', ...
        'surface_sparsity'};
end

function metric_labels = basic_metric_labels()
    metric_labels = { ...
        'Mean firing rate across ME surface', ...
        'Raw ME surface amplitude', ...
        'Additive component amplitude', ...
        'Residual amplitude', ...
        'Max |residual|', ...
        'Residual variance fraction', ...
        'Additive / residual amplitude ratio', ...
        'Residual / additive amplitude ratio', ...
        'Mean |residual| on congruent+anti axes', ...
        'Mean |residual| off congruent+anti axes', ...
        'Relational/off-axis |residual| ratio', ...
        'Relational minus off-axis |residual|', ...
        'Residual energy fraction on congruent+anti axes', ...
        'Hand marginal amplitude', ...
        'Eye marginal amplitude', ...
        'Hand vector strength', ...
        'Eye vector strength', ...
        'Surface sparsity'};
end

function C = compute_basic_metric_comparisons(insep, sep, nperm)
    metric_names = basic_metric_names();
    metric_labels = basic_metric_labels();

    C = struct();
    C.metric_names = metric_names;
    C.metric_labels = metric_labels;
    C.nperm = nperm;

    for i = 1:numel(metric_names)
        f = metric_names{i};

        x1 = insep.(f);
        x2 = sep.(f);

        x1 = x1(:);
        x2 = x2(:);

        x1 = x1(isfinite(x1));
        x2 = x2(isfinite(x2));

        out = struct();
        out.name = f;
        out.label = metric_labels{i};
        out.n_inseparable = numel(x1);
        out.n_separable = numel(x2);

        out.mean_inseparable = mean(x1, 'omitnan');
        out.mean_separable = mean(x2, 'omitnan');
        out.median_inseparable = median(x1, 'omitnan');
        out.median_separable = median(x2, 'omitnan');
        out.sem_inseparable = std(x1, 'omitnan') / sqrt(max(1, numel(x1)));
        out.sem_separable = std(x2, 'omitnan') / sqrt(max(1, numel(x2)));

        out.delta_mean = out.mean_inseparable - out.mean_separable;
        out.delta_median = out.median_inseparable - out.median_separable;

        out.p_perm_mean = permutation_test_difference_of_means(x1, x2, nperm);
        out.p_ranksum = ranksum_if_available(x1, x2);

        C.(f) = out;
    end
end

function p = permutation_test_difference_of_means(x1, x2, nperm)
    x1 = x1(:);
    x2 = x2(:);
    x1 = x1(isfinite(x1));
    x2 = x2(isfinite(x2));

    if numel(x1) < 2 || numel(x2) < 2
        p = NaN;
        return;
    end

    obs = mean(x1) - mean(x2);
    x = [x1; x2];
    n1 = numel(x1);
    n = numel(x);

    count = 0;
    for k = 1:nperm
        idx = randperm(n);
        xp1 = x(idx(1:n1));
        xp2 = x(idx(n1+1:end));
        d = mean(xp1) - mean(xp2);
        if abs(d) >= abs(obs)
            count = count + 1;
        end
    end

    p = (count + 1) / (nperm + 1);
end

function p = ranksum_if_available(x1, x2)
    x1 = x1(:);
    x2 = x2(:);
    x1 = x1(isfinite(x1));
    x2 = x2(isfinite(x2));

    if numel(x1) < 2 || numel(x2) < 2
        p = NaN;
        return;
    end

    if exist('ranksum', 'file') == 2
        p = ranksum(x1, x2);
    else
        p = NaN;
    end
end

function q = iqr_local(x)
    x = x(:);
    x = x(isfinite(x));

    if isempty(x)
        q = NaN;
        return;
    end

    q = prctile(x, 75) - prctile(x, 25);
end

% ============================================================
% 26-panel gradient figure
% ============================================================
function plot_gradient_axis_distribution_and_6_examples_full(insep, nbins, area, random_example_seed, num_random_examples)
    % This figure keeps the same 4-column arrangement for every example:
    %   col 1: overall F heatmap/quiver
    %   col 2: overall F gradient-axis polar plot
    %   col 3: residual R heatmap/quiver
    %   col 4: residual R gradient-axis polar plot
    %
    % It now shows random globally inseparable example neurons.
    % The random seed is controlled by RandomExampleSeed at the top of the
    % main function. This makes the random selection reproducible.

    n_examples_available = numel(insep.session);
    n_examples_to_show = min(num_random_examples, n_examples_available);

    rng(random_example_seed, 'twister');

    if n_examples_to_show > 0
        example_idx = randperm(n_examples_available, n_examples_to_show).';
    else
        example_idx = [];
    end

    example_type = repmat({'Random inseparable neuron'}, numel(example_idx), 1);

    n_rows = 1 + max(1, num_random_examples);
    t = tiledlayout(n_rows,4,'TileSpacing','loose','Padding','loose');

    ax1 = nexttile(t, [1 2]);
    plot_gradient_axis_histogram(ax1, insep.axis_angle);
    title(ax1, sprintf('%s inseparable: overall/raw gradient axis distribution', area), ...
        'Interpreter','none');

    ax2 = nexttile(t, [1 2]);
    plot_gradient_axis_histogram(ax2, insep.resid_axis_angle);
    title(ax2, sprintf('%s inseparable: residual gradient axis distribution', area), ...
        'Interpreter','none');

    for k = 1:num_random_examples
        ax_Fq = nexttile(t);
        ax_Fp = nexttile(t);
        ax_Rq = nexttile(t);
        ax_Rp = nexttile(t);

        if k <= numel(example_idx)
            idx = example_idx(k);
            label = example_type{k};

            plot_surface_gradient_quiver(ax_Fq, insep.F_surface{idx}, nbins, false);
            title(ax_Fq, sprintf('%d. Overall F quiver\n%s | seed=%d | s%d n%d', ...
                k, label, random_example_seed, insep.session(idx), insep.neuron(idx)), ...
                'Interpreter','none');

            plot_gradient_axis_polar_cartesian(ax_Fp, insep.axis_angle(idx), insep.axis_strength(idx));
            title(ax_Fp, sprintf('Overall F axis\n%.1f deg, R=%.2f', ...
                rad2deg(insep.axis_angle(idx)), insep.axis_strength(idx)), ...
                'Interpreter','none');

            plot_surface_gradient_quiver(ax_Rq, insep.R_surface{idx}, nbins, true);
            title(ax_Rq, sprintf('%d. Residual R quiver\n%s | seed=%d | s%d n%d', ...
                k, label, random_example_seed, insep.session(idx), insep.neuron(idx)), ...
                'Interpreter','none');

            plot_gradient_axis_polar_cartesian(ax_Rp, insep.resid_axis_angle(idx), insep.resid_axis_strength(idx));
            title(ax_Rp, sprintf('Residual R axis\n%.1f deg, R=%.2f', ...
                rad2deg(insep.resid_axis_angle(idx)), insep.resid_axis_strength(idx)), ...
                'Interpreter','none');
        else
            draw_empty_axis(ax_Fq, 'No random example available');
            draw_empty_axis(ax_Fp, 'No overall polar available');
            draw_empty_axis(ax_Rq, 'No residual example available');
            draw_empty_axis(ax_Rp, 'No residual polar available');
        end
    end
end

function idx_out = choose_top_examples(idx_in, score, nmax)
    idx_in = idx_in(:);
    idx_in = idx_in(isfinite(idx_in));

    if isempty(idx_in)
        idx_out = [];
        return;
    end

    sc = score(idx_in);
    sc(~isfinite(sc)) = -Inf;

    [~, ord] = sort(sc, 'descend');
    idx_out = idx_in(ord(1:min(nmax, numel(ord))));
end

function draw_empty_axis(ax, msg)
    cla(ax);
    set(ax, 'Visible','off');
    text(ax, 0.5, 0.5, msg, ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end


function plot_surface_gradient_quiver(ax, F, nbins, is_residual)
    cla(ax);

    dir_labels = {'-150','-90','-30','30','90','150'};

    imagesc(ax, F);
    axis(ax, 'xy');
    axis(ax, 'square');

    if is_residual
        maxabs = max(abs(F(:)), [], 'omitnan');
        if isfinite(maxabs) && maxabs > 0
            caxis(ax, [-maxabs maxabs]);
        end
    end

    colormap(ax, parula(256));
    cb = colorbar(ax);
    if is_residual
        ylabel(cb, 'Residual');
    else
        ylabel(cb, 'Firing rate');
    end

    hold(ax, 'on');

    [X, Y] = meshgrid(1:nbins, 1:nbins);

    dH = 0.5 * (circshift(F, [-1 0]) - circshift(F, [1 0]));
    dE = 0.5 * (circshift(F, [0 -1]) - circshift(F, [0 1]));

    % Heatmap convention:
    % x-axis = Eye direction
    % y-axis = Hand direction
    U = dE;
    V = dH;

    mag = sqrt(U.^2 + V.^2);
    med_mag = median(mag(:), 'omitnan');

    if isfinite(med_mag) && med_mag > 0
        U = U / med_mag * 0.35;
        V = V / med_mag * 0.35;
    end

    quiver(ax, X, Y, U, V, 0, ...
        'k', ...
        'LineWidth', 1.0, ...
        'MaxHeadSize', 0.7);

    for h = 1:nbins
        e = h;
        rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
            'EdgeColor', 'k', ...
            'LineWidth', 1.0, ...
            'LineStyle', '-');
    end

    if mod(nbins,2) == 0
        half = nbins / 2;
        for h = 1:nbins
            e = mod(h - 1 + half, nbins) + 1;
            rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
                'EdgeColor', [0.65 0 0.85], ...
                'LineWidth', 1.0, ...
                'LineStyle', '--');
        end
    end

    xlabel(ax, 'Eye direction (deg)');
    ylabel(ax, 'Hand direction (deg)');

    set(ax, ...
        'XTick', 1:nbins, ...
        'YTick', 1:nbins, ...
        'XTickLabel', dir_labels, ...
        'YTickLabel', dir_labels, ...
        'FontSize', 8);
end


function plot_gradient_axis_polar_cartesian(ax, axis_angle, axis_strength)
    cla(ax);
    hold(ax, 'on');

    set(ax, ...
        'Visible','off', ...
        'XLim',[-1.35 1.35], ...
        'YLim',[-1.35 1.35], ...
        'DataAspectRatio',[1 1 1], ...
        'PlotBoxAspectRatio',[1 1 1]);

    if ~isfinite(axis_angle)
        text(ax, 0, 0, 'No finite axis', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle');
        return;
    end

    if ~isfinite(axis_strength)
        axis_strength = 1;
    end

    axis_angle = wrapToPiHalf(axis_angle);
    r = max(0.15, min(1.0, axis_strength));

    th = linspace(0, 2*pi, 361);
    plot(ax, cos(th), sin(th), 'k-', 'LineWidth', 0.9);

    plot(ax, [-1 1], [0 0], 'k:', 'LineWidth', 0.8);
    plot(ax, [0 0], [-1 1], 'k:', 'LineWidth', 0.8);

    plot(ax, [0 cos(pi/4)], [0 sin(pi/4)], '--', ...
        'Color', [0.65 0 0.85], 'LineWidth', 0.9);

    plot(ax, [0 cos(-pi/4)], [0 sin(-pi/4)], '--', ...
        'Color', [0.65 0 0.85], 'LineWidth', 0.9);

    x1 = r * cos(axis_angle);
    y1 = r * sin(axis_angle);
    x2 = -r * cos(axis_angle);
    y2 = -r * sin(axis_angle);

    plot(ax, [x2 x1], [y2 y1], 'r-', 'LineWidth', 2.5);
    plot(ax, [x1 x2], [y1 y2], 'ro', ...
        'MarkerFaceColor','r', ...
        'MarkerSize',4);

    % Consistent with heatmap:
    % x-axis = Eye, y-axis = Hand
    text(ax, 1.08, 0, 'Eye', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'FontSize',8);

    text(ax, 0, 1.08, 'Hand', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',8);

    text(ax, 0.76, 0.76, '+45°', ...
        'HorizontalAlignment','center', ...
        'FontSize',7, ...
        'Color', [0.45 0 0.65]);

    text(ax, 0.76, -0.76, '-45°', ...
        'HorizontalAlignment','center', ...
        'FontSize',7, ...
        'Color', [0.45 0 0.65]);
end


function plot_gradient_axis_histogram(ax, axis_angles)
    cla(ax);

    axis_angles = axis_angles(:);
    axis_angles = axis_angles(isfinite(axis_angles));
    axis_angles = arrayfun(@wrapToPiHalf, axis_angles);

    edges = linspace(-pi/2, pi/2, 13);
    counts = histcounts(axis_angles, edges);
    centers = edges(1:end-1) + diff(edges)/2;

    bar(ax, rad2deg(centers), counts, 1.0, ...
        'FaceColor', [0.35 0.35 0.35], ...
        'EdgeColor', 'k');

    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    xline(ax, 0, 'k-', 'LineWidth', 1.2);
    xline(ax, 45, '--', 'Color', [0.65 0 0.85], 'LineWidth', 1.2);
    xline(ax, -45, '--', 'Color', [0.65 0 0.85], 'LineWidth', 1.2);
    xline(ax, 90, 'k:', 'LineWidth', 1.0);
    xline(ax, -90, 'k:', 'LineWidth', 1.0);

    xlabel(ax, 'Gradient resultant axis (deg)');
    ylabel(ax, 'Neuron count');
    xlim(ax, [-90 90]);
        
    set(ax, 'XTick', [-90 -45 0 45 90], ...
        'XTickLabel', {'Hand axis', '-45°', 'Eye axis', '+45°', 'Hand axis'}, ...
        'FontSize', 10);

    text(ax, 0.03, 0.95, sprintf('n = %d', numel(axis_angles)), ...
        'Units','normalized', ...
        'FontSize', 9, ...
        'FontWeight','bold', ...
        'BackgroundColor','w', ...
        'Margin', 2);
end

% ============================================================
% Standard plots
% ============================================================
function plot_summary_fraction_bars(Stats, nbins, area)
    labels = {'Congruent','Anti-parallel','Incongruent'};

    metrics = {'raw_peak','resid_peak','PD_pair'};
    metric_names = {'Raw peak','Residual peak','PD pair'};

    tiledlayout(3,1,'TileSpacing','loose','Padding','loose');

    for i = 1:numel(metrics)
        ax = nexttile;
        m = metrics{i};

        in = Stats.inseparable.(m);
        sp = Stats.separable.(m);

        vals = [
            in.frac_congruent, in.frac_antiparallel, in.frac_incongruent;
            sp.frac_congruent, sp.frac_antiparallel, sp.frac_incongruent
        ] * 100;

        errs = [
            in.sem_congruent, in.sem_antiparallel, in.sem_incongruent;
            sp.sem_congruent, sp.sem_antiparallel, sp.sem_incongruent
        ] * 100;

        b = bar(vals');
        hold on; box on; grid on;

        for j = 1:numel(b)
            errorbar(b(j).XEndPoints, vals(j,:), errs(j,:), ...
                'k.', 'LineWidth', 1.2, 'CapSize', 8);
        end

        chance_diag = 100 / nbins;
        chance_incongruent = 100 * (1 - 2/nbins);

        yline(chance_diag, 'k--', 'LineWidth', 1.0);
        text(1.25, chance_diag, sprintf(' chance congruent/anti = %.1f%%', chance_diag), ...
            'VerticalAlignment','bottom', ...
            'FontSize', 8, ...
            'BackgroundColor','w');

        yline(chance_incongruent, 'k:', 'LineWidth', 1.2);
        text(2.15, chance_incongruent, sprintf(' chance incongruent = %.1f%%', chance_incongruent), ...
            'VerticalAlignment','bottom', ...
            'FontSize', 8, ...
            'BackgroundColor','w');

        xticks(1:3);
        xticklabels(labels);
        xtickangle(15);
        ylabel('% neurons');
        ylim([0 max(100, 1.15 * max(vals(:) + errs(:)))]);
        title(sprintf('%s: %s', area, metric_names{i}), 'Interpreter','none');

        if i == 1
            legend({'Inseparable','Separable'}, ...
                'Location','northeast');
        end

        set(ax, 'FontSize', 10);
    end
end

function plot_basic_metric_comparison_bars(insep, sep, comparisons, area)
    metric_names = basic_metric_names();
    metric_labels = basic_metric_labels();

    tiledlayout(numel(metric_names), 1, ...
        'TileSpacing','loose', ...
        'Padding','loose');

    for i = 1:numel(metric_names)
        ax = nexttile;
        f = metric_names{i};

        x1 = insep.(f);
        x2 = sep.(f);

        x1 = x1(:);
        x2 = x2(:);

        x1 = x1(isfinite(x1));
        x2 = x2(isfinite(x2));

        vals = [mean(x1, 'omitnan'), mean(x2, 'omitnan')];
        errs = [
            std(x1, 'omitnan') / sqrt(max(1, numel(x1))), ...
            std(x2, 'omitnan') / sqrt(max(1, numel(x2)))
        ];

        b = bar(1:2, vals, 0.65);
        hold on; box on; grid on;

        errorbar(1:2, vals, errs, ...
            'k.', 'LineWidth', 1.2, 'CapSize', 8);

        jitter_width = 0.11;
        rng(1 + i);

        xj1 = 1 + jitter_width * randn(size(x1));
        xj2 = 2 + jitter_width * randn(size(x2));

        scatter(xj1, x1, 15, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.35, ...
            'MarkerEdgeAlpha', 0.35);

        scatter(xj2, x2, 15, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.20, ...
            'MarkerEdgeAlpha', 0.20);

        set(b, 'FaceAlpha', 0.65);

        xticks(1:2);
        xticklabels({'Inseparable','Separable'});
        ylabel(metric_labels{i}, 'Interpreter','none');

        c = comparisons.(f);
        ttl = sprintf('%s: %s | p_{perm}=%.3g', area, metric_labels{i}, c.p_perm_mean);

        if isfinite(c.p_ranksum)
            ttl = sprintf('%s, p_{rank}=%.3g', ttl, c.p_ranksum);
        end

        title(ttl, 'Interpreter','tex');

        ylo = min([x1; x2; vals(:)-errs(:)], [], 'omitnan');
        yhi = max([x1; x2; vals(:)+errs(:)], [], 'omitnan');

        if ~isfinite(ylo) || ~isfinite(yhi) || ylo == yhi
            ylo = 0;
            yhi = 1;
        end

        pad = 0.10 * (yhi - ylo);
        ylim([ylo - pad, yhi + pad]);

        set(ax, 'FontSize', 9);
    end
end

function plot_angular_difference_distribution(ax, diff_insep_deg, diff_sep_deg, nbins)
    max_diff_bins = floor(nbins / 2);
    bin_deg = 360 / nbins;

    edges = (-0.5:1:(max_diff_bins + 0.5)) * bin_deg;
    centers = (0:max_diff_bins) * bin_deg;

    c_in = histcounts(diff_insep_deg, edges);
    c_sp = histcounts(diff_sep_deg, edges);

    p_in = c_in / max(1, sum(c_in));
    p_sp = c_sp / max(1, sum(c_sp));

    se_in = sqrt(p_in .* (1 - p_in) / max(1, sum(c_in)));
    se_sp = sqrt(p_sp .* (1 - p_sp) / max(1, sum(c_sp)));

    vals = [p_in(:), p_sp(:)] * 100;
    errs = [se_in(:), se_sp(:)] * 100;

    b = bar(ax, centers, vals, 'grouped');
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    for j = 1:numel(b)
        errorbar(ax, b(j).XEndPoints, vals(:,j), errs(:,j), ...
            'k.', 'LineWidth', 1.2, 'CapSize', 8);
    end

    xlabel(ax, '|Hand peak - Eye peak| circular angular difference');
    ylabel(ax, '% neurons');
    xticks(ax, centers);
    xticklabels(ax, arrayfun(@(x)sprintf('%d°', x), centers, 'UniformOutput', false));
    legend(ax, {'Inseparable','Separable'}, 'Location','northeast');
    set(ax, 'FontSize', 12);
end

function plot_count_heatmap(ax, C, nbins, ttl, xlab, ylab)
    imagesc(ax, C);
    axis(ax, 'xy');
    axis(ax, 'square');
    colormap(ax, parula(256));

    cb = colorbar(ax);
    ylabel(cb, 'Neuron count');

    xlabel(ax, xlab);
    ylabel(ax, ylab);
    title(ax, ttl, 'Interpreter','none');

    set(ax, 'XTick', 1:nbins, ...
            'YTick', 1:nbins, ...
            'FontSize', 11);

    hold(ax, 'on');

    for h = 1:nbins
        e = h;
        rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
            'EdgeColor', 'k', ...
            'LineWidth', 1.2, ...
            'LineStyle', '-');
    end

    if mod(nbins,2) == 0
        half = nbins / 2;
        for h = 1:nbins
            e = mod(h - 1 + half, nbins) + 1;
            rectangle(ax, 'Position', [e-0.5, h-0.5, 1, 1], ...
                'EdgeColor', [0.65 0 0.85], ...
                'LineWidth', 1.2, ...
                'LineStyle', '--');
        end
    end

    max_count = max(C(:));

    for h = 1:nbins
        for e = 1:nbins
            if C(h,e) > 0
                if max_count > 0 && C(h,e) > 0.55 * max_count
                    txt_color = 'k';
                else
                    txt_color = 'w';
                end

                text(ax, e, h, sprintf('%d', C(h,e)), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'FontSize', 10, ...
                    'FontWeight','bold', ...
                    'Color', txt_color);
            end
        end
    end
end

% ============================================================
% Math helpers
% ============================================================
function q = bh_fdr(p)
    p = p(:);
    q = nan(size(p));
    good = isfinite(p) & p >= 0 & p <= 1;

    pg = p(good);
    mg = numel(pg);

    if mg == 0
        return;
    end

    [ps, idx] = sort(pg, 'ascend');
    ranks = (1:mg)';

    qs = ps .* (mg ./ ranks);

    for i = mg-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end

    qs = min(qs, 1);

    qg = nan(size(pg));
    qg(idx) = qs;

    q(good) = qg;
end

function S2 = smooth_circular_surface(S)
    k = [0.25 0.5 0.25];

    S2 = zeros(size(S));

    for dh = -1:1
        for de = -1:1
            wh = k(dh+2);
            we = k(de+2);
            S2 = S2 + wh * we * circshift(S, [dh de]);
        end
    end
end


function [axis_angle, resultant_strength] = gradient_resultant_axis(S)
    nbins = size(S,1);
    dtheta = 2*pi / nbins;

    dH = (circshift(S, [-1 0]) - circshift(S, [1 0])) / (2*dtheta);
    dE = (circshift(S, [0 -1]) - circshift(S, [0 1])) / (2*dtheta);

    mag = sqrt(dH.^2 + dE.^2);

    valid = isfinite(dH) & isfinite(dE) & isfinite(mag) & mag > 0;

    if ~any(valid(:))
        axis_angle = NaN;
        resultant_strength = NaN;
        return;
    end

    % Consistent with heatmap/quiver:
    % x-axis = Eye direction
    % y-axis = Hand direction
    %
    % Therefore angle = atan2(y-component, x-component)
    %                = atan2(dF/dHand, dF/dEye)
    phi = atan2(dH(valid), dE(valid));

    w = mag(valid);

    % Axis, not signed vector, hence doubled angle.
    C = sum(w .* cos(2*phi));
    S2 = sum(w .* sin(2*phi));

    axis_angle = 0.5 * atan2(S2, C);
    axis_angle = wrapToPiHalf(axis_angle);

    resultant_strength = sqrt(C.^2 + S2.^2) / sum(w);
end


function a = wrapToPiHalf(a)
    a = mod(a + pi/2, pi) - pi/2;
end

% ============================================================
% Console output
% ============================================================
function print_stats_to_console(area, Stats, nbins)
    fprintf('\n[%s] Extended separable/inseparable summary:\n', area);

    print_one_group('Inseparable', Stats.inseparable, nbins);
    print_one_group('Separable', Stats.separable, nbins);

    print_basic_metric_comparisons(area, Stats.basic_metric_comparisons);
end

function print_one_group(name, S, nbins)
    fprintf('\n  %s neurons: n=%d\n', name, S.n);

    print_metric_line('Raw peak', S.raw_peak, nbins);
    print_metric_line('Residual peak', S.resid_peak, nbins);
    print_metric_line('PD pair', S.PD_pair, nbins);

    if isfield(S, 'manuscript_numbers')
        print_manuscript_numbers(name, S.manuscript_numbers, nbins);
    end
end

function print_metric_line(name, C, nbins)
    fprintf('    %s:\n', name);
    fprintf('      Congruent:     %d/%d = %.1f%% +/- %.1f%%\n', ...
        C.n_congruent, C.n, 100*C.frac_congruent, 100*C.sem_congruent);
    fprintf('      Anti-parallel: %d/%d = %.1f%% +/- %.1f%%\n', ...
        C.n_antiparallel, C.n, 100*C.frac_antiparallel, 100*C.sem_antiparallel);
    fprintf('      Incongruent:   %d/%d = %.1f%% +/- %.1f%%\n', ...
        C.n_incongruent, C.n, 100*C.frac_incongruent, 100*C.sem_incongruent);
    fprintf('      Chance congruent/anti class: %.1f%%\n', 100/nbins);
    fprintf('      Chance incongruent class: %.1f%%\n', 100*(1 - 2/nbins));
end


function print_manuscript_numbers(name, M, nbins)
    fprintf('\n    Manuscript-ready %s numbers:\n', name);
    fprintf('      Component amplitudes:\n');
    fprintf('        Additive amplitude: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.add_amp.mean, M.add_amp.sem, M.add_amp.median, M.add_amp.n);
    fprintf('        Residual amplitude: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.resid_amp.mean, M.resid_amp.sem, M.resid_amp.median, M.resid_amp.n);
    fprintf('        Additive/residual amplitude ratio: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.add_to_resid_amp_ratio.mean, M.add_to_resid_amp_ratio.sem, ...
        M.add_to_resid_amp_ratio.median, M.add_to_resid_amp_ratio.n);
    fprintf('        Residual/additive amplitude ratio: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.resid_to_add_amp_ratio.mean, M.resid_to_add_amp_ratio.sem, ...
        M.resid_to_add_amp_ratio.median, M.resid_to_add_amp_ratio.n);
    fprintf('        Residual variance fraction: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.resid_var_frac.mean, M.resid_var_frac.sem, M.resid_var_frac.median, ...
        M.resid_var_frac.n);

    fprintf('      Residual relational peak locations using max |R|:\n');
    fprintf('        Congruent axis: %d/%d = %.1f%% +/- %.1f%% SEM; chance=%.1f%%; p_right=%.4g\n', ...
        M.resid_peak_congruent.k, M.resid_peak_congruent.n, ...
        M.resid_peak_congruent.percent, M.resid_peak_congruent.sem_percent, ...
        100*M.chance_congruent_peak_frac, M.resid_peak_congruent.p_binomial_right_tail);
    fprintf('        Anti-congruent axis: %d/%d = %.1f%% +/- %.1f%% SEM; chance=%.1f%%; p_right=%.4g\n', ...
        M.resid_peak_antiparallel.k, M.resid_peak_antiparallel.n, ...
        M.resid_peak_antiparallel.percent, M.resid_peak_antiparallel.sem_percent, ...
        100*M.chance_antiparallel_peak_frac, M.resid_peak_antiparallel.p_binomial_right_tail);
    fprintf('        Congruent OR anti-congruent axes: %d/%d = %.1f%% +/- %.1f%% SEM; chance=%.1f%%; p_right=%.4g\n', ...
        M.resid_peak_relational.k, M.resid_peak_relational.n, ...
        M.resid_peak_relational.percent, M.resid_peak_relational.sem_percent, ...
        100*M.chance_relational_peak_frac, M.resid_peak_relational.p_binomial_right_tail);

    fprintf('      Residual amplitude on relational axes:\n');
    fprintf('        Mean |R| on congruent+anti axes: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.resid_abs_mean_relational.mean, M.resid_abs_mean_relational.sem, ...
        M.resid_abs_mean_relational.median, M.resid_abs_mean_relational.n);
    fprintf('        Mean |R| off these axes: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.resid_abs_mean_offaxis.mean, M.resid_abs_mean_offaxis.sem, ...
        M.resid_abs_mean_offaxis.median, M.resid_abs_mean_offaxis.n);
    fprintf('        Relational/off-axis |R| ratio: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d\n', ...
        M.resid_abs_relational_to_offaxis_ratio.mean, ...
        M.resid_abs_relational_to_offaxis_ratio.sem, ...
        M.resid_abs_relational_to_offaxis_ratio.median, ...
        M.resid_abs_relational_to_offaxis_ratio.n);
    fprintf('        Residual energy fraction on congruent+anti axes: mean=%.4g +/- %.4g SEM, median=%.4g, n=%d; uniform chance=%.4g\n', ...
        M.resid_energy_relational_frac.mean, M.resid_energy_relational_frac.sem, ...
        M.resid_energy_relational_frac.median, M.resid_energy_relational_frac.n, ...
        M.chance_relational_peak_frac);
end


function print_basic_metric_comparisons(area, C)
    fprintf('\n[%s] Basic scalar metric comparisons: inseparable vs separable\n', area);

    metric_names = C.metric_names;

    for i = 1:numel(metric_names)
        f = metric_names{i};
        x = C.(f);

        fprintf('  %s:\n', x.label);
        fprintf('    Inseparable: n=%d, mean=%.4g, median=%.4g, SEM=%.4g\n', ...
            x.n_inseparable, x.mean_inseparable, x.median_inseparable, x.sem_inseparable);
        fprintf('    Separable:   n=%d, mean=%.4g, median=%.4g, SEM=%.4g\n', ...
            x.n_separable, x.mean_separable, x.median_separable, x.sem_separable);
        fprintf('    Delta mean: %.4g, p_perm=%.4g, p_ranksum=%.4g\n', ...
            x.delta_mean, x.p_perm_mean, x.p_ranksum);
    end
end
