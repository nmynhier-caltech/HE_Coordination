function check_FR_gaussianity(area, savepath, varargin)
% Per-neuron diagnostics focused on two ANOVA assumptions:
%   (A) Error normality (pooled, per-cell studentized residuals; Lilliefors)
%   (B) Mean–variance coupling via within-neuron spread–level slope
%       (log IQR vs log median across that neuron's (hand, eye) cells)
%
% For each neuron n:
%   • Build within-(hand,eye)-cell residuals: r_ij(t) = FR - mean_cell.
%   • Studentize by that cell's SD: z_ij(t) = r_ij(t) / sd_cell.
%   • Normality: Lilliefors on pooled z_ij(t) for the neuron.
%   • Spread–level slope: robust regression of log10(IQR_cell) on log10(median_cell).
%
% Plots (one 8.5×11 PDF per neuron; 2 panels):
%   (1) Spread–level within neuron: log10(IQR) vs log10(median), robust slope,
%       with Poisson-like reference slope 0.5.
%   (2) Q–Q plot of pooled studentized residuals (vs N(0,1)) with Lilliefors p.
%
% Saves:
%   models/<area>_ANOVA_Assumption_Diagnostics_perNeuron.mat   (struct 'Diag')
%   plots/anova_diag_single/neuron_###_anova_diag.pdf           (one file per neuron)
%
% Args (name/value):
%   'NBins'            (default 6)      : number of angle bins per effector
%   'EdgeShift'        (default 0.1)    : circular bin shift
%   'MinTrialsPerCell' (default 3)      : minimum trials per (H,E) cell to include
%   'MakePlots'        (default true)   : generate per-neuron PDFs
%   'Neurons'          (default 'all')  : vector of neuron indices to run/plot

    p = inputParser;
    addParameter(p, 'NBins', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=3);
    addParameter(p, 'EdgeShift', 0.1, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'MinTrialsPerCell', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(p, 'MakePlots', true, @(x)islogical(x)||ismember(x,[0 1]));
    addParameter(p, 'Neurons', 'all');
    parse(p, varargin{:});

    nbins   = p.Results.NBins;
    edge_shift = p.Results.EdgeShift;
    min_trials_per_cell = p.Results.MinTrialsPerCell;
    make_plots = logical(p.Results.MakePlots);

    % ----------------------------
    % Load data & subset to ME trials
    % ----------------------------
    S = load(processed_data_file(savepath, area));
    D = S.data;

    isME = D.TPi2 ~= 7 & D.TPi1 ~= 7;
    FR      = D.FR(isME, :);      % [T x N]
    TP1     = D.TP1(isME, :);     % hand (x,y)
    TP2     = D.TP2(isME, :);     % eye  (x,y)
    [T, N]  = size(FR);

    if isnumeric(p.Results.Neurons)
        neur_list = p.Results.Neurons(:)'; 
        neur_list = neur_list(neur_list>=1 & neur_list<=N);
        if isempty(neur_list), neur_list = 1:N; end
    else
        neur_list = 1:N;
    end

    % ----------------------------
    % Angle binning (fixed across neurons)
    % ----------------------------
    edges_theta = linspace(-pi, pi, nbins+1) + edge_shift;
    thetaH = atan2(TP1(:,2), TP1(:,1));
    thetaE = atan2(TP2(:,2), TP2(:,1));
    Hb     = discretize(thetaH, edges_theta);
    Eb     = discretize(thetaE, edges_theta);

    % Per-cell trial indices
    cell_idx = cell(nbins, nbins);
    cell_id_of_trial = nan(T,1);
    cid = 0;
    for i = 1:nbins
        for j = 1:nbins
            cid = cid + 1;
            idx = find(Hb==i & Eb==j);
            cell_idx{i,j} = idx;
            cell_id_of_trial(idx) = cid;
        end
    end

    % ----------------------------
    % Outputs (per neuron)
    % ----------------------------
    Diag = struct();
    Diag.area       = area;
    Diag.nbins      = nbins;
    Diag.edge_shift = edge_shift;

    Diag.p_Lillie         = nan(N,1);    % normality p-value (pooled studentized residuals)
    Diag.n_residuals      = nan(N,1);
    Diag.slope_SP_within  = nan(N,1);    % robust slope: log(IQR) ~ log(median) within neuron
    Diag.n_cells_SP       = nan(N,1);    % cells in slope fit
    Diag.median_trials_per_cell = nan(N,1);

    % Store per-neuron cell summaries (support plotting & downstream)
    Diag.cell_mean   = cell(N,1);
    Diag.cell_median = cell(N,1);
    Diag.cell_sd     = cell(N,1);
    Diag.cell_iqr    = cell(N,1);
    Diag.cell_trials = cell(N,1);
    Diag.cell_ids    = cell(N,1);

    % Plot directory
    plotdir = fullfile(savepath, 'plots', 'anova_diag_single');
    if make_plots && ~exist(plotdir,'dir'), mkdir(plotdir); end

    % ----------------------------
    % Loop over requested neurons
    % ----------------------------
    for n = neur_list
        % Collect per-cell quantities for this neuron
        c_mean = nan(nbins*nbins,1);
        c_median = nan(nbins*nbins,1);
        c_sd   = nan(nbins*nbins,1);
        c_iqr  = nan(nbins*nbins,1);
        c_n    = zeros(nbins*nbins,1);
        c_ids  = (1:(nbins*nbins))';

        z_pool = [];     % pooled studentized residuals for this neuron

        cid = 0;
        for i = 1:nbins
            for j = 1:nbins
                cid = cid + 1;
                idx = cell_idx{i,j};
                if numel(idx) < min_trials_per_cell, continue; end

                x = FR(idx, n);
                x = x(isfinite(x));
                if numel(x) < min_trials_per_cell, continue; end

                mu  = mean(x, 'omitnan');
                sd  = std(x, 0, 'omitnan');
                med = median(x, 'omitnan');
                q1  = prctile(x, 25);
                q3  = prctile(x, 75);
                iq  = q3 - q1;

                c_mean(cid)   = mu;
                c_median(cid) = med;
                c_sd(cid)     = sd;
                c_iqr(cid)    = iq;
                c_n(cid)      = numel(x);

                if isfinite(sd) && sd > 0
                    z_pool = [z_pool; (x - mu) / sd]; %#ok<AGROW>
                end
            end
        end

        % Save per-cell summaries
        keep = isfinite(c_mean) & isfinite(c_sd) & isfinite(c_iqr) & c_n > 0;
        Diag.cell_mean{n}   = c_mean(keep);
        Diag.cell_median{n} = c_median(keep);
        Diag.cell_sd{n}     = c_sd(keep);
        Diag.cell_iqr{n}    = c_iqr(keep);
        Diag.cell_trials{n} = c_n(keep);
        Diag.cell_ids{n}    = c_ids(keep);

        % Per-neuron median trials per cell
        if any(keep)
            Diag.median_trials_per_cell(n) = median(c_n(keep));
        end

        % --- Normality: Lilliefors on pooled, studentized residuals
        z_pool = z_pool(isfinite(z_pool));
        Diag.n_residuals(n) = numel(z_pool);
        if numel(z_pool) >= 8
            try
                [~, pL] = lillietest(z_pool, 'Distribution','norm');
            catch
                % fallback if Statistics Toolbox is unavailable
                if std(z_pool) > 0
                    [~, pL] = kstest( (z_pool - mean(z_pool))/std(z_pool) );
                else
                    pL = NaN;
                end
            end
            Diag.p_Lillie(n) = pL;
        end

        % --- Within-neuron spread–level slope (robust)
        x_med = Diag.cell_median{n};
        y_iqr = Diag.cell_iqr{n};
        mask_sl = isfinite(x_med) & isfinite(y_iqr) & x_med>0 & y_iqr>0;
        Diag.n_cells_SP(n) = sum(mask_sl);
        slope_val = NaN;
        if sum(mask_sl) >= 3
            lx = log10(x_med(mask_sl));
            ly = log10(y_iqr(mask_sl));
            try
                b = robustfit(lx, ly);  % [intercept; slope]
                slope_val = b(2);
            catch
                X = [ones(size(lx)), lx];
                bb = X \ ly;
                slope_val = bb(2);
            end
            Diag.slope_SP_within(n) = slope_val;
        end

        % --- Plots: one page per neuron (2 panels)
        if make_plots
            fig = figure('Visible','off','PaperUnits','inches');
            set(fig,'PaperSize',[8.5 11],'PaperPosition',[0.6 0.6 7.3 9.8], 'PaperPositionMode','manual');
            set(gca, 'FontSize', 20);
            tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact','Padding','compact');

            % (1) Spread–level within neuron
            ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
            if sum(mask_sl) >= 3
                lx = log10(x_med(mask_sl));
                ly = log10(y_iqr(mask_sl));
                try
                    bsl = robustfit(lx, ly);  % [intercept; slope]
                    xs = linspace(min(lx), max(lx), 200);
                    ys = bsl(1) + bsl(2)*xs;
                    slope_val = bsl(2);
                catch
                    xs = linspace(min(lx), max(lx), 200);
                    bb = [ones(numel(lx),1) lx] \ ly;
                    ys = bb(1) + bb(2)*xs;
                    slope_val = bb(2);
                end
                scatter(ax1, lx, ly, 36, 'k', 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5);
                plot(ax1, xs, ys, 'r-', 'LineWidth', 2);
                % Poisson-ish reference slope 0.5 through same left point
                yref = ys(1) + 0.5*(xs - xs(1));
                plot(ax1, xs, yref, 'b--', 'LineWidth', 1.5);
                xlabel(ax1,'log_{10}(cell median FR)'); ylabel(ax1,'log_{10}(cell IQR)');
                title(ax1, sprintf('Neuron %d: spread–level (slope=%.2f); cells=%d', n, slope_val, Diag.n_cells_SP(n)));
                legend(ax1, {'cells','robust fit','Poisson ref (0.5)'}, 'Location','northwest');
            else
                text(ax1, 0.5,0.5,'Insufficient cells for slope','HorizontalAlignment','center');
            end

            % (2) Q–Q plot of pooled studentized residuals (vs N(0,1))
            ax2 = nexttile(tl); grid(ax2,'on'); box(ax2,'on'); hold(ax2,'on');
            if numel(z_pool) >= 8 && std(z_pool) > 0 && all(isfinite(z_pool))
                z_sorted = sort(z_pool(:));
                m = numel(z_sorted);
                % theoretical normal quantiles
                q_theory = sqrt(2) * erfinv( 2*(( (1:m)' - 0.5)/m) - 1 );
                plot(ax2, q_theory, z_sorted, 'k.', 'MarkerSize', 8);
                % matched limits (avoid axis equal pitfalls)
                lo = min([q_theory; z_sorted]);
                hi = max([q_theory; z_sorted]);
                pad = 0.05*(hi-lo + eps);
                lim = [lo-pad, hi+pad];
                plot(ax2, lim, lim, 'r-');
                xlim(ax2, lim); ylim(ax2, lim);
                xlabel(ax2,'Theoretical N(0,1) quantiles'); ylabel(ax2,'Sample quantiles (studentized)');
                title(ax2, sprintf('Q–Q (Lilliefors p=%.3g, residuals=%d)', Diag.p_Lillie(n), m));
            else
                text(ax2, 0.5,0.5,'Insufficient residuals for normality','HorizontalAlignment','center');
            end

            print(fig, '-dpdf', fullfile(plotdir, sprintf('neuron_%03d_anova_diag.pdf', n)), '-painters');
            close(fig);
        end
    end

    % Save all diagnostics
    outmat = fullfile(savepath, 'models', sprintf('%s_ANOVA_Assumption_Diagnostics_perNeuron.mat', area));
    save(outmat, 'Diag', '-v7.3');
    fprintf('Saved per-neuron diagnostics: %s\n', outmat);
end
