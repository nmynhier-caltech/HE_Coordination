function summarize_FR_gaussianity(area, savepaths, varargin)
% Summarize per-neuron diagnostics focused ONLY on:
%   (1) Error-term normality (Lilliefors on pooled, studentized residuals)
%   (2) Within-neuron spread–level slope: log10(IQR_cell) ~ log10(median_cell)
%
% Input files (per session):
%   models/<area>_ANOVA_Assumption_Diagnostics_perNeuron.mat
%     required fields in 'Diag':
%       p_Lillie (Nx1), n_residuals (Nx1),
%       slope_SP_within (Nx1), n_cells_SP (Nx1),
%       (optional) median_trials_per_cell (Nx1)
%
% Outputs:
%   plots/summary/gaussianity_diagnostics/<area>_ECDF_Lillie.pdf
%   plots/summary/gaussianity_diagnostics/<area>_slope_hist.pdf
%   plots/summary/gaussianity_diagnostics/<area>_per_session_rates.pdf
%   (the final output is optional)
%
% Saved MAT:
%   models/<area>_ANOVA_Assumption_Summary.mat
%
% Args (name/value):
%   'Alpha'              : significance threshold for Lilliefors (default 0.05)
%   'MinResids'          : minimum pooled residuals for Lilliefors ECDF (default 8)
%   'MinCellsSlope'      : minimum cells to include neuron in slope summary (default 3)
%   'MakePerSessionBars' : true/false (default true)
%
% Example:
%   summarize_FR_gaussianity('M1', {'/sess1','/sess2'}, 'Alpha',0.05);

    p = inputParser;
    addParameter(p, 'Alpha', 0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
    addParameter(p, 'MinResids', 8, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(p, 'MinCellsSlope', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(p, 'MakePerSessionBars', true, @(x)islogical(x)||ismember(x,[0 1]));
    parse(p, varargin{:});
    alpha         = p.Results.Alpha;
    minResids     = p.Results.MinResids;
    minCellsSlope = p.Results.MinCellsSlope;
    makeBars      = logical(p.Results.MakePerSessionBars);

    % ----------------------------
    % Accumulate across sessions
    % ----------------------------
    all_p_Lillie = [];
    all_slope    = [];
    all_nres     = [];
    all_ncellsSP = [];
    all_medtrpc  = [];
    sess_of_neuron = [];
    used_files     = strings(0,1);

    per_session = struct('N',[],'N_Lillie_used',[],'N_Slope_used',[], ...
                         'frac_Lillie_sig',[],'frac_slope_gt0',[]);

    for d = 1:numel(savepaths)
        sp = savepaths{d};
        fdiag = fullfile(sp, 'models', sprintf('%s_ANOVA_Assumption_Diagnostics_perNeuron.mat', area));
        if ~exist(fdiag, 'file')
            warning('Missing file: %s', fdiag);
            continue;
        end
        S = load(fdiag);
        if ~isfield(S, 'Diag')
            warning('Missing struct Diag in %s', fdiag);
            continue;
        end
        D = S.Diag;

        pL = get_field_vec(D, 'p_Lillie');
        sl = get_field_vec(D, 'slope_SP_within');
        nr = get_field_vec(D, 'n_residuals');
        nc = get_field_vec(D, 'n_cells_SP');
        mt = get_field_vec(D, 'median_trials_per_cell');

        nN = max([numel(pL), numel(sl), numel(nr), numel(nc), numel(mt)]);
        all_p_Lillie = [all_p_Lillie; pad_to(pL, nN)];
        all_slope    = [all_slope;    pad_to(sl, nN)];
        all_nres     = [all_nres;     pad_to(nr, nN)];
        all_ncellsSP = [all_ncellsSP; pad_to(nc, nN)];
        all_medtrpc  = [all_medtrpc;  pad_to(mt, nN)];
        sess_of_neuron = [sess_of_neuron; d*ones(nN,1)];
        used_files(end+1,1) = string(fdiag);

        % Per-session nominal rates at alpha
        maskL = isfinite(pL) & nr>=minResids;
        maskS = isfinite(sl) & nc>=minCellsSlope;
        per_session.N(end+1,1)              = nN;
        per_session.N_Lillie_used(end+1,1)  = sum(maskL);
        per_session.N_Slope_used(end+1,1)   = sum(maskS);
        per_session.frac_Lillie_sig(end+1,1)= safe_mean(pL(maskL) < alpha);
        per_session.frac_slope_gt0(end+1,1) = safe_mean(sl(maskS) > 0);   % CHANGED: >0 instead of >=0.5
    end

    if isempty(all_p_Lillie) && isempty(all_slope)
        error('No diagnostics found for area "%s" in provided savepaths.', area);
    end

    % Usability masks overall
    maskL_all = isfinite(all_p_Lillie) & all_nres>=minResids;
    maskS_all = isfinite(all_slope)    & all_ncellsSP>=minCellsSlope;

    % ----------------------------
    % Text summary
    % ----------------------------
    fprintf('=== ANOVA assumption summary (Lilliefors & slope only) for %s ===\n', area);
    fprintf('Sessions analyzed: %d\n', numel(used_files));
    fprintf('Normality (Lilliefors): usable neurons = %d, frac p<%.2f = %.1f%%\n', ...
        sum(maskL_all), alpha, 100*safe_mean(all_p_Lillie(maskL_all) < alpha));
    if any(maskS_all)
        medSlope = median(all_slope(maskS_all));
        fprintf('Spread–level slope: usable neurons = %d, median slope = %.3f; frac >0 = %.1f%%\n', ...
            sum(maskS_all), medSlope, 100*safe_mean(all_slope(maskS_all) > 0));
    end

    % ----------------------------
    % Output directory
    % ----------------------------
    plotdir = fullfile(savepaths{1}, '..', 'plots', 'summary', ...
        'gaussianity_diagnostics');
    if (~exist(plotdir, 'dir')), mkdir(plotdir); end

    % ----------------------------
    % 1) ECDF of Lilliefors p-values (overall) — SQUARE AXES
    % ----------------------------
    fig1 = figure('Visible','off','PaperUnits','inches');
    set(fig1,'PaperSize',[8.5 11],'PaperPosition',[1.75 2.0 5.0 5.0], 'PaperPositionMode','manual'); % square print area
    set(gca, 'FontSize', 20);
    hold on; grid on; box on;
    pl = all_p_Lillie(maskL_all);
    if isempty(pl), pl = NaN; end
    ecdf(pl);
    xlim([0 1]); ylim([0 1]); axis square;                         % <- equal axes
    xline(alpha,'r--','LineWidth',1.2);
    xlabel(sprintf('p (Lilliefors) – usable N=%d', sum(maskL_all)), 'Interpreter','none');
    ylabel('ECDF', 'Interpreter','none');
    title(sprintf('%s: Error normality (Lilliefors) ECDF', area), 'Interpreter','none');
    print(fig1, '-dpdf', fullfile(plotdir, sprintf('%s_ECDF_Lillie.pdf', area)), '-painters');
    close(fig1);

    % ----------------------------
    % 2) Histogram of within-neuron slopes — SQUARE AXES
    % ----------------------------
    fig2 = figure('Visible','off','PaperUnits','inches');
    set(fig2,'PaperSize',[8.5 11],'PaperPosition',[1.75 2.0 5.0 5.0], 'PaperPositionMode','manual'); % square print area
    set(gca, 'FontSize', 20);
    hold on; grid on; box on;
    sl = all_slope(maskS_all);
    if ~isempty(sl)
        histogram(sl, 30, 'BinLimits', [-0.5 1.5], 'Normalization', 'pdf');
        xline(0,'k--','LineWidth',1.2);
        xline(0.5,'b--','Poisson ref (0.5)','LabelVerticalAlignment','middle','LineWidth',1.2);
        xlabel('Within-neuron spread–level slope  (log_{10} IQR ~ log_{10} median)','Interpreter','none');
        ylabel('Density','Interpreter','none');
        title(sprintf('%s: Distribution of spread–level slopes (N=%d)', area, numel(sl)), 'Interpreter','none');
        axis square;                                              % <- square axes
    else
        text(0.5,0.5,'No usable slope estimates','HorizontalAlignment','center','Interpreter','none');
        axis([0 1 0 1]); axis square;
    end
    print(fig2, '-dpdf', fullfile(plotdir, sprintf('%s_slope_hist.pdf', area)), '-painters');
    close(fig2);

    % ----------------------------
    % 3) Per-session bars (optional): Lilliefors rejection rate & slope > 0
    % ----------------------------
    if makeBars && ~isempty(per_session.N)
        K = numel(per_session.N);
        xcats = 1:K;

        fig3 = figure('Visible','off','PaperUnits','inches');
        set(fig3,'PaperSize',[8.5 11],'PaperPosition',[0.75 1.0 7.0 9.0], 'PaperPositionMode','manual');
        set(gca, 'FontSize', 20);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        % Lilliefors rejection rate
        ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
        bar(ax1, xcats, 100*per_session.frac_Lillie_sig, 0.6, 'FaceColor',[0.4 0.4 0.9]);
        ylim(ax1,[0 100]); yticks(ax1,0:20:100);
        xlabel(ax1,'Session','Interpreter','none'); ylabel(ax1,'% neurons with p < alpha','Interpreter','none');
        title(ax1, sprintf('%s: Lilliefors normality (alpha=%.2f)', area, alpha), 'Interpreter','none');
        for k = 1:K
            text(ax1, k, 5, sprintf('N=%d', per_session.N_Lillie_used(k)), ...
                 'HorizontalAlignment','center','Interpreter','none');
        end

        % Slope > 0 rate
        ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
        bar(ax2, xcats, 100*per_session.frac_slope_gt0, 0.6, 'FaceColor',[0.3 0.7 0.3]);
        ylim(ax2,[0 100]); yticks(ax2,0:20:100);
        xlabel(ax2,'Session','Interpreter','none'); ylabel(ax2,'% neurons with slope > 0','Interpreter','none');
        title(ax2, sprintf('%s: Spread-level slope (> 0)', area), 'Interpreter','none');
        for k = 1:K
            text(ax2, k, 5, sprintf('N=%d', per_session.N_Slope_used(k)), ...
                 'HorizontalAlignment','center','Interpreter','none');
        end

        print(fig3, '-dpdf', fullfile(plotdir, sprintf('%s_per_session_rates.pdf', area)), '-painters');
        close(fig3);
    end

    % ----------------------------
    % Save pooled summary
    % ----------------------------
    Summary = struct();
    Summary.area          = area;
    Summary.alpha         = alpha;
    Summary.usedFiles     = used_files;

    Summary.p_Lillie      = all_p_Lillie(:);
    Summary.slope_SP      = all_slope(:);
    Summary.n_residuals   = all_nres(:);
    Summary.n_cells_SP    = all_ncellsSP(:);
    Summary.median_trials_per_cell = all_medtrpc(:);
    Summary.session_id    = sess_of_neuron(:);

    Summary.masks = struct('Lillie',maskL_all(:), 'Slope',maskS_all(:));
    Summary.per_session = per_session;

    outmat = fullfile(savepaths(1), '..', 'models', sprintf('%s_ANOVA_Assumption_Summary.mat', area));
    save(outmat{1}, 'Summary', '-v7.3');
    fprintf('Saved summary: %s\n', outmat{1});
end

% ============================
% Helpers
% ============================
function v = get_field_vec(S, field)
    if isfield(S, field), v = S.(field)(:);
    else, v = [];
    end
end

function u = pad_to(v, n)
    v = v(:);
    if numel(v) >= n
        u = v(1:n);
    else
        u = [v; nan(n-numel(v),1)];
    end
end

function m = safe_mean(x)
    x = x(isfinite(x));
    if isempty(x), m = NaN; else, m = mean(x); end
end
