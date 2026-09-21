function plot_cross_derivative_rms_summary(area, savepaths)
% Aggregates Hessian per-neuron results across sessions and makes:
%   (1) 2x2 heatmap of MEAN smoothed RAW RMS (Hz/rad^2) with fixed caxis [0, 1]
%   (2) 2x2 heatmap of SIGNIFICANT NEURON COUNTS and PERCENTS (BH-FDR q<0.05)
%       with fixed caxis [0, 100]
%   (3) Linear q-value histograms (0..1) with annotations of separable vs inseparable
%   (4) log10(q) histograms with fixed x-axis [-3, 0.5] and 0.1 bins + same annotations
%
% Uses HH, HE, EE only. HE mirrored for off-diagonals.

    % ----------------------------
    % Aggregate
    % ----------------------------
    HH = []; HE = []; EE = [];
    pHH = []; pHE = []; pEE = [];

    for d = 1:length(savepaths)
        sp = savepaths{d};
        hs_file = fullfile(sp, sprintf('models/%s_hessian_results.mat', area));
        if ~exist(hs_file, 'file')
            warning("Missing Hessian results in %s", sp);
            continue
        end

        S = load(hs_file);
        r = S.results(:);

        HH  = [HH;  arrayfun(@(s) s.raw_rms_HH, r)];
        HE  = [HE;  arrayfun(@(s) s.raw_rms_HE, r)];
        EE  = [EE;  arrayfun(@(s) s.raw_rms_EE, r)];

        pHH = [pHH; arrayfun(@(s) s.p_HH, r)];
        pHE = [pHE; arrayfun(@(s) s.p_HE, r)];
        pEE = [pEE; arrayfun(@(s) s.p_EE, r)];
    end

    mask = isfinite(HH) & isfinite(HE) & isfinite(EE) & ...
           isfinite(pHH) & isfinite(pHE) & isfinite(pEE) & ...
           pHH>=0 & pHH<=1 & pHE>=0 & pHE<=1 & pEE>=0 & pEE<=1;

    HH  = HH(mask);  HE  = HE(mask);  EE  = EE(mask);
    pHH = pHH(mask); pHE = pHE(mask); pEE = pEE(mask);

    n = numel(HH);
    if n == 0
        warning('[%s] No valid neurons.', area);
        return
    end

    % ----------------------------
    % BH-FDR
    % ----------------------------
    qHH = bh_fdr(pHH);
    qHE = bh_fdr(pHE);
    qEE = bh_fdr(pEE);

    alpha = 0.05;
    sigHH = qHH < alpha;
    sigHE = qHE < alpha;
    sigEE = qEE < alpha;

    fprintf('[%s] n=%d\n', area, n);
    fprintf('  HH: %d (%.2f%%)\n', sum(sigHH), 100*mean(sigHH));
    fprintf('  HE: %d (%.2f%%)\n', sum(sigHE), 100*mean(sigHE));
    fprintf('  EE: %d (%.2f%%)\n', sum(sigEE), 100*mean(sigEE));

    plotdir = fullfile(savepaths{1}, '..', 'plots', 'summary', ...
        'separability_hessian');
    if ~exist(plotdir, 'dir'), mkdir(plotdir); end

    % ============================================================
    % 1) Mean RAW RMS heatmap (fixed colorbar [0,1])
    % ============================================================
    muHH = mean(HH); seHH = std(HH)/sqrt(n);
    muHE = mean(HE); seHE = std(HE)/sqrt(n);
    muEE = mean(EE); seEE = std(EE)/sqrt(n);

    Hraw = [muHH, muHE; muHE, muEE];

    fig1 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 4.2 4]);
    ax1 = axes('Parent',fig1); hold(ax1,'on'); set(ax1,'FontSize',14);

    imagesc(ax1, [1 2], [1 2], Hraw);
    axis(ax1,'image'); axis(ax1,'ij');
    colormap(ax1, flipud(gray)); cb1 = colorbar(ax1);
    ylabel(cb1, 'Mean RAW RMS (Hz/rad^2)');
    caxis(ax1, [0 1]); % FIXED

    set(ax1,'XTick',1:2,'XTickLabel',{'H','E'});
    set(ax1,'YTick',1:2,'YTickLabel',{'H','E'});

    text(1,1,sprintf('%.3g\n±%.2g',muHH,seHH),'HorizontalAlignment','center','FontWeight','bold');
    text(2,1,sprintf('%.3g\n±%.2g',muHE,seHE),'HorizontalAlignment','center','FontWeight','bold');
    text(1,2,sprintf('%.3g\n±%.2g',muHE,seHE),'HorizontalAlignment','center','FontWeight','bold');
    text(2,2,sprintf('%.3g\n±%.2g',muEE,seEE),'HorizontalAlignment','center','FontWeight','bold');

    title(ax1, sprintf('%s: Mean smoothed RAW RMS (n=%d)', area, n));
    grid(ax1,'on'); box(ax1,'on');
    print(fig1, fullfile(plotdir, sprintf('%s_Hessian_rawRMS_2x2.pdf', area)), '-dpdf','-painters');
    close(fig1);

    % ============================================================
    % 2) Significant fraction heatmap (fixed colorbar [0,100])
    % ============================================================
    Hsig = [100*mean(sigHH), 100*mean(sigHE);
            100*mean(sigHE), 100*mean(sigEE)];

    fig2 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 4.2 4]);
    ax2 = axes('Parent',fig2); hold(ax2,'on'); set(ax2,'FontSize',14);

    imagesc(ax2, [1 2], [1 2], Hsig);
    axis(ax2,'image'); axis(ax2,'ij');
    colormap(ax2, flipud(gray)); cb2 = colorbar(ax2);
    ylabel(cb2, '% neurons significant (q<0.05)');
    caxis(ax2, [0 100]); % FIXED

    set(ax2,'XTick',1:2,'XTickLabel',{'H','E'});
    set(ax2,'YTick',1:2,'YTickLabel',{'H','E'});

    text(1,1,sprintf('%d\n(%.1f%%)',sum(sigHH),100*mean(sigHH)),'HorizontalAlignment','center','FontWeight','bold');
    text(2,1,sprintf('%d\n(%.1f%%)',sum(sigHE),100*mean(sigHE)),'HorizontalAlignment','center','FontWeight','bold');
    text(1,2,sprintf('%d\n(%.1f%%)',sum(sigHE),100*mean(sigHE)),'HorizontalAlignment','center','FontWeight','bold');
    text(2,2,sprintf('%d\n(%.1f%%)',sum(sigEE),100*mean(sigEE)),'HorizontalAlignment','center','FontWeight','bold');

    title(ax2, sprintf('%s: Significant Hessian entries (BH-FDR, n=%d)', area, n));
    grid(ax2,'on'); box(ax2,'on');
    print(fig2, fullfile(plotdir, sprintf('%s_Hessian_sigFrac_2x2.pdf', area)), '-dpdf','-painters');
    close(fig2);

    % ============================================================
    % 3) Linear q-value histograms (0..1) with annotations
    % ============================================================
    plot_q_hist_annot(qHH, alpha, sprintf('%s: HH q-values',area), ...
        fullfile(plotdir, sprintf('%s_Hessian_HH_q_hist.pdf', area)));
    plot_q_hist_annot(qHE, alpha, sprintf('%s: HE q-values',area), ...
        fullfile(plotdir, sprintf('%s_Hessian_HE_q_hist.pdf', area)));
    plot_q_hist_annot(qEE, alpha, sprintf('%s: EE q-values',area), ...
        fullfile(plotdir, sprintf('%s_Hessian_EE_q_hist.pdf', area)));

    % ============================================================
    % 4) log10(q) histograms with FIXED x-axis [-3, 0.5] and 0.1 bins
    %    + same annotations
    % ============================================================
    plot_logq_hist_fixed_annot(qHH, alpha, sprintf('%s: HH log10(q)',area), ...
        fullfile(plotdir, sprintf('%s_Hessian_HH_logq_hist.pdf', area)));
    plot_logq_hist_fixed_annot(qHE, alpha, sprintf('%s: HE log10(q)',area), ...
        fullfile(plotdir, sprintf('%s_Hessian_HE_logq_hist.pdf', area)));
    plot_logq_hist_fixed_annot(qEE, alpha, sprintf('%s: EE log10(q)',area), ...
        fullfile(plotdir, sprintf('%s_Hessian_EE_logq_hist.pdf', area)));
end

% ============================================================
% Linear q histogram with separable/inseparable annotations
% ============================================================
function plot_q_hist_annot(qvals, alpha, ttl, outpath)
    fig = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 5 4]);
    ax = axes('Parent',fig); hold(ax,'on'); set(ax,'FontSize',14);

    qvals = qvals(isfinite(qvals));
    if isempty(qvals)
        title(ax, [ttl ' (empty)']);
        print(fig,outpath,'-dpdf','-painters'); close(fig); return;
    end

    edges = linspace(0,1,41);
    histogram(ax, qvals, edges, 'Normalization','count', ...
        'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.75);

    xline(ax, alpha, 'r--', 'LineWidth',1.8);

    % Annotation counts
    n = numel(qvals);
    ninsep = sum(qvals < alpha);
    nsep   = n - ninsep;
    pinsep = 100 * ninsep / n;
    psep   = 100 * nsep   / n;

    yl = ylim(ax);

    % Place labels near top, split left/right of threshold
    x_left  = max(alpha*0.4, 0.02);
    x_right = min(alpha + (1-alpha)*0.55, 0.98);
    y_text  = yl(2) * 0.92;

    text(ax, x_left,  y_text, sprintf('Inseparable: %d (%.1f%%)', ninsep, pinsep), ...
        'HorizontalAlignment','left','FontWeight','bold');
    text(ax, x_right, y_text, sprintf('Separable: %d (%.1f%%)', nsep, psep), ...
        'HorizontalAlignment','left','FontWeight','bold');

    xlabel(ax,'q-value'); ylabel(ax,'Neuron count');
    title(ax,ttl); grid(ax,'on'); box(ax,'on');

    print(fig,outpath,'-dpdf','-painters');
    close(fig);
end

% ============================================================
% log10(q) histogram with fixed axis and annotations
% x-axis: [-3, 0.5], bins: 0.1
% ============================================================
function plot_logq_hist_fixed_annot(qvals, alpha, ttl, outpath)
    fig = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 5 4]);
    ax = axes('Parent',fig); hold(ax,'on'); set(ax,'FontSize',14);

    qvals = qvals(isfinite(qvals));
    if isempty(qvals)
        title(ax, [ttl ' (empty)']);
        print(fig,outpath,'-dpdf','-painters'); close(fig); return;
    end

    qclip = max(qvals, eps); % avoid log10(0)
    logq = log10(qclip);

    % Fixed edges
    edges = -3:0.1:0.5;
    histogram(ax, logq, edges, 'Normalization','count', ...
        'FaceColor',[0.4 0.4 0.4],'EdgeColor','k','FaceAlpha',0.75);

    % Fixed x-limits
    xlim(ax, [-3 0.5]);

    % Threshold line at log10(alpha)
    xthr = log10(alpha);
    xline(ax, xthr, 'r--', 'LineWidth',1.8);

    % Annotation counts based on q (not logq)
    n = numel(qvals);
    ninsep = sum(qvals < alpha);
    nsep   = n - ninsep;
    pinsep = 100 * ninsep / n;
    psep   = 100 * nsep   / n;

    yl = ylim(ax);
    y_text = yl(2) * 0.92;

    % Place left label to left of threshold, right label to right
    x_left  = max(-2.9, xthr - 1.2);
    x_right = min(0.45, xthr + 0.3);

    text(ax, x_left,  y_text, sprintf('Inseparable: %d (%.1f%%)', ninsep, pinsep), ...
        'HorizontalAlignment','left','FontWeight','bold');
    text(ax, x_right, y_text, sprintf('Separable: %d (%.1f%%)', nsep, psep), ...
        'HorizontalAlignment','left','FontWeight','bold');

    xlabel(ax,'log_{10}(q)'); ylabel(ax,'Neuron count');
    title(ax,ttl); grid(ax,'on'); box(ax,'on');

    print(fig,outpath,'-dpdf','-painters');
    close(fig);
end

% ============================================================
% BH-FDR q-values (Benjamini–Hochberg), returns q aligned to input order
% ============================================================
function q = bh_fdr(p)
    p = p(:);
    m = numel(p);
    [ps, idx] = sort(p);
    ranks = (1:m)';

    qs = ps .* (m ./ ranks);
    for i = m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end
    qs = min(qs, 1);

    q = zeros(size(p));
    q(idx) = qs;
end
