%% Weighted baseline calculation and visual adjustments
function analyze_additivity_surfaces(area, savepath)
% Analyze additive separability of multi-effector tuning surfaces.
% For each neuron n:
%   - ORIGINAL ME surface F (mean FR per (hand, eye) bin)
%   - ME ADDITIVE surface F_add (row/col marginals; LS under missingness)
%   - SE ADDITIVE surface F_SEadd (from SE hand/eye tunings; baseline is a
%       WEIGHTED average of SE_H and SE_E means, with weights equal to the
%       respective peak-to-trough amplitudes)
%   - Per-cell SE of the mean from trial-level variance on ME
%   - Additivity tests:
%       • ME: wild-bootstrap p-value with refitting (p_boot)
%       • SE: fixed-null χ² p-value without refitting (p_boot_SE)
%   - Correlations:
%       • r(ME,ME) = corr(F, F_add) over observed cells
%       • r(ME,SE) = corr(F, F_SEadd) over observed cells
%
% Saves (per neuron):
%   plots/additivity/neuron_###_additivity3D.pdf   [ME raw, ME additive, SE additive]
%   plots/additivity/neuron_###_H_tuning.pdf       [3-panel: both, ME-only, SE-only]
%   plots/additivity/neuron_###_E_tuning.pdf       [3-panel: both, ME-only, SE-only]
%
% Notes on plotting:
%   • Tuning-curve figures are now 3 subplots (both / ME-only / SE-only).
%   • Tuning-curve figure size has height = 3 × width, so each subplot is
%     approximately square. Text size on tuning curves is doubled.
%   • The FR range is synchronized ACROSS the surface plots (z-axis) and BOTH
%     tuning-curve figures (y-axis) by using the common [min,max] over all
%     relevant ME/SE surfaces and curves, per neuron.
%
% Saves (MAT):
%   models/%s_Additivity_Surface_Test_NoClass.mat with fields used by summary and extras:
%     p_chi2, p_boot, p_boot_SE, R2_add, eta_res,
%     corr_F_Fadd, p_corr_F_Fadd, r_ME_SE,
%     c_SE_H, c_SE_E, c_ME_HE,
%     c_SE_weighted, ...                    % <--- NEW: weighted SE baseline
%     peak_ME_HE, peak_ME_add, peak_SE_add  % peaks (max FR) metrics

    % ----------------------------
    % Config
    % ----------------------------
    nbins        = 6;
    edge_shift   = 0.1;
    B_boot       = 1000;     % wild bootstrap iterations for ME additivity
    rng_seed     = 13;       % reproducible across runs
    eps_se       = 1e-6;     % SE floor for stability
    center_target = 7;

    % ----------------------------
    % Load processed data (raw FR)
    % ----------------------------
    S = load(processed_data_file(savepath, area));
    D = S.data;
    FR   = D.FR;    % [trials x neurons]
    TP1  = D.TP1;   % hand targets (x,y)
    TP2  = D.TP2;   % eye  targets (x,y)
    TPi1 = D.TPi1;  % hand effector code
    TPi2 = D.TPi2;  % eye effector code

    % Trial masks
    multiTP = find(TPi2 ~= center_target & TPi1 ~= center_target); % ME trials
    handTP  = find(TPi2 == center_target & TPi1 ~= center_target); % SE hand-only
    eyeTP   = find(TPi1 == center_target & TPi2 ~= center_target); % SE eye-only

    % Binning (angles)
    edges_theta = linspace(-pi, pi, nbins+1) + edge_shift;
    hcent = edges_theta(1:end-1) + diff(edges_theta)/2;
    ecent = hcent;

    % Angles & bin indices (fixed across neurons)
    thetaH_multi = atan2(TP1(multiTP,2), TP1(multiTP,1));
    thetaE_multi = atan2(TP2(multiTP,2), TP2(multiTP,1));
    Hb_multi = discretize(thetaH_multi, edges_theta);
    Eb_multi = discretize(thetaE_multi, edges_theta);

    thetaH_hand = atan2(TP1(handTP,2), TP1(handTP,1));
    thetaE_eye  = atan2(TP2(eyeTP,2),  TP2(eyeTP,1));
    Hb_hand = discretize(thetaH_hand, edges_theta);
    Eb_eye  = discretize(thetaE_eye,  edges_theta);

    % Precompute ME per-cell trial indices (speed)
    cell_indices = cell(nbins, nbins);
    for i = 1:nbins
        for j = 1:nbins
            cell_indices{i,j} = multiTP(Hb_multi == i & Eb_multi == j);
        end
    end

    % Precompute SE trial-index lists per direction bin
    hand_bin_trials = cell(nbins,1);
    eye_bin_trials  = cell(nbins,1);
    for i = 1:nbins
        hand_bin_trials{i} = handTP(Hb_hand == i);
        eye_bin_trials{i}  = eyeTP (Eb_eye  == i);
    end

    % Meshgrid for surfaces (X = Eye, Y = Hand)
    [Egrid, Hgrid] = meshgrid(ecent, hcent);

    % ----------------------------
    % Output directory
    % ----------------------------
    plot_root = fullfile(savepath, 'plots', 'additivity');
    if ~exist(plot_root, 'dir'), mkdir(plot_root); end

    % ----------------------------
    % Preallocate results needed by summary + extras
    % ----------------------------
    N = size(FR,2);
    p_chi2        = nan(N,1);  % ME (approx χ² with df correction)
    p_boot        = nan(N,1);  % ME (wild bootstrap, refit)
    p_boot_SE     = nan(N,1);  % SE (fixed-null χ²; see helper)
    R2_add        = nan(N,1);
    eta_res       = nan(N,1);
    corr_F_Fadd   = nan(N,1);  % r(ME,ME)
    p_corr_F_Fadd = nan(N,1);  % p-value for r(ME,ME)
    r_ME_SE       = nan(N,1);  % r(ME,SE)
    c_SE_H        = nan(N,1);  % mean of SE hand tuning
    c_SE_E        = nan(N,1);  % mean of SE eye  tuning
    c_ME_HE       = nan(N,1);  % grand mean of ME surface
    c_SE_weighted = nan(N,1);  % NEW: weighted SE baseline

    % NEW: peak FR (max over finite cells) for each surface
    peak_ME_HE   = nan(N,1);   % ME raw surface peak
    peak_ME_add  = nan(N,1);   % ME additive surface peak
    peak_SE_add  = nan(N,1);   % SE additive surface peak

    % ----------------------------
    % Parallel loop over neurons
    % ----------------------------
    fprintf('Analyzing additivity for %d neurons (B=%d bootstrap ME)...\n', N, B_boot);
    parfor n = 1:N
        try
            % ---------- ORIGINAL ME surface F & per-cell SE ----------
            F      = nan(nbins, nbins);
            SE_mat = nan(nbins, nbins);   % SE of per-cell mean (from ME trials)
            counts = zeros(nbins, nbins);

            for i = 1:nbins
                for j = 1:nbins
                    idx = cell_indices{i,j};
                    if ~isempty(idx)
                        fr_ij = FR(idx, n);
                        fr_ij = fr_ij(~isnan(fr_ij));
                        if ~isempty(fr_ij)
                            counts(i,j) = numel(fr_ij);
                            F(i,j)      = mean(fr_ij);
                            if numel(fr_ij) >= 2
                                v_ij = var(fr_ij, 0);
                            else
                                v_ij = NaN;
                            end
                            SE_mat(i,j) = sqrt(v_ij / max(counts(i,j),1));
                        end
                    end
                end
            end
            M = ~isnan(F);  % observed cells

            % Stabilize SE where needed
            if any(M(:))
                se_med = median(SE_mat(M), 'omitnan');
                if ~isfinite(se_med) || se_med <= 0, se_med = 1; end
                SE_mat(~isfinite(SE_mat)) = se_med;
                SE_mat(SE_mat < eps_se)   = eps_se;
            end

            % ---------- ME additive fit (row/col LS under missingness) ----------
            F_add = additive_fit_from_marginals(F);

            % ---------- Residuals & stats for ME ----------
            R = F - F_add; R(~M) = NaN;
            rows_obs = find(any(M, 2));  r = numel(rows_obs);
            cols_obs = find(any(M, 1));  c = numel(cols_obs);
            m_obs = nnz(M);
            df_local = max(m_obs - (r + c - 1), 1);  % df for refit model

            Z = R ./ SE_mat;

            % ME: Chi-square approx (not displayed)
            if df_local > 0 && all(isfinite(SE_mat(M)))
                p_chi2_local = 1 - chi2cdf(nansum((R(M)./SE_mat(M)).^2), df_local);
            else
                p_chi2_local = NaN;
            end

            % ME: Wild-bootstrap p-value with refitting (displayed as p_add)
            rng(rng_seed + n, 'twister');
            p_boot_local = wild_bootstrap_pval_ME(F, F_add, R, M, SE_mat, B_boot);

            % Effect sizes on ME
            mu_ME = mean(F(M), 'omitnan');
            SSE = nansum((F(:) - F_add(:)).^2);
            SST = nansum((F(:) - mu_ME).^2);
            R2_add_local  = max(0, 1 - SSE / max(SST, eps));
            eta_res_local = SSE / max(SST, eps);

            % Correlation r(ME,ME)
            fvec = F(M); avec = F_add(M);
            if numel(fvec) >= 3 && std(fvec,'omitnan')>0 && std(avec,'omitnan')>0
                [r_ME_ME, p_ME_ME] = corr(fvec, avec, 'Type','Pearson','Rows','pairwise');
            else
                r_ME_ME = NaN; p_ME_ME = NaN;
            end

            % ---------- SE additive surface (from SE curves; WEIGHTED baseline) ----------
            % SE curves (trial means)
            SE_hand_mu = nan(nbins,1);
            SE_eye_mu  = nan(1,nbins);
            for i = 1:nbins
                tidx = hand_bin_trials{i};
                if ~isempty(tidx)
                    fr = FR(tidx, n);
                    SE_hand_mu(i) = mean(fr, 'omitnan');
                end
            end
            for j = 1:nbins
                tidx = eye_bin_trials{j};
                if ~isempty(tidx)
                    fr = FR(tidx, n);
                    SE_eye_mu(j) = mean(fr, 'omitnan');
                end
            end

            % Build SE additive grid with WEIGHTED SE-only baseline
            A = SE_hand_mu(:);
            B = SE_eye_mu(:).';
            A0 = mean(A, 'omitnan');
            B0 = mean(B, 'omitnan');

            ampH = ptp_amplitude(A);   % peak-to-trough hand
            ampE = ptp_amplitude(B);   % peak-to-trough eye

            wH = max(ampH, 0);
            wE = max(ampE, 0);

            if isfinite(wH) && wH > 0 && isfinite(wE) && wE > 0
                c_SE = (wH*A0 + wE*B0) / (wH + wE);
            elseif isfinite(wH) && wH > 0
                c_SE = A0;
            elseif isfinite(wE) && wE > 0
                c_SE = B0;
            else
                % both amplitudes zero or invalid: fall back to unweighted mean
                c_SE = mean([A0, B0], 'omitnan');
            end

            F_SEadd = (A - A0) + (B - B0) + c_SE;       % [nbins x nbins]

            % Correlation r(ME,SE) over observed ME cells
            r_ME_SE_local = NaN;
            if any(M(:))
                gvec = F_SEadd(M);
                if numel(fvec) >= 3 && std(gvec,'omitnan')>0 && std(fvec,'omitnan')>0
                    r_ME_SE_local = corr(fvec, gvec, 'Type','Pearson','Rows','pairwise');
                end
            end

            % SE: fixed-null χ² p-value (no refit; df = # observed ME cells)
            p_boot_SE_local = pval_fixed_chi2(F, F_SEadd, M, SE_mat);

            % ---------- Tuning curves (means & SEs) ----------
            % ME projected curves (row/col marginals of F) + SE raw curves
            ME_hand_mu = mean(F, 2, 'omitnan');   % [nbins x 1]
            ME_eye_mu  = mean(F, 1, 'omitnan');   % [1 x nbins]

            % SE bars for ME curves from per-cell SE (combine across row/col)
            ME_hand_se = nan(nbins,1);
            ME_eye_se  = nan(1,nbins);
            for i = 1:nbins
                cols_i = find(M(i,:)); k = numel(cols_i);
                if k>0, ME_hand_se(i) = sqrt(nansum(SE_mat(i,cols_i).^2)) / k; end
            end
            for j = 1:nbins
                rows_j = find(M(:,j)); k = numel(rows_j);
                if k>0, ME_eye_se(j) = sqrt(nansum(SE_mat(rows_j,j).^2)) / k; end
            end

            % SE raw curve SE bars (trial-level)
            SE_hand_se = nan(nbins,1);
            SE_eye_se  = nan(1,nbins);
            for i = 1:nbins
                tidx = hand_bin_trials{i};
                if ~isempty(tidx)
                    fr = FR(tidx, n);
                    SE_hand_se(i) = std(fr, 'omitnan')/sqrt(sum(isfinite(fr)));
                end
            end
            for j = 1:nbins
                tidx = eye_bin_trials{j};
                if ~isempty(tidx)
                    fr = FR(tidx, n);
                    SE_eye_se(j) = std(fr, 'omitnan')/sqrt(sum(isfinite(fr)));
                end
            end

            % ---------- COMMON FR RANGE across SURFACES and TUNING CURVES ----------
            vals_surface = [F(M); F_add(M)];
            vals_surface = [vals_surface; F_SEadd(isfinite(F_SEadd))];
            vals_curves  = [ME_hand_mu(:); ME_eye_mu(:); SE_hand_mu(:); SE_eye_mu(:)];
            fr_min = min([vals_surface; vals_curves], [], 'omitnan');
            fr_max = max([vals_surface; vals_curves], [], 'omitnan');
            if ~isfinite(fr_min), fr_min = 0; end
            if ~isfinite(fr_max), fr_max = 1; end
            if ~(fr_max > fr_min)
                pad = max(1, abs(fr_min)*0.1);
                fr_ylim = [fr_min - pad, fr_max + pad];
            else
                fr_ylim = [fr_min, fr_max];
            end

            % ---------- 3D surfaces figure (ME raw, ME additive, SE additive) ----------
            fig = figure('Visible','Off', ...
                'Units','inches','Position',[1 1 8 11], ...
                'PaperUnits','inches','PaperPosition',[0 0 8 11]);
            cm = parula;

            % (1) ME raw surface
            subplot(3,1,1);
            surf(Egrid, Hgrid, F, 'EdgeColor','none'); colormap(cm); shading interp;
            view(45,30); axis tight; pbaspect([1 1 1]);
            zlim(fr_ylim); caxis(fr_ylim);
            set(gca,'FontSize',20); xlabel('\theta_E'); ylabel('\theta_H'); zlabel('FR [Hz]');
            title(sprintf('Original ME surface F (Neuron %d)', n));

            % (2) ME additive surface
            subplot(3,1,2);
            surf(Egrid, Hgrid, F_add, 'EdgeColor','none'); colormap(cm); shading interp;
            view(45,30); axis tight; pbaspect([1 1 1]);
            zlim(fr_ylim); caxis(fr_ylim);
            set(gca,'FontSize',20); xlabel('\theta_E'); ylabel('\theta_H'); zlabel('FR [Hz]');
            title(sprintf('ME additive F_{add} | r(ME,ME)=%.2f, p_{add}=%.3g', r_ME_ME, p_boot_local));

            % (3) SE additive surface
            subplot(3,1,3);
            surf(Egrid, Hgrid, F_SEadd, 'EdgeColor','none'); colormap(cm); shading interp;
            view(45,30); axis tight; pbaspect([1 1 1]);
            zlim(fr_ylim); caxis(fr_ylim);
            set(gca,'FontSize',20); xlabel('\theta_E'); ylabel('\theta_H'); zlabel('FR [Hz]');
            title(sprintf('SE additive F_{SE} | r(ME,SE)=%.2f, p_{add,SE}=%.3g', r_ME_SE_local, p_boot_SE_local));

            sgtitle(sprintf('%s – Additivity Surfaces (Neuron %d)', area, n), 'FontSize', 16);

            outpdf = fullfile(plot_root, sprintf('neuron_%03d_additivity3D.pdf', n));
            print(fig, '-dpdf', '-painters', outpdf);
            close(fig);

            % ---------- Tuning-curve figures (3 subplots; y-lims synced to FR range) ----------
            % Hand tuning curves
            outpdfH = fullfile(plot_root, sprintf('neuron_%03d_H_tuning.pdf', n));
            plot_tuning_curves(hcent, ...
                               SE_hand_mu, SE_hand_se, ...
                               ME_hand_mu, ME_hand_se, ...
                               '\theta_H (rad)', ...
                               outpdfH, fr_ylim);

            % Eye tuning curves
            outpdfE = fullfile(plot_root, sprintf('neuron_%03d_E_tuning.pdf', n));
            plot_tuning_curves(ecent, ...
                               SE_eye_mu(:), SE_eye_se(:), ...
                               ME_eye_mu(:), ME_eye_se(:), ...
                               '\theta_E (rad)', ...
                               outpdfE, fr_ylim);

            % ---------- NEW: Peak FRs ----------
            % ME raw/add peaks over observed cells; SE-add peak over all finite cells
            peak_ME_HE_local  = NaN;
            peak_ME_add_local = NaN;
            peak_SE_add_local = NaN;

            if any(M(:))
                peak_ME_HE_local  = max(F(M));
                peak_ME_add_local = max(F_add(M));
            end
            ms = isfinite(F_SEadd);
            if any(ms(:))
                peak_SE_add_local = max(F_SEadd(ms));
            end

            % ---------- Store outputs ----------
            p_chi2(n)        = p_chi2_local;
            p_boot(n)        = p_boot_local;
            p_boot_SE(n)     = p_boot_SE_local;
            R2_add(n)        = R2_add_local;
            eta_res(n)       = eta_res_local;
            corr_F_Fadd(n)   = r_ME_ME;
            p_corr_F_Fadd(n) = p_ME_ME;
            r_ME_SE(n)       = r_ME_SE_local;
            c_SE_H(n)        = A0;
            c_SE_E(n)        = B0;
            c_ME_HE(n)       = mu_ME;
            c_SE_weighted(n) = c_SE;      % NEW: store weighted SE baseline

            % NEW peaks
            peak_ME_HE(n)   = peak_ME_HE_local;
            peak_ME_add(n)  = peak_ME_add_local;
            peak_SE_add(n)  = peak_SE_add_local;

        catch MEe
            warning('Additivity analysis failed for neuron %d: %s', n, MEe.message);
        end
    end

    % ----------------------------
    % Save results
    % ----------------------------
    Results = struct();
    Results.area         = area;
    Results.B_boot       = B_boot;
    Results.nbins        = nbins;
    Results.edge_shift   = edge_shift;

    Results.p_chi2       = p_chi2;        % ME (df-corrected)
    Results.p_boot       = p_boot;        % ME (wild bootstrap, refit)
    Results.p_boot_SE    = p_boot_SE;     % SE (fixed-null χ²)
    Results.R2_add       = R2_add;
    Results.eta_res      = eta_res;
    Results.corr_F_Fadd   = corr_F_Fadd;  % r(ME,ME)
    Results.p_corr_F_Fadd = p_corr_F_Fadd;
    Results.r_ME_SE       = r_ME_SE;

    % Constants to compare later
    Results.c_SE_H        = c_SE_H;
    Results.c_SE_E        = c_SE_E;
    Results.c_ME_HE       = c_ME_HE;
    Results.c_SE_weighted = c_SE_weighted; % NEW

    % NEW: peaks to compare later
    Results.peak_ME_HE  = peak_ME_HE;
    Results.peak_ME_add = peak_ME_add;
    Results.peak_SE_add = peak_SE_add;

    outmat = fullfile(savepath, 'models', sprintf('%s_Additivity_Surface_Test_NoClass.mat', area));
    save(outmat, 'Results', '-v7.3');
    fprintf('Saved results: %s\n', outmat);
end

% =========================
% Helpers
% =========================
function F_add = additive_fit_from_marginals(F)
% Least-squares additive fit under missingness via row/col means (centered).
% Identifiability: sum(u)=sum(v)=0; then add grand mean.
    mu   = mean(F(:), 'omitnan');
    rbar = mean(F, 2, 'omitnan');
    cbar = mean(F, 1, 'omitnan');
    a = rbar - mean(rbar, 'omitnan');
    b = cbar - mean(cbar, 'omitnan');
    F_add = a * ones(1, size(F,2)) + ones(size(F,1), 1) * b + mu;
    F_add(isnan(F)) = NaN;
end

function p_boot = wild_bootstrap_pval_ME(F, F_add, R, M, SE_mat, B)
% ME additivity: wild sign-flip bootstrap under additive null WITH refitting.
%   F* = F_add + S .* R, S_{ij} ~ Rademacher on observed cells;
%   Refit additive on F*; compute T* = sum((R*/SE)^2); p = (1+sum(T*>=Tobs))/(1+B).
    if ~any(M(:)), p_boot = NaN; return; end
    T_obs = nansum( (R(M)./SE_mat(M)).^2 );

    Tstars = zeros(B,1);
    [m, n] = size(F);
    for b = 1:B
        Sgn = nan(m,n);
        Sgn(M) = 2*(rand(nnz(M),1) > 0.5) - 1;
        Fstar = F_add;
        Fstar(M) = F_add(M) + Sgn(M) .* R(M);

        Fadd_star = additive_fit_from_marginals(Fstar);
        Rstar     = Fstar - Fadd_star;
        Rstar(~M) = NaN;

        Tstars(b) = nansum( (Rstar(M)./SE_mat(M)).^2 );
    end
    p_boot = (1 + sum(Tstars >= T_obs)) / (1 + B);
end

function p = pval_fixed_chi2(F, F0, M, SE_mat)
% Fixed-null χ² p-value (used for SE additive):
%   T = sum_{obs} ((F - F0)/SE)^2  ~  χ²_df with df = #observed cells.
% No refit (null is fixed by SE curves), so df is simply nnz(M).
    if ~any(M(:)), p = NaN; return; end
    T = nansum( ( (F(M) - F0(M)) ./ SE_mat(M) ).^2 );
    df = nnz(M);
    p = 1 - chi2cdf(T, df);
end

function plot_tuning_curves(binc, mu_SE, se_SE, mu_ME, se_ME, xlab, spath, fr_ylim)
% Three-panel raw tuning curves (means ± SE):
%   (1) Combined:  SE (black) + ME (gray)
%   (2) ME-only   : ME (gray)
%   (3) SE-only   : SE (black)
% Figure height = 3 × width; text size doubled; y-lims set to fr_ylim.
    H = 11; W = H/3; % inches: height = 3 × width (each subplot ~ square)
    fs = 20;        % doubled from 20
    fst = 20;       % title

    % Fix both on-screen and paper geometry. Subplot positions are computed
    % from Figure.Position, so setting only PaperPosition makes batch-mode
    % output depend on MATLAB's default (landscape) figure size and clips the
    % narrow portrait panels. These dimensions match the intended legacy PDF.
    h = figure('Visible','Off', ...
        'Units','inches','Position',[1 1 W H], ...
        'PaperUnits','inches','PaperPosition',[0 0 W H]);

    % ---------- (1) Combined ----------
    subplot(3,1,1);
    hold on; grid on;
    plot(binc, mu_SE, 'k-', 'LineWidth', 1.6, 'HandleVisibility','off');
    plot(binc, mu_ME, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.6, 'HandleVisibility','off');

    errorbar(binc, mu_SE, se_SE, 'ko', 'MarkerFaceColor','k', ...
             'LineWidth',1.8, 'DisplayName', 'Raw Data (Single-EFfector)');
    errorbar(binc, mu_ME, se_ME, 'o', 'Color',[0.5 0.5 0.5], ...
             'MarkerFaceColor',[0.5 0.5 0.5], 'LineWidth',1.8, ...
             'DisplayName', 'Raw Data (Multi-EFfector)');

    % Amplitude markers for ME curve
    if any(isfinite(mu_ME))
        [peak_val, peak_idx] = max(mu_ME);
        trough_val = min(mu_ME);
        if isfinite(peak_val) && isfinite(trough_val)
            peak_x = binc(peak_idx);
            amp_color = [0.3 0.3 0.3 0.3];
            xl = [-pi, pi];
            plot(xl, [peak_val peak_val], '-', 'Color', amp_color, 'LineWidth', 2, 'HandleVisibility','off');
            plot(xl, [trough_val trough_val], '-', 'Color', amp_color, 'LineWidth', 2, 'HandleVisibility','off');
            plot([peak_x peak_x], [trough_val peak_val], '-', 'Color', amp_color, 'LineWidth', 2, 'HandleVisibility','off');
        end
    end
    xlim([-pi pi]); ylim(fr_ylim);
    set(gca,'FontSize',fs);
    xlabel(xlab); ylabel('FR [Hz]');
    lg = legend('Location','south'); set(lg,'FontSize',fs*0.5);

    % ---------- (2) ME-only ----------
    subplot(3,1,2);
    hold on; grid on;
    plot(binc, mu_ME, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.8, 'DisplayName','Raw Data (Multi-EFfector)');
    errorbar(binc, mu_ME, se_ME, 'o', 'Color',[0.5 0.5 0.5], ...
             'MarkerFaceColor',[0.5 0.5 0.5], 'LineWidth',1.8, ...
             'HandleVisibility','off');

    % Amplitude markers for ME-only
    if any(isfinite(mu_ME))
        [peak_val, peak_idx] = max(mu_ME);
        trough_val = min(mu_ME);
        if isfinite(peak_val) && isfinite(trough_val)
            peak_x = binc(peak_idx);
            amp_color = [0.3 0.3 0.3 0.3];
            xl = [-pi, pi];
            plot(xl, [peak_val peak_val], '-', 'Color', amp_color, 'LineWidth', 2, 'HandleVisibility','off');
            plot(xl, [trough_val trough_val], '-', 'Color', amp_color, 'LineWidth', 2, 'HandleVisibility','off');
            plot([peak_x peak_x], [trough_val peak_val], '-', 'Color', amp_color, 'LineWidth', 2, 'HandleVisibility','off');
        end
    end
    xlim([-pi pi]); ylim(fr_ylim);
    set(gca,'FontSize',fs);
    xlabel(xlab); ylabel('FR [Hz]');

    % ---------- (3) SE-only ----------
    subplot(3,1,3);
    hold on; grid on;
    plot(binc, mu_SE, 'k-', 'LineWidth', 1.8, 'DisplayName','Raw Data (Single-EFfector)');
    errorbar(binc, mu_SE, se_SE, 'ko', 'MarkerFaceColor','k', ...
             'LineWidth',1.8, 'HandleVisibility','off');

    xlim([-pi pi]); ylim(fr_ylim);
    set(gca,'FontSize',fs);
    xlabel(xlab); ylabel('FR [Hz]');

    print(h,'-dpdf','-painters','-r300', spath);
    close(h);
end

function amp = ptp_amplitude(x)
% Peak-to-trough amplitude (ignoring NaNs). Returns NaN if all-NaN.
    if all(~isfinite(x))
        amp = NaN;
    else
        amp = max(x,[],'omitnan') - min(x,[],'omitnan');
    end
end
