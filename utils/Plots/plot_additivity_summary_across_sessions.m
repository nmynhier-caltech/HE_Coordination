function plot_additivity_summary_across_sessions(area, savepaths, varargin)
% Summary plots for analyze_additivity_surfaces() across many sessions.
%
% Aggregates Results.* from:
%   models/%s_Additivity_Surface_Test_NoClass.mat  (preferred if exists)
%   models/%s_Additivity_Surface_Test.mat          (fallback)
%
% Plots:
%   1) Histogram of corr_F_Fadd (if available)
%   2) Histogram of p-values with corrected-threshold line (FDR-BH or Bonf)
%   3) Histogram of R2_add (optional, on by default)
%   4) Histogram of Pearson r comparing SE_add (SE-only) vs ME-smoothed (ME-only CV) by class
%   5) Baseline scatter (single panel): c_SE^weighted  vs  c_ME,HE  [RAW Hz]
%   6) Peak FR scatters (two subplots): (a) peak(ME_add) vs peak(ME_raw), (b) peak(SE_add) vs peak(ME_raw)
%
% Saves:
%   plots/summary/additive_reconstruction/<area>_*.pdf
%   models/<area>_Additivity_Summary_Accum.mat  (pooled vectors + MCP outputs)
%
% Args (name/value):
%   'PField'        : 'p_boot' (default) or 'p_chi2'
%   'Method'        : 'FDR-BH' (default) | 'bonferroni'
%   'Alpha'         : scalar in (0,1) (default 0.05)
%   'ShowR2'        : true (default) | false
%   'NBins'         : histogram bins for p-values & R2 (default 30)
%   'NBinsCorr'     : histogram bins for correlations (default 15)
%   'SECVFolds'     : K-folds for SE-only CV of kappa (default 10)
%   'PerCellTrain'  : ME per-cell training target (default 6)
%   'PerCellTotal'  : ME per-cell total (default 9; unused here)
%
% Example:
%   plot_additivity_summary_across_sessions('M1', my_paths, ...
%       'PField','p_boot', 'Method','FDR-BH', 'Alpha',0.05);

    % ----------------------------
    % Parse inputs
    % ----------------------------
    p = inputParser;
    addParameter(p, 'PField', 'p_boot', @(s)ischar(s)||isstring(s));
    addParameter(p, 'Method', 'FDR-BH', @(s)ischar(s)||isstring(s));
    addParameter(p, 'Alpha', 0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
    addParameter(p, 'ShowR2', true, @(x)islogical(x)||ismember(x,[0 1]));
    addParameter(p, 'NBins', 30, @(x)isnumeric(x)&&isscalar(x)&&x>=5);
    addParameter(p, 'NBinsCorr', 15, @(x)isnumeric(x)&&isscalar(x)&&x>=5);
    addParameter(p, 'SECVFolds', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(p, 'PerCellTrain', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(p, 'PerCellTotal', 9, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    parse(p, varargin{:});
    PField         = string(p.Results.PField);
    Method         = string(p.Results.Method);
    Alpha          = p.Results.Alpha;
    ShowR2         = logical(p.Results.ShowR2);
    NBins          = p.Results.NBins;
    NBinsCorr      = p.Results.NBinsCorr;
    Kfold_SE       = p.Results.SECVFolds;
    per_cell_train = p.Results.PerCellTrain;
    per_cell_total = p.Results.PerCellTotal; %#ok<NASGU>

    % ----------------------------
    % Accumulate across sessions
    % ----------------------------
    all_p  = [];   % chosen p-values (p_boot or p_chi2)
    all_q  = [];   % adjusted p-values
    all_r  = [];   % correlation corr_F_Fadd (if present)
    all_pr = [];   % p-value for correlation (if present)
    all_R2 = [];   % R2_add (legacy variance-explained)
    all_eta= [];   % eta_res

    % Legacy pooled variance-explained (kept for compatibility)
    all_R2_SE_MEraw = [];
    all_R2_SE_MEsm  = [];

    % NEW pooled Pearson r / r^2
    all_r_SE_MEraw  = [];
    all_r2_SE_MEraw = [];
    all_r_SE_MEsm   = [];
    all_r2_SE_MEsm  = [];

    % NEW pooled baselines in RAW (loaded from Results)
    all_c_SE   = [];   % weighted SE baseline per neuron (RAW Hz)
    all_c_ME   = [];   % ME HE baseline per neuron (RAW Hz)
    all_c_rat  = [];
    all_c_dif  = [];

    % NEW: pooled peaks (RAW)
    all_peak_ME_HE  = [];
    all_peak_ME_add = [];
    all_peak_SE_add = [];

    % NEW: pooled selectivity class labels
    all_class_labels = strings(0,1);

    sess_id = [];  % session index per neuron
    used_files = strings(0,1);

    % constants taken from the decoding script
    nbins = 6;
    theta_edges   = linspace(-pi, pi, nbins+1) + 0.1; % offset to avoid boundary bleed
    theta_centers = midpoints_circ(theta_edges);
    kappa_grid = [2, 4, 8, 12, 16, 24, 32];
    shrink_tau = 0.3; % heteroscedastic shrinkage weight (fixed)

    class_names  = ["H","E","Both","None"];
    class_colors = [ 0.80 0.20 0.20;
                     0.20 0.20 0.80;
                     0.20 0.60 0.20;
                     0.50 0.50 0.50 ];

    for d = 1:numel(savepaths)
        sp = savepaths{d};
        f1 = fullfile(sp, 'models', sprintf('%s_Additivity_Surface_Test_NoClass.mat', area));
        f2 = fullfile(sp, 'models', sprintf('%s_Additivity_Surface_Test.mat', area));
        file_to_use = '';

        if exist(f1, 'file'), file_to_use = f1;
        elseif exist(f2, 'file'), file_to_use = f2;
        else
            warning('No additivity results found in %s', sp);
        end

        R = [];
        if ~isempty(file_to_use)
            S = load(file_to_use);
            if isfield(S, 'Results')
                R = S.Results;

                % Choose p-values
                if isfield(R, PField)
                    pvals = R.(PField)(:);
                else
                    warning('Field %s not in %s; skipping p-values from this session.', PField, file_to_use);
                    pvals = [];
                end

                % Correlations (optional)
                if isfield(R, 'corr_F_Fadd'), all_r = [all_r; R.corr_F_Fadd(:)]; end
                if isfield(R, 'p_corr_F_Fadd'), all_pr = [all_pr; R.p_corr_F_Fadd(:)]; end

                % R2_add (effect size; legacy)
                if isfield(R, 'R2_add'), all_R2 = [all_R2; R.R2_add(:)]; end
                if isfield(R, 'eta_res'), all_eta = [all_eta; R.eta_res(:)]; end

                % Peaks (RAW) for peak scatters
                if isfield(R,'peak_ME_HE'),  all_peak_ME_HE  = [all_peak_ME_HE;  R.peak_ME_HE(:)]; end
                if isfield(R,'peak_ME_add'), all_peak_ME_add = [all_peak_ME_add; R.peak_ME_add(:)]; end
                if isfield(R,'peak_SE_add'), all_peak_SE_add = [all_peak_SE_add; R.peak_SE_add(:)]; end

                % Baselines (WEIGHTED SE baseline + ME baseline) from Results
                if isfield(R,'c_SE_weighted') && isfield(R,'c_ME_HE')
                    all_c_SE  = [all_c_SE;  R.c_SE_weighted(:)];
                    all_c_ME  = [all_c_ME;  R.c_ME_HE(:)];
                    all_c_rat = [all_c_rat; R.c_SE_weighted(:)./R.c_ME_HE(:)];
                    all_c_dif = [all_c_dif; R.c_ME_HE(:) - R.c_SE_weighted(:)];
                else
                    warning('Missing c_SE_weighted or c_ME_HE in %s; baselines will be unavailable for this session.', file_to_use);
                end

                % Append p-values and session ids
                all_p   = [all_p; pvals];
                sess_id = [sess_id; d*ones(numel(pvals),1)];
                used_files(end+1,1) = string(file_to_use);
            else
                warning('No Results struct in %s', file_to_use);
            end
        end

        % ------------------------------------------------------------
        % Build SE-additive & ME surfaces in RAW units for Pearson r/r^2
        % (SE baseline uses per-neuron c_SE_weighted from Results;
        %  ME baseline equality is guaranteed by loading c_ME_HE above.)
        % ------------------------------------------------------------
        pdata_file = processed_data_file(sp, area);
        if ~exist(pdata_file, 'file')
            warning('Missing processed data for SE-add vs ME comparison in %s; skipping this session.', sp);
            continue;
        end

        D = load(pdata_file);
        if ~isfield(D, 'data')
            warning('Processed data missing ''data'' struct in %s; skipping.', pdata_file);
            continue;
        end
        data = D.data;

        % Indices
        idxME  = find(data.TPi2 ~= 7 & data.TPi1 ~= 7); % ME trials
        idxH   = find(data.TPi2 == 7 & data.TPi1 ~= 7); % SE hand-only
        idxE   = find(data.TPi1 == 7 & data.TPi2 ~= 7); % SE eye-only

        if isempty(idxME) || isempty(idxH) || isempty(idxE)
            warning('Insufficient ME or SE data in %s for comparisons; skipping.', sp);
            continue;
        end

        % Pull FR (RAW)
        FR_me = data.FR(idxME, :);   % [T_me x N]
        FR_h  = data.FR(idxH, :);
        FR_e  = data.FR(idxE, :);
        [T_me, N] = size(FR_me);

        % If any prior normalization existed, un-normalize back to RAW
        have_norm_ME = isfield(data,'norm') && isfield(data.norm,'mu_ME') && isfield(data.norm,'sd_ME') ...
                       && numel(data.norm.mu_ME)==N && numel(data.norm.sd_ME)==N;
        have_norm_SE = isfield(data,'norm') && isfield(data.norm,'mu_SE') && isfield(data.norm,'sd_SE') ...
                       && numel(data.norm.mu_SE)==N && numel(data.norm.sd_SE)==N;

        if have_norm_ME
            muME = data.norm.mu_ME(:)'; sdME = data.norm.sd_ME(:)'; sdME(~isfinite(sdME))=1;
            FR_me_raw = bsxfun(@plus, bsxfun(@times, FR_me, sdME), muME);
        else
            FR_me_raw = FR_me;
        end

        if have_norm_SE
            muSE = data.norm.mu_SE(:)'; sdSE = data.norm.sd_SE(:)'; sdSE(~isfinite(sdSE))=1;
            FR_h_raw = bsxfun(@plus, bsxfun(@times, FR_h, sdSE), muSE);
            FR_e_raw = bsxfun(@plus, bsxfun(@times, FR_e, sdSE), muSE);
        else
            FR_h_raw = FR_h;
            FR_e_raw = FR_e;
        end

        % Target positions
        TP1_me = data.TP1(idxME, :);
        TP2_me = data.TP2(idxME, :);
        TP1_h  = data.TP1(idxH, :);
        TP2_e  = data.TP2(idxE, :);  % <-- fixed MATLAB indexing

        % Selectivity class labels (optional)
        class_file = fullfile(sp, 'models', sprintf('%s_FR_amplitude_selectivity.mat', area));
        if exist(class_file, 'file')
            C = load(class_file);
            if isfield(C, 'class_labels')
                cls = string(C.class_labels(:));
                if numel(cls) ~= N
                    warning('class_labels length (%d) != N (%d) in %s; marking Unknown.', numel(cls), N, class_file);
                    cls = repmat("Unknown", N, 1);
                end
            else
                cls = repmat("Unknown", N, 1);
            end
        else
            cls = repmat("Unknown", N, 1);
        end

        % Angles & binning for ME
        th_me_H = atan2(TP1_me(:,2), TP1_me(:,1));
        th_me_E = atan2(TP2_me(:,2), TP2_me(:,1));
        Hb_me   = discretize(th_me_H, theta_edges);
        Eb_me   = discretize(th_me_E, theta_edges);

        % SE angles (full pools)
        th_h_all = atan2(TP1_h(:,2), TP1_h(:,1));
        th_e_all = atan2(TP2_e(:,2), TP2_e(:,1));
        Hb_h_all = discretize(th_h_all, theta_edges);
        Eb_e_all = discretize(th_e_all, theta_edges);

        % Equalize SE pools (random subsample per bin)
        perbin_h = accumarray(Hb_h_all(~isnan(Hb_h_all)), 1, [nbins,1], @sum, 0);
        perbin_e = accumarray(Eb_e_all(~isnan(Eb_e_all)), 1, [nbins,1], @sum, 0);
        nh_min = max(1, min(perbin_h));
        ne_min = max(1, min(perbin_e));

        keep_h_eq = false(length(th_h_all), 1);
        keep_e_eq = false(length(th_e_all), 1);
        rng(2);
        for ih = 1:nbins
            idx = find(Hb_h_all == ih);
            if ~isempty(idx)
                sel = idx(randperm(numel(idx), min(nh_min, numel(idx))));
                keep_h_eq(sel) = true;
            end
        end
        for ie = 1:nbins
            idx = find(Eb_e_all == ie);
            if ~isempty(idx)
                sel = idx(randperm(numel(idx), min(ne_min, numel(idx))));
                keep_e_eq(sel) = true;
            end
        end

        % RAW units for both pools
        FR_h_eq_raw = FR_h_raw(keep_h_eq, :);
        FR_e_eq_raw = FR_e_raw(keep_e_eq, :);
        th_h_eq = th_h_all(keep_h_eq);
        th_e_eq = th_e_all(keep_e_eq);

        % Balanced ME train split (match nonparametric_MAP)
        train_mask_me = false(T_me,1);
        rng(1);
        for ih = 1:nbins
            for ie = 1:nbins
                idx_cell = find(Hb_me == ih & Eb_me == ie);
                if isempty(idx_cell), continue; end
                Kt  = min(per_cell_train, numel(idx_cell));
                sel = idx_cell(randperm(numel(idx_cell)));
                train_mask_me(sel(1:Kt)) = true;
            end
        end
        FR_me_train_raw = FR_me_raw(train_mask_me, :);
        thH_tr = th_me_H(train_mask_me);
        thE_tr = th_me_E(train_mask_me);

        % Raw ME per-cell means (legacy)
        ME_raw = nan(nbins, nbins, N);
        for ih = 1:nbins
            for ie = 1:nbins
                cell_mask = (Hb_me == ih) & (Eb_me == ie);
                if any(cell_mask)
                    FR_cell_raw = FR_me_raw(cell_mask, :); % RAW
                    ME_raw(ih, ie, :) = mean(FR_cell_raw, 1, 'omitnan');
                end
            end
        end

        % Target grid (6x6 centers)
        [Egrid, Hgrid] = meshgrid(theta_centers, theta_centers);
        Hgrid = Hgrid(:); Egrid = Egrid(:);

        % Per-neuron arrays
        R2_raw_vec = nan(N,1);   % legacy
        R2_sm_vec  = nan(N,1);   % legacy
        r_raw_vec  = nan(N,1);
        r2_raw_vec = nan(N,1);
        r_sm_vec   = nan(N,1);
        r2_sm_vec  = nan(N,1);

        % Require c_SE_weighted from Results for this session
        if isempty(R) || ~isfield(R,'c_SE_weighted') || numel(R.c_SE_weighted) ~= N
            warning('c_SE_weighted missing or wrong length in %s; skipping SE/ME comparison for this session.', file_to_use);
            % still append class labels so sizes match downstream expectations
            all_class_labels = [all_class_labels; repmat("Unknown", N, 1)];
            % also append legacy placeholders for r/R2 so dimensionality is consistent
            all_R2_SE_MEraw = [all_R2_SE_MEraw; R2_raw_vec(:)];
            all_R2_SE_MEsm  = [all_R2_SE_MEsm;  R2_sm_vec(:)];
            all_r_SE_MEraw  = [all_r_SE_MEraw;  r_raw_vec(:)];
            all_r2_SE_MEraw = [all_r2_SE_MEraw; r2_raw_vec(:)];
            all_r_SE_MEsm   = [all_r_SE_MEsm;   r_sm_vec(:)];
            all_r2_SE_MEsm  = [all_r2_SE_MEsm;  r_sm_vec(:).^2];
            continue;
        end
        cSE_weighted_session = R.c_SE_weighted(:);

        % Precompute sizes for SE CV
        nH = numel(th_h_eq);  nE = numel(th_e_eq);

        % Selectivity classes appended once per session (aligned later by masks)
        all_class_labels = [all_class_labels; cls(:)];

        for n = 1:N
            % ---- SE-only CV for kappa ----
            if (nH + nE) < Kfold_SE
                kappa_star_SE = kappa_grid(min(3,length(kappa_grid)));
            else
                cv_idx_h = kfold_partition(max(nH,1), Kfold_SE, 10000+n);
                cv_idx_e = kfold_partition(max(nE,1), Kfold_SE, 20000+n);

                nll_se = zeros(numel(kappa_grid),1);
                for ik = 1:numel(kappa_grid)
                    kap = kappa_grid(ik);
                    fold_nll = 0; valid = false;

                    for k = 1:Kfold_SE
                        % hand-only fold
                        if nH>0
                            val_h = (cv_idx_h == k); tr_h = ~val_h;
                            if any(val_h) && any(tr_h)
                                muH = kernel_mean_1d(th_h_eq(tr_h), FR_h_eq_raw(tr_h,n), th_h_eq(val_h), kap);
                                [~, sH] = kernel_local_sigma_1d(th_h_eq(tr_h), FR_h_eq_raw(tr_h,n), th_h_eq(val_h), kap, shrink_tau);
                                sH = max(sH,1e-3);
                                rH = FR_h_eq_raw(val_h,n);
                                fold_nll = fold_nll + sum(0.5*log(2*pi*sH.^2) + 0.5*((rH - muH).^2)./(sH.^2));
                                valid = true;
                            end
                        end
                        % eye-only fold
                        if nE>0
                            val_e = (cv_idx_e == k); tr_e = ~val_e;
                            if any(val_e) && any(tr_e)
                                muE = kernel_mean_1d(th_e_eq(tr_e), FR_e_eq_raw(tr_e,n), th_e_eq(val_e), kap);
                                [~, sE] = kernel_local_sigma_1d(th_e_eq(tr_e), FR_e_eq_raw(tr_e,n), th_e_eq(val_e), kap, shrink_tau);
                                sE = max(sE,1e-3);
                                rE = FR_e_eq_raw(val_e,n);
                                fold_nll = fold_nll + sum(0.5*log(2*pi*sE.^2) + 0.5*((rE - muE).^2)./(sE.^2));
                                valid = true;
                            end
                        end
                    end

                    if ~valid, fold_nll = inf; end
                    nll_se(ik) = fold_nll;
                end

                [~,bestSE] = min(nll_se);
                if isinf(nll_se(bestSE))
                    kappa_star_SE = kappa_grid(min(3,length(kappa_grid)));
                else
                    kappa_star_SE = kappa_grid(bestSE);
                end
            end

            % ---- Final SE additive surface on the 6 centers (RAW assembly) ----
            A_H = kernel_mean_1d(th_h_eq, FR_h_eq_raw(:,n), theta_centers, kappa_star_SE);
            B_E = kernel_mean_1d(th_e_eq, FR_e_eq_raw(:,n), theta_centers, kappa_star_SE);
            A0  = mean(A_H, 'omitnan');
            B0  = mean(B_E, 'omitnan');

            % >>> USE WEIGHTED SE BASELINE FROM Results <<<
            c_n_SE = cSE_weighted_session(n);

            mu_SE_grid = (A_H(:) - A0) + (B_E(:).' - B0) + c_n_SE;  % 6x6, WEIGHTED baseline
            yhat = mu_SE_grid(:);  % RAW units

            % ---- ME-only CV for kappa (ME TRAINING ONLY, RAW) ----
            r_me_tr = FR_me_train_raw(:, n);
            thH_tr_n = thH_tr; thE_tr_n = thE_tr;
            if numel(r_me_tr) < 10
                kappa_star_ME = kappa_grid(min(3,length(kappa_grid)));
            else
                Kfold_ME = 10;
                cv_idx = kfold_partition(numel(r_me_tr), Kfold_ME, n);
                nll = zeros(numel(kappa_grid),1);
                for ik = 1:numel(kappa_grid)
                    kap = kappa_grid(ik);
                    fold_nll = 0;
                    for k = 1:Kfold_ME
                        val_mask = (cv_idx == k);
                        tr_mask  = ~val_mask;
                        if ~any(val_mask) || ~any(tr_mask), continue; end

                        mu_val = kernel_mean_2d(thH_tr_n(tr_mask), thE_tr_n(tr_mask), r_me_tr(tr_mask), ...
                                                thH_tr_n(val_mask),  thE_tr_n(val_mask), kap);
                        [~, sigma_val] = kernel_local_sigma_2d(thH_tr_n(tr_mask), thE_tr_n(tr_mask), r_me_tr(tr_mask), ...
                                                               thH_tr_n(val_mask),  thE_tr_n(val_mask), kap, shrink_tau);
                        sigma_val = max(sigma_val, 1e-3);
                        rv = r_me_tr(val_mask);
                        fold_nll = fold_nll + sum( 0.5*log(2*pi*sigma_val.^2) + 0.5*((rv - mu_val).^2)./(sigma_val.^2) );
                    end
                    nll(ik) = fold_nll;
                end
                [~,best] = min(nll);
                if isempty(best) || ~isfinite(nll(best))
                    kappa_star_ME = kappa_grid(min(3,length(kappa_grid)));
                else
                    kappa_star_ME = kappa_grid(best);
                end
            end

            % ---- Final ME-smoothed grid (centers) using ME TRAINING ONLY (RAW) ----
            mu_ME_grid = kernel_mean_2d(thH_tr_n, thE_tr_n, r_me_tr, Hgrid, Egrid, kappa_star_ME);
            y_me_sm = mu_ME_grid(:);  % RAW units

            % Raw ME surface vector (per-cell means)
            y_raw = reshape(ME_raw(:,:,n), [], 1);

            % ---- Legacy R2 computations ----
            m_raw = isfinite(y_raw) & isfinite(yhat);
            if any(m_raw)
                SSE = sum((y_raw(m_raw)-yhat(m_raw)).^2);
                SST = sum((y_raw(m_raw)-mean(y_raw(m_raw))).^2);
                R2_raw_vec(n) = 1 - SSE/max(SST,eps);
            end

            m_sm = isfinite(y_me_sm) & isfinite(yhat);
            if any(m_sm)
                SSEs = sum((y_me_sm(m_sm)-yhat(m_sm)).^2);
                SSTs = sum((y_me_sm(m_sm)-mean(y_me_sm(m_sm))).^2);
                R2_sm_vec(n) = 1 - SSEs/max(SSTs,eps);
            end

            % ---- NEW: Pearson r / r^2 (shape-only) ----
            r_raw = nan_pearsonr(y_raw, yhat);
            r_sm  = nan_pearsonr(y_me_sm, yhat);
            r_raw_vec(n) = r_raw;
            r2_raw_vec(n)= r_raw.^2;
            r_sm_vec(n)  = r_sm;
            r2_sm_vec(n) = r_sm.^2;
        end

        % Append pooled (PER-NEURON)
        all_R2_SE_MEraw = [all_R2_SE_MEraw; R2_raw_vec(:)];
        all_R2_SE_MEsm  = [all_R2_SE_MEsm;  R2_sm_vec(:)];
        all_r_SE_MEraw  = [all_r_SE_MEraw;  r_raw_vec(:)];
        all_r2_SE_MEraw = [all_r2_SE_MEraw; r2_raw_vec(:)];
        all_r_SE_MEsm   = [all_r_SE_MEsm;   r_sm_vec(:)];
        all_r2_SE_MEsm  = [all_r2_SE_MEsm;  r_sm_vec(:).^2];
    end

    if isempty(all_p)
        error('No p-values aggregated. Check that your per-session scripts were run and saved.');
    end

    % Drop NaNs defensively
    mask_p = isfinite(all_p);
    p_vec  = all_p(mask_p);

    % ----------------------------
    % Multiple-comparisons correction (global)
    % ----------------------------
    switch lower(Method)
        case {'fdr-bh','fdr','bh'}
            [sig_bh, q_vec, p_star] = fdr_bh_threshold(p_vec, Alpha);
            thr_label = sprintf('BH cutoff p* = %.2g', p_star);
            all_q = nan(size(all_p)); all_q(mask_p) = q_vec;
            sig_global = false(size(all_p)); sig_global(mask_p) = sig_bh;
            xline_pos = p_star;
        case 'bonferroni'
            m = sum(mask_p);
            p_line = Alpha / max(m,1);
            thr_label = sprintf('Bonferroni p_{thr} = %.2g (m=%d)', p_line, m);
            all_q = nan(size(all_p));
            all_q(mask_p) = min(1, p_vec * m);
            sig_global = false(size(all_p));
            sig_global(mask_p) = p_vec <= p_line;
            xline_pos = p_line;
        otherwise
            error('Unknown Method: %s', Method);
    end

    % ----------------------------
    % Output directories
    % ----------------------------
    plotdir = fullfile(savepaths{1}, '..', 'plots', 'summary', ...
        'additive_reconstruction');
    if ~exist(plotdir, 'dir'), mkdir(plotdir); end

    % ----------------------------
    % 1) Histogram of correlations
    % ----------------------------
    if ~isempty(all_r)
        rvals = all_r(isfinite(all_r));
        if ~isempty(rvals)
            mu_r = mean(rvals, 'omitnan');
            fprintf('[%s] corr(F, F_add): mean r = %.4f  (N=%d)\n', area, mu_r, numel(rvals));

            figure('Visible','off');
            set(gca, 'FontSize', 20); hold on;
            histogram(rvals, NBinsCorr, 'BinLimits', [-1 1]);
            xlim([-1 1]);
            xline(mu_r, 'k-', 'LineWidth', 2);
            yl = ylim;
            text(mu_r, yl(2)*0.95, sprintf('mean = %.3f', mu_r), ...
                 'Color','k', 'Rotation',90, 'HorizontalAlignment','right', 'VerticalAlignment','top');
            xlabel('corr(F, F_{add})'); ylabel('Count');
            title(sprintf('%s: Correlation (N=%d)', area, numel(rvals)));
            print(fullfile(plotdir, sprintf('%s_corr_hist.pdf', area)), '-dpdf','-painters');
            close;
        end
    end

    % ----------------------------
    % 2) Histogram of p-values + threshold
    % ----------------------------
    figure('Visible','off');
    set(gca, 'FontSize', 20); hold on;
    histogram(p_vec, NBins, 'BinLimits',[0 1]);
    xlim([0 1]);
    xlabel(sprintf('%s p-values (%s)', PField, Method));
    ylabel('Count');
    title(sprintf('%s: p-value distribution (N=%d)', area, numel(p_vec)));
    if isfinite(xline_pos) && xline_pos > 0 && xline_pos < 1
        xline(xline_pos, 'r--');
        yl = ylim;
        text(xline_pos, yl(2)*0.95, thr_label, 'Color','r', 'Rotation',90, ...
             'HorizontalAlignment','right', 'VerticalAlignment','top');
    end
    print(fullfile(plotdir, sprintf('%s_pvalues_%s_%s.pdf', area, PField, lower(Method))), '-dpdf','-painters');
    close;

    % ----------------------------
    % 3) Histogram of R2_add (optional)  << ADDED MEAN PRINT + MEAN LINE >>
    % ----------------------------
    if ShowR2 && ~isempty(all_R2)
        r2vals = all_R2(isfinite(all_R2));
        if ~isempty(r2vals)
            mu_R2add = mean(r2vals, 'omitnan');
            fprintf('[%s] R2_add: mean = %.4f  (N=%d)\n', area, mu_R2add, numel(r2vals));

            figure('Visible','off');
            set(gca, 'FontSize', 20); hold on;
            histogram(r2vals, NBins, 'BinLimits',[0 1]);
            xlim([0 1]);
            xline(mu_R2add, 'k-', 'LineWidth', 2);
            yl = ylim;
            text(mu_R2add, yl(2)*0.95, sprintf('mean = %.3f', mu_R2add), ...
                 'Color','k', 'Rotation',90, 'HorizontalAlignment','right', 'VerticalAlignment','top');
            xlabel('R^2_{add}'); ylabel('Count');
            title(sprintf('%s: Additive R^2 (N=%d)', area, numel(r2vals)));
            print(fullfile(plotdir, sprintf('%s_R2_add.pdf', area)), '-dpdf','-painters');
            close;
        end
    end

    % ----------------------------
    % 4) Histogram of Pearson r by class  << ADDED OVERALL MEAN PRINT + MEAN LINE >>
    % ----------------------------
    if ~isempty(all_r_SE_MEsm)
        mask_r = isfinite(all_r_SE_MEsm);
        rvals3 = all_r_SE_MEsm(mask_r);
        cls_all = all_class_labels(mask_r);
        if ~isempty(rvals3)
            mu_r_SE = mean(rvals3, 'omitnan');
            fprintf('[%s] Pearson r (SE_add^{weighted} vs ME-smoothed): mean = %.4f  (N=%d)\n', ...
                    area, mu_r_SE, numel(rvals3));

            figure('Visible','off');
            set(gca, 'FontSize', 20); hold on;
            edges = linspace(-1, 1, NBins+1);
            legends = strings(0,1);
            for ci = 1:numel(class_names)
                cname = class_names(ci);
                vals  = rvals3(cls_all == cname);
                if ~isempty(vals)
                    histogram(vals, 'BinEdges', edges);
                    legends(end+1) = sprintf('%s (N=%d)', cname, numel(vals));
                end
            end
            xlim([-1 1]);

            % overall mean line (single number to reference in text)
            xline(mu_r_SE, 'k-', 'LineWidth', 2);
            yl = ylim;
            text(mu_r_SE, yl(2)*0.95, sprintf('mean = %.3f', mu_r_SE), ...
                 'Color','k', 'Rotation',90, 'HorizontalAlignment','right', 'VerticalAlignment','top');

            xlabel('Pearson r (SE_{add}^{weighted} vs ME-smoothed)'); ylabel('Count');
            title(sprintf('%s: r by class', area));
            if ~isempty(legends), legend(legends, 'Location','northwest'); end
            print(fullfile(plotdir, sprintf('%s_r_SEaddWeighted_vs_MEsmooth_byClass.pdf', area)), '-dpdf','-painters');
            close;
        end

        % Manuscript Figure 4: surface similarity for mixed-selective
        % (Both-class) neurons only. Keep this separate from the general
        % four-class diagnostic above so that neither plot changes meaning.
        both_mask = mask_r & all_class_labels == "Both";
        both_r = all_r_SE_MEsm(both_mask);
        if strcmpi(area, 'SPL') && ~isempty(both_r)
            both_color = class_colors(class_names == "Both", :);
            h_both = figure('Visible','off','PaperUnits','inches', ...
                'PaperSize',[5 4], 'PaperPosition',[0.25 0.25 4.5 3.5], ...
                'PaperPositionMode','manual');
            ax_both = axes('Parent', h_both); hold(ax_both, 'on');
            set(ax_both, 'FontSize', 14);
            histogram(ax_both, both_r, 'BinEdges', linspace(-1, 1, NBins+1), ...
                'FaceColor', both_color, 'EdgeColor', 'none');
            xlim(ax_both, [-1 1]);
            xlabel(ax_both, 'Correlation (r)');
            ylabel(ax_both, 'Neuron count');
            title(ax_both, sprintf('(N=%d)', numel(both_r)));
            legend(ax_both, 'Mixed Selective Neurons', 'Location', 'northwest');
            grid(ax_both, 'on'); box(ax_both, 'on');
            print(h_both, fullfile(plotdir, sprintf( ...
                '%s_r_SEaddWeighted_vs_MEsmooth_Both.pdf', area)), ...
                '-dpdf', '-painters');
            close(h_both);
        end
    end

    % -------------------------
    % 5) Baseline scatter: c_SE^weighted vs c_ME
    % -------------------------
    if ~isempty(all_c_ME) && ~isempty(all_c_SE)
        maskc = isfinite(all_c_ME) & isfinite(all_c_SE);
        figure('Visible','off');
        set(gca, 'FontSize', 20); hold on;
        if any(maskc)
            y = all_c_SE(maskc);
            x = all_c_ME(maskc);
            scatter(x, y, "filled"); hold on;
            lim_max = max([x; y]); if ~isfinite(lim_max) || lim_max<=0, lim_max = 1; end
            xlim([0 lim_max]); ylim([0 lim_max]);
            plot([0 lim_max], [0 lim_max], 'k--');
            ab = [x ones(numel(x),1)] \ y; a = ab(1); b = ab(2);
            xs = linspace(0, lim_max, 200);
            plot(xs, a*xs + b, 'r-');
            xlabel('c_{ME,HE}  [Hz]');
            ylabel('c_{SE}^{weighted}  [Hz]');
            title(sprintf('%s: Baseline (s=%.2f, b=%.2f)', area, a, b));
            grid on; box on;
        else
            text(0.5,0.5,'No overlapping finite pairs','HorizontalAlignment','center');
        end
        print(fullfile(plotdir, sprintf('%s_cME_vs_cSEweighted_scatter.pdf', area)), '-dpdf','-painters');
        close;
    end

    % -------------------------
    % 6) Peak FR scatters (two subplots)
    % -------------------------
    if ~isempty(all_peak_ME_HE) && (~isempty(all_peak_ME_add) || ~isempty(all_peak_SE_add))
        maskA = isfinite(all_peak_ME_HE) & isfinite(all_peak_ME_add);
        maskS = isfinite(all_peak_ME_HE) & isfinite(all_peak_SE_add);

        % Fix the hidden figure geometry so the two legacy panels retain the
        % portrait, near-square aspect used for the manuscript export.
        h_peak = figure('Visible','off','Units','pixels', ...
            'Position',[100 100 500 1000], ...
            'PaperUnits','inches','PaperPosition',[0 0 5 10]);
        % (a) ME_add vs ME_raw
        subplot(2,1,1);
        set(gca, 'FontSize', 20); hold on;
        if any(maskA)
            x = all_peak_ME_HE(maskA);
            y = all_peak_ME_add(maskA);
            scatter(x, y, "filled"); hold on;
            lim_max = max([x; y]); if ~isfinite(lim_max) || lim_max<=0, lim_max=1; end
            xlim([0 lim_max]); ylim([0 lim_max]);
            plot([0 lim_max],[0 lim_max],'k--');
            ab = [x ones(numel(x),1)] \ y; a=ab(1); b=ab(2);
            xs = linspace(0, lim_max, 200);
            plot(xs, a*xs + b, 'r-');
            xlabel('peak(F_{ME,raw})  [Hz]');
            ylabel('peak(F_{ME,add})  [Hz]');
            title(sprintf('%s: ME_{add} vs ME_{raw}', area));
            grid on; box on;
        else
            text(0.5,0.5,'No finite pairs','HorizontalAlignment','center');
        end

        % (b) SE_add vs ME_raw
        subplot(2,1,2);
        set(gca, 'FontSize', 20); hold on;
        if any(maskS)
            x = all_peak_ME_HE(maskS);
            y = all_peak_SE_add(maskS);
            scatter(x, y, "filled"); hold on;
            lim_max = max([x; y]); if ~isfinite(lim_max) || lim_max<=0, lim_max=1; end
            xlim([0 lim_max]); ylim([0 lim_max]);
            plot([0 lim_max],[0 lim_max],'k--');
            ab = [x ones(numel(x),1)] \ y; a=ab(1); b=ab(2);
            xs = linspace(0, lim_max, 200);
            plot(xs, a*xs + b, 'r-');
            xlabel('peak(F_{ME,raw})  [Hz]');
            ylabel('peak(F_{SE,add})  [Hz]');
            title(sprintf('%s: SE_{add} vs ME_{raw}', area));
            grid on; box on;
        else
            text(0.5,0.5,'No finite pairs','HorizontalAlignment','center');
        end

        print(fullfile(plotdir, sprintf('%s_peak_MEadd_and_SEadd_vs_MEraw_scatters.pdf', area)), '-dpdf','-painters');
        close(h_peak);
    end

    % ----------------------------
    % Save pooled vectors + MCP outputs
    % ----------------------------
    Summary = struct();
    Summary.area      = area;
    Summary.usedFiles = used_files;
    Summary.method    = Method;
    Summary.alpha     = Alpha;
    Summary.PField    = PField;

    Summary.p_raw     = all_p(:);
    Summary.p_mask    = mask_p(:);
    Summary.p_used    = p_vec(:);
    Summary.q_adj     = all_q(:);
    Summary.sig       = sig_global(:);
    Summary.threshold = xline_pos;

    if ~isempty(all_r),      Summary.corr_F_Fadd   = all_r(:);  end
    if ~isempty(all_pr),     Summary.p_corr_F_Fadd = all_pr(:); end
    if ~isempty(all_R2),     Summary.R2_add        = all_R2(:); end
    if ~isempty(all_eta),    Summary.eta_res       = all_eta(:);end

    % Legacy (kept)
    if ~isempty(all_R2_SE_MEraw), Summary.R2_SE_MEraw = all_R2_SE_MEraw(:); end
    if ~isempty(all_R2_SE_MEsm),  Summary.R2_SE_MEsm  = all_R2_SE_MEsm(:);  end

    % NEW Pearson metrics
    if ~isempty(all_r_SE_MEraw),  Summary.r_SE_MEraw  = all_r_SE_MEraw(:);  end
    if ~isempty(all_r2_SE_MEraw), Summary.r2_SE_MEraw = all_r2_SE_MEraw(:); end
    if ~isempty(all_r_SE_MEsm),   Summary.r_SE_MEsm   = all_r_SE_MEsm(:);   end
    if ~isempty(all_r2_SE_MEsm),  Summary.r2_SE_MEsm  = all_r2_SE_MEsm(:);  end

    % NEW baselines & peaks (weighted SE baseline + loaded ME baseline)
    if ~isempty(all_c_SE),        Summary.c_SE_weighted = all_c_SE(:); Summary.c_SE = all_c_SE(:); end
    if ~isempty(all_c_ME),        Summary.c_ME_HE       = all_c_ME(:); Summary.c_ME = all_c_ME(:); end
    if ~isempty(all_c_rat),       Summary.c_ratio       = all_c_rat(:); end
    if ~isempty(all_c_dif),       Summary.c_diff        = all_c_dif(:); end
    if ~isempty(all_peak_ME_HE),  Summary.peak_ME_HE    = all_peak_ME_HE(:);  end
    if ~isempty(all_peak_ME_add), Summary.peak_ME_add   = all_peak_ME_add(:); end
    if ~isempty(all_peak_SE_add), Summary.peak_SE_add   = all_peak_SE_add(:); end

    Summary.session_id = sess_id(:);

    outmat = fullfile(savepaths{1}, '..', 'models', sprintf('%s_Additivity_Summary_Accum.mat', area));
    save(outmat, 'Summary', '-v7.3');
    fprintf('Saved pooled summary: %s\n', outmat);
end

% =============================
% Helpers (stats/math only)
% =============================
function [sig, q, p_star] = fdr_bh_threshold(p, alpha)
    p = p(:);
    [ps, idx] = sort(p);
    m = numel(ps);
    ranks = (1:m)';

    thresh = (ranks/m)*alpha;
    pass = ps <= thresh;
    k = find(pass, 1, 'last');

    sig = false(m,1);
    if ~isempty(k), sig(1:k) = true; end

    q = m./ranks .* ps;
    q = flipud(cummin(flipud(q)));
    q = min(q, 1);

    unsig = false(size(p));
    unsig(idx) = sig;
    uq = nan(size(p));
    uq(idx) = q;

    sig = unsig;
    q   = uq;
    p_star = ~isempty(k) * ps(max(k,1));
    if isempty(k), p_star = NaN; end
end

function mu = midpoints_circ(edges)
    mu = edges(1:end-1) + (edges(2)-edges(1))/2;
    mu = wrapToPi(mu);
end

function y = wrapToPi(x)
    y = mod(x+pi, 2*pi) - pi;
end

function w = vm_weight(delta, kappa)
    w = exp(kappa * cos(delta));
end

function mu = kernel_mean_1d(theta_train, r_train, theta_eval, kappa)
    mu = nan(size(theta_eval));
    if isempty(theta_train) || isempty(r_train)
        return;
    end
    for i = 1:numel(theta_eval)
        del = wrapToPi(theta_train - theta_eval(i));
        w = vm_weight(del, kappa);
        s = sum(w);
        if s > 0
            mu(i) = sum(w .* r_train) / s;
        else
            mu(i) = mean(r_train, 'omitnan');
        end
    end
end

function mu = kernel_mean_2d(thH_tr, thE_tr, r_tr, thH_eval, thE_eval, kappa)
    mu = nan(size(thH_eval));
    if isempty(thH_tr) || isempty(thE_tr) || isempty(r_tr)
        return;
    end
    for i = 1:numel(thH_eval)
        dH = wrapToPi(thH_tr - thH_eval(i));
        dE = wrapToPi(thE_tr - thE_eval(i));
        w  = vm_weight(dH, kappa) .* vm_weight(dE, kappa);
        s  = sum(w);
        if s > 0
            mu(i) = sum(w .* r_tr) / s;
        else
            mu(i) = mean(r_tr, 'omitnan');
        end
    end
end

function [mu_loc, sigma_loc] = kernel_local_sigma_2d(thH_tr, thE_tr, r_tr, thH_eval, thE_eval, kappa, tau)
    mu_loc    = nan(size(thH_eval));
    sigma_loc = nan(size(thH_eval));
    if isempty(r_tr)
        sigma_loc(:) = 1; return;
    end
    vg = var(r_tr, 'omitnan');  if ~isfinite(vg) || vg<=0, vg = 1; end
    for i = 1:numel(thH_eval)
        dH = wrapToPi(thH_tr - thH_eval(i));
        dE = wrapToPi(thE_tr - thE_eval(i));
        w  = vm_weight(dH, kappa) .* vm_weight(dE, kappa);
        s  = sum(w);
        if s > 0
            mu_i = sum(w .* r_tr) / s;
            res  = r_tr - mu_i;
            vloc = sum(w .* (res.^2)) / s;
        else
            mu_i = mean(r_tr, 'omitnan');
            vloc = vg;
        end
        mu_loc(i)    = mu_i;
        sigma_loc(i) = sqrt( tau*vg + (1-tau)*max(vloc,1e-6) );
    end
end

function [mu_loc, sigma_loc] = kernel_local_sigma_1d(theta_tr, r_tr, theta_eval, kappa, tau)
    mu_loc    = nan(size(theta_eval));
    sigma_loc = nan(size(theta_eval));
    if isempty(r_tr)
        sigma_loc(:) = 1; return;
    end
    vg = var(r_tr, 'omitnan');  if ~isfinite(vg) || vg<=0, vg = 1; end
    for i = 1:numel(theta_eval)
        del = wrapToPi(theta_tr - theta_eval(i));
        w   = vm_weight(del, kappa);
        s   = sum(w);
        if s > 0
            mu_i = sum(w .* r_tr) / s;
            res  = r_tr - mu_i;
            vloc = sum(w .* (res.^2)) / s;
        else
            mu_i = mean(r_tr, 'omitnan');
            vloc = vg;
        end
        mu_loc(i)    = mu_i;
        sigma_loc(i) = sqrt( tau*vg + (1-tau)*max(vloc,1e-6) );
    end
end

function idx = kfold_partition(n, K, seed)
    rng(seed);
    p = randperm(max(n,1));
    base = floor(n/K);
    rem  = mod(n, K);
    idx = zeros(n,1);
    b = 0;
    for k = 1:K
        nk = base + (k <= rem);
        if nk>0
            idx(p(b+1:b+nk)) = k;
        end
        b = b + nk;
    end
end

function r = nan_pearsonr(a, b)
% Safe Pearson r with NaN handling; returns NaN if <2 overlapping finite samples
    a = a(:); b = b(:);
    m = isfinite(a) & isfinite(b);
    if sum(m) < 2
        r = NaN; return;
    end
    aa = a(m); bb = b(m);
    sa = std(aa,0,1); sb = std(bb,0,1);
    if ~isfinite(sa) || ~isfinite(sb) || sa==0 || sb==0
        r = NaN; return;
    end
    r = ( (aa-mean(aa))' * (bb-mean(bb)) ) / ( (numel(aa)-1) * sa * sb );
end
