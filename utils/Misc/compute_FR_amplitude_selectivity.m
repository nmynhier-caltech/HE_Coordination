%% Submission Version
function compute_FR_amplitude_selectivity(area, savepath)
% Computes direction-tuning amplitude selectivity for hand (H) and eye (E),
% for both multi-effector (HE) and single-effector trials, using:
%   - Shuffle-based nulls of amplitude (max-min across direction bins)
%   - z-scores of observed amplitude vs its null (per neuron)
%   - Trial bootstrap to obtain robust per-neuron medians and PDs
%   - Multiple comparisons correction (BH-FDR) on empirical p-values
%
% Assumes preprocessing saved RAW firing rates (Hz) in:
%   <publicDataRoot>/<session_id>_<area>_Processed_data.mat
%
% Outputs (saved to <savepath>/models/<area>_FR_amplitude_selectivity.mat):
%   AMP, AMP_null, AMP_single, AMP_null_single
%   thresh, sig, class_labels, colors
%   thresh_single, sig_single, class_labels_single, colors_single
%   MAXDIR_multi, MAXDIR_single                     (robust PDs from boot)
%   BOOT, CI, MCP                                   (full provenance)
%   z_med_multi, z_med_single                       (median z per neuron)
%
% Notes:
%   • All computations use RAW Hz; the “z-scores” here are *null-standardized
%     effect sizes* (observed amplitude z-scored against its shuffle-null).
%   • Families are corrected separately: {Multi-H, Multi-E, Single-H, Single-E}.

    % ----------------------------
    % Multiple-comparisons config
    % ----------------------------
    MCP.method = 'FDR-BH';   % 'FDR-BH' | 'bonferroni' | 'none'
    MCP.alpha  = 0.05;

    % ----------------------------
    % Load preprocessed RAW data
    % ----------------------------
    data = load(processed_data_file(savepath, area)).data;

    % ============================
    % Multi-effector (HE) trials
    % ============================
    multiTP = find(data.TPi2 ~= 7 & data.TPi1 ~= 7);
    FR      = data.FR(multiTP, :);              % [trials x neurons], raw Hz
    pos_H   = data.TP1(multiTP, :);
    pos_E   = data.TP2(multiTP, :);

    num_neurons = size(FR, 2);
    nbins  = 6;
    nshuff = 500;

    % Directions and binning
    thetaH    = atan2(pos_H(:,2), pos_H(:,1));
    thetaE    = atan2(pos_E(:,2), pos_E(:,1));
    bin_edges = linspace(-pi, pi, nbins+1) + 0.1;     % conventional offset
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

    thetaH_bin = discretize(thetaH, bin_edges);
    thetaE_bin = discretize(thetaE, bin_edges);

    % Bootstrap config
    BOOT_nboot = 500;            % number of trial resamples
    BOOT_alpha = 0.05;           % CI level (two-sided)
    BOOT_store_full_FR = true;   % set false if storage is a concern

    % Preallocate (parfor-friendly)
    FR_by_dir_H = nan(num_neurons, nbins);
    FR_by_dir_E = nan(num_neurons, nbins);
    AMP_H       = nan(num_neurons, 1);
    AMP_E       = nan(num_neurons, 1);
    AMP_null_H  = nan(num_neurons, nshuff);
    AMP_null_E  = nan(num_neurons, nshuff);

    CI_FR_by_dir_H = nan(num_neurons, nbins, 2); % [lo hi]
    CI_FR_by_dir_E = nan(num_neurons, nbins, 2);
    BOOT_AMP_H       = nan(num_neurons, BOOT_nboot);
    BOOT_AMP_E       = nan(num_neurons, BOOT_nboot);
    BOOT_THRESH_multi_H = nan(num_neurons, BOOT_nboot); % per-boot 95th z-null
    BOOT_THRESH_multi_E = nan(num_neurons, BOOT_nboot);
    BOOT_MAXDIR_multi_H = nan(num_neurons, BOOT_nboot);
    BOOT_MAXDIR_multi_E = nan(num_neurons, BOOT_nboot);
    BOOT_P_multi_H   = nan(num_neurons, BOOT_nboot);
    BOOT_P_multi_E   = nan(num_neurons, BOOT_nboot);

    if BOOT_store_full_FR
        BOOT_FR_by_dir_H = nan(num_neurons, nbins, BOOT_nboot);
        BOOT_FR_by_dir_E = nan(num_neurons, nbins, BOOT_nboot);
        % (Single-effector containers are created later)
    end

    % Broadcast indices for bootstrap
    validH_all = ~isnan(thetaH_bin);
    validE_all = ~isnan(thetaE_bin);
    idx_all_H  = find(validH_all);
    idx_all_E  = find(validE_all);

    % ======= parfor over neurons: Multi-effector =======
    parfor n = 1:num_neurons
        rng(n, 'twister');  % deterministic per neuron

        fr = FR(:, n);  % raw Hz

        % ----- Hand (multi) point estimates + shuffle-null -----
        FRH_row = nan(1, nbins);
        for b = 1:nbins
            FRH_row(b) = mean(fr(thetaH_bin == b), 'omitnan');
        end
        amp_obs_H = max(FRH_row) - min(FRH_row);

        amp_null_H_vec = nan(1, nshuff);
        for s = 1:nshuff
            shuff_bin = thetaH_bin(randperm(length(thetaH_bin)));
            fr_by_bin = accumarray(shuff_bin(~isnan(shuff_bin)), fr(~isnan(shuff_bin)), ...
                                   [nbins, 1], @mean, NaN);
            amp_null_H_vec(s) = max(fr_by_bin) - min(fr_by_bin);
        end
        [z_obs_H, z_null_H] = zscore_against_null(amp_obs_H, amp_null_H_vec);
        AMP_H(n)         = z_obs_H;
        AMP_null_H(n, :) = z_null_H;

        % ----- Eye (multi) point estimates + shuffle-null -----
        FRE_row = nan(1, nbins);
        for b = 1:nbins
            FRE_row(b) = mean(fr(thetaE_bin == b), 'omitnan');
        end
        amp_obs_E = max(FRE_row) - min(FRE_row);

        amp_null_E_vec = nan(1, nshuff);
        for s = 1:nshuff
            shuff_bin = thetaE_bin(randperm(length(thetaE_bin)));
            fr_by_bin = accumarray(shuff_bin(~isnan(shuff_bin)), fr(~isnan(shuff_bin)), ...
                                   [nbins, 1], @mean, NaN);
            amp_null_E_vec(s) = max(fr_by_bin) - min(fr_by_bin);
        end
        [z_obs_E, z_null_E] = zscore_against_null(amp_obs_E, amp_null_E_vec);
        AMP_E(n)         = z_obs_E;
        AMP_null_E(n, :) = z_null_E;

        FR_by_dir_H(n, :) = FRH_row;
        FR_by_dir_E(n, :) = FRE_row;

        % ----- Bootstrap (multi) -----
        % Hand
        FR_draws_H = nan(BOOT_nboot, nbins);
        for bb = 1:BOOT_nboot
            idxb = idx_all_H(randi(numel(idx_all_H), numel(idx_all_H), 1));
            sh_bin = thetaH_bin(idxb);
            frb    = fr(idxb);
            fr_by_bin_b = accumarray(sh_bin, frb, [nbins,1], @mean, NaN);
            amp_obs_b = max(fr_by_bin_b) - min(fr_by_bin_b);

            amp_null_b = nan(1, nshuff);
            rb = numel(sh_bin);
            for s = 1:nshuff
                shuff = sh_bin(randperm(rb));
                fr_by_bin_s = accumarray(shuff, frb, [nbins,1], @mean, NaN);
                amp_null_b(s) = max(fr_by_bin_s) - min(fr_by_bin_s);
            end
            [z_obs_b, z_null_b] = zscore_against_null(amp_obs_b, amp_null_b);

            % per-boot empirical p with +1 correction
            valid = isfinite(z_null_b);
            m = sum(valid);
            if m > 0 && isfinite(z_obs_b)
                BOOT_P_multi_H(n, bb) = (1 + sum(z_null_b(valid) >= z_obs_b)) / (m + 1);
            else
                BOOT_P_multi_H(n, bb) = NaN;
            end

            BOOT_AMP_H(n, bb)        = z_obs_b;
            BOOT_THRESH_multi_H(n,bb)= prctile(z_null_b, 95);
            FR_draws_H(bb, :)        = fr_by_bin_b';
            % PD per-boot (vector-sum)
            BOOT_MAXDIR_multi_H(n, bb) = angle(sum(fr_by_bin_b(:) .* exp(1i*bin_centers(:))));
        end
        P = prctile(FR_draws_H, [100*BOOT_alpha/2, 100*(1-BOOT_alpha/2)], 1);
        CI_FR_by_dir_H(n, :, :) = permute(P, [3 2 1]);
        if BOOT_store_full_FR
            BOOT_FR_by_dir_H(n, :, :) = reshape(FR_draws_H.', [1, nbins, BOOT_nboot]);
        end

        % Eye
        FR_draws_E = nan(BOOT_nboot, nbins);
        for bb = 1:BOOT_nboot
            idxb = idx_all_E(randi(numel(idx_all_E), numel(idx_all_E), 1));
            sh_bin = thetaE_bin(idxb);
            frb    = fr(idxb);
            fr_by_bin_b = accumarray(sh_bin, frb, [nbins,1], @mean, NaN);
            amp_obs_b = max(fr_by_bin_b) - min(fr_by_bin_b);

            amp_null_b = nan(1, nshuff);
            rb = numel(sh_bin);
            for s = 1:nshuff
                shuff = sh_bin(randperm(rb));
                fr_by_bin_s = accumarray(shuff, frb, [nbins,1], @mean, NaN);
                amp_null_b(s) = max(fr_by_bin_s) - min(fr_by_bin_s);
            end
            [z_obs_b, z_null_b] = zscore_against_null(amp_obs_b, amp_null_b);

            valid = isfinite(z_null_b);
            m = sum(valid);
            if m > 0 && isfinite(z_obs_b)
                BOOT_P_multi_E(n, bb) = (1 + sum(z_null_b(valid) >= z_obs_b)) / (m + 1);
            else
                BOOT_P_multi_E(n, bb) = NaN;
            end

            BOOT_AMP_E(n, bb)        = z_obs_b;
            BOOT_THRESH_multi_E(n,bb)= prctile(z_null_b, 95);
            FR_draws_E(bb, :)        = fr_by_bin_b';
            BOOT_MAXDIR_multi_E(n, bb) = angle(sum(fr_by_bin_b(:) .* exp(1i*bin_centers(:))));
        end
        P = prctile(FR_draws_E, [100*BOOT_alpha/2, 100*(1-BOOT_alpha/2)], 1);
        CI_FR_by_dir_E(n, :, :) = permute(P, [3 2 1]);
        if BOOT_store_full_FR
            BOOT_FR_by_dir_E(n, :, :) = reshape(FR_draws_E.', [1, nbins, BOOT_nboot]);
        end
    end

    % Robust PD for multi (circular mean over boot-angles)
    MAXDIR_multi.H = robust_pd_from_boot(BOOT_MAXDIR_multi_H);
    MAXDIR_multi.E = robust_pd_from_boot(BOOT_MAXDIR_multi_E);

    % ============================
    % Single-effector trials
    % ============================
    onlyTP1 = find(data.TPi2 == 7 & data.TPi1 ~= 7);  % hand-only
    onlyTP2 = find(data.TPi1 == 7 & data.TPi2 ~= 7);  % eye-only

    pos_H_single = data.TP1(onlyTP1, :);
    pos_E_single = data.TP2(onlyTP2, :);
    FR_H_single  = data.FR(onlyTP1, :);               % raw Hz
    FR_E_single  = data.FR(onlyTP2, :);               % raw Hz

    thetaH_single      = atan2(pos_H_single(:,2), pos_H_single(:,1));
    thetaE_single      = atan2(pos_E_single(:,2), pos_E_single(:,1));
    thetaH_bin_single  = discretize(thetaH_single, bin_edges);
    thetaE_bin_single  = discretize(thetaE_single, bin_edges);

    validH_single = ~isnan(thetaH_bin_single);
    validE_single = ~isnan(thetaE_bin_single);
    idxH_all = find(validH_single);
    idxE_all = find(validE_single);

    AMP_single_H   = nan(num_neurons, 1);
    AMP_single_E   = nan(num_neurons, 1);
    MAXDIR_single_H= nan(num_neurons, 1);
    MAXDIR_single_E= nan(num_neurons, 1);
    AMP_null_single_H = nan(num_neurons, nshuff);
    AMP_null_single_E = nan(num_neurons, nshuff);

    % placeholders; filled with medians over boot thresholds
    thresh_single.H = nan(num_neurons, 1);
    thresh_single.E = nan(num_neurons, 1);
    sig_single.H    = false(num_neurons, 1);
    sig_single.E    = false(num_neurons, 1);

    CI_FR_by_dir_single_H = nan(num_neurons, nbins, 2);
    CI_FR_by_dir_single_E = nan(num_neurons, nbins, 2);
    BOOT_AMP_single_H     = nan(num_neurons, BOOT_nboot);
    BOOT_AMP_single_E     = nan(num_neurons, BOOT_nboot);
    BOOT_THRESH_single_H  = nan(num_neurons, BOOT_nboot);
    BOOT_THRESH_single_E  = nan(num_neurons, BOOT_nboot);
    BOOT_MAXDIR_single_H  = nan(num_neurons, BOOT_nboot);
    BOOT_MAXDIR_single_E  = nan(num_neurons, BOOT_nboot);
    BOOT_P_single_H       = nan(num_neurons, BOOT_nboot);
    BOOT_P_single_E       = nan(num_neurons, BOOT_nboot);

    if BOOT_store_full_FR
        BOOT_FR_by_dir_single_H = nan(num_neurons, nbins, BOOT_nboot);
        BOOT_FR_by_dir_single_E = nan(num_neurons, nbins, BOOT_nboot);
    end

    % ======= parfor over neurons: Single-effector =======
    parfor n = 1:num_neurons
        rng(n, 'twister');

        % ----- Hand-only -----
        frH = FR_H_single(:, n);  % raw Hz
        fr_by_bin_H = accumarray(thetaH_bin_single(validH_single), frH(validH_single), [nbins,1], @mean, NaN);

        amp_obs = max(fr_by_bin_H) - min(fr_by_bin_H);
        MAXDIR_single_H(n) = angle(sum(fr_by_bin_H(:) .* exp(1i*bin_centers(:))));

        amp_null = nan(1, nshuff);
        thH = thetaH_bin_single;
        vmaskH = ~isnan(thH);
        for s = 1:nshuff
            shuff_bin = thH(vmaskH);
            shuff_bin = shuff_bin(randperm(numel(shuff_bin)));
            fr_by_bin = accumarray(shuff_bin, frH(vmaskH), [nbins,1], @mean, NaN);
            amp_null(s) = max(fr_by_bin) - min(fr_by_bin);
        end
        [z_obs, z_null] = zscore_against_null(amp_obs, amp_null);
        AMP_single_H(n)        = z_obs;
        AMP_null_single_H(n,:) = z_null;

        % Bootstrap (hand-only)
        FR_draws = nan(BOOT_nboot, nbins);
        z_draws  = nan(1, BOOT_nboot);
        zthr_draws = nan(1, BOOT_nboot);
        dir_draws= nan(1, BOOT_nboot);
        p_draws  = nan(1, BOOT_nboot);
        for bb = 1:BOOT_nboot
            ib = idxH_all(randi(numel(idxH_all), numel(idxH_all), 1));
            sh_bin = thetaH_bin_single(ib);
            frb    = frH(ib);
            fr_by_bin_b = accumarray(sh_bin, frb, [nbins,1], @mean, NaN);
            amp_obs_b = max(fr_by_bin_b) - min(fr_by_bin_b);

            amp_null_b = nan(1, nshuff);
            rb = numel(sh_bin);
            for s = 1:nshuff
                shuff = sh_bin(randperm(rb));
                fr_by_bin_s = accumarray(shuff, frb, [nbins,1], @mean, NaN);
                amp_null_b(s) = max(fr_by_bin_s) - min(fr_by_bin_s);
            end
            [z_obs_b, z_null_b] = zscore_against_null(amp_obs_b, amp_null_b);

            valid = isfinite(z_null_b);
            m = sum(valid);
            if m > 0 && isfinite(z_obs_b)
                p_draws(bb) = (1 + sum(z_null_b(valid) >= z_obs_b)) / (m + 1);
            else
                p_draws(bb) = NaN;
            end

            z_draws(bb)    = z_obs_b;
            zthr_draws(bb) = prctile(z_null_b, 95);
            dir_draws(bb)  = angle(sum(fr_by_bin_b(:) .* exp(1i*bin_centers(:))));
            FR_draws(bb, :) = fr_by_bin_b';
        end
        BOOT_AMP_single_H(n, :)    = z_draws;
        BOOT_THRESH_single_H(n, :) = zthr_draws;
        BOOT_MAXDIR_single_H(n, :) = dir_draws;
        BOOT_P_single_H(n, :)      = p_draws;

        P = prctile(FR_draws, [100*BOOT_alpha/2, 100*(1-BOOT_alpha/2)], 1);
        CI_FR_by_dir_single_H(n, :, :) = permute(P, [3 2 1]);

        % Optionally stash full FR-by-dir draws
        if BOOT_store_full_FR
            BOOT_FR_by_dir_single_H(n, :, :) = reshape(FR_draws.', [1, nbins, BOOT_nboot]);
        end

        % ----- Eye-only -----
        frE = FR_E_single(:, n);  % raw Hz
        fr_by_bin_E = accumarray(thetaE_bin_single(validE_single), frE(validE_single), [nbins,1], @mean, NaN);

        amp_obs = max(fr_by_bin_E) - min(fr_by_bin_E);
        MAXDIR_single_E(n) = angle(sum(fr_by_bin_E(:) .* exp(1i*bin_centers(:))));

        amp_null = nan(1, nshuff);
        thE = thetaE_bin_single;
        vmaskE = ~isnan(thE);
        for s = 1:nshuff
            shuff_bin = thE(vmaskE);
            shuff_bin = shuff_bin(randperm(numel(shuff_bin)));
            fr_by_bin = accumarray(shuff_bin, frE(vmaskE), [nbins,1], @mean, NaN);
            amp_null(s) = max(fr_by_bin) - min(fr_by_bin);
        end
        [z_obs, z_null] = zscore_against_null(amp_obs, amp_null);
        AMP_single_E(n)        = z_obs;
        AMP_null_single_E(n,:) = z_null;

        % Bootstrap (eye-only)
        FR_draws = nan(BOOT_nboot, nbins);
        z_draws  = nan(1, BOOT_nboot);
        zthr_draws = nan(1, BOOT_nboot);
        dir_draws= nan(1, BOOT_nboot);
        p_draws  = nan(1, BOOT_nboot);
        for bb = 1:BOOT_nboot
            ib = idxE_all(randi(numel(idxE_all), numel(idxE_all), 1));
            sh_bin = thetaE_bin_single(ib);
            frb    = frE(ib);
            fr_by_bin_b = accumarray(sh_bin, frb, [nbins,1], @mean, NaN);
            amp_obs_b = max(fr_by_bin_b) - min(fr_by_bin_b);

            amp_null_b = nan(1, nshuff);
            rb = numel(sh_bin);
            for s = 1:nshuff
                shuff = sh_bin(randperm(rb));
                fr_by_bin_s = accumarray(shuff, frb, [nbins,1], @mean, NaN);
                amp_null_b(s) = max(fr_by_bin_s) - min(fr_by_bin_s);
            end
            [z_obs_b, z_null_b] = zscore_against_null(amp_obs_b, amp_null_b);

            valid = isfinite(z_null_b);
            m = sum(valid);
            if m > 0 && isfinite(z_obs_b)
                p_draws(bb) = (1 + sum(z_null_b(valid) >= z_obs_b)) / (m + 1);
            else
                p_draws(bb) = NaN;
            end

            z_draws(bb)    = z_obs_b;
            zthr_draws(bb) = prctile(z_null_b, 95);
            dir_draws(bb)  = angle(sum(fr_by_bin_b(:) .* exp(1i*bin_centers(:))));
            FR_draws(bb, :) = fr_by_bin_b';
        end
        BOOT_AMP_single_E(n, :)    = z_draws;
        BOOT_THRESH_single_E(n, :) = zthr_draws;
        BOOT_MAXDIR_single_E(n, :) = dir_draws;
        BOOT_P_single_E(n, :)      = p_draws;

        P = prctile(FR_draws, [100*BOOT_alpha/2, 100*(1-BOOT_alpha/2)], 1);
        CI_FR_by_dir_single_E(n, :, :) = permute(P, [3 2 1]);

        if BOOT_store_full_FR
            BOOT_FR_by_dir_single_E(n, :, :) = reshape(FR_draws.', [1, nbins, BOOT_nboot]);
        end
    end

    % Robust PD for single (circular mean over boot-angles)
    MAXDIR_single.H = robust_pd_from_boot(BOOT_MAXDIR_single_H);
    MAXDIR_single.E = robust_pd_from_boot(BOOT_MAXDIR_single_E);

    % ============================
    % Legacy-shaped outputs
    % ============================
    AMP.H      = AMP_H;
    AMP.E      = AMP_E;
    AMP_null.H = AMP_null_H;
    AMP_null.E = AMP_null_E;

    AMP_single.H      = AMP_single_H;
    AMP_single.E      = AMP_single_E;
    AMP_null_single.H = AMP_null_single_H;
    AMP_null_single.E = AMP_null_single_E;

    % Median z across bootstraps (for plotting)
    z_med_multi.H   = nanmedian(BOOT_AMP_H,        2);
    z_med_multi.E   = nanmedian(BOOT_AMP_E,        2);
    z_med_single.H  = nanmedian(BOOT_AMP_single_H, 2);
    z_med_single.E  = nanmedian(BOOT_AMP_single_E, 2);

    % ============================
    % Thresholds, empirical p, MCP
    % ============================
    % Multi thresholds: median per-boot 95th percentile of z-null
    thr_med_multi.H = nanmedian(BOOT_THRESH_multi_H, 2);
    thr_med_multi.E = nanmedian(BOOT_THRESH_multi_E, 2);

    % Fallback to shuffle-only if needed
    thr0.H = prctile(AMP_null.H, 95, 2);
    thr0.E = prctile(AMP_null.E, 95, 2);
    thr_med_multi.H(~isfinite(thr_med_multi.H)) = thr0.H(~isfinite(thr_med_multi.H));
    thr_med_multi.E(~isfinite(thr_med_multi.E)) = thr0.E(~isfinite(thr_med_multi.E));

    % Median empirical p-values across bootstraps (multi)
    p_med_multi.H = nanmedian(BOOT_P_multi_H, 2);
    p_med_multi.E = nanmedian(BOOT_P_multi_E, 2);

    % Apply MCP within each family (multi)
    [sig_mcp.H, q_multi.H] = mcp_apply(p_med_multi.H, MCP.method, MCP.alpha);
    [sig_mcp.E, q_multi.E] = mcp_apply(p_med_multi.E, MCP.method, MCP.alpha);

    sig.H = sig_mcp.H;
    sig.E = sig_mcp.E;
    thresh.H = thr_med_multi.H;
    thresh.E = thr_med_multi.E;

    % Class labels (multi)
    class_labels = strings(num_neurons,1);
    class_labels(~sig.H & ~sig.E) = "None";
    class_labels( sig.H & ~sig.E) = "H";
    class_labels(~sig.H &  sig.E) = "E";
    class_labels( sig.H &  sig.E) = "Both";

    colors = zeros(num_neurons, 3);
    colors(class_labels == "H",    :) = repmat([0.8, 0.2, 0.2], sum(class_labels == "H"),    1);
    colors(class_labels == "E",    :) = repmat([0.2, 0.2, 0.8], sum(class_labels == "E"),    1);
    colors(class_labels == "Both", :) = repmat([0.2, 0.6, 0.2], sum(class_labels == "Both"), 1);
    colors(class_labels == "None", :) = repmat([0.5, 0.5, 0.5], sum(class_labels == "None"), 1);

    % Single thresholds: median per-boot 95th percentile of z-null
    thr_med_single.H = nanmedian(BOOT_THRESH_single_H, 2);
    thr_med_single.E = nanmedian(BOOT_THRESH_single_E, 2);

    % Fallback if needed
    thr0s.H = prctile(AMP_null_single.H, 95, 2);
    thr0s.E = prctile(AMP_null_single.E, 95, 2);
    thr_med_single.H(~isfinite(thr_med_single.H)) = thr0s.H(~isfinite(thr_med_single.H));
    thr_med_single.E(~isfinite(thr_med_single.E)) = thr0s.E(~isfinite(thr_med_single.E));

    % Median empirical p-values (single)
    p_med_single.H = nanmedian(BOOT_P_single_H, 2);
    p_med_single.E = nanmedian(BOOT_P_single_E, 2);

    % Apply MCP within each family (single)
    [sig_mcp_single.H, q_single.H] = mcp_apply(p_med_single.H, MCP.method, MCP.alpha);
    [sig_mcp_single.E, q_single.E] = mcp_apply(p_med_single.E, MCP.method, MCP.alpha);

    sig_single.H = sig_mcp_single.H;
    sig_single.E = sig_mcp_single.E;
    thresh_single.H = thr_med_single.H;
    thresh_single.E = thr_med_single.E;

    % Class labels (single)
    class_labels_single = strings(num_neurons,1);
    class_labels_single(~sig_single.H & ~sig_single.E) = "None";
    class_labels_single( sig_single.H & ~sig_single.E) = "H";
    class_labels_single(~sig_single.H &  sig_single.E) = "E";
    class_labels_single( sig_single.H &  sig_single.E) = "Both";

    colors_single = zeros(num_neurons, 3);
    colors_single(class_labels_single == "H",    :) = repmat([0.8, 0.2, 0.2], sum(class_labels_single == "H"),    1);
    colors_single(class_labels_single == "E",    :) = repmat([0.2, 0.2, 0.8], sum(class_labels_single == "E"),    1);
    colors_single(class_labels_single == "Both", :) = repmat([0.2, 0.6, 0.2], sum(class_labels_single == "Both"), 1);
    colors_single(class_labels_single == "None", :) = repmat([0.5, 0.5, 0.5], sum(class_labels_single == "None"), 1);

    % ============================
    % Compatibility alias
    % ============================
    AMP_multi = AMP;

    % ============================
    % Pack BOOT/CI/MCP structs
    % ============================
    BOOT.nboot = BOOT_nboot;
    BOOT.alpha = BOOT_alpha;
    BOOT.store_full_FR = BOOT_store_full_FR;
    BOOT.units = 'Hz';  % FR_by_dir draws are in raw Hz

    BOOT.AMP.H = BOOT_AMP_H;                   BOOT.AMP.E = BOOT_AMP_E;
    BOOT.AMP_single.H = BOOT_AMP_single_H;     BOOT.AMP_single.E = BOOT_AMP_single_E;

    BOOT.THRESH_multi.H = BOOT_THRESH_multi_H; BOOT.THRESH_multi.E = BOOT_THRESH_multi_E;
    BOOT.THRESH_single.H = BOOT_THRESH_single_H; BOOT.THRESH_single.E = BOOT_THRESH_single_E;

    BOOT.MAXDIR_multi.H = BOOT_MAXDIR_multi_H; BOOT.MAXDIR_multi.E = BOOT_MAXDIR_multi_E;
    BOOT.MAXDIR_single.H= BOOT_MAXDIR_single_H;BOOT.MAXDIR_single.E= BOOT_MAXDIR_single_E;

    BOOT.P_multi.H   = BOOT_P_multi_H;         BOOT.P_multi.E   = BOOT_P_multi_E;
    BOOT.P_single.H  = BOOT_P_single_H;        BOOT.P_single.E  = BOOT_P_single_E;

    CI.FR_by_dir.H = CI_FR_by_dir_H;                     CI.FR_by_dir.E = CI_FR_by_dir_E;
    CI.FR_by_dir_single.H = CI_FR_by_dir_single_H;       CI.FR_by_dir_single.E = CI_FR_by_dir_single_E;

    if BOOT_store_full_FR
        BOOT.FR_by_dir.H = BOOT_FR_by_dir_H;             
        BOOT.FR_by_dir.E = BOOT_FR_by_dir_E;
        BOOT.FR_by_dir_single.H = BOOT_FR_by_dir_single_H;
        BOOT.FR_by_dir_single.E = BOOT_FR_by_dir_single_E;
    end

    % Bootstrap probability of significance (per-draw thresholds)
    BOOT.p_sig.H = mean( BOOT_AMP_H > BOOT_THRESH_multi_H , 2, 'omitnan');
    BOOT.p_sig.E = mean( BOOT_AMP_E > BOOT_THRESH_multi_E , 2, 'omitnan');
    BOOT.p_sig_single.H = mean( BOOT_AMP_single_H > BOOT_THRESH_single_H , 2, 'omitnan');
    BOOT.p_sig_single.E = mean( BOOT_AMP_single_E > BOOT_THRESH_single_E , 2, 'omitnan');

    % MCP bookkeeping
    MCP.p_med.multi.H = p_med_multi.H;   MCP.p_med.multi.E = p_med_multi.E;
    MCP.p_med.single.H= p_med_single.H;  MCP.p_med.single.E= p_med_single.E;

    MCP.p_adj.multi.H = q_multi.H;       MCP.p_adj.multi.E = q_multi.E;
    MCP.p_adj.single.H= q_single.H;      MCP.p_adj.single.E= q_single.E;

    MCP.sig_only.multi.H = sig_mcp.H;    MCP.sig_only.multi.E = sig_mcp.E;
    MCP.sig_only.single.H= sig_mcp_single.H; MCP.sig_only.single.E = sig_mcp_single.E;

    % ----------------------------
    % Save
    % ----------------------------
    outdir = fullfile(savepath, 'models');
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    save(fullfile(outdir, sprintf('%s_FR_amplitude_selectivity.mat', area)), ...
        'AMP', 'AMP_null', 'thresh', 'sig', ...
        'class_labels', 'colors', ...
        'AMP_multi', 'AMP_single', ...
        'AMP_null_single', 'thresh_single', 'sig_single', ...
        'class_labels_single', 'colors_single', ...
        'MAXDIR_multi', 'MAXDIR_single', ...
        'BOOT', 'CI', 'MCP', ...
        'z_med_multi', 'z_med_single', ...
        '-v7.3');
end

% =========================
% Helpers
% =========================
function [z_obs, z_null] = zscore_against_null(obs, null_vec)
% z-score an observed scalar against its shuffle-null distribution.
% Returns NaNs if the null SD is 0 or not finite.
    mu = mean(null_vec, 'omitnan');
    sd = std(null_vec,  'omitnan');
    if ~isfinite(sd) || sd == 0
        z_obs  = NaN;
        z_null = NaN(size(null_vec));
        return;
    end
    z_obs  = (obs - mu) ./ sd;
    z_null = (null_vec - mu) ./ sd;
end

function th = robust_pd_from_boot(TH)
% Circular mean over bootstrap PD angles (rows = neurons).
    [N, ~] = size(TH);
    th = nan(N,1);
    for n = 1:N
        thn = TH(n,:);
        good = isfinite(thn);
        if any(good)
            C = mean(cos(thn(good)));
            S = mean(sin(thn(good)));
            th(n) = atan2(S, C);
        end
    end
end

% -------------------------------
% Multiple-comparisons helpers
% -------------------------------
function [sig, p_adj] = mcp_apply(p, method, alpha)
% p: Nx1 raw p-values (may contain NaNs).
% Returns logical significance after correction and adjusted p-values.
    p = p(:);
    p_adj = nan(size(p));
    sig   = false(size(p));

    finite_idx = isfinite(p);
    p_fin = p(finite_idx);
    m = numel(p_fin);
    if m == 0, return; end

    switch lower(method)
        case {'fdr-bh','fdr','bh'}
            [q, s] = fdr_bh_internal(p_fin, alpha);
            p_adj(finite_idx) = q;
            sig(finite_idx)   = s;
        case 'bonferroni'
            padj = min(1, p_fin * m);
            p_adj(finite_idx) = padj;
            sig(finite_idx)   = padj <= alpha;
        case 'none'
            p_adj(finite_idx) = p_fin;
            sig(finite_idx)   = p_fin <= alpha;
        otherwise
            error('Unknown MCP.method "%s".', method);
    end
end

function [q, sig] = fdr_bh_internal(p, alpha)
% Benjamini–Hochberg FDR control for finite p-values.
    [p_sorted, idx] = sort(p(:));
    m = numel(p_sorted);
    ranks = (1:m)';

    thresh_vec = (ranks / m) * alpha;
    pass = p_sorted <= thresh_vec;
    k = find(pass, 1, 'last');

    sig_sorted = false(m,1);
    if ~isempty(k), sig_sorted(1:k) = true; end

    q_sorted = m ./ ranks .* p_sorted;
    q_sorted = cummin(flipud(q_sorted));
    q_sorted = flipud(q_sorted);
    q_sorted = min(q_sorted, 1);

    q   = nan(size(p));
    sig = false(size(p));
    q(idx)   = q_sorted;
    sig(idx) = sig_sorted;
end
