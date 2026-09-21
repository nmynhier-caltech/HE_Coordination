%% added weighted baseline to decoding 
function nonparametric_MAP(area, savepath, varargin)
% Nonparametric MAP decoding with matched regularization for ME vs SE-additive.
% - ME: 2D von Mises×von Mises kernel regression from ME training trials.
% - SE: additive A(theta_H)+B(theta_E)+c from SE-only trials, with its OWN kappa (SE-only CV).
% - Baseline c for SE is a WEIGHTED average of A0 and B0 using peak-to-trough amplitudes.
% - Noise: BOTH ME and SE use heteroscedastic local predictive variance with the same shrinkage.
% - Likelihoods: Gaussian predictive; uniform prior over 6x6 angle grid.
%
% CHANGE:
%   ME train/test split is now fractional within each ME condition:
%       ~2/3 train, ~1/3 test.
%   This allows sessions with 6 trials per condition to contribute:
%       4 train, 2 test.
%
% Added analyses:
%   1) Accuracy vs #neurons: dynamic cumulative groups (up to 10; fewer if N<10).
%   2) Accuracy vs #training trials: percentage growth (independent pools).

% ---------- Parameters ----------
ip = inputParser;
addParameter(ip, 'TrainFraction', 2/3, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(ip, 'NumGroups', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
parse(ip, varargin{:});

train_fraction = ip.Results.TrainFraction;
num_groups     = ip.Results.NumGroups;

% ---------- Load ----------
D = load(processed_data_file(savepath, area));
data = D.data;

% Multi-effector trials (H+E)
idxME  = find(data.TPi2 ~= 7 & data.TPi1 ~= 7);
FR_me  = data.FR(idxME, :);   % [T_me x N]
TP1_me = data.TP1(idxME, :);  % hand positions
TP2_me = data.TP2(idxME, :);  % eye positions

% Single-effector trials (full pools; equalized used for BASE SE surface)
idxH  = find(data.TPi2 == 7 & data.TPi1 ~= 7);  % hand only
idxE  = find(data.TPi1 == 7 & data.TPi2 ~= 7);  % eye only
FR_h  = data.FR(idxH, :);
FR_e  = data.FR(idxE, :);
TP1_h = data.TP1(idxH, :);
TP2_e = data.TP2(idxE, :);

[T_me, N] = size(FR_me);

% ---------- Binning conventions ----------
nbins = 6;
theta_edges   = linspace(-pi, pi, nbins+1) + 0.1;    % offset to avoid boundary bleed
theta_centers = midpoints_circ(theta_edges);

% Angles
th_me_H = atan2(TP1_me(:,2), TP1_me(:,1));
th_me_E = atan2(TP2_me(:,2), TP2_me(:,1));
th_h_all= atan2(TP1_h(:,2),  TP1_h(:,1));
th_e_all= atan2(TP2_e(:,2),  TP2_e(:,1));

% Discretize ME angles to 6x6 cells
Hb_me = discretize(th_me_H, theta_edges);
Eb_me = discretize(th_me_E, theta_edges);

% ---------- UN-ZSCORE to RAW if preprocessing saved norms ----------
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

% ---------- Balanced ME train/test split ----------
% Uses approximately 2/3 train and 1/3 test within each ME condition.
% For 9 trials/cell: 6 train, 3 test.
% For 6 trials/cell: 4 train, 2 test.
train_mask_me = false(T_me,1);
test_mask_me  = false(T_me,1);

rng(1);
for ih = 1:nbins
    for ie = 1:nbins
        idx_cell = find(Hb_me == ih & Eb_me == ie);
        ncell = numel(idx_cell);

        if ncell < 2
            continue;
        end

        sel = idx_cell(randperm(ncell));

        Kt = floor(train_fraction * ncell);
        Kt = max(1, min(Kt, ncell-1)); % at least 1 train and 1 test

        train_mask_me(sel(1:Kt)) = true;
        test_mask_me(sel(Kt+1:end)) = true;
    end
end

% Extract ME splits in RAW
FR_me_train_raw = FR_me_raw(train_mask_me, :);
FR_me_test_raw  = FR_me_raw(test_mask_me,  :);
th_me_H_tr  = th_me_H(train_mask_me);
th_me_E_tr  = th_me_E(train_mask_me);
th_me_H_te  = th_me_H(test_mask_me);
th_me_E_te  = th_me_E(test_mask_me);
Hb_me_te    = Hb_me(test_mask_me);
Eb_me_te    = Eb_me(test_mask_me);
nTest = sum(test_mask_me);

fprintf('[%s] ME split: %d train trials, %d test trials, train_fraction=%.3f\n', ...
    area, sum(train_mask_me), nTest, train_fraction);

% ---------- TRAIN-ONLY STANDARDIZATION for ME ----------
mu_ME_train = mean(FR_me_train_raw, 1, 'omitnan');
sd_ME_train = std( FR_me_train_raw, 0, 1, 'omitnan');
sd_ME_train(~isfinite(sd_ME_train) | sd_ME_train==0) = 1;

FR_me_train = bsxfun(@rdivide, bsxfun(@minus, FR_me_train_raw, mu_ME_train), sd_ME_train);
FR_me_test  = bsxfun(@rdivide, bsxfun(@minus, FR_me_test_raw,  mu_ME_train), sd_ME_train);

% ---------- Build equalized SE pools ----------
Hb_h_all = discretize(th_h_all, theta_edges);
Eb_e_all = discretize(th_e_all, theta_edges);

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

FR_h_eq_raw = FR_h_raw(keep_h_eq, :);
FR_e_eq_raw = FR_e_raw(keep_e_eq, :);
th_h_eq = th_h_all(keep_h_eq);
th_e_eq = th_e_all(keep_e_eq);

% ---------- Candidate kappas and CV ----------
kappa_grid = [2, 4, 8, 12, 16, 24, 32];
Kfold = num_groups;
shrink_tau = 0.3;

% Precompute target grid
[Egrid, Hgrid] = meshgrid(theta_centers, theta_centers);
Hgrid = Hgrid(:);
Egrid = Egrid(:);
nGrid = numel(Hgrid);

% ---------- Per-neuron log-likelihood contributions ----------
llME_cell = cell(N,1);
llSE_cell = cell(N,1);
kappa_ME_vec = nan(N,1);
kappa_SE_vec = nan(N,1);

cME_train_raw = nan(N,1);
cSE_train_wt_raw = nan(N,1);
cSE_train_wt_z   = nan(N,1);

parfor n = 1:N
    % ----------------- ME CV -----------------
    r_me_tr = FR_me_train(:, n);
    thH_tr  = th_me_H_tr;
    thE_tr  = th_me_E_tr;

    if numel(r_me_tr) < Kfold
        kappa_star_ME = kappa_grid(min(3,length(kappa_grid)));
    else
        cv_idx = kfold_partition(numel(r_me_tr), Kfold, n);
        nll = zeros(numel(kappa_grid),1);

        for ik = 1:numel(kappa_grid)
            kap = kappa_grid(ik);
            fold_nll = 0;

            for k = 1:Kfold
                val_mask = (cv_idx == k);
                tr_mask  = ~val_mask;
                if ~any(val_mask) || ~any(tr_mask), continue; end

                mu_val = kernel_mean_2d(thH_tr(tr_mask), thE_tr(tr_mask), r_me_tr(tr_mask), ...
                                        thH_tr(val_mask),  thE_tr(val_mask), kap);
                [~, sigma_val] = kernel_local_sigma_2d(thH_tr(tr_mask), thE_tr(tr_mask), r_me_tr(tr_mask), ...
                                                       thH_tr(val_mask),  thE_tr(val_mask), kap, shrink_tau);
                rv = r_me_tr(val_mask);
                sigma_val = max(sigma_val, 1e-3);
                fold_nll = fold_nll + sum(0.5*log(2*pi*sigma_val.^2) + ...
                    0.5*((rv - mu_val).^2)./(sigma_val.^2));
            end

            nll(ik) = fold_nll;
        end

        [~,best] = min(nll);
        kappa_star_ME = kappa_grid(best);
    end

    % ----------------- SE CV -----------------
    pool_n = [FR_h_eq_raw(:,n); FR_e_eq_raw(:,n)];
    mu_SE_eq = mean(pool_n, 'omitnan');
    sd_SE_eq = std(pool_n, 0, 'omitnan');
    if ~isfinite(sd_SE_eq) || sd_SE_eq==0, sd_SE_eq = 1; end

    rH_SE = (FR_h_eq_raw(:,n) - mu_SE_eq) ./ sd_SE_eq;
    rE_SE = (FR_e_eq_raw(:,n) - mu_SE_eq) ./ sd_SE_eq;
    thH_SE = th_h_eq;
    thE_SE = th_e_eq;
    nH = numel(thH_SE);
    nE = numel(thE_SE);

    if (nH+nE) < Kfold
        kappa_star_SE = kappa_grid(min(3,length(kappa_grid)));
    else
        cv_idx_h = kfold_partition(max(nH,1), Kfold, 10000+n);
        cv_idx_e = kfold_partition(max(nE,1), Kfold, 20000+n);

        nll_se = zeros(numel(kappa_grid),1);

        for ik = 1:numel(kappa_grid)
            kap = kappa_grid(ik);
            fold_nll = 0;
            valid = false;

            for k = 1:Kfold
                if nH > 0
                    val_h = (cv_idx_h == k);
                    tr_h  = ~val_h;
                    if any(val_h) && any(tr_h)
                        muH = kernel_mean_1d(thH_SE(tr_h), rH_SE(tr_h), thH_SE(val_h), kap);
                        [~, sH] = kernel_local_sigma_1d(thH_SE(tr_h), rH_SE(tr_h), thH_SE(val_h), kap, shrink_tau);
                        sH = max(sH,1e-3);
                        rH = rH_SE(val_h);
                        fold_nll = fold_nll + sum(0.5*log(2*pi*sH.^2) + 0.5*((rH - muH).^2)./(sH.^2));
                        valid = true;
                    end
                end

                if nE > 0
                    val_e = (cv_idx_e == k);
                    tr_e  = ~val_e;
                    if any(val_e) && any(tr_e)
                        muE = kernel_mean_1d(thE_SE(tr_e), rE_SE(tr_e), thE_SE(val_e), kap);
                        [~, sE] = kernel_local_sigma_1d(thE_SE(tr_e), rE_SE(tr_e), thE_SE(val_e), kap, shrink_tau);
                        sE = max(sE,1e-3);
                        rE = rE_SE(val_e);
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

    % ---------- Final ME model ----------
    mu_ME_grid   = kernel_mean_2d(thH_tr, thE_tr, r_me_tr, Hgrid, Egrid, kappa_star_ME);
    [~, sg_ME_g] = kernel_local_sigma_2d(thH_tr, thE_tr, r_me_tr, Hgrid, Egrid, kappa_star_ME, shrink_tau);
    sg_ME_g = max(sg_ME_g, 1e-3);

    % ---------- Final SE-additive model ----------
    A_H = kernel_mean_1d(thH_SE, rH_SE, theta_centers, kappa_star_SE);
    B_E = kernel_mean_1d(thE_SE, rE_SE, theta_centers, kappa_star_SE);
    A0  = mean(A_H, 'omitnan');
    B0  = mean(B_E, 'omitnan');
    ampH = ptp_amplitude(A_H);
    ampE = ptp_amplitude(B_E);
    c_wt_z = weighted_baseline(A0, B0, ampH, ampE);

    mu_SE_grid = (A_H(:) - A0) + (B_E(:).' - B0) + c_wt_z;
    mu_SE_grid = mu_SE_grid(:);

    [~, sH_cent] = kernel_local_sigma_1d(thH_SE, rH_SE, theta_centers, kappa_star_SE, shrink_tau);
    [~, sE_cent] = kernel_local_sigma_1d(thE_SE, rE_SE, theta_centers, kappa_star_SE, shrink_tau);
    sH_cent = max(sH_cent,1e-3);
    sE_cent = max(sE_cent,1e-3);

    [sEgrid, sHgrid] = meshgrid(sE_cent, sH_cent);
    sg_SE_g = sqrt(sHgrid.^2 + sEgrid.^2);
    sg_SE_g = max(sg_SE_g(:),1e-3);
    sg_SE_g = reshape(sg_SE_g, [], 1);

    % ---------- Baselines ----------
    cME_train_raw(n) = mean(FR_me_train_raw(:,n), 'omitnan');

    A_H_raw = kernel_mean_1d(th_h_eq, FR_h_eq_raw(:,n), theta_centers, kappa_star_SE);
    B_E_raw = kernel_mean_1d(th_e_eq, FR_e_eq_raw(:,n), theta_centers, kappa_star_SE);
    A0_raw  = mean(A_H_raw,'omitnan');
    B0_raw  = mean(B_E_raw,'omitnan');
    ampH_raw = ptp_amplitude(A_H_raw);
    ampE_raw = ptp_amplitude(B_E_raw);
    c_wt_raw = weighted_baseline(A0_raw, B0_raw, ampH_raw, ampE_raw);

    cSE_train_wt_raw(n) = c_wt_raw;
    cSE_train_wt_z(n)   = c_wt_z;

    % ---------- Log-likelihoods on ME test set ----------
    r_test = FR_me_test(:, n);
    Rmat = repmat(r_test, 1, nGrid);

    ll_me = -0.5*log(2*pi) - log(repmat(sg_ME_g.', nTest, 1)) ...
            - 0.5*((Rmat - repmat(mu_ME_grid.', nTest, 1))./repmat(sg_ME_g.', nTest, 1)).^2;

    ll_se = -0.5*log(2*pi) - log(repmat(sg_SE_g.', nTest, 1)) ...
            - 0.5*((Rmat - repmat(mu_SE_grid.', nTest, 1))./repmat(sg_SE_g.', nTest, 1)).^2;

    llME_cell{n} = ll_me;
    llSE_cell{n} = ll_se;
    kappa_ME_vec(n) = kappa_star_ME;
    kappa_SE_vec(n) = kappa_star_SE;
end

% Sum per-neuron contributions
loglik_ME = zeros(nTest, nGrid);
loglik_SE = zeros(nTest, nGrid);
for n = 1:N
    if ~isempty(llME_cell{n}), loglik_ME = loglik_ME + llME_cell{n}; end
    if ~isempty(llSE_cell{n}), loglik_SE = loglik_SE + llSE_cell{n}; end
end

% ---------- Decode ----------
[acc_ME_H0, acc_ME_E0, C_ME_H_row, C_ME_E_row] = decode_from_LL(loglik_ME, Hb_me_te, Eb_me_te, nbins);
[acc_SE_H0, acc_SE_E0, C_SE_H_row, C_SE_E_row] = decode_from_LL(loglik_SE, Hb_me_te, Eb_me_te, nbins);

[~, idx_hat_ME] = max(loglik_ME, [], 2);
[~, idx_hat_SE] = max(loglik_SE, [], 2);
true_flat = sub2ind([nbins nbins], Hb_me_te(:), Eb_me_te(:));

% ---------- Minimal Visualization ----------
plotdir = fullfile(savepath, 'plots');
decdir  = fullfile(plotdir, 'decoding');
if ~exist(decdir, 'dir'), mkdir(decdir); end

fHS = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 8]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile; imagesc(C_ME_H_row,[0 1]); axis square; colorbar;
title(sprintf('%s: H Confusion — ME', area)); xlabel('Pred H'); ylabel('True H'); set(gca,'FontSize',13);

nexttile; imagesc(C_SE_H_row,[0 1]); axis square; colorbar;
title(sprintf('%s: H Confusion — SE-add', area)); xlabel('Pred H'); ylabel('True H'); set(gca,'FontSize',13);

nexttile; imagesc(C_ME_E_row,[0 1]); axis square; colorbar;
title(sprintf('%s: E Confusion — ME', area)); xlabel('Pred E'); ylabel('True E'); set(gca,'FontSize',13);

nexttile; imagesc(C_SE_E_row,[0 1]); axis square; colorbar;
title(sprintf('%s: E Confusion — SE-add', area)); xlabel('Pred E'); ylabel('True E'); set(gca,'FontSize',13);

print(fHS, fullfile(decdir, sprintf('%s_MAP_Confusions_H_E_Smoothed.pdf', area)), '-dpdf', '-painters');
close(fHS);

% =====================================================================
% ANALYSIS 1: Accuracy vs #NEURONS
% =====================================================================
max_groups = min(num_groups, max(1,N));
neu_levels = unique(max(1, round(N * (1:max_groups)/max_groups)), 'stable');
rng(3);
perm_neu = randperm(max(N,1));

acc_ME_H = nan(numel(neu_levels),1);
acc_ME_E = nan(numel(neu_levels),1);
acc_SE_H = nan(numel(neu_levels),1);
acc_SE_E = nan(numel(neu_levels),1);

for ii = 1:numel(neu_levels)
    k = neu_levels(ii);
    pick = perm_neu(1:min(k,N));

    LL_me = zeros(nTest, nGrid);
    LL_se = zeros(nTest, nGrid);
    for nn = pick
        if ~isempty(llME_cell{nn}), LL_me = LL_me + llME_cell{nn}; end
        if ~isempty(llSE_cell{nn}), LL_se = LL_se + llSE_cell{nn}; end
    end

    [acc_ME_H(ii), acc_ME_E(ii)] = accuracy_from_LL(LL_me, Hb_me_te, Eb_me_te, nbins);
    [acc_SE_H(ii), acc_SE_E(ii)] = accuracy_from_LL(LL_se, Hb_me_te, Eb_me_te, nbins);
end

fN = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 6]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on; box on; set(gca,'FontSize',13);
h1 = plot(neu_levels, acc_ME_H, 'ro-','LineWidth',1.5);
h2 = plot(neu_levels, acc_SE_H, 'bs-','LineWidth',1.5);
[xfit,yfit] = fit_saturating(neu_levels, acc_ME_H); h3 = []; if ~isempty(xfit), h3 = plot(xfit,yfit,'r.','LineWidth',1.2); end
[xfit,yfit] = fit_saturating(neu_levels, acc_SE_H); h4 = []; if ~isempty(xfit), h4 = plot(xfit,yfit,'b.','LineWidth',1.2); end
xlabel('# Neurons'); ylabel('Accuracy (H)'); title(sprintf('%s: Accuracy vs #Neurons', area));
legend_nonempty([h1 h2 h3 h4], {'ME-H','SE-H','ME-H fit','SE-H fit'}, 'southeast');

nexttile; hold on; grid on; box on; set(gca,'FontSize',13);
plot(neu_levels, acc_ME_E, 'ro-','LineWidth',1.5);
plot(neu_levels, acc_SE_E, 'bs-','LineWidth',1.5);
[xfit,yfit] = fit_saturating(neu_levels, acc_ME_E); if ~isempty(xfit), plot(xfit,yfit,'r.','LineWidth',1.2); end
[xfit,yfit] = fit_saturating(neu_levels, acc_SE_E); if ~isempty(xfit), plot(xfit,yfit,'b.','LineWidth',1.2); end
xlabel('# Neurons'); ylabel('Accuracy (E)');

print(fN, fullfile(decdir, sprintf('%s_Accuracy_vs_Neurons.pdf', area)), '-dpdf','-painters');
close(fN);

% =====================================================================
% ANALYSIS 2: Accuracy vs #TRAINING TRIALS
% =====================================================================
lvls = min(num_groups, 6);
pct_levels = linspace(max(1/lvls, 0.05), 1.0, lvls);

acc_ME_H_tr = nan(numel(pct_levels),1);
acc_ME_E_tr = nan(numel(pct_levels),1);
acc_SE_H_tr = nan(numel(pct_levels),1);
acc_SE_E_tr = nan(numel(pct_levels),1);
x_me_trials = nan(numel(pct_levels),1);
x_se_trials = nan(numel(pct_levels),1);

idx_cell_me = cell(nbins,nbins);
for ih = 1:nbins
    for ie = 1:nbins
        idx_cell_me{ih,ie} = find(Hb_me(train_mask_me)==ih & Eb_me(train_mask_me)==ie);
    end
end

idx_h_bins_all = cell(nbins,1);
idx_e_bins_all = cell(nbins,1);
for ih = 1:nbins, idx_h_bins_all{ih} = find(Hb_h_all==ih); end
for ie = 1:nbins, idx_e_bins_all{ie} = find(Eb_e_all==ie); end

rng(4);

for jj = 1:numel(pct_levels)
    p = pct_levels(jj);

    sel_me = false(size(FR_me_train,1),1);
    me_count = 0;
    for ih = 1:nbins
        for ie = 1:nbins
            idxc = idx_cell_me{ih,ie};
            if isempty(idxc), continue; end
            k_here = numel(idxc);
            k_use = max(1, round(p * k_here));
            take = idxc(randperm(k_here, k_use));
            sel_me(take) = true;
            me_count = me_count + numel(take);
        end
    end

    thH_tr_p = th_me_H_tr(sel_me);
    thE_tr_p = th_me_E_tr(sel_me);
    FR_me_tr_p = FR_me_train(sel_me, :);
    x_me_trials(jj) = me_count;

    sel_h = false(size(FR_h_raw,1),1);
    sel_e = false(size(FR_e_raw,1),1);
    se_count = 0;

    for ih = 1:nbins
        idxc = idx_h_bins_all{ih};
        if isempty(idxc), continue; end
        k_here = numel(idxc);
        k_use  = max(1, round(p * k_here));
        take = idxc(randperm(k_here, k_use));
        sel_h(take) = true;
        se_count = se_count + numel(take);
    end

    for ie = 1:nbins
        idxc = idx_e_bins_all{ie};
        if isempty(idxc), continue; end
        k_here = numel(idxc);
        k_use  = max(1, round(p * k_here));
        take = idxc(randperm(k_here, k_use));
        sel_e(take) = true;
        se_count = se_count + numel(take);
    end

    th_h_p = th_h_all(sel_h);
    FR_h_p_raw = FR_h_raw(sel_h, :);
    th_e_p = th_e_all(sel_e);
    FR_e_p_raw = FR_e_raw(sel_e, :);
    x_se_trials(jj) = se_count;

    llME_p = zeros(nTest, nGrid);
    llSE_p = zeros(nTest, nGrid);

    parfor n = 1:N
        kapME = kappa_ME_vec(n); if ~isfinite(kapME) || kapME<=0, kapME = 8; end
        kapSE = kappa_SE_vec(n); if ~isfinite(kapSE) || kapSE<=0, kapSE = 8; end

        if ~isempty(FR_me_tr_p)
            mu_ME = kernel_mean_2d(thH_tr_p, thE_tr_p, FR_me_tr_p(:,n), Hgrid, Egrid, kapME);
            [~, sg_ME] = kernel_local_sigma_2d(thH_tr_p, thE_tr_p, FR_me_tr_p(:,n), Hgrid, Egrid, kapME, shrink_tau);
            sg_ME = max(sg_ME, 1e-3);
        else
            mu_ME = mean(FR_me_train(:,n),'omitnan')*ones(nGrid,1);
            sg_ME = std(FR_me_train(:,n),[],'omitnan')*ones(nGrid,1);
            sg_ME = max(sg_ME,1e-3);
        end

        pool_p_n = [FR_h_p_raw(:,n); FR_e_p_raw(:,n)];
        mu_SE_p = mean(pool_p_n, 'omitnan');
        sd_SE_p = std(pool_p_n, 0, 'omitnan');
        if ~isfinite(sd_SE_p) || sd_SE_p==0, sd_SE_p = 1; end

        rH_p = (FR_h_p_raw(:,n) - mu_SE_p) ./ sd_SE_p;
        rE_p = (FR_e_p_raw(:,n) - mu_SE_p) ./ sd_SE_p;

        A_H = kernel_mean_1d(th_h_p, rH_p, theta_centers, kapSE);
        B_E = kernel_mean_1d(th_e_p, rE_p, theta_centers, kapSE);
        A0 = mean(A_H,'omitnan');
        B0 = mean(B_E,'omitnan');
        ampH = ptp_amplitude(A_H);
        ampE = ptp_amplitude(B_E);
        c_wt = weighted_baseline(A0, B0, ampH, ampE);

        mu_SE = (A_H(:)-A0) + (B_E(:).'-B0) + c_wt;
        mu_SE = mu_SE(:);

        [~, sH_p] = kernel_local_sigma_1d(th_h_p, rH_p, theta_centers, kapSE, shrink_tau);
        [~, sE_p] = kernel_local_sigma_1d(th_e_p, rE_p, theta_centers, kapSE, shrink_tau);
        sH_p = max(sH_p,1e-3);
        sE_p = max(sE_p,1e-3);

        [sEgrid, sHgrid] = meshgrid(sE_p, sH_p);
        sg_SE = sqrt(sHgrid.^2 + sEgrid.^2);
        sg_SE = reshape(sg_SE, [], 1);

        r_test = FR_me_test(:, n);
        Rmat = repmat(r_test, 1, nGrid);

        llME_n = -0.5*log(2*pi) - log(repmat(sg_ME.', nTest, 1)) ...
                 - 0.5*((Rmat - repmat(mu_ME.', nTest, 1))./repmat(sg_ME.', nTest, 1)).^2;

        llSE_n = -0.5*log(2*pi) - log(repmat(sg_SE.', nTest, 1)) ...
                 - 0.5*((Rmat - repmat(mu_SE.', nTest, 1))./repmat(sg_SE.', nTest, 1)).^2;

        llME_p = llME_p + llME_n;
        llSE_p = llSE_p + llSE_n;
    end

    [acc_ME_H_tr(jj), acc_ME_E_tr(jj)] = accuracy_from_LL(llME_p, Hb_me_te, Eb_me_te, nbins);
    [acc_SE_H_tr(jj), acc_SE_E_tr(jj)] = accuracy_from_LL(llSE_p, Hb_me_te, Eb_me_te, nbins);
end

fT = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 6]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on; box on; set(gca,'FontSize',13);
h1 = plot(x_me_trials, acc_ME_H_tr,'ro-','LineWidth',1.5);
h2 = plot(x_se_trials, acc_SE_H_tr,'bs-','LineWidth',1.5);
[xf,yf] = fit_saturating(x_me_trials, acc_ME_H_tr); h3 = []; if ~isempty(xf), h3 = plot(xf,yf,'r.','LineWidth',1.2); end
[xf,yf] = fit_saturating(x_se_trials, acc_SE_H_tr); h4 = []; if ~isempty(xf), h4 = plot(xf,yf,'b.','LineWidth',1.2); end
xlabel('Training trials'); ylabel('Accuracy (H)');
title(sprintf('%s: Accuracy vs Training Trials', area));
legend_nonempty([h1 h2 h3 h4], {'ME-H','SE-H','ME-H fit','SE-H fit'}, 'southeast');

nexttile; hold on; grid on; box on; set(gca,'FontSize',13);
plot(x_me_trials, acc_ME_E_tr,'ro-','LineWidth',1.5);
plot(x_se_trials, acc_SE_E_tr,'bs-','LineWidth',1.5);
[xf,yf] = fit_saturating(x_me_trials, acc_ME_E_tr); if ~isempty(xf), plot(xf,yf,'r.','LineWidth',1.2); end
[xf,yf] = fit_saturating(x_se_trials, acc_SE_E_tr); if ~isempty(xf), plot(xf,yf,'b.','LineWidth',1.2); end
xlabel('Training trials'); ylabel('Accuracy (E)');

print(fT, fullfile(decdir, sprintf('%s_Accuracy_vs_TrainTrials.pdf', area)), '-dpdf','-painters');
close(fT);

% ---------- Save ----------
outmat = fullfile(savepath, 'models', sprintf('%s_MAP_Smoothed.mat', area));
Results = struct();

Results.params.train_fraction = train_fraction;
Results.params.split_mode = 'within_cell_fractional_train_test';
Results.params.num_groups = num_groups;
Results.params.nTest = nTest;
Results.params.nTrain = sum(train_mask_me);

Results.kappa_grid   = kappa_grid;
Results.shrink_tau   = shrink_tau;
Results.theta_centers= theta_centers;

Results.kappa_ME = kappa_ME_vec;
Results.kappa_SE = kappa_SE_vec;

Results.C_ME_H_row = C_ME_H_row;
Results.C_SE_H_row = C_SE_H_row;
Results.C_ME_E_row = C_ME_E_row;
Results.C_SE_E_row = C_SE_E_row;

Results.acc_base = struct('ME_H',acc_ME_H0,'ME_E',acc_ME_E0,'SE_H',acc_SE_H0,'SE_E',acc_SE_E0);

Results.acc_vs_neurons = struct( ...
    'neu_levels', neu_levels, ...
    'ME_H', acc_ME_H, ...
    'ME_E', acc_ME_E, ...
    'SE_H', acc_SE_H, ...
    'SE_E', acc_SE_E);

Results.acc_vs_trials = struct( ...
    'x_me_trials', x_me_trials, ...
    'x_se_trials', x_se_trials, ...
    'ME_H', acc_ME_H_tr, ...
    'ME_E', acc_ME_E_tr, ...
    'SE_H', acc_SE_H_tr, ...
    'SE_E', acc_SE_E_tr);

Results.baseline = struct( ...
    'ME_raw_train', cME_train_raw, ...
    'SE_weighted_raw_train', cSE_train_wt_raw, ...
    'SE_weighted_z_train', cSE_train_wt_z);

% Joint condition accuracy fields
Results.true_flat  = true_flat;
Results.idx_hat_ME = idx_hat_ME;
Results.idx_hat_SE = idx_hat_SE;
Results.Hb_me_te   = Hb_me_te(:);
Results.Eb_me_te   = Eb_me_te(:);
Results.nTest      = nTest;

save(outmat, 'Results', '-v7.3');
fprintf('Saved: %s\n', outmat);

end

% ================================
% ---- Helper functions below ----
% ================================

function R = rownorm(C)
    s = sum(C,2);
    s(s==0)=1;
    R = C ./ s;
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
        sigma_loc(:) = 1;
        return;
    end

    vg = var(r_tr, 'omitnan');
    if ~isfinite(vg) || vg<=0, vg = 1; end

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
        sigma_loc(i) = sqrt(tau*vg + (1-tau)*max(vloc,1e-6));
    end
end

function [mu_loc, sigma_loc] = kernel_local_sigma_1d(theta_tr, r_tr, theta_eval, kappa, tau)
    mu_loc    = nan(size(theta_eval));
    sigma_loc = nan(size(theta_eval));

    if isempty(r_tr)
        sigma_loc(:) = 1;
        return;
    end

    vg = var(r_tr, 'omitnan');
    if ~isfinite(vg) || vg<=0, vg = 1; end

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
        sigma_loc(i) = sqrt(tau*vg + (1-tau)*max(vloc,1e-6));
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
        if nk > 0
            idx(p(b+1:b+nk)) = k;
        end
        b = b + nk;
    end
end

function [accH, accE, C_H_row, C_E_row] = decode_from_LL(LL, Hb_true, Eb_true, nbins)
    if isempty(LL) || isempty(Hb_true) || isempty(Eb_true)
        accH = NaN;
        accE = NaN;
        C_H_row = zeros(nbins);
        C_E_row = zeros(nbins);
        return;
    end

    [~, idx_hat] = max(LL, [], 2);
    [Hhat, Ehat] = ind2sub([nbins nbins], idx_hat);

    accH = mean(Hhat(:) == Hb_true(:), 'omitnan');
    accE = mean(Ehat(:) == Eb_true(:), 'omitnan');

    C_H = accumarray([Hb_true, Hhat], 1, [nbins nbins], @sum, 0);
    C_E = accumarray([Eb_true, Ehat], 1, [nbins nbins], @sum, 0);

    C_H_row = rownorm(C_H);
    C_E_row = rownorm(C_E);
end

function [accH, accE] = accuracy_from_LL(LL, Hb_true, Eb_true, nbins)
    if isempty(LL) || isempty(Hb_true) || isempty(Eb_true)
        accH = NaN;
        accE = NaN;
        return;
    end

    [~, idx_hat] = max(LL, [], 2);
    [Hhat, Ehat] = ind2sub([nbins nbins], idx_hat);

    accH = mean(Hhat(:) == Hb_true(:), 'omitnan');
    accE = mean(Ehat(:) == Eb_true(:), 'omitnan');
end

function [xfit, yfit] = fit_saturating(x, y)
    x = x(:);
    y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    if isempty(x)
        xfit = [];
        yfit = [];
        return;
    elseif numel(x) == 1
        xfit = x;
        yfit = y;
        return;
    elseif numel(x) == 2
        xfit = linspace(x(1), x(2), 100);
        yfit = interp1(x, y, xfit, 'linear', 'extrap');
        return;
    end

    y0 = y(1);
    L0 = max(y);
    k0 = 1./max(x);

    obj = @(p) sum(((p(1) - (p(1) - p(2))*exp(-p(3)*x)) - y).^2, 'omitnan');
    p = fminsearch(obj, [L0, y0, k0], optimset('Display','off'));

    xfit = linspace(min(x), max(x), 200);
    yfit = p(1) - (p(1) - p(2))*exp(-p(3)*xfit);
end

function amp = ptp_amplitude(x)
    if all(~isfinite(x))
        amp = NaN;
    else
        amp = max(x,[],'omitnan') - min(x,[],'omitnan');
    end
end

function c = weighted_baseline(A0, B0, ampH, ampE)
    wH = max(ampH, 0);
    wE = max(ampE, 0);

    if isfinite(wH) && wH>0 && isfinite(wE) && wE>0
        c = (wH*A0 + wE*B0) / (wH + wE);
    elseif isfinite(wH) && wH>0
        c = A0;
    elseif isfinite(wE) && wE>0
        c = B0;
    else
        c = mean([A0, B0], 'omitnan');
    end
end

function legend_nonempty(handles, labels, loc)
    good = arrayfun(@(h) ~isempty(h) && isgraphics(h), handles);
    handles = handles(good);
    labels = labels(good);

    if ~isempty(handles)
        legend(handles, labels, 'Location', loc);
    end
end
