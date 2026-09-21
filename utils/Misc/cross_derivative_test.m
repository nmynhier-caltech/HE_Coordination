% single neuron test statistic with raw RMS
function cross_derivative_test(area, savepath)
% Per-neuron Hessian derivative tests on ME trials only.
%
% Grid: nbins x nbins in (thetaH, thetaE). Build mean-FR map R.
% Apply 2-pass NaN-aware periodic separable smoothing (default).
% Compute periodic central finite-difference second derivatives:
%   HH = d2/dH2, HE = d2/(dH dE), EE = d2/dE2
%
% Per-neuron statistic (RAW, smoothed):
%   raw_rms_* = sqrt(mean(F*(m).^2)) over a validity stencil mask m.
%
% Null:
%   Shuffle FR across trials (angles fixed), rebuild map, smooth, derivatives,
%   compute same RMS => null_vec_* of length nPerm.
%
% Per-neuron p-value (right-tailed on RAW RMS):
%   p = (#{null >= real} + 1) / (nPerm + 1)
%
% Saves: models/<area>_hessian_results.mat with struct array "results":
%   neuron_id, range_FR
%   raw_rms_HH, raw_rms_HE, raw_rms_EE
%   p_HH, p_HE, p_EE
%   null_vec_HH, null_vec_HE, null_vec_EE
% Also saves the legacy HE-only compatibility output used by
% plot_relative_position_summary:
%   models/<area>_cross_deriv_results.mat (struct array "cross_results")

    % ----------------------------
    % Load data
    % ----------------------------
    S = load(processed_data_file(savepath, area));
    data = S.data;

    % ME trial indices and variables
    multiTP = find(data.TPi2 ~= 7 & data.TPi1 ~= 7);
    thetaH  = atan2(data.TP1(multiTP,2), data.TP1(multiTP,1));
    thetaE  = atan2(data.TP2(multiTP,2), data.TP2(multiTP,1));
    FR      = data.FR(multiTP, :);   % [trials x neurons], raw Hz
    [~, N]  = size(FR);

    % ----------------------------
    % Binning setup
    % ----------------------------
    nbins       = 6;
    edges_theta = linspace(-pi, pi, nbins + 1) + 0.1;
    dth         = 2*pi/nbins;  % rad/bin

    thetaH_bin = discretize(thetaH, edges_theta);
    thetaE_bin = discretize(thetaE, edges_theta);

    valid_all = isfinite(thetaH_bin) & isfinite(thetaE_bin);
    h_b = thetaH_bin(valid_all);
    e_b = thetaE_bin(valid_all);

    % ----------------------------
    % Params
    % ----------------------------
    nPerm     = 1000;
    min_count = 1;

    doSmoothing  = true;
    smoothK      = [0.25 0.50 0.25];
    smoothPasses = 2;

    % ----------------------------
    % Results structs
    % ----------------------------
    template = struct( ...
        'neuron_id', 0, ...
        'range_FR', NaN, ...
        'raw_rms_HH', NaN, 'raw_rms_HE', NaN, 'raw_rms_EE', NaN, ...
        'p_HH', NaN, 'p_HE', NaN, 'p_EE', NaN, ...
        'null_vec_HH', [], 'null_vec_HE', [], 'null_vec_EE', [] ...
    );
    results = repmat(template, N, 1);

    cross_template = struct( ...
        'neuron_id', 0, ...
        'rms_dydx', NaN, ...
        'raw_rms_dydx', NaN, ...
        'rms_dydx_z', NaN, ...
        'range_FR', NaN, ...
        'p_value', NaN, ...
        'null_mean', NaN, ...
        'null_sd', NaN, ...
        'null_rms', NaN, ...
        'null_rms_z', NaN ...
    );
    cross_results = repmat(cross_template, N, 1);

    % ----------------------------
    % Main loop
    % ----------------------------
    for n = 1:N
        rng(n, 'twister');
        FR_n = FR(valid_all, n);

        % Build mean-FR map
        R = nan(nbins, nbins);
        for i = 1:nbins
            for j = 1:nbins
                idx = (h_b == i) & (e_b == j);
                if nnz(idx) >= min_count
                    R(i,j) = mean(FR_n(idx));
                end
            end
        end

        % Validity mask for second-derivative stencils (8-neighborhood)
        isfin = isfinite(R);
        f = @(di,dj) circshift(isfin, [di, dj]);
        m = isfin & f(1,0) & f(-1,0) & f(0,1) & f(0,-1) & ...
                   f(1,1) & f(1,-1) & f(-1,1) & f(-1,-1);

        if ~any(m(:))
            continue
        end

        % Range from UNSMOOTHED map
        results(n).range_FR = max(R(:),[],'omitnan') - min(R(:),[],'omitnan');
        results(n).neuron_id = n;

        % Smooth real map
        Ruse = R;
        if doSmoothing
            for sp = 1:smoothPasses
                Ruse = smooth_periodic_2d(Ruse, smoothK);
            end
        end

        % Derivatives (periodic central)
        FHH = (circshift(Ruse,[-1,0]) - 2*Ruse + circshift(Ruse,[1,0])) / (dth*dth);
        FEE = (circshift(Ruse,[0,-1]) - 2*Ruse + circshift(Ruse,[0,1])) / (dth*dth);
        FHE = (circshift(Ruse,[-1,-1]) - circshift(Ruse,[-1,1]) ...
             - circshift(Ruse,[ 1,-1]) + circshift(Ruse,[ 1, 1])) / (4*dth*dth);

        real_HH = sqrt(mean(FHH(m).^2, 'omitnan'));
        real_HE = sqrt(mean(FHE(m).^2, 'omitnan'));
        real_EE = sqrt(mean(FEE(m).^2, 'omitnan'));

        results(n).raw_rms_HH = real_HH;
        results(n).raw_rms_HE = real_HE;
        results(n).raw_rms_EE = real_EE;

        % Null vectors
        null_HH = nan(1, nPerm);
        null_HE = nan(1, nPerm);
        null_EE = nan(1, nPerm);

        for p = 1:nPerm
            shuff = FR_n(randperm(numel(FR_n)));

            Rn = nan(nbins, nbins);
            for i = 1:nbins
                for j = 1:nbins
                    idx = (h_b == i) & (e_b == j);
                    if nnz(idx) >= min_count
                        Rn(i,j) = mean(shuff(idx));
                    end
                end
            end

            Rn_use = Rn;
            if doSmoothing
                for sp = 1:smoothPasses
                    Rn_use = smooth_periodic_2d(Rn_use, smoothK);
                end
            end

            FHH_n = (circshift(Rn_use,[-1,0]) - 2*Rn_use + circshift(Rn_use,[1,0])) / (dth*dth);
            FEE_n = (circshift(Rn_use,[0,-1]) - 2*Rn_use + circshift(Rn_use,[0,1])) / (dth*dth);
            FHE_n = (circshift(Rn_use,[-1,-1]) - circshift(Rn_use,[-1,1]) ...
                   - circshift(Rn_use,[ 1,-1]) + circshift(Rn_use,[ 1, 1])) / (4*dth*dth);

            null_HH(p) = sqrt(mean(FHH_n(m).^2, 'omitnan'));
            null_HE(p) = sqrt(mean(FHE_n(m).^2, 'omitnan'));
            null_EE(p) = sqrt(mean(FEE_n(m).^2, 'omitnan'));
        end

        results(n).null_vec_HH = null_HH;
        results(n).null_vec_HE = null_HE;
        results(n).null_vec_EE = null_EE;

        % Empirical right-tailed p-values on RAW RMS
        results(n).p_HH = (sum(null_HH >= real_HH) + 1) / (nPerm + 1);
        results(n).p_HE = (sum(null_HE >= real_HE) + 1) / (nPerm + 1);
        results(n).p_EE = (sum(null_EE >= real_EE) + 1) / (nPerm + 1);

        % Preserve the legacy HE-only result consumed by the relative-position
        % summary. The two draws reproduce the former HH-then-HE sampling
        % order; the first draw is intentionally discarded.
        null_mean_HE = mean(null_HE, 'omitnan');
        null_sd_HE = std(null_HE, 'omitnan');
        randi(nPerm);
        representative_HE = null_HE(randi(nPerm));

        cross_results(n).neuron_id = n;
        cross_results(n).raw_rms_dydx = real_HE;
        cross_results(n).rms_dydx = real_HE - null_mean_HE;
        cross_results(n).range_FR = results(n).range_FR;
        cross_results(n).p_value = results(n).p_HE;
        cross_results(n).null_mean = null_mean_HE;
        cross_results(n).null_sd = null_sd_HE;
        cross_results(n).null_rms = representative_HE - null_mean_HE;
        if isfinite(null_sd_HE) && null_sd_HE > 0
            cross_results(n).rms_dydx_z = ...
                (real_HE - null_mean_HE) / null_sd_HE;
            cross_results(n).null_rms_z = ...
                (representative_HE - null_mean_HE) / null_sd_HE;
        end
    end

    save(fullfile(savepath, 'models', sprintf('%s_hessian_results.mat', area)), ...
         'results', '-v7.3');
    save(fullfile(savepath, 'models', sprintf('%s_cross_deriv_results.mat', area)), ...
         'cross_results', '-v7.3');
end

% =====================================================================
% Periodic 2D smoothing with separable kernel (NaN-aware)
% =====================================================================
function Rs = smooth_periodic_2d(R, k)
    if all(~isfinite(R(:)))
        Rs = R;
        return
    end

    mask = isfinite(R);
    R0 = R; R0(~mask) = 0;

    % Smooth along rows
    num1 = k(1)*circshift(R0,[-1,0]) + k(2)*R0 + k(3)*circshift(R0,[1,0]);
    den1 = k(1)*circshift(mask,[-1,0]) + k(2)*mask + k(3)*circshift(mask,[1,0]);
    R1 = num1 ./ max(den1, eps);
    R1(den1 == 0) = NaN;

    % Smooth along cols
    mask1 = isfinite(R1);
    R1z = R1; R1z(~mask1) = 0;

    num2 = k(1)*circshift(R1z,[0,-1]) + k(2)*R1z + k(3)*circshift(R1z,[0,1]);
    den2 = k(1)*circshift(mask1,[0,-1]) + k(2)*mask1 + k(3)*circshift(mask1,[0,1]);
    Rs = num2 ./ max(den2, eps);
    Rs(den2 == 0) = NaN;
end
