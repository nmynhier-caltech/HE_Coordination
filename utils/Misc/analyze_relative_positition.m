function analyze_relative_positition(area, savepath)
% Analyze gradient orientation of HE tuning surfaces (6x6) and summarize.
% For each neuron n:
%   1) Build 6x6 HE surface F (mean FR per (hand, eye) bin) from ME trials
%   2) If any cells are missing, impute with additive LS surface
%   3) Compute circular (torus) central-difference gradients dE, dH on F
%   4) Summarize per-neuron orientation by the **dominant axis** (weighted PCA on unit directions)
%   5) Plot (single PDF per neuron):
%       (i)   Heatmap (FR) + quiver (2D; all grid centers), square axes, no whitespace
%       (ii)  Unit-circle compass with gradient points and **dominant axis line only**
%       (iii) Centered cross-derivative RMS: null histogram, blue line at centered real,
%             red dashed line at Bonferroni-corrected cutoff quantile (1 - 0.05/N)
%
% Outputs:
%   plots/gradient/neuron_###_gradient_combo.pdf
%   models/%s_Gradient_HE.mat  (includes axis angle/strength; theta_mean/rho kept but unused in plots)

    % ----------------------------
    % Config
    % ----------------------------
    nbins       = 6;
    edge_shift  = 0.1;           % consistent with your other scripts
    font_sz     = 18;            % figure font size
    quiv_autoSF = 0.9;           % quiver AutoScaleFactor
    nPerm_hist  = 400;           % on-the-fly null draws per neuron for histogram
    min_count   = 1;             % >=1 trial per cell
    dth         = 2*pi/nbins;    % rad/bin for cross-derivative units

    % ----------------------------
    % Load processed data (raw FR)
    % ----------------------------
    S = load(processed_data_file(savepath, area));
    D = S.data;
    FR   = D.FR;    % [trials x neurons]
    TP1  = D.TP1;   % hand targets (x,y)
    TP2  = D.TP2;   % eye  targets (x,y)
    TPi1 = D.TPi1;  % hand effector code
    TPi2 = D.TPi2;  % eye  effector code

    % ME trials where both vary
    multiTP = find(TPi2 ~= 7 & TPi1 ~= 7);

    % ----------------------------
    % Binning (angles) and centers
    % ----------------------------
    edges_theta = linspace(-pi, pi, nbins+1) + edge_shift;
    hcent = edges_theta(1:end-1) + diff(edges_theta)/2;
    ecent = hcent;
    dEdge = 2*pi/nbins;  % bin width

    % Angles & bin indices
    thetaH = atan2(TP1(multiTP,2), TP1(multiTP,1));
    thetaE = atan2(TP2(multiTP,2), TP2(multiTP,1));
    Hb     = discretize(thetaH, edges_theta);
    Eb     = discretize(thetaE, edges_theta);
    valid  = isfinite(Hb) & isfinite(Eb);
    Hb     = Hb(valid);
    Eb     = Eb(valid);
    FR_me  = FR(multiTP(valid), :);   % trials x neurons

    % Pre-list trials per cell (indices into FR_me rows)
    cell_trials = cell(nbins, nbins);
    for i = 1:nbins
        for j = 1:nbins
            cell_trials{i,j} = find(Hb == i & Eb == j);
        end
    end

    % Plot grid
    [Egrid, Hgrid] = meshgrid(ecent, hcent);

    % Output directories
    plot_root = fullfile(savepath, 'plots', 'gradient');
    if ~exist(plot_root, 'dir'), mkdir(plot_root); end

    % ----------------------------
    % Results containers
    % ----------------------------
    N = size(FR_me,2);                 % number of neurons
    alpha_bonf = 0.05 / max(N,1);      % Bonferroni-corrected alpha

    % (We keep theta_mean/rho in the MAT for completeness, but they are not plotted.)
    theta_mean = nan(N,1);             % radians in [0,2pi)
    rho        = nan(N,1);             % 0..1
    n_grad     = nan(N,1);             % number of gradient vectors (36 here)
    theta_all  = cell(N,1);            % per-neuron list of angles (radians)
    w_all      = cell(N,1);            % per-neuron list of magnitudes
    axis_phi   = nan(N,1);             % dominant axis angle (radians)
    axis_alpha = nan(N,1);             % axis strength (anisotropy index)

    fprintf('Analyzing gradient orientation for %d neurons...\n', N);

    % ----------------------------
    % Loop over neurons
    % ----------------------------
    parfor n = 1:N
        try
            % ---------- Build HE surface F (means) ----------
            F = nan(nbins, nbins);
            for i = 1:nbins
                for j = 1:nbins
                    idx = cell_trials{i,j};
                    if numel(idx) >= min_count
                        F(i,j) = mean(FR_me(idx, n), 'omitnan');
                    end
                end
            end
            if any(~isfinite(F(:)))
                F = fill_additive(F);
            end

            % ---------- Gradients on torus ----------
            [dE, dH, mag, theta] = gradients_torus(F);   % radians
            th_vec = mod(theta(:), 2*pi);
            w_vec  = mag(:);

            % ---------- (kept but unused in plots) mean resultant ----------
            [th_mean, rho_val] = resultant_direction(th_vec, w_vec);
            theta_mean(n) = th_mean;
            rho(n)        = rho_val;

            % ---------- Dominant axis via weighted PCA ----------
            [phi_axis, alpha_axis] = dominant_axis_from_angles(th_vec, w_vec);

            n_grad(n)     = numel(th_vec);
            theta_all{n}  = th_vec(:);
            w_all{n}      = w_vec(:);
            axis_phi(n)   = phi_axis;
            axis_alpha(n) = alpha_axis;

            % ---------- Per-neuron cross-derivative: real & null (hist, centered) ----------
            real_rms = cross_deriv_rms_from_surface(F, dth);

            null_rms = nan(nPerm_hist,1);
            rng(1000 + n, 'twister');   % reproducible per neuron
            for p = 1:nPerm_hist
                shuff = FR_me(randperm(size(FR_me,1)), n); % permute trial FR for this neuron
                R_null = nan(nbins, nbins);
                for i = 1:nbins
                    for j = 1:nbins
                        idx = cell_trials{i,j};
                        if numel(idx) >= min_count
                            R_null(i,j) = mean(shuff(idx), 'omitnan');
                        end
                    end
                end
                if any(~isfinite(R_null(:)))
                    R_null = fill_additive(R_null);
                end
                null_rms(p) = cross_deriv_rms_from_surface(R_null, dth);
            end

            % ---------- Figure ----------
            outpdf = fullfile(plot_root, sprintf('neuron_%03d_gradient_combo.pdf', n));
            plot_combo_heatmap_compass_hist(Egrid, Hgrid, F, dE, dH, ...
                                            th_vec, w_vec, phi_axis, alpha_axis, ...
                                            real_rms, null_rms, ...
                                            area, n, font_sz, quiv_autoSF, ...
                                            ecent, hcent, dEdge, ...
                                            alpha_bonf, outpdf);

        catch MEe
            warning('Gradient analysis failed for neuron %d: %s', n, MEe.message);
        end
    end

    % ----------------------------
    % Save results
    % ----------------------------
    Results = struct();
    Results.area        = area;
    Results.nbins       = nbins;
    Results.edge_shift  = edge_shift;
    Results.theta_mean  = theta_mean;   % kept for completeness (not plotted)
    Results.rho         = rho;          % kept for completeness (not plotted)
    Results.n_grad      = n_grad;
    Results.theta_all   = theta_all;
    Results.w_all       = w_all;
    Results.axis_phi    = axis_phi;     % dominant axis angle (radians)
    Results.axis_alpha  = axis_alpha;   % axis strength

    outmat = fullfile(savepath, 'models', sprintf('%s_Gradient_HE.mat', area));
    save(outmat, 'Results', '-v7.3');
    fprintf('Saved results: %s\n', outmat);
end

% =========================
% Helpers
% =========================

function F_add = additive_fit_from_marginals(F)
    mu   = mean(F(:), 'omitnan');
    rbar = mean(F, 2, 'omitnan');
    cbar = mean(F, 1, 'omitnan');
    a = rbar - mean(rbar, 'omitnan');
    b = cbar - mean(cbar, 'omitnan');
    F_add = a * ones(1, size(F,2)) + ones(size(F,1), 1) * b + mu;
end

function F_filled = fill_additive(F)
    F_filled = F;
    F_add = additive_fit_from_marginals(F);
    mask = ~isfinite(F);
    F_filled(mask) = F_add(mask);
end

function [dE, dH, mag, theta] = gradients_torus(F)
% Central-difference gradients with circular wrap on both axes (unit step).
% dE ~ dF/dE (unit/bin), dH ~ dF/dH (unit/bin); theta in radians.
    [r,c] = size(F);
    ip = @(ii) (ii < r) .* (ii+1) + (ii == r) .* 1;
    im = @(ii) (ii > 1) .* (ii-1) + (ii == 1) .* r;
    jp = @(jj) (jj < c) .* (jj+1) + (jj == c) .* 1;
    jm = @(jj) (jj > 1) .* (jj-1) + (jj == 1) .* c;

    dE = nan(r,c); dH = nan(r,c);
    for i = 1:r
        for j = 1:c
            dE(i,j) = (F(i, jp(j)) - F(i, jm(j))) / 2;
            dH(i,j) = (F(ip(i), j) - F(im(i), j)) / 2;
        end
    end

    mag   = hypot(dE, dH);
    theta = atan2(dH, dE);
end

function rms_val = cross_deriv_rms_from_surface(R, dth)
% Mixed partial via finite differences:
% 1) d/dE via gradient(...,1)
% 2) d/dH of that result
% 3) scale by (dth^2) to get Hz/rad^2.
    [~, R_dE_unit]   = gradient(R, 1);      % dF/dE (unit/bin)
    [R_dHdE_unit, ~] = gradient(R_dE_unit); % d/dH (unit/bin^2)
    R_dHdE = R_dHdE_unit / (dth * dth);
    rms_val = sqrt(mean(R_dHdE(:).^2, 'omitnan')); % raw RMS
end

function [theta_mean, rho] = resultant_direction(theta, w)
% (Kept for completeness, not used in plots.)
    if isempty(theta), theta_mean = NaN; rho = NaN; return; end
    if nargin < 2 || isempty(w), w = ones(size(theta)); end
    z = sum(w(:) .* exp(1i*theta(:)));
    if z == 0, theta_mean = NaN; rho = 0; return; end
    theta_mean = angle(z);
    if theta_mean < 0, theta_mean = theta_mean + 2*pi; end
    rho = abs(z) / sum(w(:));
end

function [phi, alpha] = dominant_axis_from_angles(theta, w)
% Weighted PCA on unit direction vectors to get dominant axis (bidirectional).
% Returns:
%   phi   : axis angle (radians) of the principal eigenvector (mod 2π)
%   alpha : (lambda1 - lambda2) / (lambda1 + lambda2) in [0,1] (anisotropy)
    if isempty(theta)
        phi = NaN; alpha = NaN; return;
    end
    if nargin < 2 || isempty(w), w = ones(size(theta)); end
    x = cos(theta(:)); y = sin(theta(:));
    w = w(:);
    % weighted second-moment / covariance-like matrix
    Sxx = sum(w .* x .* x);
    Syy = sum(w .* y .* y);
    Sxy = sum(w .* x .* y);
    S = [Sxx Sxy; Sxy Syy];
    [V, D] = eig(S);
    [lam, idx] = sort(diag(D), 'descend');
    v1 = V(:, idx(1));
    phi = atan2(v1(2), v1(1));
    if phi < 0, phi = phi + 2*pi; end
    if lam(1)+lam(2) > 0
        alpha = (lam(1) - lam(2)) / (lam(1) + lam(2));
    else
        alpha = NaN;
    end
end

function plot_combo_heatmap_compass_hist(Egrid, Hgrid, F, dE, dH, ...
                                         theta_vec, w_vec, phi_axis, alpha_axis, ...
                                         real_rms, null_rms, ...
                                         area, n, font_sz, quiv_autoSF, ...
                                         ecent, hcent, dEdge, alpha_bonf, outpdf)
% 3-row figure:
%   Row 1: heatmap (FR) + quiver (arrows at all 6x6 centers), square axes, no whitespace
%   Row 2: unit-circle compass with gradient points and **dominant axis line only**
%   Row 3: centered null histogram with Bonferroni cutoff and centered real line

    fig = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 8 11]);

    % ---- Row 1: Heatmap + quiver (square, no whitespace) ----
    subplot(3,1,1);
    hold on; box on;
    imagesc(Egrid(1,:), Hgrid(:,1), F); axis xy;
    colormap(parula); colorbar;
    set(gca,'FontSize',font_sz);
    xlabel('\theta_E (rad)'); ylabel('\theta_H (rad)');
    title(sprintf('%s – Neuron %03d: HE heatmap + gradient', area, n), 'FontSize', font_sz);

    % Quiver at all centers (6x6)
    Xq = Egrid(:);  Yq = Hgrid(:);
    Uq = dE(:);     Vq = dH(:);
    quiver(Xq, Yq, Uq, Vq, 'AutoScale','on','AutoScaleFactor',quiv_autoSF, ...
           'Color','k','LineWidth',1.2);

    % Tight square axes based on bin edges (centers ± Δ/2)
    xlim([min(ecent)-dEdge/2, max(ecent)+dEdge/2]);
    ylim([min(hcent)-dEdge/2, max(hcent)+dEdge/2]);
    axis equal tight;

    % ---- Row 2: Unit-circle compass with **dominant axis only** ----
    subplot(3,1,2);
    draw_unit_circle_axis(theta_vec, w_vec, phi_axis, font_sz);
    title(sprintf('Dominant axis (Neuron %03d): angle=%.1f°  axis\\_strength=%.2f', ...
          n, phi_axis*180/pi, alpha_axis), 'FontSize', font_sz);

    % ---- Row 3: Centered NULL histogram vs centered REAL line (+Bonferroni cutoff) ----
    subplot(3,1,3);
    hold on; grid on; box on; set(gca,'FontSize',font_sz);
    null_vals = null_rms(isfinite(null_rms));
    if isempty(null_vals), null_vals = 0; end
    mu_null  = mean(null_vals, 'omitnan');
    null_ctr = null_vals - mu_null;
    real_ctr = real_rms - mu_null;

    lo = min([null_ctr(:); real_ctr], [], 'omitnan');
    hi = max([null_ctr(:); real_ctr], [], 'omitnan');
    if ~isfinite(lo) || ~isfinite(hi) || lo==hi
        lo = min(0, real_ctr); hi = max(1, real_ctr);
    end
    NBINS = max(15, ceil(sqrt(numel(null_ctr))));
    edges = linspace(lo, hi, NBINS+1);

    histogram(null_ctr, edges, 'Normalization','count', ...
              'FaceColor',[0.7 0.7 0.7], 'EdgeColor','k', 'FaceAlpha',0.75, ...
              'DisplayName','Centered null RMS');

    % Bonferroni-corrected cutoff at quantile (1 - alpha_bonf)
    qCut = 1 - alpha_bonf;
    if ~isempty(null_ctr)
        cBonf = quantile(null_ctr, qCut);
        xline(cBonf, 'r--', 'LineWidth', 1.8, ...
              'DisplayName', sprintf('Bonf cutoff q=%.4f', qCut));
    end

    % Centered real RMS
    xline(real_ctr, 'b-', 'LineWidth', 1.8, 'DisplayName','Centered real RMS');

    xlabel('Centered Cross-Derivative RMS (Hz/rad^2)');
    ylabel('Count');
    title(sprintf('%s – Neuron %03d: Centered cross-deriv RMS (null vs real)', area, n));
    legend('Location','best');

    print(fig, outpdf, '-dpdf', '-painters');
    close(fig);
end

function draw_unit_circle_axis(theta_vec, w_vec, phi_axis, font_sz)
% Robust Cartesian compass:
% - Draw unit circle and axes
% - Label +E (0°), +H (90°), -E (180°), -H (270°)
% - Plot gradient directions as points on the unit circle (size ~ |∇F|)
% - Draw the **dominant axis** as a bidirectional line through the origin

    hold on; axis equal; box on;
    R = 1.0;
    t = linspace(0, 2*pi, 400);
    plot(R*cos(t), R*sin(t), 'k-', 'LineWidth', 1.2); % unit circle

    % Axes
    plot([-R R], [0 0], 'k:', 'LineWidth', 1); % E-axis
    plot([0 0], [-R R], 'k:', 'LineWidth', 1); % H-axis

    % Cardinal labels
    text( R*1.08, 0,    '+E', 'HorizontalAlignment','center', 'FontSize', font_sz);
    text(-R*1.08, 0,    '-E', 'HorizontalAlignment','center', 'FontSize', font_sz);
    text(0,  R*1.08,    '+H', 'HorizontalAlignment','center', 'FontSize', font_sz);
    text(0, -R*1.08,    '-H', 'HorizontalAlignment','center', 'FontSize', font_sz);

    % Points for individual gradient directions
    if ~isempty(theta_vec)
        xv = cos(theta_vec(:));
        yv = sin(theta_vec(:));
        if isempty(w_vec)
            w_vec = ones(size(theta_vec(:)));
        else
            w_vec = w_vec(:);
        end
        % Normalize weights for marker sizes
        w_norm = w_vec / max(w_vec + eps);
        ms = 20 + 60 * w_norm;  % between 20 and 80
        scatter(xv, yv, ms, 'MarkerFaceColor', [0.5 0.7 1.0], ...
                'MarkerEdgeColor','k', 'MarkerFaceAlpha', 0.7);
    end

    % Dominant axis (bidirectional)
    if isfinite(phi_axis)
        xd = cos(phi_axis);  yd = sin(phi_axis);
        plot([-xd, xd], [-yd, yd], 'k-', 'LineWidth', 2.0);
    end

    % Limits and ticks
    xlim([-1.25 1.25]); ylim([-1.25 1.25]);
    set(gca,'XTick',[],'YTick',[],'FontSize',font_sz);
end
