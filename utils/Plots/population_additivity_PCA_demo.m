function population_additivity_PCA_demo(area, savepath)
% Population additivity demo using a Demixed Additive Basis (DAB)
% (updated: principled corner selection from α/β, sign convention,
%  top-2 α and top-2 β neurons, letter-sized figures)
%
% Outputs under plots/pop_demo_dab/:
%   top_neuron_tuning_alpha_beta.pdf   (4 neurons × {hand, eye} = 8 panels)
%   dab_tuning_curves.pdf              ({alpha, beta} × {hand, eye} = 4 panels)
%   dab_parallelogram.pdf              (α vs β plane)
% Also saves:
%   <area>_DABmap.mat with fields: alpha, beta, mu, sigma, z_rowmean, bins

    % ----------------------------
    % Config
    % ----------------------------
    nbins      = 6;
    edge_shift = 0.1;

    % ----------------------------
    % Load processed data
    % ----------------------------
    S = load(processed_data_file(savepath, area));
    D = S.data;
    FR   = D.FR;    % [trials x neurons]
    TP1  = D.TP1;   % hand (x,y)
    TP2  = D.TP2;   % eye  (x,y)
    TPi1 = D.TPi1;  % hand effector code
    TPi2 = D.TPi2;  % eye  effector code

    % ME trials (both engaged)
    multiTP = find(TPi2 ~= 7 & TPi1 ~= 7);

    % ----------------------------
    % Angle binning
    % ----------------------------
    edges_theta = linspace(-pi, pi, nbins+1) + edge_shift;
    hcent = edges_theta(1:end-1) + diff(edges_theta)/2;
    ecent = hcent;

    thetaH_multi = atan2(TP1(multiTP,2), TP1(multiTP,1));
    thetaE_multi = atan2(TP2(multiTP,2), TP2(multiTP,1));
    Hb_multi = discretize(thetaH_multi, edges_theta);
    Eb_multi = discretize(thetaE_multi, edges_theta);

    % (i,j) trial indices
    cell_indices = cell(nbins, nbins);
    for i = 1:nbins
        for j = 1:nbins
            cell_indices{i,j} = multiTP(Hb_multi == i & Eb_multi == j);
        end
    end

    % ----------------------------
    % Per-neuron ME surfaces  F_all(i,j,n)
    % ----------------------------
    N = size(FR,2);
    F_all = nan(nbins, nbins, N);

    for n = 1:N
        for i = 1:nbins
            for j = 1:nbins
                idx = cell_indices{i,j};
                if ~isempty(idx)
                    fr_ij = FR(idx, n);
                    fr_ij = fr_ij(~isnan(fr_ij));
                    if ~isempty(fr_ij)
                        F_all(i,j,n) = mean(fr_ij);
                    end
                end
            end
        end
    end

    % ----------------------------
    % Arrange condition-by-neuron matrix + z-score
    % ----------------------------
    H = nbins; E = nbins; C = H*E;
    X = reshape(F_all, [C, N]);    % C x N (condition means)

    mu    = nanmean(X, 1);         % per-neuron mean
    sigma = nanstd(X, 0, 1);       % per-neuron std
    sigma(~isfinite(sigma) | sigma==0) = 1;

    Z = (X - mu) ./ sigma;         % C x N (NaNs allowed)
    Z(isnan(Z)) = 0;               % fill at 0 (z-mean) to define axes

    % Row-centering (affine invariance)
    z_rowmean = mean(Z, 1);        % 1 x N
    Zc = Z - z_rowmean;            % C x N

    % Reshape to H x E x N tensor for marginals
    S3 = reshape(Zc, [H, E, N]);   % H x E x N

    % ----------------------------
    % Demixed additive axes via marginals + SVD
    % ----------------------------
    RH = squeeze(mean(S3,2));  RH = RH - mean(RH,1);  RH = RH';   % N x H (hand)
    RE = squeeze(mean(S3,1));  RE = RE - mean(RE,1);  RE = RE';   % N x E (eye)

    [U_H,~,~] = svds(RH, 1);  alpha = U_H(:,1);       % N x 1
    [U_E,~,~] = svds(RE, 1);  beta  = U_E(:,1);       % N x 1

    % Orthonormalize
    alpha = alpha / norm(alpha);
    beta  = beta - alpha*(alpha'*beta);
    if norm(beta) > 0
        beta = beta / norm(beta);
    else
        v = zeros(N,1); v(1) = 1;
        beta = v - alpha*(alpha'*v); beta = beta / norm(beta);
    end

    % Store the linear map
    DABmap.alpha     = alpha;
    DABmap.beta      = beta;
    DABmap.mu        = mu;
    DABmap.sigma     = sigma;
    DABmap.z_rowmean = z_rowmean;
    DABmap.hcent     = hcent;
    DABmap.ecent     = ecent;

    out_root = fullfile(savepath, 'plots', 'pop_demo_dab');
    if ~exist(out_root, 'dir'), mkdir(out_root); end
    save(fullfile(out_root, sprintf('%s_DABmap.mat', area)), 'DABmap');

    % ----------------------------
    % Project each condition into the DAB plane
    % ----------------------------
    AlphaBeta = [alpha, beta];        % N x 2
    AB_scores = (Zc) * AlphaBeta;     % C x 2   (s = col1, t = col2)
    S_coord = reshape(AB_scores(:,1), H, E);  % α coordinate
    T_coord = reshape(AB_scores(:,2), H, E);  % β coordinate

    % DAB “tuning curves” (marginals)
    A_H = mean(S_coord, 2, 'omitnan');  % α vs hand
    A_E = mean(S_coord, 1, 'omitnan');  % α vs eye
    B_H = mean(T_coord, 2, 'omitnan');  % β vs hand
    B_E = mean(T_coord, 1, 'omitnan');  % β vs eye

    % ----------------------------
    % Sign convention:
    %   flip α so max(A_H) > 0; flip β so max(B_E) > 0
    % ----------------------------
    if any(isfinite(A_H))
        if max(A_H) < abs(min(A_H))
            alpha   = -alpha;
            S_coord = -S_coord;
            A_H     = -A_H;  A_E = -A_E;
        end
    end
    if any(isfinite(B_E))
        if max(B_E) < abs(min(B_E))
            beta    = -beta;
            T_coord = -T_coord;
            B_H     = -B_H;  B_E = -B_E;
        end
    end

    % (repack map in case of flips)
    DABmap.alpha = alpha; DABmap.beta = beta;
    save(fullfile(out_root, sprintf('%s_DABmap.mat', area)), 'DABmap');

    % ----------------------------
    % Corner selection (principled demixed):
    %   x1,x2 from α vs hand (A_H)
    %   y1,y2 from β vs eye  (B_E)
    % ----------------------------
    [~, ixH_peak]   = max(A_H);  [~, ixH_trough] = min(A_H);
    [~, iyE_peak]   = max(B_E);  [~, iyE_trough] = min(B_E);

    x1 = hcent(ixH_peak);   x2 = hcent(ixH_trough);
    y1 = ecent(iyE_peak);   y2 = ecent(iyE_trough);

    % ----------------------------
    % Top neurons: 2 by |α| and 2 by |β| (prefer non-overlap)
    % ----------------------------
    [~, ordA] = sort(abs(alpha), 'descend');
    [~, ordB] = sort(abs(beta ), 'descend');

    top2_alpha = unique(ordA(1:min(2,N)),'stable');
    top2_beta  = []; k = 1;
    for ii = 1:min( N, max(4,2*numel(top2_alpha)) )
        if numel(top2_beta) >= 2, break; end
        cand = ordB(ii);
        if ~ismember(cand, top2_alpha)
            top2_beta(end+1) = cand; %#ok<AGROW>
        end
    end
    if numel(top2_beta) < 2
        % allow overlaps if necessary
        need = 2 - numel(top2_beta);
        top2_beta = [top2_beta, ordB(1:need)];
    end

%     show_neurons = [top2_alpha(:); top2_beta(:)];
    show_neurons = [top2_alpha(1); top2_beta(1)]; % only top neuron
    show_tags    = [repmat("α-top",numel(top2_alpha),1); repmat("β-top",numel(top2_beta),1)];

    % ----------------------------
    % FIGURE 1: Top α/β neurons – hand & eye curves (8 panels)
    % ----------------------------
    fig1 = figure('Visible','off','PaperUnits','inches');
    set(fig1,'PaperSize',[8.5 11],'PaperPosition',[0.25 0.25 8 10.5], ...
        'PaperPositionMode','manual');
    t = tiledlayout(fig1, numel(show_neurons), 2, 'TileSpacing','compact','Padding','compact');

    for k = 1:numel(show_neurons)
        n = show_neurons(k);
        F = F_all(:,:,n);
        hand_mu = mean(F, 2, 'omitnan');
        eye_mu  = mean(F, 1, 'omitnan');

        % get ylims
        ylim_high = max(max(hand_mu), max(eye_mu));
        ylim_low = min(min(hand_mu), min(eye_mu));
        
        % plot hand curve
        nexttile(t);
        set(gca, 'FontSize', 20); hold on;
        plot(hcent, hand_mu, 'k-o','LineWidth',1.4); grid on;
        xlim([-pi, pi]); xlabel('\theta_H (rad)'); ylim([ylim_low, ylim_high]); ylabel('FR [Hz]');
        title(sprintf('%s | Neuron %d (|\\alpha|=%.3f, |\\beta|=%.3f): Hand', ...
              show_tags(k), n, abs(alpha(n)), abs(beta(n))));
        hold off;

        % plot eye curve
        nexttile(t);
        set(gca, 'FontSize', 20); hold on;
        plot(ecent, eye_mu, 'k-o','LineWidth',1.4); grid on;
        xlim([-pi, pi]); xlabel('\theta_E (rad)'); ylim([ylim_low, ylim_high]); ylabel('FR [Hz]');
        title(sprintf('%s | Neuron %d (|\\alpha|=%.3f, |\\beta|=%.3f): Eye', ...
              show_tags(k), n, abs(alpha(n)), abs(beta(n))));
        hold off;
        
    end
    print(fig1, '-dpdf', fullfile(out_root, 'top_neuron_tuning_alpha_beta.pdf'), '-painters');
    close(fig1);

    % ----------------------------
    % DAB tuning curves
    % ----------------------------
    fig2 = figure('Visible','off','PaperUnits','inches');
    set(fig2,'PaperSize',[8.5 11],'PaperPosition',[0.5 0.75 7.5 9.5], ...
        'PaperPositionMode','manual');
    t2 = tiledlayout(fig2, 2, 2, 'TileSpacing','compact','Padding','compact');

    vline = @(xv) xline(xv, 'r--', 'LineWidth', 1.2);

    % α vs Hand
    nexttile(t2);
    set(gca, 'FontSize', 20); hold on;
    plot(hcent, A_H, 'k-o','LineWidth',1.6); grid on;
    vline(x1); vline(x2);
    xlim([-pi,pi]); xlabel('\theta_H (rad)'); ylabel('\alpha-score');
    title('\alpha vs Hand');
    hold off;

    % α vs Eye
    nexttile(t2);
    set(gca, 'FontSize', 20); hold on;
    plot(ecent, A_E, 'k-o','LineWidth',1.6); grid on;
    vline(y1); vline(y2);  % tag β-chosen eye corners here too for reference
    xlim([-pi,pi]); xlabel('\theta_E (rad)'); ylabel('\alpha-score');
    title('\alpha vs Eye');
    hold off;

    % β vs Hand
    nexttile(t2);
    set(gca, 'FontSize', 20); hold on;
    plot(hcent, B_H, 'k-o','LineWidth',1.6); grid on;
    vline(x1); vline(x2);  % tag α-chosen hand corners
    xlim([-pi,pi]); xlabel('\theta_H (rad)'); ylabel('\beta-score');
    title('\beta vs Hand');
    hold off;

    % β vs Eye
    nexttile(t2);
    set(gca, 'FontSize', 20); hold on;
    plot(ecent, B_E, 'k-o','LineWidth',1.6); grid on;
    vline(y1); vline(y2);
    xlim([-pi,pi]); xlabel('\theta_E (rad)'); ylabel('\beta-score');
    title('\beta vs Eye');
    hold off;

    sgtitle(sprintf('%s – DAB tuning curves (x1=%.2f, x2=%.2f, y1=%.2f, y2=%.2f)', ...
                    area, x1, x2, y1, y2), 'FontSize', 14);
    print(fig2, '-dpdf', fullfile(out_root, 'dab_tuning_curves.pdf'), '-painters');
    close(fig2);

    % ----------------------------
    % Parallelogram in the DAB plane from four corner conditions
    %   corners: (x1,x2) from α–hand ; (y1,y2) from β–eye
    % ----------------------------
    i_x1 = ixH_peak;   i_x2 = ixH_trough;
    j_y1 = iyE_peak;   j_y2 = iyE_trough;

    pts = [
        S_coord(i_x1, j_y1), T_coord(i_x1, j_y1);  % (x1,y1)
        S_coord(i_x1, j_y2), T_coord(i_x1, j_y2);  % (x1,y2)
        S_coord(i_x2, j_y2), T_coord(i_x2, j_y2);  % (x2,y2)
        S_coord(i_x2, j_y1), T_coord(i_x2, j_y1)   % (x2,y1)
    ];

    fig3 = figure('Visible','off','PaperUnits','inches');
    set(fig3,'PaperSize',[8.5 11],'PaperPosition',[1 2 6.5 6.5], ...
        'PaperPositionMode','manual');
    hold on; grid on; axis equal;
    set(gca, 'FontSize', 20);
    if any(~isfinite(pts), 'all')
        warning('Some parallelogram points are NaN; check coverage/missingness.');
    end
    plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor','k','LineWidth',1.5);
    cyc = [1 2 3 4 1];
    plot(pts(cyc,1), pts(cyc,2), 'k-','LineWidth',1.5);

    labels = {'(x1,y1)','(x1,y2)','(x2,y2)','(x2,y1)'};
    for m = 1:4
        text(pts(m,1), pts(m,2), ['  ' labels{m}], 'FontSize',10);
    end
    xlabel('\alpha-score'); ylabel('\beta-score');
    title(sprintf('%s – Parallelogram in DAB plane', area));
    print(fig3, '-dpdf', fullfile(out_root, 'dab_parallelogram.pdf'), '-painters');
    close(fig3);

    fprintf('Saved:\n  %s\n  %s\n  %s\n', ...
        fullfile(out_root,'top_neuron_tuning_alpha_beta.pdf'), ...
        fullfile(out_root,'dab_tuning_curves.pdf'), ...
        fullfile(out_root,'dab_parallelogram.pdf'));
end
