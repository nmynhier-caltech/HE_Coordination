function plot_relative_position_summary(area, savepaths)
% Class-conditioned summary for analyze_relative_positition outputs.
% For each class in {'H','E','Both'}, generate three PDFs:
%   1) Dominant axis angle histogram (0–180°, modulo π)
%   2) Axis strength histogram (alpha in [0,1])
%   3) Overlay histogram: real vs null cross-derivative RMS (centered), with 95% cutoff
%
% Inputs:
%   area:      area tag used in filenames (e.g., 'MC', 'PPC')
%   savepaths: cellstr of session directories
%
% Outputs (PDFs):
%   ../plots/summary/relative_position/
%     <area>_<CLASS>_dominant_axis_angle_hist.pdf
%     <area>_<CLASS>_dominant_axis_strength_hist.pdf
%     <area>_<CLASS>_RMS_hist_real_vs_null.pdf

    classes = {'H','E','Both'};

    % -------- Accumulators per class --------
    ACC = struct();
    for k = 1:numel(classes)
        c = classes{k};
        ACC.(c).axis_phi   = [];  % radians (dominant axis)
        ACC.(c).axis_alpha = [];  % [0,1]
        ACC.(c).real_rms   = [];  % centered real RMS
        ACC.(c).null_rms   = [];  % one centered null RMS
    end

    % -------- Pool across sessions --------
    for d = 1:numel(savepaths)
        sp = savepaths{d};
        try
            % Selectivity class labels
            amp_file = fullfile(sp, 'models', sprintf('%s_FR_amplitude_selectivity.mat', area));
            if ~exist(amp_file, 'file')
                warning('Missing class file: %s', amp_file);
                continue;
            end
            A = load(amp_file);
            assert(isfield(A, 'class_labels'), 'Missing class_labels in %s', amp_file);
            labels = string(A.class_labels(:));

            % Gradient (dominant axis)
            grad_file = fullfile(sp, 'models', sprintf('%s_Gradient_HE.mat', area));
            if ~exist(grad_file, 'file')
                warning('Missing gradient file: %s', grad_file);
                continue;
            end
            G = load(grad_file);
            assert(isfield(G, 'Results'), 'Missing Results in %s', grad_file);
            axis_phi   = G.Results.axis_phi(:);     % radians
            axis_alpha = G.Results.axis_alpha(:);   % [0,1]

            % Cross-derivative per-neuron (centered)
            cd_file = fullfile(sp, 'models', sprintf('%s_cross_deriv_results.mat', area));
            if ~exist(cd_file, 'file')
                warning('Missing cross-derivative file: %s', cd_file);
                continue;
            end
            C = load(cd_file);
            if isfield(C, 'results')
                rr = C.results;
            elseif isfield(C, 'cross_results')
                rr = C.cross_results;
            else
                error('Missing results/cross_results in %s', cd_file);
            end
            real_vec = arrayfun(@(s) s.rms_dydx, rr);  % centered real RMS
            null_vec = arrayfun(@(s) s.null_rms, rr);  % one centered null RMS

            % ---- Align lengths robustly (truncate to min) ----
            n1 = numel(labels);
            n2 = numel(axis_phi);
            n3 = numel(axis_alpha);
            n4 = numel(real_vec);
            n5 = numel(null_vec);
            n  = min([n1 n2 n3 n4 n5]);
            if n == 0
                warning('No overlap for session %s (n=0 after alignment).', sp);
                continue;
            end
            labels    = labels(1:n);
            axis_phi  = axis_phi(1:n);
            axis_alpha= axis_alpha(1:n);
            real_vec  = real_vec(1:n);
            null_vec  = null_vec(1:n);

            % ---- Append per class ----
            for k = 1:numel(classes)
                cls = classes{k};
                mask = labels == string(cls);
                if any(mask)
                    ACC.(cls).axis_phi   = [ACC.(cls).axis_phi;   axis_phi(mask)];   
                    ACC.(cls).axis_alpha = [ACC.(cls).axis_alpha; axis_alpha(mask)]; 
                    ACC.(cls).real_rms   = [ACC.(cls).real_rms;   real_vec(mask)];   
                    ACC.(cls).null_rms   = [ACC.(cls).null_rms;   null_vec(mask)];  
                end
            end

        catch ME
            warning('Error in %s: %s', sp, ME.message);
        end
    end

    % -------- Save directory --------
    plotdir = fullfile(savepaths{1}, '..', 'plots', 'summary', ...
        'relative_position');
    if ~exist(plotdir, 'dir'), mkdir(plotdir); end

    % -------- Colors (consistent with your palette) --------
    cmap = containers.Map( ...
        {'H','E','Both'}, ...
        {[0.8 0.2 0.2], [0.2 0.2 0.8], [0.2 0.6 0.2]} ...
    );

    % -------- Generate 3 plots per class --------
    for k = 1:numel(classes)
        cls = classes{k};
        baseColor = cmap(cls);

        % --- Angles: modulo π → [0, 180) deg ---
        phi = ACC.(cls).axis_phi;
        phi = phi(isfinite(phi));
        phi_deg = mod(phi, pi) * 180/pi;  % axis is bidirectional

        fig1 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6 4.5]);
        ax1 = axes('Parent',fig1); hold(ax1,'on'); set(ax1,'FontSize',16); box(ax1,'on'); grid(ax1,'on');
        if ~isempty(phi_deg)
            NBINS = 18;  % 10° bins
            edges = linspace(0, 180, NBINS+1);
            histogram(ax1, phi_deg, edges, 'FaceColor', baseColor, ...
                      'EdgeColor', 'k', 'FaceAlpha', 0.85, 'Normalization','count');
        else
            warning('No axis angles to plot for class %s.', cls);
        end
        set(gca, 'FontSize', 20);
        xlabel(ax1, sprintf('%s Dominant axis angle (deg, mod 180)', cls));
        ylabel(ax1, 'Neuron count');
        title(ax1, sprintf('%s – %s: Dominant axis distribution', area, cls));
        print(fig1, fullfile(plotdir, sprintf('%s_%s_dominant_axis_angle_hist.pdf', area, cls)), '-dpdf', '-painters');
        close(fig1);

        % --- Strengths: alpha in [0,1] ---
        alpha_vals = ACC.(cls).axis_alpha;
        alpha_vals = alpha_vals(isfinite(alpha_vals) & alpha_vals >= 0 & alpha_vals <= 1);

        fig2 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6 4.5]);
        ax2 = axes('Parent',fig2); hold(ax2,'on'); set(ax2,'FontSize',16); box(ax2,'on'); grid(ax2,'on');
        if ~isempty(alpha_vals)
            NBINS = 20;
            edges = linspace(0, 1, NBINS+1);
            histogram(ax2, alpha_vals, edges, 'FaceColor', baseColor, ...
                      'EdgeColor', 'k', 'FaceAlpha', 0.85, 'Normalization','count');
        else
            warning('No axis strengths to plot for class %s.', cls);
        end
        set(gca, 'FontSize', 20);
        xlabel(ax2, sprintf('%s Axis strength \\alpha', cls));
        ylabel(ax2, 'Neuron count');
        title(ax2, sprintf('%s – %s: Dominant axis strength', area, cls));
        print(fig2, fullfile(plotdir, sprintf('%s_%s_dominant_axis_strength_hist.pdf', area, cls)), '-dpdf', '-painters');
        close(fig2);

        % --- Cross-derivative RMS: real vs null (centered values) ---
        real_rms = ACC.(cls).real_rms;  real_rms = real_rms(isfinite(real_rms));
        null_rms = ACC.(cls).null_rms;  null_rms = null_rms(isfinite(null_rms));

        fig3 = figure('Visible','Off','PaperUnits','inches','PaperPosition',[0 0 6 4.5]);
        ax3 = axes('Parent',fig3); hold(ax3,'on'); set(ax3,'FontSize',16); box(ax3,'on'); grid(ax3,'on');

        if isempty(real_rms)
            warning('No real RMS values to plot for class %s.', cls);
        end
        if isempty(null_rms)
            warning('No null RMS values to plot for class %s.', cls);
        end
        if isempty(real_rms) && isempty(null_rms)
            close(fig3);
        else
            lo = min([real_rms; null_rms], [], 'omitnan');
            hi = max([real_rms; null_rms], [], 'omitnan');
            if ~isfinite(lo) || ~isfinite(hi) || lo==hi
                lo = 0; hi = 1;
            end
            NB = 40;
            edges = linspace(lo, hi, NB+1);

            if ~isempty(null_rms)
                histogram(ax3, null_rms, edges, 'Normalization','count', ...
                          'FaceColor',[0.7 0.7 0.7], 'EdgeColor','k', 'FaceAlpha',0.6, ...
                          'DisplayName','Null RMS (per neuron)');
                % 95% cutoff (like your example)
                alpha = 0.05;
                null_thr = quantile(null_rms, 1 - alpha);
                xline(ax3, null_thr, 'r--', 'LineWidth', 1.8, ...
                      'DisplayName', sprintf('Null %.0f%% cutoff', 100*(1-alpha)));
            end

            if ~isempty(real_rms)
                histogram(ax3, real_rms, edges, 'Normalization','count', ...
                          'FaceColor', baseColor, 'EdgeColor','k', 'FaceAlpha',0.7, ...
                          'DisplayName','Real RMS (per neuron)');
                mu_real = mean(real_rms, 'omitnan');
                xline(ax3, mu_real, 'b-.', 'LineWidth', 1.6, 'DisplayName','Mean Real RMS');
            end
            set(gca, 'FontSize', 20);
            xlabel(ax3, 'Centered Cross-Derivative RMS (Hz/rad^2)');
            ylabel(ax3, 'Neuron count');
            title(ax3, sprintf('%s – %s: Real vs Null RMS', area, cls));
            legend(ax3, 'Location', 'best');

            print(fig3, fullfile(plotdir, sprintf('%s_%s_RMS_hist_real_vs_null.pdf', area, cls)), '-dpdf', '-painters');
            close(fig3);
        end
    end
end
