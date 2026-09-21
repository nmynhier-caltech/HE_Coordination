function manifest = Collect_Paper_Figure_Plots(outputRoot)
%COLLECT_PAPER_FIGURE_PLOTS Copy source PDFs into paper-figure folders.
%
% This is the final, non-analytical step in the public pipeline. It does
% not regenerate, edit, or rename any producer output. Instead, it copies
% only the PDFs used to assemble each manuscript or supplemental figure to
% plots/paper_figures/<figure_name>/.
%
% The public script intentionally uses pseudonymized session IDs only.

repoRoot = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(outputRoot)
    outputRoot = fullfile(fileparts(fileparts(repoRoot)), 'results', 'public_analysis');
end

figureRoot = fullfile(outputRoot, 'plots', 'paper_figures');

% Columns: destination folder, source path relative to outputRoot, copied name.
spec = {
    % Figure 1: session_01 SPL, Both-selective neuron 003.
    'figure_1', 'session_01_SPL_results/plots/additivity/neuron_003_additivity3D.pdf', 'session_01_SPL__neuron_003_additivity3D.pdf';
    'figure_1', 'session_01_SPL_results/plots/additivity/neuron_003_H_tuning.pdf',    'session_01_SPL__neuron_003_H_tuning.pdf';
    'figure_1', 'session_01_SPL_results/plots/additivity/neuron_003_E_tuning.pdf',    'session_01_SPL__neuron_003_E_tuning.pdf';

    % Figure 2: M1 hand-selective neuron 017, SPL Both-selective neuron 003,
    % and the matching population separability/additivity summaries.
    'figure_2', 'session_02_M1_results/plots/additivity/neuron_017_additivity3D.pdf', 'session_02_M1__neuron_017_additivity3D.pdf';
    'figure_2', 'session_02_M1_results/plots/additivity/neuron_017_H_tuning.pdf',    'session_02_M1__neuron_017_H_tuning.pdf';
    'figure_2', 'session_02_M1_results/plots/additivity/neuron_017_E_tuning.pdf',    'session_02_M1__neuron_017_E_tuning.pdf';
    'figure_2', 'session_01_SPL_results/plots/additivity/neuron_003_additivity3D.pdf', 'session_01_SPL__neuron_003_additivity3D.pdf';
    'figure_2', 'session_01_SPL_results/plots/additivity/neuron_003_H_tuning.pdf',    'session_01_SPL__neuron_003_H_tuning.pdf';
    'figure_2', 'session_01_SPL_results/plots/additivity/neuron_003_E_tuning.pdf',    'session_01_SPL__neuron_003_E_tuning.pdf';
    'figure_2', 'plots/summary/separability_hessian/M1_Hessian_HE_logq_hist.pdf',      'M1_Hessian_HE_logq_hist.pdf';
    'figure_2', 'plots/summary/separability_hessian/M1_Hessian_rawRMS_2x2.pdf',       'M1_Hessian_rawRMS_2x2.pdf';
    'figure_2', 'plots/summary/additive_reconstruction/M1_R2_add.pdf',                'M1_R2_add.pdf';
    'figure_2', 'plots/summary/additive_reconstruction/M1_corr_hist.pdf',             'M1_corr_hist.pdf';
    'figure_2', 'plots/summary/separability_hessian/SPL_Hessian_HE_logq_hist.pdf',     'SPL_Hessian_HE_logq_hist.pdf';
    'figure_2', 'plots/summary/separability_hessian/SPL_Hessian_rawRMS_2x2.pdf',      'SPL_Hessian_rawRMS_2x2.pdf';
    'figure_2', 'plots/summary/additive_reconstruction/SPL_R2_add.pdf',               'SPL_R2_add.pdf';
    'figure_2', 'plots/summary/additive_reconstruction/SPL_corr_hist.pdf',            'SPL_corr_hist.pdf';

    % Figure 3: the same two example neurons and population selectivity /
    % generalizability summaries.
    'figure_3', 'session_02_M1_results/plots/additivity/neuron_017_H_tuning.pdf', 'session_02_M1__neuron_017_H_tuning.pdf';
    'figure_3', 'session_02_M1_results/plots/additivity/neuron_017_E_tuning.pdf', 'session_02_M1__neuron_017_E_tuning.pdf';
    'figure_3', 'session_01_SPL_results/plots/additivity/neuron_003_H_tuning.pdf','session_01_SPL__neuron_003_H_tuning.pdf';
    'figure_3', 'session_01_SPL_results/plots/additivity/neuron_003_E_tuning.pdf','session_01_SPL__neuron_003_E_tuning.pdf';
    'figure_3', 'plots/summary/selectivity_generalizability/M1_Amp_Selectivity_Across_Sessions_MULTI_vs_SINGLE.pdf', 'M1_Amp_Selectivity_Across_Sessions_MULTI_vs_SINGLE.pdf';
    'figure_3', 'plots/summary/selectivity_generalizability/M1_Single_vs_Multi_Amplitude_HandEye_byClass_4col.pdf',  'M1_Single_vs_Multi_Amplitude_HandEye_byClass_4col.pdf';
    'figure_3', 'plots/summary/selectivity_generalizability/M1_Single_vs_Multi_TuningDir_HandEye_byClass_4col.pdf',  'M1_Single_vs_Multi_TuningDir_HandEye_byClass_4col.pdf';
    'figure_3', 'plots/summary/selectivity_generalizability/SPL_Amp_Selectivity_Across_Sessions_MULTI_vs_SINGLE.pdf','SPL_Amp_Selectivity_Across_Sessions_MULTI_vs_SINGLE.pdf';
    'figure_3', 'plots/summary/selectivity_generalizability/SPL_Single_vs_Multi_Amplitude_HandEye_byClass_4col.pdf', 'SPL_Single_vs_Multi_Amplitude_HandEye_byClass_4col.pdf';
    'figure_3', 'plots/summary/selectivity_generalizability/SPL_Single_vs_Multi_TuningDir_HandEye_byClass_4col.pdf', 'SPL_Single_vs_Multi_TuningDir_HandEye_byClass_4col.pdf';

    % Figure 4: session_01 SPL neuron 003 plus compositional reconstruction
    % and decoding summaries.
    'figure_4', 'session_01_SPL_results/plots/additivity/neuron_003_additivity3D.pdf', 'session_01_SPL__neuron_003_additivity3D.pdf';
    'figure_4', 'plots/summary/additive_reconstruction/SPL_r_SEaddWeighted_vs_MEsmooth_Both.pdf', 'SPL_r_SEaddWeighted_vs_MEsmooth_Both.pdf';
    'figure_4', 'plots/summary/additive_reconstruction/SPL_cME_vs_cSEweighted_scatter.pdf', 'SPL_cME_vs_cSEweighted_scatter.pdf';
    'figure_4', 'plots/summary/additive_reconstruction/SPL_peak_MEadd_and_SEadd_vs_MEraw_scatters.pdf', 'SPL_peak_MEadd_and_SEadd_vs_MEraw_scatters.pdf';
    'figure_4', 'plots/summary/decoding/SPL_AllPoints_Accuracy_vs_Neurons_ME.pdf', 'SPL_AllPoints_Accuracy_vs_Neurons_ME.pdf';
    'figure_4', 'plots/summary/decoding/SPL_AllPoints_Accuracy_vs_Neurons_SE.pdf', 'SPL_AllPoints_Accuracy_vs_Neurons_SE.pdf';
    'figure_4', 'plots/summary/decoding/SPL_MAP_Confusions_BestSession_HE.pdf', 'SPL_MAP_Confusions_BestSession_HE.pdf';

    % Supplemental Figure 1. The three example neurons are embedded in the
    % generated summary PDF rather than imported as separate neuron plots.
    'figure_s1', 'plots/summary/inseparable_neurons/SPL_SeparableVsInseparable_BasicMetrics.pdf', 'SPL_SeparableVsInseparable_BasicMetrics.pdf';
    'figure_s1', 'plots/summary/inseparable_neurons/SPL_PeakClassFractions_Inseparable_vs_Separable.pdf', 'SPL_PeakClassFractions_Inseparable_vs_Separable.pdf';
    'figure_s1', 'plots/summary/inseparable_neurons/SPL_Inseparable_GradientAxis.pdf', 'SPL_Inseparable_GradientAxis.pdf';

    % Supplemental Figure 2: session_04 SPL neurons 005 and 010 are selected
    % within these population-DAB outputs.
    'figure_s2', 'session_04_SPL_results/plots/pop_demo_dab/top_neuron_tuning_alpha_beta.pdf', 'session_04_SPL__top_neuron_tuning_alpha_beta.pdf';
    'figure_s2', 'session_04_SPL_results/plots/pop_demo_dab/dab_tuning_curves.pdf', 'session_04_SPL__dab_tuning_curves.pdf';
    'figure_s2', 'session_04_SPL_results/plots/pop_demo_dab/dab_parallelogram.pdf', 'session_04_SPL__dab_parallelogram.pdf';

    % Supplemental Figure 3: session_01 M1 neuron 028 and SPL neuron 003.
    'figure_s3', 'session_01_M1_results/plots/anova_diag_single/neuron_028_anova_diag.pdf', 'session_01_M1__neuron_028_anova_diag.pdf';
    'figure_s3', 'session_01_SPL_results/plots/anova_diag_single/neuron_003_anova_diag.pdf', 'session_01_SPL__neuron_003_anova_diag.pdf';
    'figure_s3', 'plots/summary/gaussianity_diagnostics/M1_ECDF_Lillie.pdf', 'M1_ECDF_Lillie.pdf';
    'figure_s3', 'plots/summary/gaussianity_diagnostics/M1_slope_hist.pdf', 'M1_slope_hist.pdf';
    'figure_s3', 'plots/summary/gaussianity_diagnostics/SPL_ECDF_Lillie.pdf', 'SPL_ECDF_Lillie.pdf';
    'figure_s3', 'plots/summary/gaussianity_diagnostics/SPL_slope_hist.pdf', 'SPL_slope_hist.pdf';

    % Supplemental Figure 4: population encoding and condition-level
    % decoding plots; no individual-neuron PDF is used.
    'figure_s4', 'plots/summary/selectivity_generalizability/M1_Amp_Tuning_Classification_Bar_MULTI_vs_SINGLE.pdf', 'M1_Amp_Tuning_Classification_Bar_MULTI_vs_SINGLE.pdf';
    'figure_s4', 'plots/summary/selectivity_generalizability/M1_PreferredDirection_Distribution_SignificantOnly.pdf', 'M1_PreferredDirection_Distribution_SignificantOnly.pdf';
    'figure_s4', 'plots/summary/decoding/M1_ME_NeuronMeanFR_Congruent_vs_Incongruent.pdf', 'M1_ME_NeuronMeanFR_Congruent_vs_Incongruent.pdf';
    'figure_s4', 'plots/summary/decoding/M1_MAP_ConditionError_6x6_AcrossSessions.pdf', 'M1_MAP_ConditionError_6x6_AcrossSessions.pdf';
    'figure_s4', 'plots/summary/decoding/M1_MAP_Congruent_Incongruent_Error_Bar_Split_ME_SE.pdf', 'M1_MAP_Congruent_Incongruent_Error_Bar_Split_ME_SE.pdf';
    'figure_s4', 'plots/summary/selectivity_generalizability/SPL_Amp_Tuning_Classification_Bar_MULTI_vs_SINGLE.pdf', 'SPL_Amp_Tuning_Classification_Bar_MULTI_vs_SINGLE.pdf';
    'figure_s4', 'plots/summary/selectivity_generalizability/SPL_PreferredDirection_Distribution_SignificantOnly.pdf', 'SPL_PreferredDirection_Distribution_SignificantOnly.pdf';
    'figure_s4', 'plots/summary/decoding/SPL_ME_NeuronMeanFR_Congruent_vs_Incongruent.pdf', 'SPL_ME_NeuronMeanFR_Congruent_vs_Incongruent.pdf';
    'figure_s4', 'plots/summary/decoding/SPL_MAP_ConditionError_6x6_AcrossSessions.pdf', 'SPL_MAP_ConditionError_6x6_AcrossSessions.pdf';
    'figure_s4', 'plots/summary/decoding/SPL_MAP_Congruent_Incongruent_Error_Bar_Split_ME_SE.pdf', 'SPL_MAP_Congruent_Incongruent_Error_Bar_Split_ME_SE.pdf';

    % Supplemental Figure 5: ME decoding curves and best-session confusion
    % maps; no individual-neuron PDF is used.
    'figure_s5', 'plots/summary/decoding/all_AllPoints_Accuracy_vs_Neurons_ME.pdf', 'all_AllPoints_Accuracy_vs_Neurons_ME.pdf';
    'figure_s5', 'plots/summary/decoding/all_MAP_Confusions_BestSession_HE.pdf', 'all_MAP_Confusions_BestSession_HE.pdf';
    'figure_s5', 'plots/summary/decoding/M1_AllPoints_Accuracy_vs_Neurons_ME.pdf', 'M1_AllPoints_Accuracy_vs_Neurons_ME.pdf';
    'figure_s5', 'plots/summary/decoding/M1_MAP_Confusions_BestSession_HE.pdf', 'M1_MAP_Confusions_BestSession_HE.pdf';
    'figure_s5', 'plots/summary/decoding/SPL_AllPoints_Accuracy_vs_Neurons_ME.pdf', 'SPL_AllPoints_Accuracy_vs_Neurons_ME.pdf';
    'figure_s5', 'plots/summary/decoding/SPL_MAP_Confusions_BestSession_HE.pdf', 'SPL_MAP_Confusions_BestSession_HE.pdf';
    };

sourceFiles = cellfun(@(p) fullfile(outputRoot, strrep(p, '/', filesep)), ...
    spec(:, 2), 'UniformOutput', false);
missing = ~cellfun(@isfile, sourceFiles);
if any(missing)
    error('CollectPaperFigurePlots:MissingSource', ...
        'Missing %d paper-figure source file(s):\n%s', sum(missing), ...
        strjoin(sourceFiles(missing), newline));
end

figureNames = unique(spec(:, 1), 'stable');
for f = 1:numel(figureNames)
    figureName = figureNames{f};
    outputFolder = fullfile(figureRoot, figureName);
    if ~isfolder(outputFolder), mkdir(outputFolder); end

    rows = strcmp(spec(:, 1), figureName);
    expectedNames = spec(rows, 3);
    oldPdfs = dir(fullfile(outputFolder, '*.pdf'));
    for p = 1:numel(oldPdfs)
        if ~ismember(oldPdfs(p).name, expectedNames)
            delete(fullfile(oldPdfs(p).folder, oldPdfs(p).name));
        end
    end

    rowIndices = find(rows);
    for r = reshape(rowIndices, 1, [])
        destination = fullfile(outputFolder, spec{r, 3});
        [ok, message] = copyfile(sourceFiles{r}, destination, 'f');
        if ~ok
            error('CollectPaperFigurePlots:CopyFailed', ...
                'Could not copy %s to %s: %s', sourceFiles{r}, ...
                destination, message);
        end
    end
end

manifest = cell2table(spec, 'VariableNames', ...
    {'FigureFolder', 'SourceRelativeToOutputRoot', 'CopiedFile'});
if ~isfolder(figureRoot), mkdir(figureRoot); end
writetable(manifest, fullfile(figureRoot, 'plot_manifest.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

fprintf('Copied %d source PDFs into %d paper-figure folders under:\n  %s\n', ...
    height(manifest), numel(figureNames), figureRoot);
end
