%% Repository and run configuration
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'utils')));

% Load the editable pipeline configuration. Values defined by the caller
% take precedence because pipeline_config.m only fills missing variables.
run(fullfile(repoRoot, 'config', 'pipeline_config.m'));

% The downloaded processed data are inputs only. All newly generated files
% go beneath outputRoot, which must be separate from publicDataRoot.
publicDataRoot = char(java.io.File(publicDataRoot).getCanonicalPath());
outputRoot = char(java.io.File(outputRoot).getCanonicalPath());
if ~isfolder(publicDataRoot)
    error('AnalysisMacro:MissingPublicDataRoot', ...
        'The configured publicDataRoot does not exist:\n  %s', publicDataRoot);
end
inputPrefix = [lower(publicDataRoot) filesep];
outputPrefix = [lower(outputRoot) filesep];
if strcmpi(publicDataRoot, outputRoot) || startsWith(outputPrefix, inputPrefix)
    error('AnalysisMacro:OverlappingRoots', ...
        ['outputRoot must be separate from, and not nested inside, the ' ...
         'read-only publicDataRoot.\npublicDataRoot: %s\noutputRoot: %s'], ...
        publicDataRoot, outputRoot);
end
configure_public_data_root(publicDataRoot);

% Fix the global random stream for any stochastic MATLAB routine that does
% not set a more specific per-analysis/per-neuron seed. Individual analyses
% retain their existing local seeds, so their numerical outputs are unchanged.
rng(randomSeed, 'twister');

% Load the public manifest. It intentionally has no dates or raw filenames.
sessions = public_session_manifest();

% Shared destinations used by the across-session summary functions.
sharedOutputDirectories = {fullfile(outputRoot, 'models'), ...
    fullfile(outputRoot, 'plots')};
for d = 1:numel(sharedOutputDirectories)
    if ~isfolder(sharedOutputDirectories{d}), mkdir(sharedOutputDirectories{d}); end
end

% Preallocate one output-path list per area.
savepaths_M1  = cell(sum(strcmp(sessions(:, 2), 'M1')), 1);
savepaths_SPL = cell(sum(strcmp(sessions(:, 2), 'SPL')), 1);
savepaths_all = cell(sum(strcmp(sessions(:, 2), 'all')), 1);
areaIndex = struct('M1', 0, 'SPL', 0, 'all', 0);

% Loop over sessions
for s = 1:size(sessions, 1)

    sessionID = sessions{s, 1};
    area = sessions{s, 2};
    savepath = fullfile(outputRoot, sprintf('%s_%s_results', sessionID, area));
    outputDirectories = {savepath, ...
        fullfile(savepath, 'models'), fullfile(savepath, 'plots')};
    for d = 1:numel(outputDirectories)
        if ~isfolder(outputDirectories{d}), mkdir(outputDirectories{d}); end
    end

    areaIndex.(area) = areaIndex.(area) + 1;
    idx = areaIndex.(area);
    switch area
        case 'M1'
            savepaths_M1{idx} = savepath;
        case 'SPL'
            savepaths_SPL{idx} = savepath;
        case 'all'
            savepaths_all{idx} = savepath;
        otherwise
            error('AnalysisMacro:UnknownArea', 'Unknown area "%s".', area);
    end

    if runSessionAnalyses
        processedFile = processed_data_file(savepath, area);
        if ~isfile(processedFile)
            error('AnalysisMacro:MissingProcessedData', ...
                ['Missing public pipeline input: %s\nInstall the processed ' ...
                 'data release as documented in README.md.'], ...
                processedFile);
        end
        analyze_additivity_surfaces(area, savepath);
        population_additivity_PCA_demo(area, savepath);
        cross_derivative_test(area, savepath);
        analyze_relative_positition(area, savepath);
        compute_FR_amplitude_selectivity(area, savepath);
        nonparametric_MAP(area, savepath);
        check_FR_gaussianity(area, savepath);
    end

    fprintf('Prepared session %d/%d: %s %s\n', ...
        s, size(sessions, 1), sessionID, area);
end

%% Graphical abstract
if runGraphicalAbstract
    Graphical_Abstract_3DPlot(fullfile(outputRoot, 'plots', ...
        'graphical_abstract_components'));
end


%% Summary Figures
if runSummaryFigures
    %% Fig 2
    % Plot surface-additivity summaries.
    plot_additivity_summary_across_sessions('M1', savepaths_M1);
    plot_additivity_summary_across_sessions('SPL', savepaths_SPL);
    plot_additivity_summary_across_sessions('all', savepaths_all);

    % Plot cross-derivative RMS distributions.
    plot_cross_derivative_rms_summary('M1', savepaths_M1);
    plot_cross_derivative_rms_summary('SPL', savepaths_SPL);
    plot_cross_derivative_rms_summary('all', savepaths_all);

    %% Fig 3
    % Plot selectivity-classification differences.
    plot_FR_amplitude_summary_across_sessions('M1', savepaths_M1);
    plot_FR_amplitude_summary_across_sessions('SPL', savepaths_SPL);
    plot_FR_amplitude_summary_across_sessions('all', savepaths_all);

    %% Fig 4
    % Plot nonparametric decoding results.
    plot_MAP_confusion_summary_across_sessions('M1', savepaths_M1);
    plot_MAP_confusion_summary_across_sessions('SPL', savepaths_SPL);
    plot_MAP_confusion_summary_across_sessions('all', savepaths_all);

    %% Supplemental Figure: Congruent/Anisotropic
    plot_MAP_condition_accuracy_6x6_summary('M1', savepaths_M1);
    plot_MAP_condition_accuracy_6x6_summary('SPL', savepaths_SPL);
    plot_MAP_condition_accuracy_6x6_summary('all', savepaths_all);

    %% Supplemental Figure: Inseparable Neurons
    plot_inseparable_neuron_summary('M1', savepaths_M1);
    plot_inseparable_neuron_summary('SPL', savepaths_SPL);
    plot_inseparable_neuron_summary('all', savepaths_all);

    %% Supplemental Figure: Relative Position
    plot_relative_position_summary('M1', savepaths_M1);
    plot_relative_position_summary('SPL', savepaths_SPL);

    %% Supplemental Figure: Check Gaussianity
    summarize_FR_gaussianity('M1', savepaths_M1);
    summarize_FR_gaussianity('SPL', savepaths_SPL);
end

%% Collect the exact source PDFs used in the paper figures
if runPaperFigureCollection
    Collect_Paper_Figure_Plots(outputRoot);
end






