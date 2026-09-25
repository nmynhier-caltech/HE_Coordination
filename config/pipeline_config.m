%% HE Coordination public pipeline configuration
% Edit the two paths below, then run Analysis_Macro.m. Variables set
% before the macro is called take precedence over these defaults.

% file locations
% publicDataRoot = [FILL];
% outputRoot = [FILL];

% Extract the CaltechDATA release here. The pipeline treats it as read-only.
if ~exist('publicDataRoot', 'var') || isempty(publicDataRoot)
    publicDataRoot = fullfile(fileparts(repoRoot), 'data', 'public_release');
end

% All generated models and plots are written here. This must not be the same
% directory as publicDataRoot or a directory inside publicDataRoot.
if ~exist('outputRoot', 'var') || isempty(outputRoot)
    outputRoot = fullfile(fileparts(repoRoot), 'results', 'public_analysis');
end

% A full clean reproduction uses true for all four stages. Change individual
% flags to false for faster reruns from outputs that already exist.
if ~exist('runSessionAnalyses', 'var'), runSessionAnalyses = true; end
if ~exist('runSummaryFigures', 'var'), runSummaryFigures = true; end
if ~exist('runGraphicalAbstract', 'var'), runGraphicalAbstract = true; end
if ~exist('runPaperFigureCollection', 'var')
    runPaperFigureCollection = true;
end

if ~exist('randomSeed', 'var') || isempty(randomSeed), randomSeed = 1; end
