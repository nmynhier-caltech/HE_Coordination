# HE Coordination: reproducible public analysis

This repository reproduces the downstream analyses and MATLAB-generated
source plots for the HE Coordination paper from pseudonymized, minimally
processed neural data. It contains no acquisition dates, raw file names,
continuous eye-tracking samples, or private lab code.

## What this repository reproduces

The pipeline performs all per-session analyses, across-session summaries,
decoding analyses, graphical-abstract plotting, and collection of the source
PDFs used in Figures 1-4 and Supplemental Figures 1-5.

It does not reproduce final Illustrator assembly, schematic artwork, panel
lettering, or cosmetic edits applied after the MATLAB plots were exported.

## Requirements

The code was tested with:

- MATLAB R2021a and R2024b
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox
- Windows and Mac OS

The detailed machine-readable dependency record is
[`dependencies.yaml`](dependencies.yaml). Generated models and plots are
written beneath an external data directory, not into the Git repository.

## 1. Obtain and arrange the processed data

Download the processed-data release from CaltechDATA and extract it into a
directory of your choice. This read-only directory is called `publicDataRoot`
below. Choose a separate writable `outputRoot` for generated results.

The release contains all files in one flat directory. The expected names are
listed in [`config/data_manifest.csv`](config/data_manifest.csv). There are 11
pseudonymized sessions for each of three area labels (`M1`, `SPL`, and
combined `all`), for a total of 33 input MAT files:

```text
<publicDataRoot>/
    session_01_M1_Processed_data.mat
    session_01_SPL_Processed_data.mat
    session_01_all_Processed_data.mat
    ...
    session_11_M1_Processed_data.mat
    session_11_SPL_Processed_data.mat
    session_11_all_Processed_data.mat
```

Each MAT file contains a scalar struct named `data` with fields `FR`, `TP1`,
`TP2`, `TPi1`, and `TPi2`. See the complete
[`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) for matrix dimensions, units,
condition codes, and trial-alignment details. `TP2` and `TPi2` are task
target/condition labels, not continuous eye-tracking measurements.

## 2. Configure the two pipeline locations

Open [`config/pipeline_config.m`](config/pipeline_config.m) and set:

- `publicDataRoot` to the extracted CaltechDATA directory; and
- `outputRoot` to a separate writable directory for generated results.

The supplied defaults use `data/public_release/` and
`results/public_analysis/` beside the code. These folders are ignored by Git.
The macro rejects an `outputRoot` that is equal to or nested inside
`publicDataRoot`, preventing results from being mixed into the downloaded
processed-data release.

All analysis stages default to enabled, so a new user normally edits only
these two paths.

## 3. Run the complete pipeline

Open a terminal in the code-repository root (PowerShell on Windows, or a
shell such as zsh/bash on macOS or Linux). If MATLAB's `bin` directory is on
your shell's `PATH`, run:

```powershell
matlab -batch "run('Analysis_Macro.m')"
```

On macOS, MATLAB may be installed as an application without its command
being on `PATH`. Invoke the executable inside the application directly,
adjusting the release name to match your installation:

```sh
/Applications/MATLAB_R2024a.app/bin/matlab -batch "run('Analysis_Macro.m')"
```

Alternatively, override the configuration without editing it:

```powershell
matlab -batch "publicDataRoot='D:/downloaded_caltechdata'; outputRoot='D:/HE_coord_results'; run('Analysis_Macro.m')"
```

Forward slashes are recommended inside MATLAB path strings.
The override examples below use Windows paths; on macOS or Linux, replace
them with absolute paths appropriate to your system, such as
`/Users/yourname/data` or `/home/yourname/data`.

The command runs these steps in order:

1. loads the fixed pseudonymized session list from
   [`public_session_manifest.m`](public_session_manifest.m);
2. validates each processed input and runs all per-session analyses;
3. regenerates the across-session summary models and plots;
4. generates graphical-abstract components; and
5. copies the exact manuscript source PDFs into figure-specific folders.

## Run modes

[`Analysis_Macro.m`](Analysis_Macro.m) accepts variables set before `run(...)`:

| Variable | Default | Purpose |
|---|---:|---|
| `publicDataRoot` | repository-adjacent `data/public_release` | Root of the extracted CaltechDATA input tree |
| `outputRoot` | repository-adjacent `results/public_analysis` | Root for generated models and plots |
| `runSessionAnalyses` | `true` | Recompute per-session analyses and models |
| `runSummaryFigures` | `true` | Recompute across-session summaries |
| `runGraphicalAbstract` | `true` | Generate graphical-abstract source plots |
| `runPaperFigureCollection` | `true` | Copy manuscript sources into figure folders |
| `randomSeed` | `1` | Global MATLAB `twister` seed |

For a full clean reproduction, retain the four default `true` values.

### Regenerate summaries from existing session results

```powershell
matlab -batch "publicDataRoot='D:/downloaded_caltechdata'; outputRoot='D:/HE_coord_results'; runSessionAnalyses=false; runSummaryFigures=true; runGraphicalAbstract=true; runPaperFigureCollection=true; run('Analysis_Macro.m')"
```

### Refresh only the paper-figure folders

After all source PDFs exist:

```powershell
matlab -batch "outputRoot='D:/HE_coord_results'; Collect_Paper_Figure_Plots(outputRoot)"
```

[`Collect_Paper_Figure_Plots.m`](Collect_Paper_Figure_Plots.m) is copy-only
with respect to analytical outputs: it leaves each producer PDF in place and
copies the manuscript sources into managed figure folders. It also removes
obsolete PDFs from those managed folders so their contents remain exact.

## Analysis map

The entry point calls the following analysis blocks:

| Analysis | Per-session producer | Across-session/output producer |
|---|---|---|
| Additive reconstruction | `analyze_additivity_surfaces.m` | `plot_additivity_summary_across_sessions.m` |
| Population additivity geometry | `population_additivity_PCA_demo.m` | Per-session PDFs |
| Hessian/cross-derivative separability | `cross_derivative_test.m` | `plot_cross_derivative_rms_summary.m` |
| Effector selectivity | `compute_FR_amplitude_selectivity.m` | `plot_FR_amplitude_summary_across_sessions.m` |
| Nonparametric MAP decoding | `nonparametric_MAP.m` | `plot_MAP_confusion_summary_across_sessions.m`, `plot_MAP_condition_accuracy_6x6_summary.m` |
| Gaussianity diagnostics | `check_FR_gaussianity.m` | `summarize_FR_gaussianity.m` |
| Relative position | `analyze_relative_positition.m` | `plot_relative_position_summary.m` |
| Inseparable-neuron summaries | Uses outputs above | `plot_inseparable_neuron_summary.m` |

All implementation files are under [`utils/`](utils/).

## Where outputs are written

Every generated file is beneath `outputRoot`; the CaltechDATA download stays
unchanged beneath `publicDataRoot`:

```text
<outputRoot>/
    session_XX_<area>_results/
        models/     per-session numerical results
        plots/      per-session and individual-neuron PDFs
    models/         across-session pooled MAT files
    plots/
        graphical_abstract_components/
        summary/
            additive_reconstruction/
            separability_hessian/
            selectivity_generalizability/
            decoding/
            inseparable_neurons/
            relative_position/
            gaussianity_diagnostics/
        paper_figures/
            figure_1/
            figure_2/
            figure_3/
            figure_4/
            figure_s1/
            figure_s2/
            figure_s3/
            figure_s4/
            figure_s5/
            plot_manifest.tsv
```

The `summary` folders contain all population-level producer outputs. The
`paper_figures` folders contain only the source PDFs used in the paper.

## Finding the plots used in each paper figure

| Paper item | Output folder | Principal producing code |
|---|---|---|
| Figure 1 | `plots/paper_figures/figure_1/` | `analyze_additivity_surfaces.m` |
| Figure 2 | `plots/paper_figures/figure_2/` | `analyze_additivity_surfaces.m`, `plot_additivity_summary_across_sessions.m`, `plot_cross_derivative_rms_summary.m` |
| Figure 3 | `plots/paper_figures/figure_3/` | `compute_FR_amplitude_selectivity.m`, `plot_FR_amplitude_summary_across_sessions.m` |
| Figure 4 | `plots/paper_figures/figure_4/` | `plot_additivity_summary_across_sessions.m`, `plot_MAP_confusion_summary_across_sessions.m` |
| Supplemental Figure 1 | `plots/paper_figures/figure_s1/` | `plot_inseparable_neuron_summary.m` |
| Supplemental Figure 2 | `plots/paper_figures/figure_s2/` | `population_additivity_PCA_demo.m` |
| Supplemental Figure 3 | `plots/paper_figures/figure_s3/` | `check_FR_gaussianity.m`, `summarize_FR_gaussianity.m` |
| Supplemental Figure 4 | `plots/paper_figures/figure_s4/` | `plot_FR_amplitude_summary_across_sessions.m`, `plot_MAP_condition_accuracy_6x6_summary.m` |
| Supplemental Figure 5 | `plots/paper_figures/figure_s5/` | `plot_MAP_confusion_summary_across_sessions.m` |

The authoritative file-by-file mapping is the `spec` table in
[`Collect_Paper_Figure_Plots.m`](Collect_Paper_Figure_Plots.m). After a run,
the same mapping is available without MATLAB code in
`<outputRoot>/plots/paper_figures/plot_manifest.tsv`. The selected example
session, area, and neuron number are encoded in the copied filenames and
documented in the collector comments.

## Reproducibility details

- The pipeline initializes MATLAB's `twister` random-number generator with
  seed `1` unless `randomSeed` is explicitly overridden.
- Stochastic routines use fixed local or per-neuron seeds so their results do
  not depend on execution order.
- The session order and input names are fixed by `public_session_manifest.m`
  and `config/data_manifest.csv`.
- The tested software versions and required fields are recorded in
  `dependencies.yaml`.

## Troubleshooting

- **MATLAB cannot find an input:** place all 33 MAT files directly under
  `publicDataRoot` and compare their names with `config/data_manifest.csv`.
- **The collector reports a missing PDF:** run the complete pipeline with
  `runSessionAnalyses=true` and `runSummaryFigures=true` first.
- **An old output remains in a summary folder:** summary folders retain
  producer outputs from previous runs. The managed `paper_figures` folders
  are the definitive minimal source sets for the manuscript.
- **The shell cannot find `matlab` (`command not found` or not recognized):**
  add MATLAB's `bin` directory to the shell's `PATH` or invoke the executable
  by its full path (`matlab.exe` on Windows; see the macOS example above).
- **A different output location is used:** always set an absolute `outputRoot`.

## Data availability, citation, and license

Replace these placeholders when the records are public:

- **Processed data DOI:** `https://doi.org/10.22002/p50pc-49s61`
- **Processed data record:** `https://data.caltech.edu/records/p50pc-49s61`
- **Public code repository:** `https://github.com/nmynhier-caltech/HE_Coordination`

The code is released under the terms in [`LICENSE`](LICENSE).
