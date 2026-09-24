# Processed public data dictionary

The CaltechDATA release contains 33 MATLAB files in one flat directory. File
names follow this pattern:

```text
session_XX_<area>_Processed_data.mat
```

`session_XX` is a pseudonymized session identifier. `<area>` is `M1` (motor
cortex), `SPL` (superior parietal lobule), or `all`; `all` contains the
combined recorded population for that session.
Every file contains one scalar MATLAB struct named `data` with the five fields
below. All fields are numeric `double` arrays, and their rows are aligned.

| Field | Dimensions | Meaning |
|---|---|---|
| `FR` | `T x N` | Firing rate in Hz. Each row is one analyzed trial and each column is one retained single neuron. `N` varies across sessions and areas. Rates are calculated for the configured analysis window during the selected task phase. |
| `TP1` | `T x 2` | Two-dimensional hand/cursor target position for each trial, in task-coordinate units. The two columns are the planar coordinates used to calculate hand direction with `atan2(TP1(:,2), TP1(:,1))`. |
| `TP2` | `T x 2` | Two-dimensional gaze target position for each trial, in task-coordinate units. The two columns are the planar coordinates used to calculate eye direction with `atan2(TP2(:,2), TP2(:,1))`. These are task target labels, not continuous eye-tracker samples. |
| `TPi1` | `T x 1` | Integer hand/cursor target-condition index. Values `1` through `6` identify the six movement directions; `7` indicates that the hand/cursor effector is inactive or at the center condition. |
| `TPi2` | `T x 1` | Integer gaze target-condition index. Values `1` through `6` identify the six movement directions; `7` indicates that the gaze effector is inactive or at the center condition. |

Here, `T` is the number of analyzed trials in that file and `N` is the number
of retained neurons. `T` and `N` therefore vary between files.

The condition-index fields identify the trial type:

| Trial type | Condition rule |
|---|---|
| Coordinated hand-eye movement | `TPi1 ~= 7` and `TPi2 ~= 7` |
| Hand-only movement | `TPi1 ~= 7` and `TPi2 == 7` |
| Eye-only movement | `TPi1 == 7` and `TPi2 ~= 7` |

The public files do not contain acquisition dates, raw file names, electrode
identifiers, continuous eye-tracking samples, or the indices of units removed
during preprocessing. Private provenance remains in the lab-only Stage 1
output.
