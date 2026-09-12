# LipiRich Changelog

> LipiRich is licensed under the [GNU Affero General Public License v3.0](LICENSE) (AGPL-3.0).

---

## [Unreleased]

---

## [0.9.0] — 2026-09-12

### Changed

#### "ISTD" renamed to "IS" throughout the app and README
Every user-facing label, help/validation/warning message, and the PCA plot's legend (labels, colours, caption) now say "IS" instead of "ISTD," matching the terminology already used for the `[IS]` metabolite-name tag. Internal-only identifiers (input IDs, output IDs, the internal-standard CSV's own dataframe column) are unchanged, since users never see these. The internal-standard amount/units CSV now accepts either an `IS` or `ISTD` column header, so existing CSV files (including `Demo data_IS_amounts.csv`, itself renamed from `Demo data ISTD_amounts.csv`) keep working without modification.

While tracing the rename, found that the MS-DIAL processing tutorial instructs users to name internal-standard-only injections `IS_pos`, but the sample-detection function only matched an `istd`/`itsd` substring — a sample named exactly per the tutorial's own instructions wouldn't have been recognised as IS-only. Extended `is_istd_sample()` (and the separate group-label exclusion regex) to also match a delimiter-bounded `IS` token (`IS_pos`, `IS_1`, bare `IS`), alongside the existing `ISTD`/`ITSD` patterns. Tested against the demo dataset's real sample names plus deliberately adversarial ones (`Fish_1_pos`, `Island_5`, `This_sample`) to confirm no false positives.

---

## [0.8.0] — 2026-09-12

### Added

#### Ion mode detection now falls back to adduct sign when sample names lack a `_pos`/`_neg` suffix
Previously, sample columns were only recognised if their header ended in `_pos` or `_neg` — the naming convention documented on the landing page and in the README. If none of a file's sample columns carried that suffix (a single-polarity upload, or sample names that simply didn't follow the convention), the suffix-based selection matched nothing and every sample column was silently dropped from `data_clean()`, with no data reaching any downstream tab.

`data_clean()` now falls back, when no suffixed columns are found, to identifying sample columns by elimination against MS-DIAL 5's fixed, documented alignment-result metadata columns (Alignment ID, Ontology, INCHIKEY, SMILES, spectrum references, and so on — 34 columns in total). Each value's ion mode is then resolved from its own row's Adduct type sign (already used elsewhere to set `ion.mode`) instead of from the column name. An early version of this fallback used a numeric-content heuristic instead of the metadata denylist; tested against this project's demo dataset, it incorrectly captured several genuinely numeric MS-DIAL metadata columns (Alignment ID, Fill %, Reference RT, Reference m/z, Annotation tag) as fake samples, so it was replaced with the denylist approach. A second edge case surfaced in the same testing: MS-DIAL's optional "Average/Stdev by class" summary columns reuse bare group names (e.g. `chow`) and would otherwise be indistinguishable from a genuine unsuffixed sample column; these are now excluded via the `...N` suffix `readr` applies when it deduplicates an exact-duplicate column name within one uploaded file — a pattern no real sample name ever produces, since sample names are unique to begin with.

Regression-tested against this project's demo dataset (both files, real suffixes intact) and a synthetic copy with every `_pos`/`_neg` suffix stripped: the fallback path selected exactly the same 29 sample columns as the suffix-based path (58 suffixed columns collapsing to 29 base names), and the two paths' final wide-format tables were confirmed identical across all 906 features × 29 samples, with zero MS-DIAL metadata columns leaking into either.

### Fixed

#### "Download long-format CSV" on the Export tab had no server-side handler at all
The button was wired up in the UI but `output$download_export_long` was never defined anywhere in the server logic, so clicking it did nothing. Fixed by factoring the tab's filtering/grouping logic (blank exclusion, iQC exclusion, class selection, grouping, percent-of-class conversion) out of `export_wide_data()` into a new shared reactive, `export_prepped_data()`, and adding `export_long_data()` and the missing `download_export_long` handler alongside it — following the project's existing convention of downstream tabs sharing reactives rather than duplicating logic.

#### "Include group column in wide export header" checkbox had no effect
`input$export_include_metadata` was read nowhere in the server code, so the checkbox did nothing regardless of its state. It now appends the group name to each sample's column header when checked (e.g. `sample1` → `sample1 (chow)`), since a wide table has one column per sample rather than one row per sample-group pair. The long-format CSV always includes group as its own column, checkbox notwithstanding.

#### Several Synthesis Pathway scores were matching against a lipid class that doesn't exist in the data
Five pathway scores (PG synthesis, TG storage, DGAT activity, PC synthesis, PE synthesis) used the literal string `"DAG"` as a denominator class to match against `plot_class`, but the app's actual class label for diacylglycerol is `"DG"` (confirmed elsewhere in the app's own default class lists). Since `"DAG"` never matched any real class, these denominators silently resolved to zero; given the score-computation code's zero-handling logic, this produced near-infinite/garbage ratio values rather than sensible numbers or a clean `NA`. Fixed by correcting the matching values, not just display labels, in `pathway_scores_defs()`.

### Changed

#### "Metabolite"/"metabolite name" replaced with lipid-specific terminology throughout the UI
MS-DIAL's own column is literally named "Metabolite name," but LipiRich is lipid-specific, so every place the app itself displays or labels this concept — plot titles, axis labels, dropdown labels, validation/warning messages, download filenames, and (found during this pass) a "metabolite" column header that was leaking directly into the Statistics results and post-hoc results tables — now reads "Lipid species" or "lipid" as appropriate. MS-DIAL's own "Metabolite name" field (the actual column read from the uploaded file, referenced 100+ times internally) and code comments describing it are unchanged, since renaming those would risk breaking data parsing for no user-visible benefit.

#### Tab renames
- "Export - wide" → "Export", since the tab has always held both the wide- and long-format CSV downloads (the landing page's Step 6 description already called it "Export").
- "Heatmap - significant" → "Heatmap" (the landing page's Step 9c description already called it "Heatmap").

#### "TAG"/"DAG" standardised to "TG"/"DG" in Synthesis Pathway labels
"TAG storage" → "TG storage"; every "…/DAG" label → "…/DG", matching the class-label convention used throughout the rest of the app. "CDP-DAG" is unchanged in three pathway names, since that's the standard biochemical name for that specific metabolic intermediate, not a class label.

### Docs

#### README: clarified ether-species pathway-score nomenclature
Added a note to the Synthesis Pathways section explaining that the plasmanyl/plasmenyl score notation (`PE-O XX:0`, `PE-O XX:≥1`) refers to the double-bond count on the ether-linked chain specifically, from MS-DIAL's molecular-species (chain-resolved) identification — not sum composition, which cannot distinguish plasmanyl from plasmenyl since it only reports the total double bonds across both chains. Also corrected a stale line claiming all scores use class-level totals with chain-resolved scoring "planned for a future release" — the plasmanyl/plasmenyl scores already require chain-resolved identification and exclude sum-composition-only species from their numerator/denominator.

Also updated the sample-naming section to describe the new adduct-based fallback (see Added, above), while still recommending the `_pos`/`_neg` suffix convention as the more robust choice when merging two files.

---

## [0.7.0] — 2026-09-01

### Added

#### Plot download buttons on three tabs that were missing them
The **IS Plots**, **Plot single lipid**, and **Plot lipids by class** tabs had no way to export their plots — only the interactive Plotly view, plus (for the latter two) a summary-table CSV download. All three now have the same **Export plot** panel (Width/Height/DPI/Scale fraction/Base font size numeric inputs with PNG and SVG buttons) used throughout the rest of the app. Each tab's plot-building code was factored out of its `renderPlotly` block into a standalone function (`.build_is_plot()`, `.build_met_plot()`, `.build_class_all_plot()`) shared by the on-screen render and both new download handlers, so the exported image always matches what's on screen — the same shared-function pattern already used for PCA, Volcano Plot, Class Bar Plots, and the other export-enabled tabs (and the same divergence bug that pattern was adopted to prevent, per the v0.6.2 fix below). Every other tab with a plot — Outlier Detection, PCA, Statistics, Volcano Plot, Class Bar Plots, Heatmap, Enrichment, Correlation Network, and Synthesis Pathways — already had PNG/SVG export and is unaffected.

### Fixed

#### IS Plots didn't follow Group Order — and lost its "individual samples" view along the way
IS Plots had its own independent, older grouping mechanism — a "How to group samples" control (`grouping_method_is`: By sample / Delimiter-based / Regex capture group, with its own delimiter/token/regex inputs) that resolved group labels completely separately from every other tab, via `resolve_group_labels()`. Since the app's Group Order feature is defined in terms of the *central* grouping reactive's group names (`.current_group_choices()`, which calls `get_active_grouping()`, driven by the Group Preview tab), IS Plots' independently-resolved group labels almost never matched an entry in that list, so `get_active_group_order()` had nothing to reorder and the plot silently fell back to alphabetical order regardless of a CSV `order` column or the Group Preview drag list.

Fixed by rewiring `is_plot_data()` to assign `group` via `get_active_grouping()`, the same central reactive used by `met_plot_data()`, `class_all_plot_data()`, PCA, Statistics, and every other tab — the same "dead parallel grouping UI" pattern documented as fixed for Synthesis Pathways in v0.0.2/v0.5.0, which IS Plots had never been migrated onto.

That rewire initially removed the tab's original "By sample" option entirely — every IS plot became a grouped (e.g. chow/fishoil/lard) view with no way back to one bar per sample. That behaviour (each individual sample as its own x-axis category, useful for QC — e.g. spotting a single bad injection) was legitimate functionality, not part of the disconnected-grouping bug; it had just been entangled with the same code. Restored as an explicit **"How to display samples"** choice on the IS Plots tab: **Individual samples** (the default, matching the tab's original behaviour) or **Grouped (from Group Preview)**. `is_plot_data()` sets `group` to the raw sample name in individual mode, or resolves it via `get_active_grouping()` in grouped mode. Group Order (in `.build_is_plot()` and the `isPlotSummary` table below) is applied only in grouped mode, since it has no group names to act on when every row is its own sample. The Group Preview group-selection note/mirror is shown only while grouped mode is selected.

Regression-tested against this project's demo dataset with `shiny::testServer()`, exercising the real reactive pipeline rather than reasoning from the code alone: the final `ggplotly()` widget's `layout.xaxis.categoryorder`/`categoryarray` (what plotly.js actually uses to decide bar order at render time) was confirmed to match the intended order under both a manually reordered Group Preview drag list and a grouping CSV's `order` column; individual mode correctly lists every sample (`chow_1`, `chow_13`, …) as its own category; and switching between individual and grouped mode reverts cleanly in both directions.

Individual mode initially left sample order plain alphabetical, since Group Order is defined over group names and individual mode's x-axis categories are sample IDs, not group names. `is_plot_data()` now also resolves and carries each sample's central Group Preview group as `central_group` (on both the points-level data and the per-sample summary), even in individual mode, purely for ordering purposes — the x-axis still shows one bar per sample either way. `.build_is_plot()` and `isPlotSummary` share a new `.is_group_levels()` helper: in grouped mode it applies Group Order to the group names directly as before; in individual mode it ranks each sample by its `central_group`'s position in Group Order (primary key) and alphabetically by sample name within a group (secondary key), so setting Group Order to lard → fishoil → chow clusters all `lard_*` bars together, then all `fishoil_*` bars, then all `chow_*` bars, instead of interleaving them alphabetically. Verified with `shiny::testServer()`: after setting a manual lard → fishoil → chow order, individual-mode plot and summary-table category order both group into exactly that sequence (with any sample group absent from the order list, e.g. `blank`, appended afterwards — the same "unseen groups" fallback `get_active_group_order()` already used), samples within each group remain alphabetically ordered, every sample still appears as its own bar (not collapsed), and grouped mode's own ordering is unaffected.

#### IS Plots summary table was completely blank
`DTOutput("isPlotSummary")` (the "Summary table (group-level IS statistics)" box under the IS plot) has been in the UI since IS Plots was first built, but no matching `output$isPlotSummary` render was ever defined server-side — a dangling output, same class of bug as the `sigHeatmapInfo` case fixed in v0.0.3. The box has always rendered empty regardless of data or settings. Added the missing render, sourced from the same `is_plot_data()$summary` used by the plot, with Group Order applied to its row order (in grouped mode) so the table and the plot above it agree.

---

## [0.6.2] — 2026-08-28

### Fixed

#### Synthesis Pathways score-statistics export didn't match the on-screen plot
The score-statistics lollipop plot ("Higher in [group]") had two independent implementations — one for on-screen display, one for PNG/SVG export. The v0.6.0 Group Order fix was applied only to the on-screen copy, so the exported image fell back to ggplot's default alphabetical ordering for the group legend/colours while the screen showed the correct Group Order. Rather than patch both copies, they're now consolidated into one shared function (`.build_path_plot()`); `pathStatsPlot` calls it directly, so there is exactly one implementation to keep in sync going forward. The pathway heatmap's on-screen/export pair was left as separate implementations, since its export path uses `pheatmap(..., silent = TRUE)` for manual device control — consolidating risked breaking the on-screen render for no benefit, since heatmap column order is always clustering-determined in both versions and was never affected by this bug.

---

## [0.6.1] — 2026-08-28

### Fixed

#### Synthesis Pathways stats plot crashed after the Group Order fix
`path_score_stats()` (feeding the Synthesis Pathways stats plot, `pathStatsPlot`) called `nzchar(group)` to drop unassigned samples. Since v0.6.0 made `group` an ordered factor (to support Group Order), `nzchar()` — which only accepts plain character vectors — threw `'nzchar()' requires a character vector` and crashed the plot. Fixed by coercing to character just for that check (`nzchar(as.character(group))`), leaving the factor and its ordering intact everywhere else in the pipeline. Also hardened the ANOVA branch's `direction` value to be explicitly character (it was inheriting factor-ness from a `group_by()`/`summarise()` result, while the t-test branch's `direction` was plain character — mixing the two across rows in `bind_rows()` is fragile depending on dplyr version). Audited every other consumer of `path_scores_long()`'s group column (heatmap pivot, table display, CSV download); none assume character type, so no other regressions from the same cause.

---

## [0.6.0] — 2026-08-28

### Added

#### Correlation network image export
The Correlation Network tab gains two ways to export a publication-quality image, alongside the existing edge-list CSV download. A one-click **"Export network image (PNG)"** button on the interactive widget itself captures the current on-screen arrangement (including any manual node dragging) as a screenshot, browser-resolution PNG. A new **"Export publication image"** panel in the sidebar renders a separate static image via `ggraph`/`igraph`, with the same width/height/DPI/scale controls and PNG/SVG output used elsewhere in the app (default 300 DPI, 2400×1800px). A new **"Match on-screen layout"** button reads the exact node positions from the interactive view and feeds them into the static export, so the two can be made to match node-for-node rather than relying on a different auto-layout algorithm; a status line confirms when the export is matched. New dependencies: `igraph`, `tidygraph`, `ggraph` (loaded on demand, consistent with `visNetwork`'s existing pattern).

### Fixed

#### Group Order not applied on several tabs
Four tabs resolved grouping independently from `get_active_grouping()` rather than through the shared reactives touched in v0.5.0, so the Group Order setting silently had no effect on them:
- **Class Bar Plots** — group bar order and colour-legend order now follow Group Order (`cbp_data()` / `.build_cbp()`).
- **Enrichment** — the "Group(s) to analyze" list, its default selection, and the facet order on the enrichment plot now follow Group Order (`enrich_base_df()`, `available_groups_for_enrich()`, `lsea_plot_obj()`).
- **Synthesis Pathways** — group order in the score statistics (t-test/ANOVA group labelling, fold-change direction), the "Higher in" legend, and group-mean heatmap columns now follow Group Order (`path_scores_long()`, `path_score_stats()`, `pathStatsPlot`).
- **PCA, Statistics, single-metabolite plot, class-all plot, and Correlation Network group selectors** — all five were populated by a single shared helper (`.refresh_available_groups()`) that sorted alphabetically regardless of Group Order; now resolved through `get_active_group_order()`.

#### FGSEA results not reproducible between identical runs
`fgsea::fgseaMultilevel()` estimates its null distribution via internal Monte Carlo permutation and was never seeded, so re-running an unchanged one-vs-rest comparison (e.g. triggered by toggling an unrelated group in the selector) produced slightly different NES/p-values each time, even though the enrichment score (ES) itself was always identical and correct. A fixed seed is now set immediately before the FGSEA call. ORA was unaffected (closed-form hypergeometric test, no randomness) and was confirmed bit-for-bit identical across runs.

#### Enrichment plot ranking sets from a global rather than per-group rank
`reorder(label, NES)` (and the ORA equivalent) ranked each lipid set **once globally**, using the mean value across every group sharing that set name — so with multiple groups faceted together, a group's panel could show a ranking blended with other groups' values rather than its own. Sets are now ranked independently within each group's own facet, via a label+group composite key (`forcats::fct_reorder`), with the group suffix stripped back off for display. Table and CSV outputs were never affected, only the plot's visual ranking.

### Changed

#### Group Order UI safety net
The drag-to-reorder list in the Group Preview tab's Group Order section is now wrapped in `tryCatch()`, so a rendering failure surfaces a visible error message instead of leaving the section silently blank.

---

## [0.5.0] — 2026-08-20

> **Note:** this project's `app.R` was already at v0.4.0 when this entry was written, with no corresponding CHANGELOG entry — a second undocumented gap alongside the known 0.1.0→0.2.0 gap. This entry covers only the changes made in the session that produced v0.5.0; the 0.3.0→0.4.0 changes are not retroactively documented here.

### Added

#### Manual group plotting order
Groups can now be manually reordered for plotting, via a new **Group Order** section in the Group Preview tab. Two methods are supported: a numeric `order` column in the grouping CSV (the smallest value found among a group's samples sets that group's position), or a drag-to-reorder list for when no grouping CSV is in use. Priority is CSV `order` column → drag list → default (first-appearance) order. The resolved order is applied consistently to group axes in bar plots (single-metabolite and class bar plots), group-mean heatmap columns, and post-hoc significance annotation matching. Individual sample order within a group (e.g. sample-level heatmap columns) is unaffected and remains alphabetical.

### Changed

#### Statistical plotting steps regrouped under "Plotting"
The Volcano Plot, Class Bar Plots, and Heatmap steps are renumbered from Step 8a / 9 / 9b to **Step 9a / 9b / 9c** respectively, and grouped under a new "Plotting" sub-heading on the landing page workflow diagram and as a label above each tab's header. Step 8 (Statistics) and Steps 10–12 (Enrichment, Correlation Network, Synthesis Pathways) are unchanged.

---

## [0.3.0] — 2026-08-01

### Added

#### Cross-mode identity reconciliation
For the glycerophospholipid classes (PC, PE, PG, PI, PS, PA) and cardiolipin (CL), LipiRich can now pair each species across ion modes by sum composition and retention time, taking quantitation from the adduct-matched positive-mode feature and the acyl-resolved identity from the negative-mode feature. A species is quantified once rather than double-counted across modes. Unmatched features are never dropped: a positive-only feature is shown at the depth positive mode can determine (sum composition, or the two combined halves for CL), a negative-only feature keeps its full acyl identity, and each is quantified from its own-mode internal standard. The behaviour, the reconciled class list, and the retention-time tolerance are configurable in the **Settings** tab (enabled by default); turning it off reverts those classes to prefer-negative handling. Internal standards are never reconciled.

#### Cross-mode audit tab
A new **Cross-mode audit** tab lists every reconciliation decision — positive features re-identified from negative mode, positive features retained at sum composition (or positive substructure for CL), and negative-only features retained and negative-quantified — with the shorthand, matched retention times and ΔRT, the final identity carried downstream, and the ion mode each abundance is quantified from. Fully traceable and reloads with the data.

#### Plot single lipid — ion-mode selector and cascading filters
The Plot single lipid tab gains an **Ion mode** control (All / Positive / Negative). The ion-mode, adduct, and metabolite-name selectors now cascade: the adduct list reflects the chosen class and mode, and the name list only offers species that have data under the current mode and adduct. A selection can therefore no longer produce an empty plot.

#### Demo dataset
Added a full demonstration dataset to `demo_data/`, generated from mouse liver samples (n = 8 per group: chow, lard-based high-fat diet, 90% fish oil high-fat diet):

- `Demo_data_pos.txt` / `Demo_data_neg.txt` — MS-DIAL 5 aligned output, positive and negative ion mode, ready for direct upload
- `Demo_data_metadata.xlsx` — combined reference workbook (sample order/analytical run, protein content, internal standard amounts)
- `Demo_data_protein.csv` / `Demo_data_ISTD_amounts.csv` — protein and internal standard files, formatted for direct upload to LipiRich
- `Demo_data_MSDIAL_processing_parameters_pos.mdparameter` / `..._neg.mdparameter` — MS-DIAL 5 parameters used to generate the aligned output, for reproducing the alignment from raw data

The original LC-MS `.raw` files are archived separately on Zenodo (DOI: 10.5281/zenodo.21448733) due to size. README updated with a full description of the demo dataset and repository structure.

### Changed

#### IS Plots default to both ion modes, with standard names on hover
The Internal Standards plot now defaults to showing positive and negative standards side by side (**Both (separate)**) rather than a single mode, and each point's hover leads with the exact internal-standard name — useful for identifying which standard is responsible for an anomalous point, particularly in the pooled **Unknown** class.

#### Landing page
The front page has been redesigned (modern card layout and typography) and its documentation updated to describe cross-mode reconciliation, the audit tab, the single-lipid cascade, and the IS-plot changes.

#### Environment
Now developed and tested with **R 4.6.1** and **Bioconductor 3.23** (previously R 4.5.2 / Bioconductor 3.22). `packages.R` and the Dockerfile base image updated accordingly.

### Fixed

- **Cross-mode retention-time matching** no longer errors when `Average Rt(min)` is imported as a character column; it is coerced to numeric for the RT arithmetic while the source column type is left unchanged.
- **`case_when()` deprecation (dplyr 1.2.0):** the Class Bar Plots point fill and outline, which selected a vector palette from a scalar mode setting, were rewritten as `if`/`else` to avoid the deprecated scalar-LHS / vector-RHS pattern.
- **`.data[[ ]]` in tidyselect deprecation (tidyselect 1.2.0):** the Correlation Network `select()` and `pivot_wider()` now use `all_of()` instead of the `.data` pronoun in tidyselect position.
- **plotly `layout(height=)` deprecation:** the single-lipid plot height is now passed to `ggplotly()` directly.
- **Unknown-column warning:** `norm_units` is guarded in the single-lipid plotting path so a filtered slice or non-normalised value type no longer triggers an "unknown or uninitialised column" warning.
- **Correlation Network:** zero-variance (constant) features are dropped before `cor()`, avoiding a "standard deviation is zero" warning and spurious `NA` columns.

---

## [0.2.2] — 2026-07-29

### Added

#### Outlier Detection (Step 3b) — plot export
The Sample-level PCA Hotelling's T² and iQC replicate deviation plots now have the same standardised Width/Height/DPI/Scale/Base font size export controls with PNG/SVG buttons used throughout the rest of the app. Both plots share one set of export controls; each build was factored into its own `.build_outlier_sample_plot()` / `.build_outlier_iqc_plot()` function (parameterised by font size) so the on-screen render and both download handlers draw from a single source, matching the pattern already used for PCA and Class Bar Plots.

#### Synthesis Pathways — direct plasmanyl-PC synthesis score
Added `plasmanyl-PC synthesis, direct (PC-O XX:0/DG-O)` alongside the existing headgroup-conversion score `plasmanyl-PC synthesis (PC-O XX:0/PE-O XX:0)`. Mammalian plasmanyl-PC can be synthesised either directly from the shared DG-O ether precursor via the CDP-choline branch, or via headgroup conversion from plasmanyl-PE — the two scores now let both routes be inspected separately. Plasmenyl-PC is unaffected and remains scored only via headgroup conversion from plasmenyl-PE, since the mammalian plasmalogen desaturase (PEDS1/TMEM189) is PE-selective and there is no direct DG-O route to plasmenyl-PC.

#### Synthesis Pathways — interpretation note
Added a note under the Step 12 header flagging that these scores are based on established mammalian synthesis pathways, some of which involve enzymes with overlapping or tissue-dependent substrate preferences, and should be read as pathway-activity indicators rather than direct single-step flux measurements.

---

## [0.2.1] — 2026-07-27

### Added

#### Idle-session timeout
Sessions now track user activity (mouse movement, clicks, keyboard, scroll, touch) client-side. After 15 minutes of inactivity a warning modal appears asking the user to confirm they're still there; if no further activity follows, the session closes automatically at 20 minutes via `session$close()`, freeing the server process's memory held by an abandoned tab. Any tracked activity — including dismissing the warning — resets both timers. Implemented entirely within `app.R` (`idle_timeout_js` in the UI, paired `observeEvent` handlers server-side); `shiny-server.conf`'s `app_idle_timeout` is unchanged and continues to govern process-level cleanup only (killing the R process once zero sessions remain connected), which is separate from this per-session inactivity check.

---

## [0.1.0] — 2026-07-20

### New features

#### Protein Match tab (Step 3a)
A new diagnostic tab checks whether every imported MS-DIAL sample has a corresponding row in the uploaded protein content CSV, matched on the same normalised sample name used by the actual protein-normalisation join (`_pos`/`_neg` suffix stripped, trimmed, lower-cased). Reports which samples matched, which are missing a protein value (and so will *not* be protein-normalised), and which protein CSV rows don't correspond to any imported sample (a common sign of a typo). ISTD, Blank, and iQC/QC samples are excluded from the check on both sides, since they are not expected to have protein measurements.

#### Outlier Detection tab (Step 3b)
A combined outlier detection and exclusion workflow, positioned upstream of Plot single lipid so that exclusions apply to every downstream analysis (Plot single lipid, Class Bar Plots, Export, PCA, Statistics, Volcano, Heatmap, Enrichment, Correlation Network, Synthesis Pathways):

- **Sample-level — PCA Hotelling's T²**: flags samples beyond a user-selected confidence threshold (95/97.5/99%) on a self-contained PCA (log-transformed, unit-variance scaled, IS-normalised, non-blank/non-iQC/non-ISTD samples).
- **Sample-level — iQC replicate deviation**: for datasets with ≥ 3 iQC replicates, flags iQC samples whose median relative deviation from the cross-replicate median exceeds a MAD-based threshold. Informational only — iQC/ISTD/Blank samples are never excluded from downstream analyses regardless of this flag.
- **Feature-level — per lipid × group**: flags individual sample values within a Metabolite name × group combination using either a modified Z-score (MAD, robust to skew) or the classic IQR rule, both with adjustable thresholds.
- **Bulk exclude toggles** (off by default) apply all auto-flagged samples/points to the shared dataset used by every other tab. Excluded feature points are set to missing rather than deleted; excluded samples are removed from the shared reactive but never from the underlying file.
- **Individual review & manual override**: a selectable "Sample review" table (T² score per sample) and a selectable feature-level candidates table let specific samples/points be individually excluded or kept, independently of (and taking precedence over) the bulk toggles.
- ISTD, Blank, and iQC/QC samples are excluded from the entire outlier detection and exclusion workflow — they are never flagged, never excluded, and never contribute to the candidate pool.

#### PCA — iQC and ISTD handled as supplementary (projected) individuals
The PCA fit (axes, loadings, and centring/scaling statistics) is now computed from biological samples only. iQC and ISTD-only samples, if included, are fitted as FactoMineR supplementary individuals — projected into the fitted space afterward without being able to distort the ordination. Since neither type is typically present in the protein CSV, they are scaled by the *median* protein content of the biological samples (an assumed, not measured, value) purely so their projection lands somewhere visually comparable. Two independent toggles control inclusion: **Include iQC samples** (on by default) and **Include ISTD-only samples** (off by default, since an ISTD-only injection has no biological matrix and a very different lipid profile to a real sample). The plot legend distinguishes real groups, iQC, and ISTD with distinct colours; supplementary points are shown as triangles.

#### PCA — plot export and sample ID labels
The PCA tab now has the same standardised Width/Height/DPI/Scale/Base font size export controls with PNG/SVG buttons used throughout the rest of the app, plus a **Show sample ID labels** toggle (off by default, using ggrepel to avoid overlap).

### Bug fixes

#### iQC and ISTD-named samples were being conflated throughout the app
The internal `is_iqc_sample()` helper matched `iqc`, `qc`, `istd`, and `itsd` all as equivalent. A sample named ISTD (an internal-standard-only injection with no biological matrix) was therefore being treated identically to a pooled iQC replicate everywhere in the app, including the PCA tab — where its very different lipid profile appeared as a spurious dominant point and inflated apparent variance. Split into `is_iqc_sample()` (iQC/QC only), a new `is_istd_sample()` (ISTD/ITSD only), and `is_qc_type_sample()` (their union, used wherever the existing broad "not a real biological sample" behaviour was correct, e.g. `filter_iqc()`, outlier detection, and the Protein Match tab).

#### PCA plot crash: "number of active individuals is different from the length of the factor habillage"
`factoextra::fviz_pca_ind()` requires the `habillage` grouping factor to be sized to active individuals only; supplementary individuals must be styled separately. The PCA plot no longer uses `fviz_pca_ind()` — it is built directly from `pca_result()$ind$coord` and `pca_result()$ind.sup$coord`, which also removes an unexplained extra point that had been appearing on the plot (traced to the iQC/ISTD conflation above).

#### Backtick-quoted column name with a `\u` unicode escape
A column name written as `` `T² score` `` failed to parse — R does not support `\u` unicode escapes inside backtick-quoted names, only inside string literals. Replaced with the literal `²` character.

### Improvements

#### Step numbering
Outlier Detection is now Step 3b (moved from its original position after PCA Analysis to upstream of Plot single lipid, so its exclusions can propagate to every downstream tab). Protein Match remains Step 3a. All other step numbers are unchanged.

---

## [0.0.4] — 2026-06-21

### New features

#### Volcano plot — pairwise comparisons for ANOVA data
When the Statistics tab has run a one-way ANOVA (≥ 3 groups), a **Comparison** selector now appears in the Volcano Plot sidebar. Available pairwise comparisons are populated from `posthoc_for_plot()`. Selecting a comparison computes the pairwise log2FC from group means and the pairwise adjusted p-value via `.compute_pairwise`, then renders a standard volcano for that pair. The plot subtitle shows the active comparison. When "Overall (ANOVA)" is selected, the behaviour reverts to the existing overall F-test volcano (for datasets where log2FC is available).

#### Volcano plot — axis capping for outlier control
A new **Cap axes to reduce outlier distortion** checkbox (off by default) allows the -log10(p) and log2FC axes to be clamped to user-defined limits. Features beyond the cap are plotted *at* the cap value as filled triangles (▲) rather than being hidden or causing axis compression that makes the remainder of the plot unreadable. Hover text for capped points is annotated with `[axis-capped for display]` so the true values remain accessible interactively. The shape legend (circle = normal, triangle = capped) only appears when capping is active.

#### Volcano plot — improved label rendering
Labels for significant features are now rendered differently for the interactive and export contexts:

- **Interactive (plotly):** ggrepel is not supported by plotly's ggplot converter and was silently dropped. Labels are now added as native plotly annotations via `add_annotations()`, with arrow offsets fanned radially so labels near plot edges are not pushed off-screen.
- **PNG/SVG export:** ggrepel labels are retained for static export with `max.overlaps = Inf` (labels are never silently dropped), expanded `xlim`/`ylim` boundaries to give the repulsion algorithm room to place labels near plot edges, and increased `force` and `force_pull` values. When axis capping is active, ggrepel boundaries respect the cap limits.

#### Statistics tab — per-class detection and significance summary panel
A collapsible `<details>` panel has been added to the Statistics tab showing a summary table of results broken down by lipid class. For each class the table reports: number of species detected, number tested, number significant, and percentage significant. The full table is downloadable as CSV.

---

### Bug fixes

#### Volcano plot — t-test data crashing with `p_adj` not found
The `volcano_ggplot` reactive was setting `p_col = "p_adj"` then immediately renaming the column to `p.adj`, leaving `p_col` pointing to a non-existent column. The rename has been removed; the t-test path now uses `p_adj` consistently throughout.

#### Volcano plot — pairwise view showed only significant features
`.volcano_pairwise_data` was sourcing p-values from `posthoc_for_plot()`, which internally filters to `p < alpha` before returning. Features that did not pass the significance threshold had no matching row to join to and were silently dropped. Fixed by adding a `filter_alpha` parameter to `.compute_pairwise` (default `= alpha` preserves existing behaviour for all other callers; pass `Inf` to bypass the filter). `.volcano_pairwise_data` now calls `.compute_pairwise` directly on `stats_input_long()` with `filter_alpha = Inf`, independently of `posthoc_for_plot()`.

---

### Improvements

#### Tab step headings corrected
All tab `h4` step headings updated to reflect the current step numbering: Volcano Plot is now labelled **Step 8a** and Correlation Network is labelled **Step 11**.

#### Landing page updated with Step 8a and Step 9b cards
The landing page feature cards now include entries for **Step 8a — Volcano Plot** and **Step 9b — Heatmap (significant features)**, which were previously missing.

---

### Other

#### Copyright year updated to 2025–2026
File header, AGPL licence block, and citation text updated to reflect the 2026 release year.

---

## [0.0.3] — 2026-06-06

### New features

#### Class Bar Plots tab (new, Step 9 — between Volcano Plot and Heatmap)
A dedicated tab for visualising all lipid species within a selected class as a single grouped bar plot with dodged bars per group. Features include:

- **Class selector** — populated from the loaded data as soon as data is uploaded; shows all detected lipid classes.
- **Value type** — IS-normalised or background-subtracted; units selector for y-axis label.
- **Abundance range filter** — percentage slider filters species by mean abundance relative to the class maximum, enabling focus on high-, mid-, or low-abundance species independently.
- **Significant species only toggle** — now correctly uses `stats_results_all` (all tested features) so the checkbox meaningfully restricts to species that passed the significance threshold, rather than the pre-filtered result set.
- **Significance brackets** — pairwise comparisons drawn as annotated bracket segments directly on the plot. For t-test (2 groups): a single bracket spanning the two bars. For ANOVA (≥3 groups): one stacked bracket per significant Tukey HSD pairwise comparison. Brackets show `*`, `**`, or `***` symbols.
- **Point overlay** — individual data points dodged to align with their group bar; user-configurable fill (match palette / black / white) and outline (white / black / none).
- **Colour palettes** — 18 options across qualitative (Okabe-Ito default, Dark2, Set1, Set2, Paired, Accent, Tableau 10), sequential (Viridis, Plasma, Inferno, Blues, Greens, Purples), and diverging (RdBu, PuOr, BrBG) schemes.
- **Bar opacity** and **point size** numeric controls.
- **Species ordering** — alphabetical (default) or by descending mean abundance.
- **PNG and SVG export** with the standard Width/Height/DPI/Scale/Font size controls.

#### Shared pairwise comparison helper `.compute_pairwise`
Extracted a shared internal function used by both the Statistics tab bar plots (`posthoc_for_plot`) and the Class Bar Plots tab (`cbp_data`). Eliminates code duplication, ensures consistent pairwise test logic (Welch t-test for 2 groups; Tukey HSD for ≥3 groups via rstatix with base R fallback), and guarantees the Class Bar Plots tab always shows brackets regardless of what posthoc settings are active on the Statistics tab.

#### Heatmap — t-test (2-group) data enabled
Column clustering is now automatically disabled when the heatmap matrix has fewer than 3 columns (i.e. t-test / 2-group results), preventing a pheatmap clustering error. A note is displayed in the heatmap panel when this applies.

#### Volcano plot — label repulsion
Significant feature labels on the volcano plot now use `ggrepel::geom_text_repel` instead of `ggplot2::geom_text`, automatically repositioning overlapping labels with connector lines. `coord_cartesian(clip = "off")` ensures labels pushed outside the panel boundary are not clipped.

#### Colour palette options — heatmap tabs
Both the Heatmap and Synthesis Pathways → Scores heatmap tabs now include a Colour palette dropdown with 10 options: Viridis (default), Magma, Plasma, Inferno, Cividis, Rocket, Mako, Turbo, Blue–White–Red (diverging), Green–White–Purple (diverging).

#### Row/column label toggles — Heatmap tab
Two checkboxes — **Show row labels** and **Show column labels** — allow labels to be independently hidden.

---

### Bug fixes

#### Synthesis Pathways — TAG storage score corrected
The "TAG storage" score was computing DAG/TG (DAG in numerator). Corrected to TG/DAG (TG in numerator) — a higher score now correctly reflects greater TAG storage relative to the DAG precursor pool. The score name is updated to `TAG storage (TG/DAG)`.

#### Downstream tabs respect Statistics group selection
The Enrichment (`enrich_base_df`) and Synthesis Pathways (`path_totals_extended`) tabs now correctly filter to only the groups selected on the Statistics tab, matching the behaviour of the Statistics bar plots, heatmap, and volcano.

#### Correlation Network group dropdown
Wired into the same central `net_available_groups` reactiveVal pattern used by PCA and Statistics, and restricted to groups included in the current statistics run. Resolves a race condition where the dropdown populated before group tokens were configured.

#### Volcano plot — all features shown
The volcano plot previously only showed features that passed the significance filter. Now reads from `stats_results_all` (always full unfiltered results), so all tested features appear regardless of the Statistics tab's significance filter setting.

#### Volcano plot — PNG/SVG export producing HTML files
Export handlers were calling a stale `volcano_plot_obj()` reactive that returned a plotly object. Fixed by extracting a `volcano_ggplot` reactive (pure ggplot2) used by both the interactive display and download handlers.

#### Statistics tab group tokens not inherited by downstream tabs
Synthesis Pathways and Correlation Network tabs now correctly reflect the active token selection from the Group Preview tab. The Synthesis Pathways dead parallel grouping UI (method selector, delimiter, token inputs) has been removed.

#### sigHeatmapInfo dangling output
Fixed a textOutput in the Heatmap tab that had no corresponding server render.

---

### Improvements

#### Tab renamed: "Lipid Network" → "Correlation Network"
Clarifies that the network reflects Pearson correlation, not metabolic pathway connectivity.

#### Landing page step numbering updated
Step 9 is now Class Bar Plots; previous steps 9–11 renumbered to 10–12.

#### Plot export controls standardised
All visualisation tabs (Statistics bar plots, Volcano, Heatmap, Enrichment, Synthesis Pathways, Class Bar Plots) now use consistent Width (px) / Height (px) / DPI / Scale fraction / Base font size controls with PNG + SVG buttons.

---

## [0.0.2] — 2026-06-03

### Summary
Major UI, statistical, and correctness improvements across the export, visualisation, and analysis tabs. This release replaces the underpowered class-level ORA in the Synthesis Pathways tab with a direct score-level group comparison, fixes the volcano plot to show all tested features, corrects group-filtering across downstream tabs, and standardises plot export controls throughout the app.

---

### Licence change
Licence changed from MIT to **GNU Affero General Public License v3.0 (AGPL-3.0)**. The AGPL-3.0 extends copyleft protection to network services: anyone who deploys a modified version of LipiRich as a hosted web application must make the corresponding source code available to users. This better reflects the web-app delivery model of LipiRich and protects against commercial rebranding without attribution.

### New features

#### Standardised plot export controls (all visualisation tabs)
All tabs with downloadable plots now use a consistent set of export controls matching the MetaboRich style: **Width (px)**, **Height (px)**, **DPI**, **Scale fraction**, and **Base font size** numeric inputs, with side-by-side **PNG** and **SVG** download buttons. Affected tabs: Statistics (bar plots), Volcano Plot, Heatmap, Enrichment, Synthesis Pathways. The Statistics tab retains an additional batch **PDF** export for all significant features.

#### Colour palette selector (heatmap tabs)
Both the **Heatmap — significant features** tab and the **Synthesis Pathways → Scores heatmap** subtab now include a **Colour palette** dropdown with 10 options: Viridis, Magma, Plasma, Inferno, Cividis, Rocket, Mako, Turbo, Blue–White–Red (diverging), and Green–White–Purple (diverging). The selected palette is applied to both the live plot and all exported images.

#### Row/column label toggles (Heatmap tab)
Two new checkboxes — **Show row labels (features)** and **Show column labels (samples)** — allow labels to be hidden independently. Both default to on.

#### Synthesis Pathways — Score Statistics subtab (replaces ORA)
The **ORA enrichment** subtab has been replaced with a **Score Statistics** subtab. Rather than running Fisher's exact test over lipid class sets (which was underpowered due to the small universe size and produced adjusted p-values of 1 for most datasets), the new subtab runs a **t-test** (2 groups) or **one-way ANOVA** (≥ 3 groups) directly on each pathway score's sample values, with BH correction applied across all scores. Results are displayed as:
- A **lollipop plot** of −log10(p_adj) per score, coloured by which group has the higher mean, with a dashed significance threshold line.
- A **results table** showing score, test type, groups, test statistic, log2FC (t-test only), direction, p, and p_adj, with significant rows highlighted.

#### Volcano plot — all tested features now shown
Previously the volcano plot only displayed features that passed the significance filter set in the Statistics tab (controlled by the **Show all rows** toggle). The volcano now always shows all tested features regardless of this setting, by drawing from a new `stats_results_all` reactiveVal that is always populated with `show_all_rows_flag = TRUE`. The bar plots, heatmap, and enrichment tabs continue to respect the significance filter as before.

---

### Bug fixes

#### Volcano plot exported as HTML (fixed)
The PNG and SVG download buttons on the Volcano Plot tab were calling a stale `volcano_plot_obj()` reactive from an earlier version of the app, which returned a plotly/htmlwidget object. Passing this to `ggsave` caused it to serialise as an HTML file. Fixed by extracting a dedicated `volcano_ggplot` reactive (pure ggplot2, no `ggplotly()`) shared by both the `renderPlotly` display and the download handlers.

#### Downstream tabs computing results for unselected groups (fixed)
When a subset of groups was selected on the Statistics tab, the **Enrichment** and **Synthesis Pathways** tabs were still computing results for all groups in the dataset. Fixed by filtering `enrich_base_df` and `path_totals_extended` (the base data reactives for those tabs) to only include samples belonging to `stats_selected_groups`.

#### Correlation Network group dropdown not reflecting token selection (fixed)
The **Compute correlations within group** dropdown was being populated by an observer that fired at data load — before group tokens were set in the Group Preview tab — and was not reliably updated when tokens changed. Fixed by wiring `net_available_groups` into the same central `observe({ .refresh_available_groups() })` pattern used by the PCA and Statistics group selectors, which correctly waits for both data and token configuration. The dropdown now also restricts to groups included in the current statistics run.

#### Synthesis Pathways showing all groups despite Statistics group selection (fixed)
`path_class_totals` had the group filter applied but `path_totals_extended` — the reactive actually used by `path_scores_long` and all downstream outputs — did not. The `stats_selected_groups` filter is now applied in both reactives.

---

### Improvements

#### Tab renamed: "Lipid Network" → "Correlation Network"
The tab name and landing page description have been updated to clarify that the network is based purely on Pearson correlation among significant lipid species, not on metabolic pathway connectivity.

#### Synthesis Pathways grouping UI simplified
The Synthesis Pathways tab previously had its own parallel grouping UI (method selector, delimiter, token count, regex inputs) that was never actually read by any reactive — grouping was always sourced from `get_active_grouping()`. The dead UI has been removed and replaced with the same read-only token mirror (`uiOutput`) used by the PCA and Statistics tabs, directing users to the Group Preview tab to configure grouping.

#### Correlation Network grouping mirror added
The Correlation Network tab now shows the same active token summary as the PCA, Statistics, and Synthesis Pathways tabs.

#### Heatmap colour scheme updated to Viridis
The previous blue–yellow–red (`colorRampPalette`) colour scheme has been replaced with Viridis as the default, with additional palette options available via the new selector (see above).

#### Landing page updated
- Step 10 title corrected from "Lipid Network" to "Correlation Network".
- Step 11 description updated to reflect the score statistics approach rather than ORA.

---

### Internal / code quality

- Added `stats_results_all` reactiveVal alongside `stats_results_val`; the volcano plot uses the former (always full results) and all other downstream outputs use the latter (respects significance filter).
- Added shared `.heatmap_palette(palette_id, n)` helper function used by all five pheatmap colour calls, eliminating duplicated palette construction code.
- Removed dead `class_ovr_stats` function and `path_ora_results`, `pathOraPlot`, `pathOraTable` server blocks.
- All plot export dimension helpers (`.stats_export_dims`, `.hm_export_dims`, `.lsea_export_dims`, `.volc_export_dims`, `.path_export_dims`) use a consistent px → inches conversion: `w = (px / dpi) * scale`.
- Version strings now maintained in four consistent locations: file header comment, `# Version:` line, `APP_VERSION` variable, and citation text on the landing page.
