# LipiRich Changelog

> LipiRich is licensed under the [GNU Affero General Public License v3.0](LICENSE) (AGPL-3.0).

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