# LipiRich Changelog

---

## [0.0.2] — 2026-06-03

### Summary
Major UI, statistical, and correctness improvements across the export, visualisation, and analysis tabs. This release replaces the underpowered class-level ORA in the Synthesis Pathways tab with a direct score-level group comparison, fixes the volcano plot to show all tested features, corrects group-filtering across downstream tabs, and standardises plot export controls throughout the app.

---

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
