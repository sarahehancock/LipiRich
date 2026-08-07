# LipiRich <img src="https://img.shields.io/badge/version-0.4.0-blue" alt="v0.4.0"/> <img src="https://img.shields.io/badge/license-AGPL--3.0-blue" alt="AGPL-3.0"/> <img src="https://img.shields.io/badge/R-%3E%3D4.6.1-informational" alt="R 4.6.1"/> <img src="https://img.shields.io/badge/live%20app-lipirich.sarahehancock.com-brightgreen" alt="Live App"/>

**LipiRich** is an open-source, browser-based Shiny application for the normalisation, statistical analysis, and visualisation of untargeted lipidomics data exported from [MS-DIAL 5](https://systemsomicslab.github.io/compms/msdial/main.html). It requires no programming knowledge and runs entirely in a web browser.

> Developed and tested with **MS-DIAL 5.5.251021**, R 4.6.1, and Bioconductor 3.23.

> ⚠️ **Pre-publication software (v0.4.0):** LipiRich is under active development. A citable preprint and demonstration dataset will be released alongside v1.0.0. Please check the [GitHub repository](https://github.com/sarahehancock/LipiRich) for the latest updates and to report issues.

---

## Features

| Module | Description |
|---|---|
| **Data import** | Upload 1–2 MS-DIAL aligned `.txt` files; positive and negative ion modes merged automatically |
| **Ion mode deduplication** | Per-class mode preference rules (user-configurable) ensure each species is represented by its most informative ion mode |
| **Cross-mode identity** | For glycerophospholipids and cardiolipin, optionally quantify from the positive-mode feature while transferring the acyl-resolved identity from the retention-time-matched negative-mode feature; unmatched species are retained at the depth their own mode can determine. A dedicated audit tab records every reassignment |
| **Internal standards** | QC plots (positive and negative shown together by default, with the exact standard name on hover) and per-ISTD quantification; IS matched by class, ion mode, and adduct type |
| **Normalisation** | Blank subtraction followed by IS-based quantitative normalisation, with optional protein normalisation |
| **Protein Match** | Diagnostic check confirming every imported sample has a matching protein CSV entry (and vice versa), before relying on protein normalisation |
| **Outlier Detection** | Sample-level (PCA Hotelling's T², iQC replicate deviation) and feature-level (modified Z-score / IQR) outlier flagging, with bulk and individual exclusion that propagates to every downstream tab |
| **Visualisation** | Interactive bar plots per lipid species and per class; the single-species tab has cascading ion-mode, adduct-type, and name filters that only offer selections with data |
| **Export** | Wide and long-format CSV export with class and unit filtering |
| **PCA** | Principal component analysis with token-based group colouring, sample selection, and iQC/ISTD projected as supplementary (non-fit-influencing) points |
| **Statistics** | Auto t-test / one-way ANOVA / two-way ANOVA with post-hoc tests, volcano plot, and significance heatmap |
| **Class Bar Plots** | Faceted bar plots of all species within a selected lipid class post-statistics, with abundance range filtering and significance highlighting |
| **Enrichment** | ORA and FGSEA across lipid class, fatty acid identity, saturation, chain length, odd-chain, and ether subclass sets |
| **Correlation Network** | Pearson correlation network of significant lipid species within a selected group |
| **Synthesis Pathways** | Curated enzyme activity proxy scores with group comparison statistics |
| **Settings** | User-editable ion mode preference lists per lipid class, and the cross-mode identity reconciliation toggle, class list, and RT tolerance |

---

## Getting Started

> **New here?** The [step-by-step tutorial](TUTORIAL.md) walks through every tab using the bundled demo dataset.

### Option 1 — Use the live web app

The easiest way to use LipiRich is via the hosted instance — no installation required:

**[https://lipirich.sarahehancock.com](https://lipirich.sarahehancock.com)**

> **Data privacy:** LipiRich does not store, transmit, or retain any data you upload. All files and results exist only within your browser session and are permanently discarded when you close the tab or the session ends. No data ever leaves the server in any persistent form.

---

### Option 2 — Run locally in R

**Requirements:** R ≥ 4.6.1

1. Clone the repository:
   ```bash
   git clone https://github.com/sarahehancock/LipiRich.git
   cd LipiRich
   ```

2. Install all required packages:
   ```r
   source("packages.R")
   ```

3. Launch the app:
   ```r
   shiny::runApp("app.R")
   ```

---

### Option 3 — Run with Docker (self-hosted)

**Requirements:** Docker

1. Build the image:
   ```bash
   docker build -t lipirich:latest .
   ```

2. Run the container:
   ```bash
   docker run -d -p 3838:3838 --name lipirich lipirich:latest
   ```

3. Open your browser at `http://localhost:3838`

> For deployment behind a reverse proxy (e.g. Cloudflare Tunnel, Nginx), ensure WebSocket support is enabled.

---

## Preparing Your Data in MS-DIAL 5

LipiRich accepts the aligned output exported from MS-DIAL 5 as a tab-delimited `.txt` file. The steps below are recommended for best results.

### Step 1 — Run alignment and curate features

Complete peak picking and alignment in MS-DIAL 5 as normal. Review your alignment result and tag correctly identified lipid species with the **✓ checkmark** in the Peak spot table. Select Tag filter by **✓ checkmark** in the Peak spot table, then select **Filter by current parameter** during alignment export to include only checked features. This is the recommended way to control which species are passed to LipiRich for analysis.

> **Note:** If you do not filter by Tag filter **✓ checkmark** in the Peak spot table and by current parameter on export, all aligned features will appear in the output file and LipiRich will process them all.

### Step 2 — Export aligned data

Go to **Export → Alignment result** and select the `.txt` format. LipiRich skips the first 4 MS-DIAL header rows automatically.

### Step 3 — Sample naming conventions

Sample column names **must end in `_pos` or `_neg`** to indicate ion mode. LipiRich uses these suffixes to identify, merge, and deduplicate ion mode measurements.

```
Sample1_pos   Sample2_pos   Blank_pos
Sample1_neg   Sample2_neg   Blank_neg
```

If a sample name does not contain `_pos` or `_neg`, it will not be detected as a sample column.

### Step 4 — Include blank samples

Include at least one sample named `Blank` (case-insensitive, e.g. `Blank_pos`, `blank_neg`) in your run. LipiRich uses blank signal for background subtraction on a per-metabolite, per-ion-mode basis.

### Step 5 — Tag internal standards

Internal standard features must contain `[IS]` in their metabolite name as assigned by MS-DIAL (e.g. `PC 15:0/18:1-d7 [IS]`). These are automatically detected and separated from analyte features. IS are matched to analytes by **class + ion mode + adduct type** and are never affected by ion mode deduplication rules.

---

## Ion Mode Deduplication

When data are acquired in both positive and negative mode, some lipid species may be detected in both. LipiRich resolves this using per-class mode preference rules, retaining the most analytically appropriate measurement:

| Preference | Default classes |
|---|---|
| **Prefer negative** | PC, PE, PG, PI, PS, PA, CL, PC-O, PE-O |
| **Prefer positive** | TG, DG, MG, CE, Cer, HexCer, SM, LPC, LPE, LPG, LPI, LPS, LPA |
| **Keep both** | All other classes — both ion mode measurements are retained as separate rows |

These defaults reflect established ionisation behaviour for each lipid class (e.g. glycerophospholipids ionise more efficiently and with greater structural information in negative mode; neutral lipids in positive mode). They can be edited per-session in the **⚙️ Settings** tab without reloading data.

IS rows are always retained from whichever ion mode they were detected in, regardless of the class preference rules.

### Cross-mode identity reconciliation

For classes where negative mode resolves fatty-acyl detail that positive mode reports only as sum composition — the glycerophospholipids (PC, PE, PG, PI, PS, PA) and cardiolipin (CL) — LipiRich can go beyond simply preferring one mode. When enabled (the default), each positive-mode feature is paired to a negative-mode feature of the same sum composition and retention time, and:

- the **quantitation** is taken from the positive-mode feature (which has an adduct-matched positive internal standard), while
- the **identity** is transferred from the negative-mode feature (which resolves the individual acyl chains — all four for CL, versus the two combined halves positive mode reports).

A species is therefore quantified once, not double-counted across modes. Unmatched features are never dropped: a positive-only feature is shown at the depth positive mode can determine (sum composition, or the two combined halves for CL), a negative-only feature keeps its full acyl identity, and each is quantified from its own-mode standard. Every decision — re-identification, sum-composition retention, or negative-only retention — is listed in the **📋 Cross-mode audit** tab, showing the final identity and the ion mode each abundance is drawn from.

This mirrors the way structural lipidomics workflows resolve identity in negative mode and quantify in positive. It requires that the reconciled classes have an adduct-matched positive internal standard; the behaviour, the class list, and the retention-time tolerance are all configurable in the **⚙️ Settings** tab, and it can be turned off to revert those classes to prefer-negative handling.

---

## Grouping Samples

LipiRich supports three approaches to assigning samples to groups for plotting and statistics. All tabs inherit grouping from a single central configuration set in the **🔍 Group Preview** tab.

### Token-based grouping (recommended)

Upload your MS-DIAL file, then use the **Group Preview** tab to preview how your sample names split into tokens based on a delimiter (e.g. `_`). Select which token position(s) define your group variable. For example:

```
KO_treated_rep1   →   token 1: KO    token 2: treated    token 3: rep1
WT_vehicle_rep2   →   token 1: WT    token 2: vehicle    token 3: rep2
```

Selecting token 1 groups samples as `KO` vs `WT`. Selecting tokens 1+2 produces `KO_treated`, `KO_vehicle`, `WT_treated`, `WT_vehicle`. A regex mode is also available for more complex naming schemes.

For two-way ANOVA, select Factor A and Factor B token positions separately in the Group Preview tab.

### Grouping CSV (overrides token-based grouping)

To assign groups explicitly — particularly useful for complex experimental designs or when sample names do not follow a consistent pattern — upload a CSV file with the following columns:

| Column | Required | Description |
|---|---|---|
| `sample` | Yes | Sample name matched case-insensitively to MS-DIAL column names |
| `group` | Yes | Group label for plots and statistics |
| `factorA` | No | Factor A label for two-way ANOVA designs |
| `factorB` | No | Factor B label for two-way ANOVA designs |

When a grouping CSV is uploaded and enabled, it takes priority over token-based grouping across all tabs. Blank, iQC, and ISTD-named samples are automatically excluded from group selectors regardless of method.

---

## Internal Standard Normalisation

### Single global ISTD amount

Enter a numeric amount (in pmol, or your chosen unit) in the sidebar. All classes will be normalised using this value.

### Per-ISTD CSV (recommended for multi-class experiments)

Upload a CSV with exactly three columns to specify different amounts and units per internal standard:

| Column | Description |
|---|---|
| `ISTD` | Must match the MS-DIAL `[IS]` metabolite name exactly |
| `amount` | Numeric amount per injection |
| `units` | Units string (e.g. `pmol`, `nmol`) |

IS are matched to analytes by class, ion mode, and adduct type. If an ISTD CSV provides entries with different adducts for the same class (e.g. `[M-H]-` and `[M+NH4]+`), each matching adduct produces a separate normalised value in the output, which can then be filtered using the adduct filter on the plot tabs.

If an analyte's IS cannot be matched via the CSV, LipiRich falls back to the single global ISTD amount entered in the sidebar.

---

## Protein Normalisation & Protein Match

Optionally, normalise IS-normalised or background-subtracted values by per-sample protein content, in addition to (not instead of) internal standard normalisation. Upload a CSV with exactly two columns:

| Column | Description |
|---|---|
| `sample` | Sample name matched (trimmed) to the MS-DIAL sample name, ignoring the `_pos`/`_neg` suffix |
| `protein` | Numeric protein content per sample |

Enable **Apply protein normalisation** in the sidebar to divide `norm`/`value_bs` by each sample's protein value. Samples with no matching protein CSV row are left un-normalised rather than dropped.

### Protein Match tab (Step 3a)

Before relying on protein normalisation, the **Protein Match** tab reports whether every imported sample has a matching protein CSV row, and flags any protein CSV rows that don't correspond to an imported sample — a common sign of a typo in one file or the other. ISTD, Blank, and iQC/QC samples are excluded from this check on both sides, since they are not expected to have protein measurements.

---

## Outlier Detection (Step 3b)

A combined sample-level and feature-level outlier workflow sits upstream of all visualisation and statistics tabs, so any exclusion applies consistently everywhere downstream.

**Sample-level:**
- **PCA Hotelling's T²** — a self-contained PCA (log-transformed, unit-variance scaled) flags samples beyond a user-selected confidence threshold (95/97.5/99%).
- **iQC replicate deviation** — for datasets with ≥ 3 iQC replicates, flags iQC samples whose deviation from the cross-replicate median exceeds a MAD-based threshold. This check is informational only; iQC samples are never excluded from downstream analyses regardless of the result.

**Feature-level:** flags individual sample values within a Metabolite name × group combination using a modified Z-score (MAD) or the classic IQR rule, both with adjustable thresholds.

**Exclusion:** bulk toggles (off by default) exclude every auto-flagged sample/point. Independently, selectable review tables allow specific samples or points to be individually excluded or kept, taking precedence over the bulk toggles. Excluded feature values are set to missing rather than deleted; excluded samples are removed only from the shared in-app dataset, never from the uploaded file.

ISTD, Blank, and iQC/QC samples are excluded from the entire workflow — none of them can be flagged or excluded, since they are not biological replicates.

Both plots (PCA Hotelling's T² and iQC replicate deviation) have standard Width/Height/DPI/Scale/Base font size export controls with PNG/SVG buttons, sharing one set of dimension controls.

---

## PCA

Principal component analysis on IS-normalised or background-subtracted values, with token-based group colouring and per-group sample selection.

### iQC and ISTD as supplementary (projected) points

The PCA fit — axes, loadings, and centring/scaling statistics — is computed from biological samples only. If included, iQC and ISTD-only samples are projected into the fitted space afterward as supplementary individuals (via FactoMineR's `ind.sup`), so neither can distort the ordination:

- **Include iQC samples** — on by default.
- **Include ISTD-only samples** — off by default. An ISTD-only injection has no biological matrix and a very different lipid profile to a real sample, so it is excluded unless explicitly requested.

Since neither type is typically present in the protein CSV, when protein normalisation is active they are scaled by the *median* protein content of the biological samples purely so the projection lands somewhere visually comparable — this is an assumed, not measured, value. The plot legend distinguishes real groups, iQC, and ISTD by colour; supplementary points are shown as triangles, active points as circles.

### Export and labels

Standard Width/Height/DPI/Scale/Base font size export controls with PNG/SVG buttons, plus a **Show sample ID labels** toggle (off by default).

---

## Statistics

The Statistics tab automatically selects the appropriate test based on the number of groups:

- **2 groups** → unpaired t-test (Welch by default)
- **≥ 3 groups** → one-way ANOVA
- **Two-way design** → two-way ANOVA using Factor A and Factor B defined in the Group Preview tab

Groups to include in the comparison are selected via the group selector on the Statistics tab, which inherits from the Group Preview tab settings. Deselect groups to exclude them from all tests. All downstream tabs (Volcano Plot, Heatmap, Enrichment, Correlation Network, Synthesis Pathways) inherit the statistics results from the most recent "Run statistics" click, and restrict their outputs to the selected groups only.

Multiple testing correction (FDR/BH, Bonferroni, or none) and significance threshold (α) are configurable. Post-hoc tests (Tukey HSD, pairwise t-tests with Holm correction) are available for ANOVA results.

A collapsible **per-class summary panel** shows the number of species detected, tested, significant, and percentage significant for each lipid class, with a CSV download.

### Volcano plot

The volcano plot displays **all tested features** regardless of the significance filter applied to the bar plots and table. Points are coloured by direction (up/down/NS) according to the configured log2FC and p-value thresholds. Significant feature labels use native plotly annotations in the interactive view and ggrepel in PNG/SVG exports. The plot is interactive (hover for feature details) and can be exported as PNG or SVG.

For **ANOVA (≥ 3 groups)** results, a **Comparison** selector appears in the sidebar — select any pairwise combination to render a volcano for that specific comparison using post-hoc p-values and group-mean log2FC. Select "Overall (ANOVA)" to revert to the overall F-test result.

An optional **axis capping** mode clamps the -log10(p) and log2FC axes to user-defined limits. Features beyond the cap are plotted at the cap value as filled triangles (▲), keeping the rest of the plot readable without hiding any data. True values remain accessible via hover text.

### Significance heatmap

Displays the top N significant features as a z-scored heatmap. Row and column clustering, row/column label visibility, and colour palette (10 options including diverging schemes) are all configurable. Column clustering is automatically disabled for 2-group (t-test) results. Exported as PNG or SVG with user-specified pixel dimensions and DPI.

---

## Class Bar Plots

The Class Bar Plots tab displays all lipid species within a selected class as a **single grouped bar plot** with one colour per group, using the data and group assignments from the most recent statistics run. All groups selected on the Statistics tab are shown.

### Controls

**Data** — select a lipid class (populated immediately on data load), the value type (IS-normalised or background-subtracted), and the units for the y-axis label.

**Filtering:**
- **Show significant species only** — restricts the plot to species that passed the significance threshold. Reads from the full set of tested features so this toggle is always meaningful.
- **Abundance range filter** — a percentage slider that filters species by mean abundance relative to the most abundant species in the class. For example, 10–100% removes low-abundance species; 0–10% shows only minor species.

**Appearance:**
- **Colour palette** — 18 options across qualitative (Okabe-Ito default, Tableau 10, Dark2, Set1, Set2, Paired, Accent), sequential (Viridis, Plasma, Inferno, Blues, Greens, Purples), and diverging (RdBu, PuOr, BrBG) schemes.
- **Bar opacity** and **point size** — numeric controls.
- **Point fill** — match palette colour, black, or white.
- **Point outline** — white, black, or none.
- **Error bars** — SEM, SD, or 95% CI.
- **Species ordering** — alphabetical (default) or by descending mean abundance.

**Significance brackets** — drawn automatically on the plot using the same pairwise test logic as the Statistics tab bar plots. For t-test (2 groups): a single bracket spanning both bars. For ANOVA (≥3 groups): one stacked bracket per significant Tukey HSD pairwise comparison. Symbols follow the standard convention (`*` p < 0.05, `**` p < 0.01, `***` p < 0.001).

---

## Correlation Network

The Correlation Network tab computes pairwise **Pearson correlations** among significant lipid species within a selected group. Edges are drawn between species whose correlation exceeds a configurable threshold, with edge width scaled to correlation strength and node colour reflecting direction of change relative to the comparison group.

> This tab is purely correlational — edges reflect co-variation in abundance, not metabolic pathway connectivity.

Group selection is inherited from the Statistics tab: only groups included in the most recent statistics run are available in the group dropdown.

---

## Synthesis Pathways

The Synthesis Pathways tab provides two views of curated **enzyme activity proxy scores** — lipid class ratios and fractions designed to reflect the relative activity of key lipid synthesis enzymes (e.g. Kennedy pathway, PEMT, sphingomyelin synthase, ether lipid synthesis).

For the ether lipid branch, plasmanyl-PC is scored two ways: **direct** (PC-O plasmanyl/DG-O), reflecting the CDP-choline route straight from the shared ether precursor, and **headgroup-conversion** (PC-O plasmanyl/PE-O plasmanyl), reflecting the alternative route via plasmanyl-PE. Plasmenyl-PC is only scored via headgroup conversion (PC-O plasmenyl/PE-O plasmenyl), since mammalian plasmalogen desaturation is PE-selective and there is no direct route from DG-O to plasmenyl-PC.

> **Interpret with care:** these scores are built from established *mammalian* synthesis pathways. Some steps involve enzymes with overlapping or tissue-dependent substrate preferences, and a given score may reflect more than one biosynthetic route contributing to the same lipid pool. Treat them as pathway-activity indicators rather than direct measurements of flux through a single enzymatic step.

### Scores heatmap

Each row is a pathway score; each column is a sample (or group mean). Values are z-scored by row to highlight relative differences across groups. Colour palette, clustering, and export dimensions are configurable.

### Score statistics

Runs a **t-test** (2 groups) or **one-way ANOVA** (≥ 3 groups) on each score's sample values across the groups selected in the Statistics tab, with BH correction applied across all scores. Results are displayed as:

- A **lollipop plot** of −log10(p_adj) per score, sorted by significance and coloured by which group has the higher mean.
- A **results table** with score name, test type, groups compared, test statistic, log2FC (t-test only), direction, p, and p_adj.

> Scores are currently based on lipid class-level totals. Acyl chain-resolved scoring is planned for a future release.

---

## Plot Export

All tabs with plots include consistent export controls:

| Control | Description |
|---|---|
| **Width / Height (px)** | Output dimensions in pixels |
| **DPI** | Resolution for raster outputs (PNG) |
| **Scale fraction** | Scales the final plot size (0.25–2.0×) |
| **Base font size** | Controls axis label and annotation text size |

PNG and SVG formats are available on all plot tabs. The Statistics tab additionally offers a **batch PDF export** that produces one bar plot per significant feature across all pages.

---

## Demo Data

A full demonstration dataset is available in [`demo_data/`](demo_data/), generated from mouse liver samples (n = 8 per group: chow, lard-based high-fat diet, and 90% fish oil high-fat diet). Animal and diet details are as previously described (Liu et al. 2014, *Scientific Reports* 4:5538). Lipid extraction and LC-MS data acquisition are described in the Zenodo repository for `.raw` files (see below).

| File | Description |
|---|---|
| `Demo data pos.txt` | MS-DIAL 5 aligned output, positive ion mode |
| `Demo data neg.txt` | MS-DIAL 5 aligned output, negative ion mode |
| `Demo data metadata.xlsx` | Combined reference workbook: sample order/analytical run order, protein content, and internal standard amounts (three sheets) |
| `Demo data protein.csv` | Protein content per sample, formatted for direct upload to LipiRich's protein normalisation step |
| `Demo data  ISTD_amounts.csv` | Internal standard amounts, formatted for direct upload to LipiRich |
| `Processing parameters pos.mdparameter` | MS-DIAL 5 processing parameters used to generate `Demo data pos.txt` |
| `Processing parameters neg.mdparameter` | MS-DIAL 5 processing parameters used to generate `Demo data neg.txt` |

To try LipiRich immediately, upload `Demo data pos.txt` and `Demo data neg.txt` directly (Step 1), then `Demo data protein.csv` and `Demo data ISTD_amounts.csv` at the relevant normalisation steps.

**Raw instrument data:** The original LC-MS `.raw` files from which the above were generated in MS-DIAL 5 are archived separately due to size:
- DOI: [10.5281/zenodo.21448733](https://doi.org/10.5281/zenodo.21448733)

The `.mdparameter` files above can be used to reproduce the exact MS-DIAL 5 alignment from these raw files. For a full walkthrough of this process, including peak curation and export, see the [MS-DIAL processing tutorial](demo_data/MSDIAL_processing_tutorial.md).

---

## Repository Structure

```
LipiRich/
├── app.R            # Main Shiny application
├── packages.R       # R package installer (run once before first use)
├── Dockerfile       # Docker image definition
├── README.md        # This file
├── TUTORIAL.md      # Step-by-step walkthrough with screenshots
├── CHANGELOG.md     # Version history
├── LICENSE          # GNU Affero General Public License v3.0
├── docs/
│   └── images/      # Screenshots for TUTORIAL.md
└── demo_data/       # Example input files
    ├── MSDIAL_processing_tutorial.md  # Raw data → MS-DIAL alignment walkthrough
    ├── images/      # Screenshots for MSDIAL_processing_tutorial.md
    ├── Demo data pos.txt
    ├── Demo data neg.txt
    ├── Demo data metadata.xlsx
    ├── Demo data protein.csv
    ├── Demo data ISTD_amounts.csv
    ├── Processing parameters pos.mdparameter
    └── Processing parameters neg.mdparameter
```

---

## Citation

If you use LipiRich in your research, please cite:

> Hancock, SE. (2026). *LipiRich: A Shiny application for normalisation, statistics, and visualisation of MS-DIAL lipidomics data* (v0.4.0). GitHub: https://github.com/sarahehancock/LipiRich. DOI: [pending]

---

## License

This project is licensed under the [GNU Affero General Public License v3.0](LICENSE) (AGPL-3.0).

If you use, modify, or deploy LipiRich as a network service (e.g. a hosted web application), you must make the complete corresponding source code available to users under the same licence. See the [AGPL-3.0 licence text](https://www.gnu.org/licenses/agpl-3.0.html) for full details.

---

## Acknowledgements

LipiRich builds on the following open-source R packages and tools:

- [MS-DIAL](https://systemsomicslab.github.io/compms/msdial/main.html) — Tsugawa et al., *Nature Communications* (2024)
- [fgsea](https://bioconductor.org/packages/fgsea/) — Korotkevich et al.
- [FactoMineR](https://cran.r-project.org/package=FactoMineR) — Lê et al., *Journal of Statistical Software* (2008)
- [Shiny](https://shiny.posit.co/) — Chang et al., Posit Software
- [visNetwork](https://cran.r-project.org/package=visNetwork) — Almende B.V.
- [plotly](https://plotly.com/r/) — Sievert et al.
- [pheatmap](https://cran.r-project.org/package=pheatmap) — Kolde et al.
- [viridis](https://cran.r-project.org/package=viridis) — Garnier et al.
