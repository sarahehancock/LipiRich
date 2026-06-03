# LipiRich <img src="https://img.shields.io/badge/version-0.0.2-blue" alt="v0.0.2"/> <img src="https://img.shields.io/badge/license-AGPL--3.0-blue" alt="AGPL-3.0"/> <img src="https://img.shields.io/badge/R-%3E%3D4.5.2-informational" alt="R 4.5.2"/> <img src="https://img.shields.io/badge/live%20app-lipirich.sarahehancock.com-brightgreen" alt="Live App"/>

**LipiRich** is an open-source, browser-based Shiny application for the normalisation, statistical analysis, and visualisation of untargeted lipidomics data exported from [MS-DIAL 5](https://systemsomicslab.github.io/compms/msdial/main.html). It requires no programming knowledge and runs entirely in a web browser.

> Developed and tested with **MS-DIAL 5.5.251021**, R 4.5.2, and Bioconductor 3.22.

> ⚠️ **Pre-publication software (v0.0.2):** LipiRich is under active development. A citable preprint and demonstration dataset will be released alongside v1.0.0. Please check the [GitHub repository](https://github.com/sarahehancock/LipiRich) for the latest updates and to report issues.

---

## Features

| Module | Description |
|---|---|
| **Data import** | Upload 1–2 MS-DIAL aligned `.txt` files; positive and negative ion modes merged automatically |
| **Ion mode deduplication** | Per-class mode preference rules (user-configurable) ensure each species is represented by its most informative ion mode |
| **Internal standards** | QC plots and per-ISTD quantification; IS matched by class, ion mode, and adduct type |
| **Normalisation** | Blank subtraction followed by IS-based quantitative normalisation |
| **Visualisation** | Interactive bar plots per lipid species and per class, with adduct-type filtering |
| **Export** | Wide and long-format CSV export with class and unit filtering |
| **PCA** | Principal component analysis with token-based group colouring and sample selection |
| **Statistics** | Auto t-test / one-way ANOVA / two-way ANOVA with post-hoc tests, volcano plot, and significance heatmap |
| **Enrichment** | ORA and FGSEA across lipid class, fatty acid identity, saturation, chain length, and ether subclass sets |
| **Correlation Network** | Pearson correlation network of significant lipid species within a selected group |
| **Synthesis Pathways** | Curated enzyme activity proxy scores with group comparison statistics |
| **Settings** | User-editable ion mode preference lists per lipid class |

---

## Getting Started

### Option 1 — Use the live web app

The easiest way to use LipiRich is via the hosted instance — no installation required:

**[https://lipirich.sarahehancock.com](https://lipirich.sarahehancock.com)**

> **Data privacy:** LipiRich does not store, transmit, or retain any data you upload. All files and results exist only within your browser session and are permanently discarded when you close the tab or the session ends. No data ever leaves the server in any persistent form.

---

### Option 2 — Run locally in R

**Requirements:** R ≥ 4.5.2

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

Complete peak picking and alignment in MS-DIAL 5 as normal. Review your alignment result and tag correctly identified lipid species with the **✓ checkmark** in MS-DIAL, then select **Filter by current parameter** during alignment export to include only checked features. This is the recommended way to control which species are passed to LipiRich for analysis.

> **Note:** If you do not filter by current parameter on export, all aligned features will appear in the output file and LipiRich will process them all.

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

When a grouping CSV is uploaded and enabled, it takes priority over token-based grouping across all tabs. Blank samples and iQC samples are automatically excluded from group selectors regardless of method.

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

## Statistics

The Statistics tab automatically selects the appropriate test based on the number of groups:

- **2 groups** → unpaired t-test (Welch by default)
- **≥ 3 groups** → one-way ANOVA
- **Two-way design** → two-way ANOVA using Factor A and Factor B defined in the Group Preview tab

Groups to include in the comparison are selected via the group selector on the Statistics tab, which inherits from the Group Preview tab settings. Deselect groups to exclude them from all tests. All downstream tabs (Volcano Plot, Heatmap, Enrichment, Correlation Network, Synthesis Pathways) inherit the statistics results from the most recent "Run statistics" click, and restrict their outputs to the selected groups only.

Multiple testing correction (FDR/BH, Bonferroni, or none) and significance threshold (α) are configurable. Post-hoc tests (Tukey HSD, pairwise t-tests with Holm correction) are available for ANOVA results.

### Volcano plot

The volcano plot displays **all tested features** regardless of the significance filter applied to the bar plots and table. Points are coloured by direction (up/down/NS) according to the configured log2FC and p-value thresholds. The plot is interactive (hover for feature details) and can be exported as PNG or SVG.

### Significance heatmap

Displays the top N significant features as a z-scored heatmap. Row and column clustering, label visibility, and colour palette (10 options including diverging schemes) are all configurable. Exported as PNG or SVG with user-specified pixel dimensions and DPI.

---

## Correlation Network

The Correlation Network tab computes pairwise **Pearson correlations** among significant lipid species within a selected group. Edges are drawn between species whose correlation exceeds a configurable threshold, with edge width scaled to correlation strength and node colour reflecting direction of change relative to the comparison group.

> This tab is purely correlational — edges reflect co-variation in abundance, not metabolic pathway connectivity.

Group selection is inherited from the Statistics tab: only groups included in the most recent statistics run are available in the group dropdown.

---

## Synthesis Pathways

The Synthesis Pathways tab provides two views of curated **enzyme activity proxy scores** — lipid class ratios and fractions designed to reflect the relative activity of key lipid synthesis enzymes (e.g. Kennedy pathway, PEMT, sphingomyelin synthase, ether lipid synthesis).

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

A demonstration dataset compatible with LipiRich is available at:

> **[Link pending]** — example MS-DIAL 5 aligned output files (.txt) for positive and negative ion mode.

The raw data used to generate the demo files is publicly available via the MS-DIAL 5 tutorial repository on Zenodo:
- DOI: [10.5281/zenodo.10616947](https://doi.org/10.5281/zenodo.10616947)

---

## Repository Structure

```
LipiRich/
├── app.R            # Main Shiny application
├── packages.R       # R package installer (run once before first use)
├── Dockerfile       # Docker image definition
├── README.md        # This file
├── CHANGELOG.md     # Version history
├── LICENSE          # GNU Affero General Public License v3.0
└── demo_data/       # Example input files (coming soon)
```

---

## Citation

If you use LipiRich in your research, please cite:

> Hancock, SE. (2025). *LipiRich: A Shiny application for normalisation, statistics, and visualisation of MS-DIAL lipidomics data* (v0.0.2). GitHub: https://github.com/sarahehancock/LipiRich. DOI: [pending]

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
