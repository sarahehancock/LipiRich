# app.R
# --------------------------
# LipiRich v0.3.0
# Copyright (C) 2025–2026 Sarah E. Hancock
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
# License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# --------------------------
# A Shiny application for normalisation, statistics, and visualisation
# of MS-DIAL lipidomics alignment output data.
#
# Author:   Hancock, SE.
# GitHub:   https://github.com/sarahehancock/LipiRich
# License:  GNU Affero General Public License v3.0 (AGPL-3.0)
#
# LipiRich — Copyright (C) 2025–2026 Sarah E. Hancock
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/agpl-3.0.html>.
# Version:  0.3.0
# Tested with: MS-DIAL 5.5.251021, R 4.6.1, Bioconductor 3.23
# --------------------------

APP_VERSION <- "0.3.0"

suppressPackageStartupMessages({
  library(shiny); library(DT); library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(stringr); library(plotly); library(FactoMineR)
  library(factoextra); library(tibble); library(rlang); library(ggpubr)
  library(rstatix); library(pheatmap); library(svglite)
  # visNetwork loaded on demand in network tab
})

options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

# --- Helper functions (GLOBAL) ---

make_group_labels_generic <- function(samples, method, delim = "_", tokens = 1, regex = "^[^_]+") {
  if (is.null(method)) method <- "By sample"
  if (method == "By sample") return(samples)
  if (method == "Delimiter-based") {
    if (is.null(delim) || delim == "") delim <- "_"
    if (is.null(tokens) || is.na(tokens) || tokens < 1) tokens <- 1
    labs <- vapply(samples, function(s) {
      parts <- strsplit(s, split = delim, fixed = TRUE)[[1]]
      paste(head(parts, tokens), collapse = delim)
    }, FUN.VALUE = character(1))
    return(labs)
  }
  if (method == "Regex capture group") {
    if (is.null(regex) || regex == "") regex <- "^[^_]+"
    m <- stringr::str_match(samples, regex)
    labs <- if (ncol(m) >= 2) ifelse(is.na(m[, 2]), samples, m[, 2]) else samples
    return(labs)
  }
  samples
}

# ── resolve_group_labels: replaces all repeated grouping if/else blocks ──
resolve_group_labels <- function(samples, method, delim, tokens, regex,
                                 group_csv_input, use_csv_flag, group_map) {
  if (identical(method, "By sample")) return(samples)
  if (!is.null(group_csv_input) && isTRUE(use_csv_flag))
    return(prefer_csv_group_labels(samples, method, delim, tokens, regex, group_map))
  make_group_labels_generic(samples, method, delim, tokens, regex)
}

# ── impute_row_median: row-median imputation for matrices ──
impute_row_median <- function(mat) {
  for (i in seq_len(nrow(mat))) {
    v <- mat[i, ]
    if (anyNA(v) && !all(is.na(v)))
      mat[i, is.na(v)] <- stats::median(v, na.rm = TRUE)
  }
  mat
}

# ── load_msdial_uploaded: moved to global (no reactive deps) ──
load_msdial_uploaded <- function(upload_df) {
  validate(need(!is.null(upload_df) && nrow(upload_df) >= 1,
                "Please upload at least one .txt file."))
  validate(need(nrow(upload_df) <= 2,
                "Please upload no more than two .txt files."))
  df_list <- lapply(seq_along(upload_df$datapath), function(i)
    readr::read_tsv(upload_df$datapath[i], skip = 4,
                    col_types = readr::cols(.default = "c")))
  out <- dplyr::bind_rows(df_list)
  attr(out, "source_files") <- basename(upload_df$name)
  out
}


# ---- Canonicalization for matching sample names to CSV keys ----
canonicalize_sample <- function(x) {
  x <- if (is.null(x)) character(0) else x
  x <- stringr::str_trim(x)
  x <- stringr::str_to_lower(x)
  
  # Strip ion-mode suffixes (e.g., "_pos", "_neg") plus trailing tokens like ".1"
  # Also remove stray trailing replication markers like "-rep1", "_rep2"
  x <- gsub("(_pos|_neg)(.*)?$", "", x, ignore.case = TRUE)
  x <- gsub("[\\._-]rep[0-9]+$", "", x, ignore.case = TRUE)
  x <- gsub("\\.[0-9]+$", "", x, ignore.case = TRUE)
  
  # Collapse multiple separators: hyphens, underscores, spaces
  x <- gsub("[ \\-_]+", "_", x)
  x
}

# ---- Detect columns in uploaded CSV (sample, group, factorA/B) ----
pick_col <- function(df, candidates) {
  nm <- names(df)
  idx <- which(tolower(nm) %in% tolower(candidates))
  if (length(idx) >= 1) nm[idx[1]] else NA_character_
}

# ---- Build mapping from uploaded CSV ----
# Returns a list: list(group = named_vector, factorA = named_vector_or_NULL, factorB = named_vector_or_NULL)
build_group_map <- function(df) {
  csample <- pick_col(df, c("sample","Sample","sample_id","SampleID","filename","file","File"))
  cgroup  <- pick_col(df, c("group","Group","condition","Condition"))
  cA      <- pick_col(df, c("factorA","FactorA","A","GroupA"))
  cB      <- pick_col(df, c("factorB","FactorB","B","GroupB"))
  cBio    <- pick_col(df, c("bio_sample","biological_sample","biol_sample",
                            "parent_sample","BioSample","ParentSample","biosample"))
  
  validate(need(!is.na(csample) && !is.na(cgroup),
                "CSV must contain at least 'sample' and 'group' columns."))
  
  key    <- canonicalize_sample(df[[csample]])
  gmap   <- setNames(as.character(df[[cgroup]]), key)
  amap   <- if (!is.na(cA))   setNames(as.character(df[[cA]]),   key) else NULL
  bmap   <- if (!is.na(cB))   setNames(as.character(df[[cB]]),   key) else NULL
  biomap <- if (!is.na(cBio)) setNames(as.character(df[[cBio]]), key) else NULL
  
  list(group = gmap, factorA = amap, factorB = bmap, bio_sample = biomap,
       colnames = list(sample = csample, group = cgroup, factorA = cA, factorB = cB, bio_sample = cBio))
}

# ---- Labelers that prefer the uploaded CSV, else fallback to your existing methods ----
labels_from_map <- function(samples, named_map) {
  if (is.null(named_map)) return(rep(NA_character_, length(samples)))
  keys <- canonicalize_sample(samples)
  labs <- unname(named_map[keys])
  labs
}

prefer_csv_group_labels <- function(samples, method, delim, tokens, regex, group_map) {
  csv_labs <- labels_from_map(samples, group_map$group)
  # Use CSV label when present (non-empty), else fallback to your original strategy
  ifelse(!is.na(csv_labs) & csv_labs != "",
         csv_labs,
         make_group_labels_generic(samples, method, delim, tokens, regex))
}

prefer_csv_factor_labels <- function(samples, which, method, delim, tokens, regex, group_map) {
  named_map <- switch(which,
                      "A" = group_map$factorA,
                      "B" = group_map$factorB,
                      NULL)
  csv_labs <- labels_from_map(samples, named_map)
  ifelse(!is.na(csv_labs) & csv_labs != "",
         csv_labs,
         make_group_labels_generic(samples, method, delim, tokens, regex))
}

# ---- Loading modal helpers ----
show_loading_modal <- function(title = "Running enrichment…", detail = NULL) {
  shiny::showModal(
    shiny::modalDialog(
      title = title,
      tags$div(
        style = "display:flex; align-items:center; gap:12px;",
        # Simple CSS spinner
        tags$div(style = "
            width:22px;height:22px;border-radius:50%;
            border:3px solid #ccc;border-top-color:#0078D4;
            animation: spin 0.8s linear infinite;
        "),
        tags$div(if (is.null(detail)) "Please wait…" else detail)
      ),
      footer = NULL,
      easyClose = FALSE
    )
  )
}

hide_loading_modal <- function() {
  shiny::removeModal()
}

# Spinner keyframes (inject once in UI head)
enrichment_spinner_css <- tags$style(HTML("
@keyframes spin { from {transform: rotate(0deg);} to {transform: rotate(360deg);} }
"))

# --- Idle-session timeout ---
# Client tracks activity (mouse/keyboard/touch/scroll). After 15 min of no
# activity, an "idle_warning" event asks the server to show a warning modal.
# After 20 min of no activity, an "idle_timeout" event asks the server to
# close the session outright (frees server memory held by abandoned tabs;
# data was only ever held in-session per the app's no-persistence design).
# Any tracked activity — including clicking the warning modal's own button —
# resets both timers via the shared document-level listeners.
idle_timeout_js <- tags$script(HTML("
(function() {
  var WARN_MS    = 15 * 60 * 1000;
  var TIMEOUT_MS = 20 * 60 * 1000;
  var warnTimer, timeoutTimer;
  var warningActive = false;

  function resetIdleTimers() {
    clearTimeout(warnTimer);
    clearTimeout(timeoutTimer);
    if (warningActive) {
      warningActive = false;
      Shiny.setInputValue('idle_dismiss_warning', Math.random(), {priority: 'event'});
    }
    warnTimer = setTimeout(function() {
      warningActive = true;
      Shiny.setInputValue('idle_warning', Math.random(), {priority: 'event'});
    }, WARN_MS);
    timeoutTimer = setTimeout(function() {
      Shiny.setInputValue('idle_timeout', Math.random(), {priority: 'event'});
    }, TIMEOUT_MS);
  }

  ['mousemove', 'mousedown', 'keydown', 'scroll', 'touchstart', 'click'].forEach(function(evt) {
    document.addEventListener(evt, resetIdleTimers, {passive: true});
  });

  $(document).on('shiny:connected', resetIdleTimers);
})();
"))


# --- iQC / QC / ISTD sample tag helpers ---
# iQC and ISTD are NOT the same thing: iQC is a pooled biological QC sample
# (real matrix, monitors extraction/injection reproducibility); a sample
# literally named ISTD/ITSD is typically an internal-standard-only injection
# with no biological matrix, and has a wildly different lipid profile. They
# are kept as separate checks so callers that need to tell them apart (e.g.
# PCA) can, while callers that just want "not a real biological sample"
# (is_qc_type_sample / filter_iqc / is_protected_sample) still catch both.
is_iqc_sample <- function(x) {
  # Matches iQC or standalone QC only — NOT ISTD/ITSD
  grepl("iqc|\\bqc\\b", x, ignore.case = TRUE, perl = TRUE)
}
is_istd_sample <- function(x) {
  # Matches ISTD/ITSD only — a whole sample injection, not the "[IS]" metabolite tag
  grepl("istd|itsd", x, ignore.case = TRUE, perl = TRUE)
}
is_qc_type_sample <- function(x) {
  is_iqc_sample(x) | is_istd_sample(x)
}
filter_iqc <- function(df, include_iqc = FALSE, sample_col = "sample") {
  if (include_iqc) return(df)
  if (is.null(df) || nrow(df) == 0) return(df)
  sc <- df[[sample_col]]
  df[!is_qc_type_sample(sc), , drop = FALSE]
}

# Samples that should never be treated as biological samples for outlier
# detection/exclusion or protein-match checking: ISTD, iQC/QC, and Blank
# (flexible match: "Blank", "Blank_1", "Blank-2", etc.)
is_protected_sample <- function(x) {
  x  <- as.character(x)
  sn <- stringr::str_trim(stringr::str_to_lower(x))
  is_qc_type_sample(x) | stringr::str_detect(sn, "^blank(?:$|[-_])")
}


# ── Tokenisation helpers (ported from MetaboDash) ────────────────────────────

# Vectorised single-token extractor
# token_index: 1-based; negative counts from end (-1 = last)
extract_token <- function(x, delimiter = "_", token_index = 1) {
  x     <- as.character(x)
  parts <- stringr::str_split(x, pattern = delimiter, simplify = TRUE)
  if (!is.matrix(parts) || ncol(parts) == 0)
    return(rep(NA_character_, length(x)))
  k   <- ncol(parts)
  idx <- if (token_index < 0) (k + token_index + 1L) else token_index
  out <- if (idx >= 1 && idx <= k) parts[, idx] else rep(NA_character_, nrow(parts))
  dplyr::na_if(trimws(out), "")
}

# Combine multiple token indices into one label
# e.g. indices c(1,2) from "Group_Rep_pos" → "Group_Rep"
extract_combined_tokens <- function(x, delimiter = "_", token_indices = 1L, sep = "_") {
  x     <- as.character(x)
  parts <- stringr::str_split(x, pattern = delimiter, simplify = TRUE)
  if (!is.matrix(parts) || ncol(parts) == 0)
    return(rep(NA_character_, length(x)))
  k       <- ncol(parts)
  indices <- as.integer(token_indices)
  indices <- indices[indices >= 1 & indices <= k]
  if (length(indices) == 0) return(rep(NA_character_, length(x)))
  if (length(indices) == 1) {
    out <- parts[, indices, drop = TRUE]
  } else {
    out <- apply(parts[, indices, drop = FALSE], 1, function(r)
      paste(r[nzchar(r) & !is.na(r)], collapse = sep))
  }
  dplyr::na_if(trimws(out), "")
}

# Build labelled token choices for checkboxGroupInput/selectInput
# Returns named character vector: "1  (e.g. WT)" → "1"
build_token_choices <- function(sample_names, delimiter = "_") {
  parts    <- stringr::str_split(as.character(sample_names),
                                 pattern = delimiter, simplify = TRUE)
  n_tokens <- ncol(parts)
  if (n_tokens == 0) return(c("1" = "1"))
  token_labels <- vapply(seq_len(n_tokens), function(i) {
    vals <- parts[, i]
    vals <- vals[nzchar(vals) & !is.na(vals)]
    top  <- if (length(vals) > 0)
      names(sort(table(vals), decreasing = TRUE))[1] else "?"
    paste0(i, "  (e.g. ", top, ")")
  }, character(1))
  setNames(as.character(seq_len(n_tokens)), token_labels)
}

# ── Landing page CSS ──────────────────────────────────────────────────────────
landing_css <- tags$style(HTML("
  .lipid-landing{
    --lr-ink:#131b20; --lr-muted:#5c6b72; --lr-surface:#ffffff;
    --lr-line:#e4eaeb; --lr-primary:#0b6b66; --lr-primary-ink:#08514d;
    --lr-primary-soft:#e9f3f2; --lr-amber-soft:#fbf3e1;
    --lr-sans:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,'Helvetica Neue',Arial,sans-serif;
    --lr-mono:ui-monospace,'SF Mono','JetBrains Mono',Menlo,Consolas,monospace;
    --lr-r:14px; --lr-r-sm:10px;
    --lr-shadow:0 1px 2px rgba(16,24,32,.05),0 12px 30px -18px rgba(16,24,32,.22);
    font-family:var(--lr-sans); font-size:16px; color:var(--lr-ink);
    max-width:1000px; margin:0 auto; padding:8px 24px 56px;
    -webkit-font-smoothing:antialiased;
  }
  .lipid-landing *{ box-sizing:border-box; }
  .lipid-landing .lr-eyebrow{ font-family:var(--lr-mono); font-size:.78em; letter-spacing:.14em; text-transform:uppercase; color:var(--lr-primary); font-weight:600; }
  .lipid-landing h1.lr-brand{ font-size:2.7em; font-weight:800; letter-spacing:-.02em; color:var(--lr-ink); margin:.15em 0 0; line-height:1.04; }
  .lipid-landing .lr-ver{ font-family:var(--lr-mono); font-size:11px; font-weight:600; letter-spacing:.06em; color:var(--lr-primary-ink); background:var(--lr-primary-soft); border:1px solid #cfe6e3; padding:3px 9px; border-radius:999px; vertical-align:middle; margin-left:12px; }
  .lipid-landing .lr-tagline{ color:var(--lr-muted); font-size:1.32em; margin:.55em 0 0; max-width:62ch; }
  .lipid-landing .lr-badges{ display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-top:10px; }
  .lipid-landing .lr-badges img{ height:20px; }
  .lipid-landing .lr-section{ margin-top:36px; }
  .lipid-landing .lr-section > h2{ font-size:1.52em; font-weight:700; color:var(--lr-ink); margin:.3em 0 .2em; letter-spacing:-.01em; }
  .lipid-landing p, .lipid-landing li{ font-size:1.18em; line-height:1.7; color:#2c3a40; }
  .lipid-landing a{ color:var(--lr-primary-ink); text-decoration:none; border-bottom:1px solid #bfe0dc; }
  .lipid-landing a:hover{ border-bottom-color:var(--lr-primary); }
  .lipid-landing .lr-grid{ display:grid; grid-template-columns:repeat(auto-fill,minmax(220px,1fr)); gap:14px; margin-top:16px; }
  .lipid-landing .lr-card{ background:var(--lr-surface); border:1px solid var(--lr-line); border-radius:var(--lr-r); padding:16px 18px; transition:transform .15s ease, box-shadow .15s ease, border-color .15s ease; }
  .lipid-landing .lr-card:hover{ transform:translateY(-2px); box-shadow:var(--lr-shadow); border-color:#cfe0df; }
  .lipid-landing .lr-step{ font-family:var(--lr-mono); font-size:.74em; font-weight:600; letter-spacing:.08em; color:var(--lr-primary); text-transform:uppercase; }
  .lipid-landing .lr-card h3{ font-size:1.16em; font-weight:700; color:var(--lr-ink); margin:6px 0 4px; }
  .lipid-landing .lr-card p{ font-size:1.08em; line-height:1.55; color:#41505a; margin:0; }
  .lipid-landing .lr-feature{ background:linear-gradient(160deg,#0b6b66,#073f3c); color:#eaf6f5; border-radius:var(--lr-r); padding:22px 24px; margin-top:18px; box-shadow:var(--lr-shadow); }
  .lipid-landing .lr-feature .lr-eyebrow{ color:#8fd8d2; }
  .lipid-landing .lr-feature h3{ color:#fff; font-size:1.34em; font-weight:700; margin:6px 0 8px; letter-spacing:-.01em; }
  .lipid-landing .lr-feature p{ color:#d6ecea; font-size:1.12em; line-height:1.65; margin:0 0 8px; }
  .lipid-landing .lr-feature strong{ color:#fff; }
  .lipid-landing .lr-chips{ display:flex; flex-wrap:wrap; gap:8px; margin-top:12px; }
  .lipid-landing .lr-chip{ font-family:var(--lr-mono); font-size:.78em; background:rgba(255,255,255,.12); color:#eafffb; border:1px solid rgba(255,255,255,.22); padding:4px 10px; border-radius:999px; }
  .lipid-landing .lr-note{ background:var(--lr-primary-soft); border:1px solid #d3e7e5; border-radius:var(--lr-r-sm); padding:14px 18px; margin-top:16px; font-size:1.1em; color:#274b48; }
  .lipid-landing .lr-note strong{ color:var(--lr-primary-ink); }
  .lipid-landing .lr-notice{ background:var(--lr-amber-soft); border:1px solid #f0d79a; border-radius:var(--lr-r-sm); padding:12px 18px; margin:18px 0; font-size:1.06em; color:#6b4e12; }
  .lipid-landing .lr-notice strong{ color:#8a5a00; }
  .lipid-landing ol{ padding-left:1.3em; }
  .lipid-landing .lr-links{ display:grid; grid-template-columns:1fr 1fr; gap:14px; margin-top:14px; }
  .lipid-landing .lr-linkcard{ background:var(--lr-surface); border:1px solid var(--lr-line); border-radius:var(--lr-r-sm); padding:14px 16px; }
  .lipid-landing .lr-linkcard .lr-k{ font-family:var(--lr-mono); font-size:.72em; text-transform:uppercase; letter-spacing:.1em; color:var(--lr-muted); }
  .lipid-landing .lr-linkcard .lr-v{ font-size:1.14em; margin-top:4px; word-break:break-word; }
  .lipid-landing .lr-linkcard .lr-v.pending{ color:#9aa7ac; font-style:italic; }
  .lipid-landing .lr-cite{ background:#0f171b; color:#c9d6da; border-radius:var(--lr-r-sm); padding:14px 16px; font-family:var(--lr-mono); font-size:.96em; line-height:1.6; white-space:pre-wrap; margin-top:12px; border:1px solid #1d2a30; }
  .lipid-landing .lr-footer{ margin-top:40px; padding-top:16px; border-top:1px solid var(--lr-line); font-size:.85em; color:#93a1a6; text-align:center; font-family:var(--lr-mono); letter-spacing:.02em; }
  @media (max-width:620px){ .lipid-landing .lr-links{ grid-template-columns:1fr; } .lipid-landing h1.lr-brand{ font-size:2.1em; } }
  @media (prefers-reduced-motion:reduce){ .lipid-landing .lr-card{ transition:none; } }
"))

# ── Landing page UI ───────────────────────────────────────────────────────────
landing_page_ui <- function() {
  div(class = "lipid-landing",

      # ── Hero ──
      div(class = "lr-hero",
          div(class = "lr-eyebrow", "MS-DIAL 5 \u00b7 untargeted lipidomics"),
          tags$h1(class = "lr-brand", "LipiRich",
                  tags$span(class = "lr-ver", paste0("v", APP_VERSION))),
          tags$p(class = "lr-tagline",
                 "Normalisation, statistics, and visualisation for MS-DIAL lipidomics data \u2014 in the browser, no code required."),
          div(class = "lr-badges",
              tags$a(href = paste0("https://github.com/sarahehancock/LipiRich/releases/tag/v", APP_VERSION), target = "_blank",
                     tags$img(src = paste0("https://img.shields.io/badge/version-", APP_VERSION, "-0b6b66"), alt = paste0("v", APP_VERSION), style = "height:20px;")),
              tags$a(href = "https://github.com/sarahehancock/LipiRich", target = "_blank",
                     tags$img(src = "https://img.shields.io/badge/GitHub-sarahehancock%2FLipiRich-181717?logo=github", alt = "GitHub", style = "height:20px;")),
              tags$img(src = "https://img.shields.io/badge/license-AGPL--3.0-0b6b66", alt = "AGPL-3.0 License", style = "height:20px;"),
              tags$img(src = "https://img.shields.io/badge/R-%3E%3D4.5.2-informational", alt = "R >= 4.5.2", style = "height:20px;")
          )
      ),

      # ── Pre-publication notice ──
      div(class = "lr-notice",
          tags$strong("\u26a0\ufe0f Pre-publication software (v", APP_VERSION, ")"), tags$br(),
          "LipiRich is under active development; features and outputs may change between versions. ",
          "A citable preprint and demonstration dataset will accompany the first stable release. See the ",
          tags$a(href = "https://github.com/sarahehancock/LipiRich", target = "_blank", "GitHub repository"),
          " for updates and to report issues."
      ),

      # ── Overview ──
      div(class = "lr-section",
          div(class = "lr-eyebrow", "Overview"),
          tags$h2("What LipiRich does"),
          tags$p("LipiRich streamlines post-processing of untargeted lipidomics data exported from MS-DIAL: from raw alignment output through internal-standard normalisation, statistical testing, and lipid set enrichment \u2014 without requiring programming knowledge."),
          tags$p("It reads aligned .txt files exported directly from MS-DIAL 5 and supports positive and negative ion mode, including combined pos/neg experiments. Developed and tested with MS-DIAL 5.5.251021.")
      ),

      # ── Cross-mode feature highlight ──
      div(class = "lr-feature",
          div(class = "lr-eyebrow", "Dual-polarity reconciliation"),
          tags$h3("Structural identity from negative mode, quantitation from positive"),
          tags$p(HTML("For glycerophospholipids and cardiolipin, LipiRich pairs each species across polarities by sum-composition and retention time, then takes the acyl-resolved identity from the negative-mode feature and the abundance from the adduct-matched positive-mode feature \u2014 the way structural lipidomics workflows resolve fatty-acyl detail. A species is quantified once, not double-counted across modes.")),
          tags$p(HTML("Unmatched species are never dropped: a positive-only feature is shown at the depth positive can determine (sum composition, or the two combined halves for CL), a negative-only feature keeps its full acyl identity, and each is quantified from its own-mode standard. Every decision is listed in the <strong>Cross-mode audit</strong> tab, and the whole behaviour is a toggle in <strong>Settings</strong>.")),
          div(class = "lr-chips",
              tags$span(class = "lr-chip", "PC"), tags$span(class = "lr-chip", "PE"),
              tags$span(class = "lr-chip", "PG"), tags$span(class = "lr-chip", "PI"),
              tags$span(class = "lr-chip", "PS"), tags$span(class = "lr-chip", "PA"),
              tags$span(class = "lr-chip", "CL")
          )
      ),

      # ── Workflow ──
      div(class = "lr-section",
          div(class = "lr-eyebrow", "Analytical workflow"),
          tags$h2("From alignment file to enrichment, in one session"),
          div(class = "lr-grid",
              div(class = "lr-card", div(class = "lr-step", "Step 1"), tags$h3("Import & cleaning"), tags$p("Upload 1\u20132 MS-DIAL .txt files. Data are reshaped, ion modes resolved (including cross-mode identity transfer), and duplicate features deduplicated.")),
              div(class = "lr-card", div(class = "lr-step", "Step 2"), tags$h3("Internal standards"), tags$p("QC standard signal by class, shown in both ion modes by default; hover any point for the exact standard name. Per-ISTD amount and units via CSV.")),
              div(class = "lr-card", div(class = "lr-step", "Step 3"), tags$h3("Normalisation"), tags$p("Blank background subtraction, IS-based quantitative normalisation, and optional unit-aware protein normalisation.")),
              div(class = "lr-card", div(class = "lr-step", "Step 4\u20135"), tags$h3("Visualisation"), tags$p("Interactive bar plots for a single species or every species in a class, with cascading ion-mode, adduct, and name filters and flexible grouping.")),
              div(class = "lr-card", div(class = "lr-step", "Step 6"), tags$h3("Export"), tags$p("Download wide or long CSV, filtered by class and value type (absolute or % of class total).")),
              div(class = "lr-card", div(class = "lr-step", "Step 7"), tags$h3("PCA"), tags$p("Principal component analysis with group colouring, iQC overlay, and a loadings table.")),
              div(class = "lr-card", div(class = "lr-step", "Step 8"), tags$h3("Statistics"), tags$p("Unpaired t-test, one- and two-way ANOVA with multiple-testing correction, post-hoc tests, per-class summary, and downloadable results.")),
              div(class = "lr-card", div(class = "lr-step", "Step 8a"), tags$h3("Volcano plot"), tags$p("Interactive volcano of all tested features, coloured by direction, with significant labels repelled automatically.")),
              div(class = "lr-card", div(class = "lr-step", "Step 9"), tags$h3("Class bar plots"), tags$p("Faceted bars of every species in a class, with abundance-range filtering, significance highlighting, and group overlays.")),
              div(class = "lr-card", div(class = "lr-step", "Step 9b"), tags$h3("Heatmap"), tags$p("Z-scored heatmap of significant features with configurable clustering, palettes, and label toggles.")),
              div(class = "lr-card", div(class = "lr-step", "Step 10"), tags$h3("Enrichment"), tags$p("Lipid set enrichment (ORA or FGSEA) across class, fatty-acid identity, saturation, chain length, and ether subclass.")),
              div(class = "lr-card", div(class = "lr-step", "Step 11"), tags$h3("Correlation network"), tags$p("Pearson correlation network of significant species; node colour reflects direction, edge width reflects correlation strength.")),
              div(class = "lr-card", div(class = "lr-step", "Step 12"), tags$h3("Synthesis pathways"), tags$p("Enzyme-activity proxy scores (class ratios) as a z-scored heatmap with per-score group statistics \u2014 interpretive indicators, not flux."))
          )
      ),

      # ── Data prep ──
      div(class = "lr-section",
          div(class = "lr-eyebrow", "Before you upload"),
          tags$h2("Preparing your data in MS-DIAL"),
          tags$p("LipiRich reads the aligned result exported from MS-DIAL 5 as a tab-delimited .txt file. In MS-DIAL:"),
          tags$ol(
            tags$li(tags$strong("Run alignment"), " \u2014 complete peak picking and alignment as normal. Review the result, tag correctly identified species with the ", tags$strong("\u2713 checkmark"), ", then choose ", tags$strong("Filter by current parameter"), " on export to pass only checked features to LipiRich."),
            tags$li(tags$strong("Export"), " \u2014 ", tags$em("Export \u2192 Alignment result"), " in ", tags$strong(".txt"), " format. LipiRich reads from row 5 (the four MS-DIAL header rows are skipped automatically)."),
            tags$li(tags$strong("Sample naming"), " \u2014 end sample columns in ", tags$code("_pos"), " or ", tags$code("_neg"), " (e.g. ", tags$code("Sample1_pos"), "). These suffixes drive ion-mode resolution and sample-name deduplication."),
            tags$li(tags$strong("Internal standards"), " \u2014 standard features must contain ", tags$code("[IS]"), " in the metabolite name, and their adduct must match the analytes they normalise (class + ion mode + adduct). These are detected automatically."),
            tags$li(tags$strong("Both polarities"), " \u2014 export each mode as its own .txt and upload both together (max 2 files); LipiRich merges them."),
            tags$li(tags$strong("Blank samples"), " \u2014 include at least one sample named ", tags$code("Blank"), " (case-insensitive) for background subtraction.")
          )
      ),

      div(class = "lr-note",
          tags$strong("Grouping: "), "group samples for plots and statistics by sample name, by a delimiter token (", tags$code("Sample1_rep1"), " \u2192 ", tags$code("Sample1"), "), or by a regex capture group \u2014 or upload a grouping CSV with ", tags$code("sample"), " and ", tags$code("group"), " columns (plus optional ", tags$code("factorA"), "/", tags$code("factorB"), " for two-way designs)."
      ),
      div(class = "lr-note",
          tags$strong("Technical replicates: "), "average technical injections per biological sample before PCA, statistics, and plots \u2014 via a ", tags$code("bio_sample"), " column in the grouping CSV, or by selecting sample-name tokens in Group Preview (including combined indices like ", tags$code("1-1"), " for bio 1, tech 1)."
      ),

      # ── Access ──
      div(class = "lr-section",
          div(class = "lr-eyebrow", "Access"),
          tags$h2("Run it, or read the source"),
          div(class = "lr-links",
              div(class = "lr-linkcard", div(class = "lr-k", "Live web app"), div(class = "lr-v", tags$a(href = "https://lipirich.sarahehancock.com", target = "_blank", "lipirich.sarahehancock.com"))),
              div(class = "lr-linkcard", div(class = "lr-k", "GitHub"), div(class = "lr-v", tags$a(href = "https://github.com/sarahehancock/LipiRich", target = "_blank", "github.com/sarahehancock/LipiRich")))
          ),
          div(class = "lr-note", style = "margin-top:14px;",
              tags$strong("\U0001F512 Data privacy: "), "LipiRich does not store, transmit, or retain uploaded data. Files and results exist only in your browser session and are discarded when the tab or session closes."
          )
      ),

      # ── Citation ──
      div(class = "lr-section",
          div(class = "lr-eyebrow", "Publication & citation"),
          tags$h2("Citing LipiRich"),
          div(class = "lr-links",
              div(class = "lr-linkcard", div(class = "lr-k", "DOI"), div(class = "lr-v pending", "Pending publication")),
              div(class = "lr-linkcard", div(class = "lr-k", "Preprint"), div(class = "lr-v pending", "Pending upload")),
              div(class = "lr-linkcard", div(class = "lr-k", "Repository"), div(class = "lr-v", tags$a(href = "https://github.com/sarahehancock/LipiRich", target = "_blank", "github.com/sarahehancock/LipiRich"))),
              div(class = "lr-linkcard", div(class = "lr-k", "License"), div(class = "lr-v", "AGPL-3.0"))
          ),
          tags$p(style = "margin-top:14px;", "If you use LipiRich in your research, please cite:"),
          div(class = "lr-cite",
              paste0("Hancock, SE. (2026). LipiRich: A Shiny application for normalisation,\nstatistics, and visualisation of MS-DIAL lipidomics data (v", APP_VERSION, ").\nGitHub: https://github.com/sarahehancock/LipiRich\nDOI: [pending]")
          )
      ),

      div(class = "lr-footer",
          paste0("LipiRich v", APP_VERSION, "  \u00b7  R Shiny  \u00b7  MS-DIAL 5.5.251021 compatible  \u00b7  AGPL-3.0")
      )
  )
}

# UI
ui <- fluidPage(
  landing_css,
  enrichment_spinner_css,
  idle_timeout_js,
  titlePanel("LipiRich — Normalisation, statistics and visualisation for MS-DIAL lipidomics data"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      
      fileInput(
        inputId = "msdial_txts",
        label   = "Upload MS-DIAL .txt files (1–2)",
        multiple = TRUE,
        accept   = c(".txt"),
        buttonLabel = "Browse...",
        placeholder = "Select 1 or 2 .txt files"
      ),
      helpText("Upload aligned data exported from MS-DIAL. Min 1 file, max 2"),
      # Show uploaded file names BEFORE processing
      textOutput("uploadedNames"),
      
      # Show processed file names AFTER loading
      textOutput("loadedFiles"),
      tags$hr(),
      
      # --- Upload ISTD amount/units map (per internal standard name) ---
      fileInput(
        inputId = "istd_map_csv",
        label   = "Upload ISTD amount/units CSV",
        multiple = FALSE,
        accept   = c(".csv"),
        buttonLabel = "Browse..."
      ),
      checkboxInput(
        "use_istd_map_csv",
        "Use uploaded ISTD CSV for normalization",
        value = TRUE
      ),
      helpText(
        "CSV must contain exactly three columns: 'ISTD' (must match MSDIAL [IS] name),",
        " 'amount' (numeric), and 'units' (e.g., pmol, nmol). If missing or mismatched,",
        " the app will fall back to the single ISTD amount input below.",
        "\n If using the single ISTD amount method Metabolite names must contain '[IS]' to be matched"
      ),
      verbatimTextOutput("istd_map_summary"),
      
      # --- ISTD amount (leave as-is) ---
      numericInput("ISTD_vol", "ISTD amount (pmol):", 100, min = 1),
      
      tags$hr(),
      
      # --- Protein normalisation (optional) ---
      tags$h4("Protein normalisation (optional)"),
      helpText(tags$small(
        "Upload a CSV with columns ", tags$code("sample"), " and ", tags$code("protein"),
        " to enable protein content normalisation. Sample names should match",
        " the MS-DIAL column names (without _pos/_neg suffix).",
        " When enabled, normalised lipid values are divided by protein content."
      )),
      fileInput(
        inputId     = "protein_csv",
        label       = "Upload protein content CSV",
        multiple    = FALSE,
        accept      = c(".csv"),
        buttonLabel = "Browse..."
      ),
      checkboxInput(
        "use_protein_norm",
        "Apply protein normalisation",
        value = FALSE
      ),
      conditionalPanel(
        condition = "input.use_protein_norm == true",
        selectInput(
          "protein_units",
          "Protein units (as uploaded in CSV):",
          choices  = c("\u00b5g" = "ug", "mg" = "mg", "g" = "g"),
          selected = "mg"
        ),
        helpText(tags$small(
          "Used only to label the y-axis correctly (e.g. \u2018pmol/mg protein\u2019).",
          " It does not rescale your uploaded protein values \u2014 make sure the CSV",
          " values are already in the unit selected here."
        ))
      ),
      verbatimTextOutput("protein_csv_summary"),
      fileInput(
        inputId = "group_csv",
        label   = "Upload grouping CSV (sample, group[, factorA, factorB, bio_sample])",
        multiple = FALSE,
        accept   = c(".csv"),
        buttonLabel = "Browse..."
      ),
      checkboxInput(
        "use_group_csv",
        "Use uploaded CSV for grouping (and factors)",
        value = TRUE
      ),
      helpText("CSV must contain columns 'sample' and 'group'; optional 'factorA', 'factorB', 'bio_sample'. Matching is case-insensitive on canonicalized sample names."),
      verbatimTextOutput("group_csv_summary"),
      tags$hr(),
      # ── Technical replicate averaging ───────────────────────────────────────
      tags$h4("Technical replicates (optional)"),
      helpText(tags$small(
        "If your dataset includes multiple technical injections per biological sample,",
        " average them before PCA, statistics, and downstream plots. Identify the parent",
        " biological sample either via a ", tags$code("bio_sample"), " column in the grouping",
        " CSV above, or by selecting sample-name token(s) in the ",
        tags$strong("Group Preview"), " tab (a bio/tech sub-delimiter option there handles",
        " composite tokens like ", tags$code("1-1"), " for bio 1, tech 1)."
      )),
      checkboxInput(
        "use_tech_rep_avg",
        "Average technical replicates before analysis",
        value = FALSE
      ),
      conditionalPanel(
        condition = "input.use_tech_rep_avg == true",
        radioButtons(
          "tech_rep_source",
          "Identify technical reps via:",
          choices  = c("Grouping CSV column ('bio_sample')" = "csv",
                       "Sample-name token(s)"                = "token"),
          selected = "token"
        ),
        helpText(tags$small(
          "iQC, ISTD, and Blank samples are never averaged \u2014 each injection is kept separate."
        ))
      ),
      tags$hr(),
      # ── Token-based grouping ──────────────────────────────────────────────────
      tags$h4("Sample name parsing"),
      helpText(tags$small(
        "Set the delimiter and grouping mode here. Use the ",
        tags$strong("Group Preview"), " tab to preview tokens and assign groups."
      )),
      radioButtons("lr_group_mode", "Grouping mode:",
                   choices = c("Delimiter-based" = "delimiter", "Regex" = "regex"),
                   selected = "delimiter", inline = TRUE),
      conditionalPanel(
        condition = "input.lr_group_mode == 'delimiter'",
        textInput("lr_delimiter", "Delimiter", value = "_",
                  placeholder = "e.g. _ or - or .")
      ),
      conditionalPanel(
        condition = "input.lr_group_mode == 'regex'",
        textInput("lr_group_regex", "Regex (first capture group = group)",
                  value = "^([^_]+)",
                  placeholder = "e.g. ^([A-Za-z]+)"),
        checkboxInput("lr_regex_ignore_case", "Ignore case", value = TRUE)
      ),
      tags$hr(),
      # Click to start processing the uploaded files
      actionButton("load_data", "Load and Process Data", icon = icon("play-circle")),
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          title = tagList(tags$span("\u2139\ufe0f About LipiRich")),
          value = "about",
          br(),
          landing_page_ui()
        ),
        
        # ── Settings tab ──────────────────────────────────────────────────────
        tabPanel(
          title = tagList(tags$span("\u2699\ufe0f Settings")),
          value = "settings",
          h4("Ion mode preference settings"),
          div(class = "info-box", style = "margin-bottom:16px;",
              tags$strong("How ion mode preference works:"),
              tags$br(),
              "When the same lipid species is detected in both positive and negative ion mode, ",
              "LipiRich keeps ", tags$strong("one measurement per species"), " by selecting the ",
              "analytically preferred ion mode. The preferred mode is chosen because it typically ",
              "provides ", tags$strong("higher structural information"), " for that class — for example, ",
              "negative mode provides fatty acid pairing information for phospholipids (e.g. PC 16:0/18:1), ",
              "whereas positive mode typically yields only sum composition (e.g. PC 34:1). ",
              "Where a species is detected in only one mode, that measurement is always retained. ",
              "For classes not listed in either box below, both measurements are kept as separate rows.",
              tags$br(), tags$br(),
              tags$strong("Adding classes:"), " you can type any MS-DIAL class name directly into the ",
              "text boxes below (one per line). Class names must match the MS-DIAL ",
              tags$code("Lipid class"), " column exactly, e.g. ", tags$code("BMP"), ", ",
              tags$code("HBMP"), ", ", tags$code("Hex2Cer"), ", ", tags$code("GM3"), " etc. ",
              "Changes take effect immediately when data is reloaded."
          ),
          fluidRow(
            column(width = 4,
                   wellPanel(
                     h5("Prefer negative mode"),
                     helpText(tags$small(
                       "When a species is detected in both modes, the negative ion measurement is kept.",
                       "The positive measurement is only retained if no negative measurement exists for that species.",
                       "One class name per line — must match MS-DIAL class names exactly."
                     )),
                     textAreaInput(
                       "pref_neg_classes",
                       label = NULL,
                       value = paste(c(
                         "PC", "PE", "PG", "PI", "PS", "PA", "CL",
                         "PC-O", "PE-O"
                       ), collapse = "\n"),
                       rows = 10,
                       width = "100%"
                     )
                   )
            ),
            column(width = 4,
                   wellPanel(
                     h5("Prefer positive mode"),
                     helpText(tags$small(
                       "When a species is detected in both modes, the positive ion measurement is kept.",
                       "The negative measurement is only retained if no positive measurement exists for that species.",
                       "One class name per line — must match MS-DIAL class names exactly."
                     )),
                     textAreaInput(
                       "pref_pos_classes",
                       label = NULL,
                       value = paste(c(
                         "TG", "DG", "MG", "CE", "Cer", "HexCer", "SM",
                         "LPC", "LPE", "LPG", "LPI", "LPS", "LPA"
                       ), collapse = "\n"),
                       rows = 10,
                       width = "100%"
                     )
                   )
            ),
            column(width = 4,
                   wellPanel(
                     h5("Keep both modes"),
                     helpText(tags$small(
                       "Any class not listed in either box above is treated as 'no preference'.",
                       "If a species is detected in both modes, ", tags$strong("both measurements are kept"),
                       " as separate rows (distinguished by adduct type in the output).",
                       "This is appropriate for classes where both ion modes provide complementary,",
                       "non-redundant information."
                     )),
                     br(),
                     h5("IS matching"),
                     helpText(tags$small(
                       "Internal standards are always retained from whichever ion mode",
                       "they were detected in, regardless of the settings above.",
                       "When an ISTD CSV is provided, IS are matched to analytes by",
                       tags$strong("class + ion mode + adduct type."),
                       "If multiple ISTD entries match different adducts, one normalised",
                       "row is produced per matched adduct."
                     )),
                     br(),
                     actionButton("settings_reset", "Reset to defaults",
                                  icon = icon("rotate-left"),
                                  class = "btn-sm btn-default")
                   )
            )
          ),
          tags$hr(),
          h4("Cross-mode identity (negative ID \u2192 positive quantitation)"),
          div(class = "info-box", style = "margin-bottom:16px;",
              tags$strong("What this does:"),
              tags$br(),
              "For the classes listed below, LipiRich takes the ",
              tags$strong("quantitation from the positive-mode feature"),
              " (which has an adduct-matched positive internal standard) and transfers the ",
              tags$strong("acyl-resolved identity from the matching negative-mode feature"),
              ", pairing them by sum-composition shorthand plus retention time. ",
              "This mirrors the LipidSearch-style workflow: negative mode resolves the fatty ",
              "acids (all four chains for CL), positive mode gives the more consistent precursor ",
              "for quantitation. ",
              tags$br(), tags$br(),
              "A positive feature with no negative match is kept at sum composition; a ",
              "negative feature with no positive match is retained and quantified from its own ",
              "negative-mode IS. When the toggle is off, these classes revert to prefer-negative ",
              "handling. Internal standards are never touched. Every decision is listed in the ",
              "audit table so the reassignment is fully traceable."
          ),
          fluidRow(
            column(width = 5,
                   wellPanel(
                     checkboxInput("cross_mode_enable",
                                   "Enable cross-mode identity transfer",
                                   value = TRUE),
                     helpText(tags$small(
                       "One class per line — must match the MS-DIAL class exactly.",
                       "Only classes with adduct-matched positive internal standards",
                       "and positive analytes should be listed here."
                     )),
                     textAreaInput(
                       "cross_mode_classes",
                       label = "Classes: quant positive, identity negative",
                       value = paste(c("PC","PE","PG","PI","PS","PA","CL"), collapse = "\n"),
                       rows = 6, width = "100%"
                     ),
                     numericInput(
                       "cross_mode_rt_tol",
                       "RT match tolerance (min)",
                       value = 0.1, min = 0, max = 2, step = 0.01
                     )
                   )
            ),
            column(width = 7,
                   div(class = "info-box", style = "margin-top:8px;",
                       tags$strong("Reconciliation audit"),
                       " has its own tab for a wider view \u2014 see the ",
                       tags$strong("\U1F4CB Cross-mode audit"), " tab. Every re-identification, ",
                       "sum-composition retention, and negative-only feature is listed there, ",
                       "and it reloads whenever the data or these settings change."
                   )
            )
          )
        ),
        
        # ── Cross-mode audit tab ──────────────────────────────────────────────
        tabPanel(
          title = tagList(tags$span("\U1F4CB Cross-mode audit")),
          value = "cross_mode_audit",
          h4("Cross-mode identity reconciliation audit"),
          div(class = "info-box", style = "margin-bottom:16px;",
              "Each row records how a designated-class feature (configured under ",
              tags$strong("\u2699\ufe0f Settings \u2192 Cross-mode identity"),
              ") was resolved. ",
              tags$strong("Final identity"), " is the label the feature carries downstream and ",
              tags$strong("Quant mode"), " is the ion mode its abundance is taken from: positive ",
              "features re-identified from the matching negative feature (",
              tags$em("pos quant / neg id"), "), positive features kept at the depth positive can ",
              "determine (", tags$em("sum composition"), ", or positive substructure for CL), and ",
              "negative-only features retained and negative-quantified. Reloads with the data."
          ),
          DTOutput("crossModeAuditTable")
        ),
        
        # ── Group Preview tab ─────────────────────────────────────────────────
        tabPanel(
          title = tagList(tags$span("\U1F50D Group Preview")),
          value = "group_preview",
          h4("Sample name tokeniser — identify group positions"),
          fluidRow(
            column(width = 4,
                   wellPanel(
                     helpText(
                       "Upload your MS-DIAL file(s) then click Preview to see how sample names",
                       "split into tokens. Select which token position(s) define your groups."
                     ),
                     actionButton("lr_preview_groups", "Preview tokens",
                                  icon = icon("eye"), class = "btn-primary"),
                     tags$hr(),
                     # ── Token group assignment ────────────────────────────────────
                     h5("Group token selection"),
                     helpText(tags$small(
                       "After previewing, select which token position(s) to use as the group label.",
                       "Select multiple to combine (e.g. tokens 1+2 from KO_treated_rep1 → KO_treated)."
                     )),
                     uiOutput("lr_token_selector_ui"),
                     textInput("lr_group_sep", "Join multiple tokens with", value = "_"),
                     tags$hr(),
                     # ── Two-way design ────────────────────────────────────────────
                     h5("Two-way design (optional)"),
                     helpText(tags$small(
                       "Optionally select token positions for Factor A and Factor B.",
                       "Only used when two-way ANOVA is selected in Statistics."
                     )),
                     uiOutput("lr_factorA_selector_ui"),
                     uiOutput("lr_factorB_selector_ui"),
                     tags$hr(),
                     # ── Technical replicate parent sample ──────────────────────────
                     h5("Technical replicate parent sample (optional)"),
                     helpText(tags$small(
                       "Select token position(s) that identify the biological/parent sample",
                       " (i.e. excluding the replicate-index token, e.g. tokens 1+2 from",
                       " KO_treated_rep1 \u2192 KO_treated). Used when \u2018Average technical",
                       " replicates\u2019 is enabled in the sidebar."
                     )),
                     uiOutput("lr_techrep_token_selector_ui"),
                     textInput(
                       "lr_techrep_subdelim",
                       "Bio/tech sub-delimiter within token (optional)",
                       value = "",
                       placeholder = "e.g. - for '1-1' meaning bio 1, tech 1"
                     ),
                     helpText(tags$small(
                       "Use this if a single token encodes both the biological and technical",
                       " replicate together, e.g. ", tags$code("1-1"), " = bio rep 1, tech rep 1.",
                       " Enter the separator (", tags$code("-"), " in that example) and only the",
                       " part ", tags$strong("before"), " it is kept as the biological replicate",
                       " index. Leave blank if your tech-rep token doesn't need splitting."
                     ))
                   )
            ),
            column(width = 8,
                   h5("Token table — each column shows one token position across all samples"),
                   DT::DTOutput("lr_group_preview_tbl"),
                   br(),
                   uiOutput("lr_token_summary_ui")
            )
          )
        ),
        tabPanel("Cleaned Wide Data",
                 h4("Step 1: Cleaned and Reshaped MS-DIAL Output"),
                 helpText(tags$small(
                   "Use the ", tags$strong("Group Preview"), " tab (above) to check",
                   "how your sample names split into tokens before running analyses."
                 )),
                 DTOutput("wideData")
        ),
        tabPanel("Internal Standards",
                 h4("Step 2: Internal Standards ([IS])"),
                 DTOutput("ISTDData")
        ),
        tabPanel("IS Plots",
                 h4("Step 2a: Internal Standard Values by Class"),
                 fluidRow(
                   column(
                     width = 4,
                     wellPanel(
                       checkboxInput("exclude_blank_is", "Exclude samples named 'Blank'", value = FALSE),
                       selectInput("is_class", "Lipid Class", choices = "Loading..."),
                       selectInput(
                         "grouping_method_is",
                         "How to group samples:",
                         choices = c("By sample", "Delimiter-based", "Regex capture group"),
                         selected = "By sample"
                       ),
                       conditionalPanel(
                         condition = "input.grouping_method_is == 'Delimiter-based'",
                         selectInput(
                           "group_delim_is", "Delimiter",
                           choices = c("_" = "_", "-" = "-", "space" = " ", "." = ".", "/" = "/"),
                           selected = "_"
                         ),
                         numericInput("group_tokens_is", "Use first N tokens as group", value = 1, min = 1, max = 5, step = 1)
                       ),
                       conditionalPanel(
                         condition = "input.grouping_method_is == 'Regex capture group'",
                         textInput("group_regex_is", "Regex with ONE capture group", value = "^([^_]+)")
                       ),
                       radioButtons(
                         "is_mode",
                         "Ion mode to display",
                         choices = c(
                           "Auto (prefer negative for CL)" = "auto",
                           "Negative only" = "neg",
                           "Positive only" = "pos",
                           "Both (separate)" = "both"
                         ),
                         selected = "both"
                       ),
                       selectInput(
                         "error_type_is", "Error bars",
                         choices = c("SEM", "SD", "95% CI"), selected = "SEM"
                       )
                     )
                   ),
                   column(
                     width = 8,
                     plotlyOutput("isPlot", height = "520px"),
                     br(),
                     h5("Summary table (group-level IS statistics):"),
                     DTOutput("isPlotSummary")
                   )
                 )
        ),
        tabPanel("Background & Normalisation",
                 h4("Step 3: Background Subtraction & Normalisation"),
                 DTOutput("bgNormTable")
        ),
        tabPanel("Protein Match",
                 h4("Step 3a: Protein Normalisation Sample Matching"),
                 helpText(tags$small(
                   "Checks whether each imported MS-DIAL sample (", tags$code("_pos"), "/",
                   tags$code("_neg"), " suffix stripped) has a matching row in the uploaded",
                   " protein content CSV, and flags any protein CSV rows with no corresponding sample.",
                   " ISTD, Blank, and iQC/QC samples are excluded from this check on both sides."
                 )),
                 verbatimTextOutput("protein_match_summary"),
                 br(),
                 DT::DTOutput("proteinMatchTable")
        ),
        tabPanel("Outlier Detection",
                 h4("Step 3b: Outlier Detection"),
                 helpText(tags$small(
                   "Flags potential outlier ", tags$strong("samples"), " (PCA Hotelling's T\u00b2 and iQC replicate deviation)",
                   " and potential outlier ", tags$strong("data points"), " within a lipid feature \u00d7 group",
                   " (modified Z-score or IQR rule), and applies any exclusions to every downstream tab",
                   " from ", tags$strong("Plot single lipid"), " onward. ISTD, Blank, and iQC/QC-named samples",
                   " are never flagged or excluded here \u2014 they're outside the scope of biological outlier",
                   " detection and are handled by their own dedicated logic elsewhere in the app."
                 )),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       h5("Sample-level \u2014 PCA"),
                       selectInput(
                         "outlier_sample_conf", "Confidence level",
                         choices  = c("95%" = "0.95", "97.5%" = "0.975", "99%" = "0.99"),
                         selected = "0.975"
                       ),
                       checkboxInput("outlier_exclude_samples",
                                     "Exclude all auto-flagged samples from downstream analyses",
                                     value = FALSE),
                       tags$hr(),
                       h5("Sample-level \u2014 iQC deviation"),
                       helpText(tags$small(
                         "Requires \u2265 3 iQC replicate samples; skipped otherwise. Informational only \u2014",
                         " iQC samples themselves are never excluded from downstream analyses."
                       )),
                       numericInput("outlier_iqc_mad", "MAD multiplier", value = 3, min = 1, max = 10, step = 0.5),
                       tags$hr(),
                       h5("Feature-level \u2014 per lipid \u00d7 group"),
                       radioButtons(
                         "outlier_feature_method", "Method",
                         choices  = c("Modified Z-score (MAD)" = "mad", "IQR rule" = "iqr"),
                         selected = "mad"
                       ),
                       numericInput("outlier_feature_mult", "Threshold multiplier", value = 3, min = 1, max = 10, step = 0.5),
                       helpText(tags$small(
                         tags$strong("Modified Z:"), " flags points where |value \u2212 median| / MAD exceeds the multiplier",
                         " (3\u20133.5 is a common default).",
                         tags$br(),
                         tags$strong("IQR:"), " flags points beyond Q1/Q3 \u00b1 multiplier \u00d7 IQR (1.5 = standard boxplot rule)."
                       )),
                       checkboxInput("outlier_exclude_features",
                                     "Exclude all auto-flagged points from downstream analyses",
                                     value = FALSE),
                       helpText(tags$small(
                         "Excluded points/samples are removed from the shared dataset used by every",
                         " other tab, not deleted from your file \u2014 turn a toggle back off, or use the",
                         " manual review tables, to restore them at any time."
                       )),
                       tags$hr(),
                       h5("Export plots"),
                       numericInput("outlier_export_width",    "Width (px)",      1200, 400, 4000, 50),
                       numericInput("outlier_export_height",   "Height (px)",      500, 300, 4000, 50),
                       numericInput("outlier_export_dpi",      "DPI",              300,  72,  600, 12),
                       numericInput("outlier_export_scale",    "Scale fraction",  1.00, 0.25, 2.00, 0.05),
                       numericInput("outlier_export_fontsize", "Base font size",    13,    6,   24,  1),
                       helpText(tags$small("Sample-level: PCA Hotelling's T\u00b2")),
                       fluidRow(
                         column(6, downloadButton("download_outlier_sample_png", "PNG", class = "btn-primary btn-sm")),
                         column(6, downloadButton("download_outlier_sample_svg", "SVG", class = "btn-sm"))
                       ),
                       br(),
                       helpText(tags$small("Sample-level: iQC replicate deviation")),
                       fluidRow(
                         column(6, downloadButton("download_outlier_iqc_png", "PNG", class = "btn-primary btn-sm")),
                         column(6, downloadButton("download_outlier_iqc_svg", "SVG", class = "btn-sm"))
                       )
                     )
                   ),
                   column(
                     width = 9,
                     h5("Sample-level: PCA Hotelling's T\u00b2"),
                     plotOutput("outlierSamplePlot", height = "320px"),
                     tags$hr(),
                     h5("Sample-level: iQC replicate deviation"),
                     plotOutput("outlierIqcPlot", height = "280px"),
                     tags$hr(),
                     h5("Sample review & individual exclusion"),
                     helpText(tags$small(
                       "Select one or more rows below, then use a button to override the automatic",
                       " decision for those specific samples \u2014 independently of the toggle above.",
                       " ISTD/Blank/iQC samples are not listed here."
                     )),
                     fluidRow(
                       column(width = 4, actionButton("sample_manual_exclude_btn", "Exclude selected",
                                                      icon = icon("ban"), class = "btn-sm btn-danger")),
                       column(width = 4, actionButton("sample_manual_keep_btn", "Keep selected (override)",
                                                      icon = icon("check"), class = "btn-sm btn-success")),
                       column(width = 4, actionButton("sample_manual_clear_btn", "Clear manual overrides",
                                                      icon = icon("rotate-left"), class = "btn-sm"))
                     ),
                     br(),
                     DTOutput("sampleReviewTable"),
                     tags$hr(),
                     h5("Feature-level outliers: review & individual exclusion"),
                     verbatimTextOutput("outlier_feature_summary"),
                     helpText(tags$small(
                       "Select one or more rows below, then use a button to override the automatic",
                       " decision for those specific points \u2014 independently of the toggle above."
                     )),
                     fluidRow(
                       column(width = 4, actionButton("feature_manual_exclude_btn", "Exclude selected",
                                                      icon = icon("ban"), class = "btn-sm btn-danger")),
                       column(width = 4, actionButton("feature_manual_keep_btn", "Keep selected (override)",
                                                      icon = icon("check"), class = "btn-sm btn-success")),
                       column(width = 4, actionButton("feature_manual_clear_btn", "Clear manual overrides",
                                                      icon = icon("rotate-left"), class = "btn-sm"))
                     ),
                     br(),
                     DTOutput("outlierFeatureTable")
                   )
                 )
        ),
        tabPanel("Plot single lipid",
                 h4("Step 4: Plot Individual Metabolites or Total by Class"),
                 fluidRow(
                   column(
                     width = 4,
                     wellPanel(
                       checkboxInput("exclude_blank_met", "Exclude samples named 'Blank'", value = TRUE),
                       radioButtons(
                         "plot_value_type_met",
                         "Data to display:",
                         choices = c(
                           "Quantified (IS-normalised)" = "norm",
                           "Background-subtracted" = "value_bs"
                         ),
                         selected = "norm"
                       ),
                       radioButtons(
                         "display_mode_met",
                         "Units:",
                         choices = c(
                           "Absolute" = "absolute",
                           "Percent of class total" = "percent"
                         ),
                         selected = "absolute"
                       ),
                       selectInput("met_class", "Lipid Class", choices = "Loading..."),
                       radioButtons("met_ion_mode", "Ion mode",
                                    choices = c("All modes" = "all",
                                                "Positive"  = "positive",
                                                "Negative"  = "negative"),
                                    selected = "all", inline = TRUE),
                       selectInput("met_name", "Metabolite name", choices = "Select a class first"),
                       selectInput("met_adduct_filter", "Adduct filter",
                                   choices = c("All" = "all"),
                                   selected = "all"),
                       helpText(tags$small(
                         "Ion mode, adduct, and name cascade: the name list only shows",
                         "species that have data under the current mode and adduct, so a",
                         "selection always plots. After cross-mode reconciliation most",
                         "species resolve to a single quant mode; use these to inspect",
                         "classes kept in both modes."
                       )),
                       tags$hr(),
                       h5("Groups to plot"),
                       helpText(tags$small(
                         "Groups are inherited from the Group Preview tab.",
                         "Deselect to exclude groups from this plot."
                       )),
                       uiOutput("met_group_select_ui"),
                       tags$hr(),
                       selectInput(
                         "error_type_met", "Error bars",
                         choices = c("SEM", "SD", "95% CI"), selected = "SEM"
                       )
                     )
                   ),
                   column(
                     width = 8,
                     plotlyOutput("metPlot", height = "520px"),
                     br(),
                     h5("Summary table (group-level statistics):"),
                     DTOutput("metPlotSummary")
                   )
                 )
        ),
        tabPanel("Plot lipids by class",
                 h4("Step 5: Plot All Metabolites in a Selected Class by Group"),
                 fluidRow(
                   column(
                     width = 4,
                     wellPanel(
                       checkboxInput("exclude_blank_all", "Exclude samples named 'Blank'", value = TRUE),
                       radioButtons(
                         "plot_value_type_all",
                         "Data to display:",
                         choices = c(
                           "Quantified (IS-normalised)" = "norm",
                           "Background-subtracted" = "value_bs"
                         ),
                         selected = "norm"
                       ),
                       radioButtons(
                         "display_mode_all",
                         "Units:",
                         choices = c(
                           "Absolute" = "absolute",
                           "Percent of class total" = "percent"
                         ),
                         selected = "absolute"
                       ),
                       selectInput("class_all", "Lipid Class", choices = "Loading..."),
                       selectInput("class_all_adduct_filter", "Adduct filter",
                                   choices = c("All" = "all"),
                                   selected = "all"),
                       helpText(tags$small(
                         "Filter to a specific adduct type. Auto-populates from the selected class."
                       )),
                       selectInput(
                         "class_all_order",
                         "Order metabolites by:",
                         choices = c(
                           "Abundance (high → low)" = "abundance_desc",
                           "Alphabetical (A → Z)"   = "alpha_asc"
                         ),
                         selected = "abundance_desc"
                       ),
                       tags$hr(),
                       h5("Group to display"),
                       helpText(tags$small(
                         "Select one group to display all lipids within the chosen class.",
                         "Groups are inherited from the Group Preview tab."
                       )),
                       uiOutput("class_all_group_select_ui"),
                       tags$hr(),
                       selectInput(
                         "error_type_all", "Error bars",
                         choices = c("SEM", "SD", "95% CI"), selected = "SEM"
                       )
                     )
                   ),
                   column(
                     width = 8,
                     plotlyOutput("classAllPlot", height = "auto"),
                     br(), br(), br(),
                     h5("Summary table (per-metabolite statistics for selected group):"),
                     DTOutput("classAllSummary")
                   )
                 )
        ),
        tabPanel("Export - wide",
                 h4("Step 6: Export wide-format data"),
                 fluidRow(
                   column(
                     width = 4,
                     wellPanel(
                       checkboxInput("exclude_blank_export", "Exclude samples named 'Blank'", value = TRUE),
                       radioButtons(
                         "export_value_type",
                         "Data to include:",
                         choices = c(
                           "Quantified (IS-normalised)" = "norm",
                           "Background-subtracted"      = "value_bs"
                         ),
                         selected = "norm"
                       ),
                       radioButtons(
                         "export_display_mode",
                         "Units:",
                         choices = c(
                           "Absolute"               = "absolute",
                           "Percent of class total" = "percent"
                         ),
                         selected = "absolute"
                       ),
                       div(
                         style = "display:flex; align-items:center; justify-content:space-between; margin-top:6px;",
                         tags$label("Select classes to include:"),
                         actionButton("toggle_classes", "Select/Deselect All", icon = icon("check-square"))
                       ),
                       checkboxGroupInput(
                         "export_classes",
                         label = NULL,
                         choices = c(),
                         selected = NULL
                       ),
                       br(),
                       downloadButton("download_export_wide", "Download wide CSV"),
                       br(), br(),
                       checkboxInput("export_include_metadata",
                                     "Include group column in wide export header",
                                     value = TRUE),
                       br(),
                       downloadButton("download_export_long", "Download long-format CSV")
                     )
                   ),
                   column(
                     width = 8,
                     DTOutput("exportWideDT")
                   )
                 )
        ),
        tabPanel("PCA Analysis",
                 h4("Step 7: Principal Component Analysis (PCA)"),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       h5("PCA Options"),
                       checkboxInput("exclude_blank_pca", "Exclude samples named 'Blank'", value = TRUE),
                       checkboxInput("pca_show_iqc", "Include iQC samples (projected)", value = TRUE),
                       checkboxInput("pca_show_istd", "Include ISTD-only samples (projected)", value = FALSE),
                       helpText(tags$small(
                         "iQC (pooled biological QC) and ISTD-only injections (no biological matrix)",
                         " are ", tags$strong("not the same thing"), " and are never treated as equivalent.",
                         " ISTD samples are excluded by default \u2014 their lipid profile is so different",
                         " from a real sample that including them was previously showing up as a",
                         " spurious dominant point."
                       )),
                       conditionalPanel(
                         condition = "input.use_protein_norm == true && (input.pca_show_iqc == true || input.pca_show_istd == true)",
                         helpText(tags$small(
                           "Both are fitted as ", tags$strong("supplementary points"), ":",
                           " the PCA axes are calculated from biological samples only, so neither can",
                           " distort the ordination. Since they usually aren't in the protein CSV,",
                           " they're projected after being scaled by the ", tags$em("median"),
                           " protein content of the biological samples \u2014 an assumed, not measured,",
                           " value \u2014 purely so they land in a comparable position for visual QC."
                         ))
                       ),
                       radioButtons(
                         "pca_measure",
                         "Data to use:",
                         choices = c(
                           "Quantified (IS-normalised)" = "norm",
                           "Background-subtracted"      = "value_bs"
                         ),
                         selected = "norm"
                       ),
                       radioButtons(
                         "pca_units",
                         "Units:",
                         choices = c(
                           "Absolute"               = "absolute",
                           "Percent of class total" = "percent"
                         ),
                         selected = "absolute"
                       ),
                       radioButtons(
                         "pca_scaling",
                         "Feature scaling:",
                         choices = c(
                           "Unit variance (UV) — equal weight to all features"     = "uv",
                           "Pareto — square-root of SD, reduces dominance of large signals" = "pareto",
                           "Log transform then UV — compresses dynamic range"       = "log_uv",
                           "None (mean-centred only)"                               = "none"
                         ),
                         selected = "uv"
                       ),
                       helpText(tags$small(
                         tags$strong("UV"), " is standard for lipidomics. ",
                         tags$strong("Pareto"), " is useful when a few abundant classes dominate the scores. ",
                         tags$strong("Log + UV"), " is recommended when signals span several orders of magnitude."
                       )),
                       tags$hr(),
                       h5("Token-based grouping"),
                       helpText(tags$small(
                         "Select which token position(s) define the group for PCA colouring.",
                         "Use the Group Preview tab to identify positions."
                       )),
                       uiOutput("pca_group_token_ui"),
                       textInput("pca_group_sep", "Join tokens with", value = "_"),
                       tags$hr(),
                       h5("Group selection"),
                       helpText(tags$small("Select which groups to include in the PCA plot.")),
                       uiOutput("pca_group_select_ui"),
                       tags$hr(),
                       checkboxInput("pca_show_labels", "Show sample ID labels", value = FALSE),
                       tags$hr(),
                       h5("Export plot"),
                       numericInput("pca_export_width",    "Width (px)",      1200, 400, 4000, 50),
                       numericInput("pca_export_height",   "Height (px)",      900, 300, 4000, 50),
                       numericInput("pca_export_dpi",      "DPI",              300,  72,  600, 12),
                       numericInput("pca_export_scale",    "Scale fraction",  1.00, 0.25, 2.00, 0.05),
                       numericInput("pca_export_fontsize", "Base font size",    14,    6,   24,  1),
                       fluidRow(
                         column(6, downloadButton("download_pca_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_pca_svg", "SVG"))
                       )
                     )
                   ),
                   column(
                     width = 9,
                     plotOutput("pcaPlot", height = "600px"),
                     br(),
                     h5("PCA Loadings Table:"),
                     DTOutput("pcaLoadingsTable")
                   )
                 )
        ),
        tabPanel("Statistics",
                 
                 h4("Step 8: Statistical tests by metabolite / class"),
                 
                 # ── Per-class summary panel ──────────────────────────────────
                 tags$details(
                   style = "margin-bottom: 12px;",
                   tags$summary(
                     style = paste0(
                       "cursor:pointer; font-weight:600; font-size:0.95em;",
                       " padding:6px 10px; background:#f0f4f8;",
                       " border:1px solid #d0d7de; border-radius:4px;"
                     ),
                     "📋  Per-class detection & significance summary"
                   ),
                   div(
                     style = "padding: 8px 4px 4px 4px;",
                     fluidRow(
                       column(6,
                              helpText(tags$small(
                                "Detected: species present in at least one sample after normalisation. ",
                                "Tested: species included in the most recent statistics run. ",
                                "Significant: passing the current α threshold."
                              ))
                       ),
                       column(6, style = "text-align:right;",
                              downloadButton("download_class_summary",
                                             "Download summary (CSV)",
                                             style = "font-size:0.85em; padding:3px 10px;")
                       )
                     ),
                     DT::DTOutput("classStatsSummary", width = "100%")
                   )
                 ),
                 
                 fluidRow(
                   column(
                     width = 8,
                     h5("Bar plot of significant features"),
                     div(
                       style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:10px;",
                       actionButton("prev_plot", "Previous"),
                       actionButton("next_plot", "Next")
                     ),
                     textOutput("plotIndexInfo"),
                     uiOutput("statsPlotUI"),
                     br(),
                     h5(textOutput("stats_header")),
                     DTOutput("statsResults"),
                     br(),
                     downloadButton("download_stats", "Download results (CSV)"),
                     br(), br(),
                     conditionalPanel(
                       condition = "input.posthoc.length > 0 && input.stats_scope == 'single'",
                       h5("Post-hoc comparisons"),
                       DTOutput("posthocResults")
                     )
                   ),
                   column(
                     width = 4,
                     wellPanel(
                       actionButton("run_stats",
                                    "Run statistics",
                                    icon = icon("play"),
                                    style = "background-color:#28a745; color:#fff; border-color:#28a745; width:100%;"),
                       br(),
                       helpText(
                         tags$ol(
                           style = "padding-left:1.2em; margin:4px 0;",
                           tags$li("Select a lipid class and grouping method below."),
                           tags$li("Choose a test (or leave Auto)."),
                           tags$li("Click Run statistics. Results appear in the table and bar plots.")
                         )
                       ),
                       hr(),
                       checkboxInput("exclude_blank_stats", "Exclude samples named 'Blank'", value = TRUE),
                       radioButtons(
                         "stats_value_type",
                         "Data to test:",
                         choices = c("Quantified (IS-normalised)" = "norm",
                                     "Background-subtracted"      = "value_bs"),
                         selected = "norm"
                       ),
                       radioButtons(
                         "stats_units",
                         "Units:",
                         choices = c("Absolute" = "absolute", "Percent of class total" = "percent"),
                         selected = "absolute"
                       ),
                       hr(),
                       selectInput("stats_class", "Lipid Class", choices = "Loading..."),
                       selectInput("stats_scope",
                                   "What to test:",
                                   choices = c("All metabolites in class" = "all",
                                               "Single metabolite"       = "single",
                                               "Class total (sum)"       = "total"),
                                   selected = "all"
                       ),
                       conditionalPanel(
                         condition = "input.stats_scope == 'single'",
                         selectInput("stats_met_name", "Metabolite", choices = "Select a class first")
                       ),
                       hr(),
                       # ── Groups inherited from Group Preview ───────────────────
                       div(
                         style = "background:#f0f4fa; border-radius:6px; padding:10px 12px; margin-bottom:10px;",
                         tags$strong(tags$small("Group selection")),
                         helpText(tags$small(
                           "Groups are defined in the ", tags$strong("Group Preview"), " tab.",
                           "Deselect groups to exclude them from all tests.",
                           "If a grouping CSV is uploaded and enabled, it overrides token-based grouping."
                         )),
                         uiOutput("stats_group_select_ui")
                       ),
                       hr(),
                       h5("Test & options"),
                       radioButtons("stats_test",
                                    "Test type:",
                                    choices = c(
                                      "Auto (2 groups → t-test; ≥3 → one-way ANOVA)" = "auto",
                                      "Unpaired t-test (2 groups)"                   = "ttest",
                                      "One-way ANOVA (≥3 groups)"                    = "oneway",
                                      "Two-way ANOVA (uses Factor A & B from Group Preview)" = "twoway"
                                    ),
                                    selected = "auto"
                       ),
                       conditionalPanel(
                         condition = "input.stats_test == 'ttest'",
                         checkboxInput("ttest_equal_var", "Assume equal variances (Student)", value = FALSE)
                       ),
                       conditionalPanel(
                         condition = "input.stats_test == 'twoway'",
                         helpText(tags$small(
                           "Two-way ANOVA uses Factor A and Factor B defined in the Group Preview tab."
                         )),
                         selectInput("twoway_effect_filter",
                                     "Report/Filter effect by:",
                                     choices = c("Any significant effect" = "any",
                                                 "Main effect: Factor A"  = "A",
                                                 "Main effect: Factor B"  = "B",
                                                 "Interaction A:B"        = "AxB"),
                                     selected = "any")
                       ),
                       selectInput("padj_method",
                                   "Multiple testing correction:",
                                   choices = c("None" = "none", "FDR (BH)" = "BH", "Bonferroni" = "bonferroni"),
                                   selected = "BH"
                       ),
                       checkboxInput("use_adj_threshold", "Apply threshold to adjusted p-values (if corrected)", value = TRUE),
                       numericInput("alpha", "Significance threshold (α):", value = 0.05, min = 0.0001, max = 0.2, step = 0.005),
                       checkboxInput("show_all_rows", "Show all rows (not only significant)", value = FALSE),
                       checkboxGroupInput("posthoc",
                                          "Post-hoc (for ANOVA / two-way):",
                                          choices = c("Tukey HSD" = "tukey", "Pairwise t-tests (Holm)" = "pairwise"),
                                          selected = NULL
                       ),
                       hr(),
                       h5("Plot options"),
                       selectInput("stats_error_type", "Error bars for plot:",
                                   choices = c("SEM", "SD", "95% CI"), selected = "SEM"
                       ),
                       sliderInput("stats_plot_topn", "Top N features to plot:", min = 5, max = 50, value = 20, step = 1),
                       selectInput("stats_plot_sort",
                                   "Sort features by:",
                                   choices = c("Adjusted p-value (asc)" = "p_adj",
                                               "Raw p-value (asc)"      = "p",
                                               "Absolute log2FC (desc, t-test only)" = "log2fc"),
                                   selected = "p_adj"
                       ),
                       checkboxInput("stats_plot_points", "Overlay individual points", value = TRUE),
                       hr(),
                       h5("Export significant plots"),
                       helpText(tags$small(
                         "Downloads one bar plot per significant feature as a multi-page PDF,",
                         "matching exactly what is shown on the statistics tab including",
                         "significance brackets, error bars, and individual points."
                       )),
                       numericInput("stats_export_width",    "Width (px)",      1200, 400, 4000, 50),
                       numericInput("stats_export_height",   "Height (px)",      700, 300, 4000, 50),
                       numericInput("stats_export_dpi",      "DPI",              300,  72,  600, 12),
                       numericInput("stats_export_scale",    "Scale fraction",  1.00, 0.25, 2.00, 0.05),
                       numericInput("stats_export_fontsize", "Base font size",    14,    6,   24,  1),
                       fluidRow(
                         column(6, downloadButton("download_stats_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_stats_svg", "SVG"))
                       ),
                       hr(),
                       h5("Batch PDF export"),
                       helpText(tags$small(
                         "Exports all significant plots as a multi-page PDF.",
                         "Uses width/height/DPI settings above."
                       )),
                       radioButtons("stats_pdf_scale", "Layout:",
                                    choices = c("1 per page"        = "single",
                                                "2 per page"        = "two",
                                                "4 per page (2×2)"  = "four"),
                                    selected = "single",
                                    inline   = TRUE),
                       downloadButton("download_stats_pdf",
                                      label = "Download PDF",
                                      icon  = icon("file-pdf"),
                                      style = "width:100%; margin-top:4px;")
                     )
                   )
                 )
        ),
        tabPanel("Volcano Plot",
                 h4("Step 8a: Volcano Plot"),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       numericInput("volcano_alpha",
                                    "Significance threshold:",
                                    value = 0.05, min = 0.0001, max = 0.5, step = 0.005),
                       numericInput("volcano_log2fc",
                                    "log2FC threshold (+/-):",
                                    value = 1, min = 0, max = 10, step = 0.25),
                       checkboxInput("volcano_use_padj",
                                     "Use adjusted p-value", value = TRUE),
                       uiOutput("volcano_comparison_ui"),
                       checkboxInput("volcano_label_sig",
                                     "Label top significant points", value = TRUE),
                       numericInput("volcano_topn_labels",
                                    "Max labels:", value = 20, min = 0, max = 100),
                       hr(),
                       checkboxInput("volcano_cap_axes",
                                     "Cap axes to reduce outlier distortion", value = FALSE),
                       conditionalPanel(
                         condition = "input.volcano_cap_axes == true",
                         numericInput("volcano_cap_y",
                                      "-log10(p) cap:",
                                      value = 10, min = 1, max = 500, step = 1),
                         numericInput("volcano_cap_x",
                                      "log2FC cap (+/-):",
                                      value = 5, min = 0.5, max = 50, step = 0.5),
                         helpText(tags$small("Features beyond the cap are plotted at the cap value as triangles \u25b2."))
                       ),
                       helpText("Run Statistics first to populate this plot."),
                       hr(),
                       h5("Export plot"),
                       numericInput("volc_export_width",  "Width (px)",  1200, 400, 4000, 50),
                       numericInput("volc_export_height", "Height (px)",  700, 300, 4000, 50),
                       numericInput("volc_export_dpi",    "DPI",          300,  72,  600, 12),
                       numericInput("volc_export_scale",  "Scale fraction", 1.00, 0.25, 2.00, 0.05),
                       fluidRow(
                         column(6, downloadButton("download_volcano_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_volcano_svg", "SVG"))
                       )
                     )
                   ),
                   column(
                     width = 9,
                     plotlyOutput("volcanoPlot", height = "580px")
                   )
                 )
        ),
        tabPanel("Class Bar Plots",
                 h4("Step 9: Class bar plots of significant features"),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       h5("Data"),
                       uiOutput("cbp_class_ui"),
                       selectInput("cbp_value_type", "Values to plot",
                                   choices = c("Quantified (IS-normalised)" = "norm",
                                               "Background-subtracted"      = "value_bs"),
                                   selected = "norm"),
                       helpText(tags$small(
                         "Units are detected automatically from your ISTD CSV and,",
                         " if enabled, your protein normalisation settings."
                       )),
                       hr(),
                       h5("Filtering"),
                       checkboxInput("cbp_sig_only", "Show significant species only", value = FALSE),
                       sliderInput("cbp_abundance_range",
                                   "Abundance range filter (% of class max)",
                                   min = 0, max = 100, value = c(0, 100), step = 1),
                       helpText(tags$small(
                         "Filter species by their mean abundance relative to the most abundant",
                         "species in the class. Use this to focus on high-, mid-, or",
                         "low-abundance species."
                       )),
                       hr(),
                       h5("Appearance"),
                       selectInput("cbp_palette", "Colour palette",
                                   choices = c(
                                     "── Qualitative ──"      = "",
                                     "Okabe-Ito (CB-safe)"   = "OkabeIto",
                                     "Dark2"                  = "Dark2",
                                     "Set1"                   = "Set1",
                                     "Set2"                   = "Set2",
                                     "Paired"                 = "Paired",
                                     "Accent"                 = "Accent",
                                     "Tableau 10"             = "Tableau10",
                                     "── Sequential ──"       = "",
                                     "Viridis"                = "viridis",
                                     "Plasma"                 = "plasma",
                                     "Inferno"                = "inferno",
                                     "Blues"                  = "Blues",
                                     "Greens"                 = "Greens",
                                     "Purples"                = "Purples",
                                     "── Diverging ──"        = "",
                                     "RdBu"                   = "RdBu",
                                     "PuOr"                   = "PuOr",
                                     "BrBG"                   = "BrBG"
                                   ),
                                   selected = "OkabeIto"),
                       numericInput("cbp_bar_alpha", "Bar opacity", value = 0.85, min = 0.2, max = 1.0, step = 0.05),
                       numericInput("cbp_point_size", "Point size", value = 2.2, min = 0.5, max = 6.0, step = 0.2),
                       selectInput("cbp_point_fill", "Point fill",
                                   choices = c("Match palette" = "palette",
                                               "Black"         = "black",
                                               "White"         = "white"),
                                   selected = "palette"),
                       selectInput("cbp_point_outline", "Point outline",
                                   choices = c("White"  = "white",
                                               "Black"  = "black",
                                               "None"   = "none"),
                                   selected = "white"),
                       selectInput("cbp_error_type", "Error bars",
                                   choices = c("SEM" = "SEM", "SD" = "SD", "95% CI" = "CI95"),
                                   selected = "SEM"),
                       checkboxInput("cbp_show_points", "Overlay individual points", value = TRUE),
                       radioButtons("cbp_sort_order", "Species order",
                                    choices  = c("Alphabetical" = "alpha",
                                                 "Abundance (descending)" = "abundance"),
                                    selected = "alpha",
                                    inline   = TRUE),
                       hr(),
                       h5("Export plot"),
                       numericInput("cbp_export_width",    "Width (px)",     1600, 400, 6000, 50),
                       numericInput("cbp_export_height",   "Height (px)",    1200, 300, 6000, 50),
                       numericInput("cbp_export_dpi",      "DPI",             300,  72,  600, 12),
                       numericInput("cbp_export_scale",    "Scale fraction", 1.00, 0.25, 2.00, 0.05),
                       numericInput("cbp_export_fontsize", "Base font size",   11,    6,   24,  1),
                       fluidRow(
                         column(6, downloadButton("download_cbp_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_cbp_svg", "SVG"))
                       )
                     )
                   ),
                   column(
                     width = 9,
                     uiOutput("cbpInfoUI"),
                     plotOutput("classBarPlot",
                                width  = "100%",
                                height = "auto")
                   )
                 )
        ),
        tabPanel("Heatmap - significant",
                 h4("Step 9b: Heatmap of statistically significant features"),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       checkboxInput("hm_use_group_means", "Use group means (collapse replicates)", value = TRUE),
                       checkboxInput("hm_use_axb", "If two-way factors present, use A×B combinations", value = TRUE),
                       selectInput("hm_scale", "Row scaling",
                                   choices = c("Z-score (mean=0, sd=1)" = "z",
                                               "Mean center only"       = "center",
                                               "None"                   = "none"),
                                   selected = "z"),
                       selectInput("hm_dist", "Distance",
                                   choices = c("Euclidean" = "euclidean",
                                               "1 - Pearson correlation" = "correlation"),
                                   selected = "correlation"),
                       selectInput("hm_linkage", "Linkage",
                                   choices = c("Complete" = "complete",
                                               "Average (UPGMA)" = "average",
                                               "Ward.D2" = "ward.D2",
                                               "Single" = "single"),
                                   selected = "complete"),
                       checkboxInput("hm_cluster_rows", "Cluster rows (features)", TRUE),
                       checkboxInput("hm_cluster_cols", "Cluster columns", TRUE),
                       checkboxInput("hm_show_rownames", "Show row labels (features)", TRUE),
                       checkboxInput("hm_show_colnames", "Show column labels (samples)", TRUE),
                       selectInput("hm_palette", "Colour palette",
                                   choices = c(
                                     "Viridis"              = "viridis",
                                     "Magma"                = "magma",
                                     "Plasma"               = "plasma",
                                     "Inferno"              = "inferno",
                                     "Cividis"              = "cividis",
                                     "Rocket"               = "rocket",
                                     "Mako"                 = "mako",
                                     "Turbo"                = "turbo",
                                     "Blue–White–Red"       = "bwr",
                                     "Green–White–Purple"   = "gwp"
                                   ),
                                   selected = "viridis"),
                       numericInput("hm_topn", "Top N features (by adjusted p)", value = 50, min = 2, max = 300, step = 1),
                       br(),
                       downloadButton("download_sig_hm_matrix", "Download matrix (CSV)"),
                       hr(),
                       h5("Export heatmap"),
                       numericInput("hm_export_width",    "Width (px)",      1100, 400, 4000, 50),
                       numericInput("hm_export_height",   "Height (px)",      900, 300, 4000, 50),
                       numericInput("hm_export_dpi",      "DPI",              300,  72,  600, 12),
                       numericInput("hm_export_scale",    "Scale fraction",  1.00, 0.25, 2.00, 0.05),
                       numericInput("hm_export_fontsize", "Base font size",    10,    6,   24,  1),
                       fluidRow(
                         column(6, downloadButton("download_hm_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_hm_svg", "SVG"))
                       )
                       
                     )
                   ),
                   column(
                     width = 9,
                     textOutput("sigHeatmapInfo"),
                     plotOutput("sigHeatmap", width = "100%",
                                height = "auto")
                   )
                 )
        ),
        tabPanel("Enrichment",
                 h4("Step 10: Lipid set enrichment"),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       h5("Fatty acid interpretation"),
                       radioButtons(
                         "fa_parsing_mode", NULL,
                         choices = c("Sum composition (e.g., PC 34:1)" = "sum",
                                     "Resolved FA list (e.g., PC 16:0_18:1)" = "resolved"),
                         selected = "sum"
                       ),
                       checkboxInput(
                         "fa_resolved_fallback_sum",
                         "If FA list is missing for a lipid, fallback to sum composition",
                         value = TRUE
                       ),
                       checkboxInput(
                         "fa_stoich_weighting",
                         "Resolved mode: weight FA contributions by stoichiometry",
                         value = TRUE
                       ),
                       hr(),
                       h5("Enrichment scope"),
                       radioButtons(
                         "lsea_scope", NULL, inline = FALSE,
                         choices = c("Overall (use current Statistics results)" = "overall",
                                     "Per-group (one-vs-rest)"                 = "group"),
                         selected = "group"
                       ),
                       # Group controls (shown only in per-group mode)
                       conditionalPanel(
                         condition = "input.lsea_scope == 'group'",
                         selectizeInput(
                           "lsea_group_select", "Group(s) to analyze",
                           choices = c(), multiple = TRUE, selected = NULL,
                           options = list(placeholder = 'Choose one or more groups…')
                         ),
                         checkboxInput("lsea_facet_groups", "Facet plot by group (if multiple selected)", TRUE),
                         selectInput(
                           "lsea_univ_scope", "Feature universe",
                           choices = c("All classes" = "all", "Current class (from Statistics)" = "current"),
                           selected = "all"
                         ),
                         # Direction for ORA hits in group-vs-rest
                         checkboxInput("ora_up_only", "ORA (group): hits must be up in the selected group (log2FC > 0)", TRUE)
                       ),
                       hr(),
                       h5("Enrichment method"),
                       radioButtons(
                         "lsea_method", NULL,
                         choices = c("ORA (Fisher over-representation)" = "ora",
                                     "FGSEA (rank-based)"               = "fgsea"),
                         selected = "ora"
                       ),
                       # FGSEA hint
                       conditionalPanel(
                         condition = "input.lsea_method == 'fgsea'",
                         helpText("Tip: FGSEA works best when 'Show all rows' is ON in Statistics, so the rank covers all tested lipids.")
                       ),
                       # ORA UX
                       conditionalPanel(
                         condition = "input.lsea_method == 'ora'",
                         checkboxInput("ora_only_overlap", "ORA: show sets with ≥1 overlapping hit", TRUE),
                         radioButtons(
                           "ora_x_axis", "ORA x-axis", inline = TRUE,
                           choices = c("Adjusted p" = "padj", "Raw p" = "p"),
                           selected = "padj"
                         )
                       ),
                       checkboxInput("lsea_per_class", "Include per-class sets (in addition to global)", TRUE),
                       numericInput("lsea_min_set", "Minimum set size", value = 3, min = 2, max = 200, step = 1),
                       numericInput("lsea_topn", "Top N rows to display", value = 30, min = 5, max = 300, step = 1),
                       br(),
                       actionButton("show_refs_lsea", "Show references", icon = icon("book-open")),
                       br(), br(),
                       downloadButton("download_lsea_table", "Download results (CSV)"),
                       hr(),
                       h5("Export enrichment plot"),
                       numericInput("lsea_export_width",    "Width (px)",      1200, 400, 4000, 50),
                       numericInput("lsea_export_height",   "Height (px)",      700, 300, 4000, 50),
                       numericInput("lsea_export_dpi",      "DPI",              300,  72,  600, 12),
                       numericInput("lsea_export_scale",    "Scale fraction",  1.00, 0.25, 2.00, 0.05),
                       numericInput("lsea_export_fontsize", "Base font size",    11,    6,   24,  1),
                       fluidRow(
                         column(6, downloadButton("download_lsea_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_lsea_svg", "SVG"))
                       )
                       
                     )
                   ),
                   column(
                     width = 9,
                     textOutput("lseaContextInfo"),
                     plotOutput("lseaDotPlot", height = "560px"),
                     br(),
                     DTOutput("lseaTable")
                   )
                 )
        ),
        tabPanel("Correlation Network",
                 h4("Step 11: Correlation Network"),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       h5("Group & data"),
                       uiOutput("net_group_token_ui"),
                       selectInput("net_group",
                                   "Compute correlations within group:",
                                   choices = c()),
                       radioButtons("net_value_type", "Data:",
                                    choices = c("Quantified (IS-normalised)" = "norm",
                                                "Background-subtracted"      = "value_bs"),
                                    selected = "norm"),
                       checkboxInput("net_exclude_blank",
                                     "Exclude blanks", value = TRUE),
                       hr(),
                       h5("Node filter"),
                       helpText(tags$small(
                         "Nodes are lipid species from Statistics results.",
                         "Node colour = direction in the selected group (coral = up, teal = down).",
                         "Edge width = correlation strength."
                       )),
                       numericInput("net_alpha",
                                    "Max adj p for node inclusion:",
                                    value = 0.05, min = 0.001, max = 0.5, step = 0.005),
                       numericInput("net_max_nodes",
                                    "Max nodes (top by adj p):",
                                    value = 50, min = 5, max = 200, step = 5),
                       hr(),
                       h5("Edge filter"),
                       numericInput("net_r_thresh",
                                    "Min |r| for edge:",
                                    value = 0.7, min = 0.2, max = 0.99, step = 0.05),
                       hr(),
                       actionButton("net_run", "Build network",
                                    icon = icon("project-diagram"),
                                    style = "background-color:#28a745; color:#fff; border-color:#28a745; width:100%;"),
                       br(), br(),
                       uiOutput("net_download_ui")
                     )
                   ),
                   column(
                     width = 9,
                     helpText(tags$small(
                       "Nodes = significant lipid species. Coral/orange = elevated in selected group; teal/green = lower.",
                       "Edge colour: green = positive correlation, red = negative. Edge width = |r|.",
                       "Drag nodes to rearrange. Hover for details. Click a node to highlight its edges."
                     )),
                     visNetwork::visNetworkOutput("net_plot", height = "600px"),
                     br(),
                     DTOutput("net_table")
                   )
                 )
        ),
        tabPanel("Synthesis Pathways",
                 h4("Step 12: Lipid class synthesis pathway scores & statistics"),
                 helpText(tags$small(
                   tags$strong("Note:"), " these scores are ratios/fractions built from established ",
                   tags$em("mammalian"), " lipid synthesis pathways (e.g. the Kennedy pathway, ether",
                   " lipid biosynthesis, PEMT-mediated PE\u2192PC conversion). Some steps involve enzymes",
                   " with overlapping or tissue-dependent substrate preferences, and a given score may",
                   " reflect more than one biosynthetic route contributing to the same lipid pool.",
                   " Treat them as pathway-activity indicators rather than direct measurements of flux",
                   " through a single enzymatic step, and interpret changes with reference to the",
                   " underlying biology of your tissue/model system."
                 )),
                 fluidRow(
                   column(
                     width = 3,
                     wellPanel(
                       checkboxInput("exclude_blank_path", "Exclude samples named 'Blank'", TRUE),
                       radioButtons(
                         "path_value_type", "Data to use:",
                         choices = c("Quantified (IS-normalised)" = "norm",
                                     "Background-subtracted" = "value_bs"),
                         selected = "norm"
                       ),
                       radioButtons(
                         "path_units", "Units:",
                         choices = c("Absolute" = "absolute",
                                     "Percent of total lipids" = "percent_total"),
                         selected = "absolute"
                       ),
                       uiOutput("path_group_token_ui"),
                       selectInput("path_hm_palette", "Colour palette",
                                   choices = c(
                                     "Viridis"              = "viridis",
                                     "Magma"                = "magma",
                                     "Plasma"               = "plasma",
                                     "Inferno"              = "inferno",
                                     "Cividis"              = "cividis",
                                     "Rocket"               = "rocket",
                                     "Mako"                 = "mako",
                                     "Turbo"                = "turbo",
                                     "Blue–White–Red"       = "bwr",
                                     "Green–White–Purple"   = "gwp"
                                   ),
                                   selected = "viridis"),
                       checkboxInput("path_use_group_means", "Use group means (collapse replicates)", TRUE),
                       selectizeInput("path_scores_select", "Scores to display", multiple = TRUE, choices = c(),
                                      options = list(placeholder = "Select scores (defaults to all)")),
                       hr(),
                       h5("Score statistics"),
                       helpText(tags$small(
                         "Runs a t-test (2 groups) or one-way ANOVA (≥3 groups) on each",
                         "score's sample values. BH correction is applied across all scores."
                       )),
                       selectInput("path_padj_method", "Multiple testing correction:",
                                   choices = c("FDR (BH)" = "BH", "Bonferroni" = "bonferroni", "None" = "none"),
                                   selected = "BH"),
                       numericInput("path_alpha", "Significance threshold (α):", value = 0.05, min = 0.0001, max = 0.2, step = 0.005),
                       checkboxInput("path_stats_equal_var", "Assume equal variance (t-test)", value = FALSE),
                       br(),
                       downloadButton("download_path_scores", "Download scores (CSV)"),
                       hr(),
                       h5("Export plots"),
                       numericInput("path_export_width",    "Width (px)",      1400, 400, 5000, 50),
                       numericInput("path_export_height",   "Height (px)",      900, 300, 5000, 50),
                       numericInput("path_export_dpi",      "DPI",              300,  72,  600, 12),
                       numericInput("path_export_scale",    "Scale fraction",  1.00, 0.25, 4.00, 0.05),
                       numericInput("path_export_fontsize", "Base font size",    11,    6,   24,  1),
                       selectInput("path_export_plot_sel", "Plot to export",
                                   choices = c("Scores heatmap"   = "heatmap",
                                               "Score statistics" = "scorestats"),
                                   selected = "heatmap"),
                       fluidRow(
                         column(6, downloadButton("download_path_png", "PNG", class = "btn-primary")),
                         column(6, downloadButton("download_path_svg", "SVG"))
                       )
                     )
                   ),
                   column(
                     width = 9,
                     tabsetPanel(
                       tabPanel("Scores heatmap",
                                plotOutput("pathHeatmap", height = "560px"),
                                br(), DT::DTOutput("pathScoresTable")),
                       tabPanel("Score statistics",
                                plotOutput("pathStatsPlot", height = "520px"),
                                br(), DT::DTOutput("pathStatsTable"))
                     )
                   )
                 )
        )
        
        
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ── Idle-session timeout ──────────────────────────────────────────────────
  # Paired with idle_timeout_js in the UI. Warns at 15 min idle, closes the
  # session at 20 min idle. Session close frees the R process's held memory
  # for that session (uploaded data, computed results) once a tab has been
  # abandoned; it does not affect active users.
  observeEvent(input$idle_warning, {
    showModal(modalDialog(
      title = "Still there?",
      "This session has been idle for 15 minutes. To free up server resources, ",
      "it will close automatically after 20 minutes of inactivity. ",
      "Click below or interact with the app to stay connected.",
      footer = modalButton("Stay connected"),
      easyClose = TRUE
    ))
  }, ignoreInit = TRUE)
  
  observeEvent(input$idle_dismiss_warning, {
    removeModal()
  }, ignoreInit = TRUE)
  
  observeEvent(input$idle_timeout, {
    showModal(modalDialog(
      title = "Session closed",
      "This session was closed after 20 minutes of inactivity. Please refresh the page to start a new session.",
      footer = NULL,
      easyClose = FALSE
    ))
    session$close()
  }, ignoreInit = TRUE)
  
  
  
  # ---------- Utility helpers ----------
  `%||%` <- function(a, b) {
    if (is.null(a) || length(a) == 0) return(b)
    if (length(a) == 1 && is.na(a)) return(b)
    a
  }
  
  # ── Settings server logic ─────────────────────────────────────────────────────
  
  # Default class lists (mirror the UI defaults)
  .default_pref_neg <- c("PC","PE","PG","PI","PS","PA","CL","PC-O","PE-O")
  .default_pref_pos <- c("TG","DG","MG","CE","Cer","HexCer","SM",
                         "LPC","LPE","LPG","LPI","LPS","LPA")
  .default_cross_mode <- c("PC","PE","PG","PI","PS","PA","CL")
  
  # Reset button
  observeEvent(input$settings_reset, {
    updateTextAreaInput(session, "pref_neg_classes",
                        value = paste(.default_pref_neg, collapse = "\n"))
    updateTextAreaInput(session, "pref_pos_classes",
                        value = paste(.default_pref_pos, collapse = "\n"))
    updateTextAreaInput(session, "cross_mode_classes",
                        value = paste(.default_cross_mode, collapse = "\n"))
    updateCheckboxInput(session, "cross_mode_enable", value = TRUE)
    updateNumericInput(session, "cross_mode_rt_tol", value = 0.1)
  })
  
  # Reactive helpers to parse current preference lists from textarea
  pref_neg_classes <- reactive({
    x <- input$pref_neg_classes %||% ""
    cls <- stringr::str_trim(unlist(strsplit(x, "\n")))
    cls[nzchar(cls)]
  })
  pref_pos_classes <- reactive({
    x <- input$pref_pos_classes %||% ""
    cls <- stringr::str_trim(unlist(strsplit(x, "\n")))
    cls[nzchar(cls)]
  })
  cross_mode_classes <- reactive({
    x <- input$cross_mode_classes %||% ""
    cls <- stringr::str_trim(unlist(strsplit(x, "\n")))
    cls[nzchar(cls)]
  })
  
  # Audit map for the cross-mode reconciliation (recomputed independently of the
  # main pipeline so the table stays available even after downstream mutates).
  cross_mode_audit <- reactive({
    if (!isTRUE(input$cross_mode_enable)) return(NULL)
    cm_classes <- cross_mode_classes()
    if (length(cm_classes) == 0) return(NULL)
    req(bg_norm_long())
    cm_rt_tol <- suppressWarnings(as.numeric(input$cross_mode_rt_tol %||% 0.1))
    if (length(cm_rt_tol) != 1 || is.na(cm_rt_tol) || cm_rt_tol < 0) cm_rt_tol <- 0.1
    df <- resolve_cross_mode(bg_norm_long(), cm_classes, rt_tol = cm_rt_tol)
    attr(df, "cross_mode_map")
  })
  
  output$crossModeAuditTable <- DT::renderDT({
    m <- cross_mode_audit()
    validate(need(!is.null(m) && nrow(m) > 0,
                  "Cross-mode identity transfer is off, or no designated-class features were found."))
    m %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4))) %>%
      dplyr::rename(
        Class = class, Shorthand = shorthand,
        `Positive feature` = pos_name, `Pos RT` = pos_rt, `Pos adduct` = pos_adduct,
        `Negative identity` = neg_identity, `Neg RT` = neg_rt, `dRT` = delta_rt,
        `Final identity` = final_identity, `Quant mode` = quant_mode, Status = status
      )
  }, options = list(pageLength = 15, scrollX = TRUE, order = list()),
     rownames = FALSE)
  
  # ── Token preview server logic ─────────────────────────────────────────────────
  
  # Reactive: build token preview table from current sample names
  lr_token_preview <- eventReactive(input$lr_preview_groups, {
    # Read directly from uploaded file(s) — does NOT require Load button to be clicked
    req(input$msdial_txts)
    validate(need(nrow(input$msdial_txts) >= 1, "Upload at least one MS-DIAL .txt file first."))
    
    # Read just the header row (row 5 in MS-DIAL = first data row after skip = 4)
    # We only need column names so read 0 rows of data
    raw_cols <- tryCatch({
      df <- readr::read_tsv(input$msdial_txts$datapath[1],
                            skip          = 4,
                            n_max         = 1,
                            col_types     = readr::cols(.default = "c"),
                            show_col_types = FALSE)
      colnames(df)
    }, error = function(e) character(0))
    
    validate(need(length(raw_cols) > 0, "Could not read column names from the uploaded file."))
    
    # Known MS-DIAL metadata columns to exclude
    meta_cols <- c("Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type",
                   "Average Rt.min.", "AverageMz", "MetaboliteName", "AdductIonName",
                   "Spectrum reference file name", "MS1 isotopic spectrum",
                   "MS/MS spectrum", "Comment", "Fill %", "S/N average",
                   "Dot product", "Reverse dot product", "Fragment presence %")
    
    # Keep only columns that look like sample columns (contain _pos or _neg)
    sample_raw <- raw_cols[grepl("(?i)(_pos|_neg)", raw_cols, perl = TRUE)]
    
    # If no pos/neg suffix columns found fall back to removing known meta cols
    if (length(sample_raw) == 0)
      sample_raw <- setdiff(raw_cols, meta_cols)
    
    # Strip ion mode suffix and replicate suffixes to get clean base sample names
    samples <- sample_raw %>%
      stringr::str_replace("(?i)(_pos|_neg)(\\.[0-9]+)?$", "") %>%
      stringr::str_replace("(?i)(_pos|_neg)_.*$", "") %>%
      unique()
    samples <- samples[nzchar(samples)]
    
    validate(need(length(samples) > 0, "No sample columns detected. Check that column names contain '_pos' or '_neg'."))
    
    delim <- if (!is.null(input$lr_delimiter) && nzchar(input$lr_delimiter))
      input$lr_delimiter else "_"
    
    tbl <- tibble::tibble(sample = samples)
    max_tok <- min(
      max(lengths(strsplit(samples, delim, fixed = TRUE)), na.rm = TRUE),
      10L
    )
    for (i in seq_len(max_tok))
      tbl[[paste0("token_", i)]] <- extract_token(samples, delimiter = delim, token_index = i)
    
    tbl
  }, ignoreInit = TRUE)
  
  output$lr_group_preview_tbl <- DT::renderDT({
    req(lr_token_preview())
    DT::datatable(
      lr_token_preview(),
      options  = list(pageLength = 30, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  # Summary: show most common value per token position
  output$lr_token_summary_ui <- renderUI({
    tbl <- tryCatch(lr_token_preview(), error = function(e) NULL)
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    tok_cols <- grep("^token_", names(tbl), value = TRUE)
    if (length(tok_cols) == 0) return(NULL)
    items <- lapply(tok_cols, function(col) {
      vals    <- na.omit(tbl[[col]])
      top3    <- names(sort(table(vals), decreasing = TRUE))[seq_len(min(3, length(unique(vals))))]
      n_uniq  <- length(unique(vals))
      tags$li(
        tags$strong(col), " — ",
        n_uniq, " unique value(s). ",
        "Examples: ", paste(top3, collapse = ", ")
      )
    })
    tagList(
      h5("Token summary"),
      helpText("Each row below describes one token position across all samples:"),
      tags$ul(items)
    )
  })
  
  # ── Available groups reactiveVals (populated when token selection changes) ────
  pca_available_groups  <- reactiveVal(character(0))
  stats_available_groups <- reactiveVal(character(0))
  
  # All available groups are updated centrally via the observe(.refresh_available_groups())
  # block below — no separate per-tab observers needed.
  
  # ── PCA group selector output ─────────────────────────────────────────────────
  output$pca_group_select_ui <- renderUI({
    grps <- pca_available_groups()
    if (length(grps) == 0)
      return(helpText(tags$small(style = "color:#888;",
                                 "Select token(s) above to populate groups, or click Preview tokens in the sidebar.")))
    cur <- isolate(input$pca_selected_groups)
    sel <- if (!is.null(cur) && length(cur) > 0 && all(cur %in% grps)) cur else grps
    tagList(
      tags$label("Groups to plot", class = "control-label"),
      tags$div(
        style = "display:flex; gap:6px; margin-bottom:4px;",
        actionButton("pca_select_all",  "All",  class = "btn-xs btn-default",
                     style = "font-size:11px; padding:1px 7px;"),
        actionButton("pca_select_none", "None", class = "btn-xs btn-default",
                     style = "font-size:11px; padding:1px 7px;")
      ),
      checkboxGroupInput("pca_selected_groups", label = NULL,
                         choices = grps, selected = sel)
    )
  })
  
  observeEvent(input$pca_select_all, {
    updateCheckboxGroupInput(session, "pca_selected_groups",
                             choices  = pca_available_groups(),
                             selected = pca_available_groups())
  })
  observeEvent(input$pca_select_none, {
    updateCheckboxGroupInput(session, "pca_selected_groups",
                             choices  = pca_available_groups(),
                             selected = character(0))
  })
  
  # ── Stats group selector output ───────────────────────────────────────────────
  output$stats_group_select_ui <- renderUI({
    grps <- stats_available_groups()
    if (length(grps) == 0)
      return(helpText(tags$small(style = "color:#888;",
                                 "Select token(s) above to populate groups, or click Preview tokens in the sidebar.")))
    cur <- isolate(input$stats_selected_groups)
    sel <- if (!is.null(cur) && length(cur) > 0 && all(cur %in% grps)) cur else grps
    tagList(
      tags$label("Groups to compare", class = "control-label"),
      tags$div(
        style = "display:flex; gap:6px; margin-bottom:4px;",
        actionButton("stats_select_all",  "All",  class = "btn-xs btn-default",
                     style = "font-size:11px; padding:1px 7px;"),
        actionButton("stats_select_none", "None", class = "btn-xs btn-default",
                     style = "font-size:11px; padding:1px 7px;")
      ),
      checkboxGroupInput("stats_selected_groups", label = NULL,
                         choices = grps, selected = sel)
    )
  })
  
  observeEvent(input$stats_select_all, {
    updateCheckboxGroupInput(session, "stats_selected_groups",
                             choices  = stats_available_groups(),
                             selected = stats_available_groups())
  })
  observeEvent(input$stats_select_none, {
    updateCheckboxGroupInput(session, "stats_selected_groups",
                             choices  = stats_available_groups(),
                             selected = character(0))
  })
  # ── Central token selector UI (rendered in Group Preview tab) ─────────────────
  
  # Helper: build token choices from preview table
  .token_choices <- reactive({
    tbl <- tryCatch(lr_token_preview(), error = function(e) NULL)
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    delim <- if (identical(input$lr_group_mode, "delimiter") &&
                 !is.null(input$lr_delimiter) && nzchar(input$lr_delimiter))
      input$lr_delimiter else "_"
    build_token_choices(tbl$sample, delimiter = delim)
  })
  
  # Group token selector (shown in Group Preview tab)
  output$lr_token_selector_ui <- renderUI({
    choices <- .token_choices()
    if (is.null(choices)) {
      return(helpText(tags$small(style = "color:#888;",
                                 "Click 'Preview tokens' to populate token positions.")))
    }
    prev <- isolate(input$lr_group_token)
    sel  <- if (!is.null(prev) && length(prev) > 0 && all(as.character(prev) %in% choices))
      as.character(prev) else as.character(choices[1])
    checkboxGroupInput("lr_group_token",
                       label    = "Group token(s) — select one or more to combine",
                       choices  = choices,
                       selected = sel,
                       inline   = TRUE)
  })
  
  # Factor A token selector (Group Preview tab, two-way design)
  output$lr_factorA_selector_ui <- renderUI({
    choices <- .token_choices()
    if (is.null(choices)) return(helpText(tags$small(style="color:#888;", "Preview tokens first.")))
    prev <- isolate(input$lr_factorA_token)
    sel  <- if (!is.null(prev) && as.character(prev) %in% choices) as.character(prev) else as.character(choices[1])
    tagList(
      tags$label("Factor A token", class = "control-label"),
      selectInput("lr_factorA_token", label = NULL, choices = choices, selected = sel)
    )
  })
  
  # Factor B token selector (Group Preview tab, two-way design)
  output$lr_factorB_selector_ui <- renderUI({
    choices <- .token_choices()
    if (is.null(choices)) return(helpText(tags$small(style="color:#888;", "Preview tokens first.")))
    prev <- isolate(input$lr_factorB_token)
    n    <- length(choices)
    def  <- if (n >= 2) as.character(choices[2]) else as.character(choices[1])
    sel  <- if (!is.null(prev) && as.character(prev) %in% choices) as.character(prev) else def
    tagList(
      tags$label("Factor B token", class = "control-label"),
      selectInput("lr_factorB_token", label = NULL, choices = choices, selected = sel)
    )
  })
  
  # Technical replicate parent-sample token selector (Group Preview tab)
  output$lr_techrep_token_selector_ui <- renderUI({
    choices <- .token_choices()
    if (is.null(choices)) {
      return(helpText(tags$small(style = "color:#888;",
                                 "Click 'Preview tokens' to populate token positions.")))
    }
    prev <- isolate(input$lr_techrep_token)
    sel  <- if (!is.null(prev) && length(prev) > 0 && all(as.character(prev) %in% choices))
      as.character(prev) else character(0)
    checkboxGroupInput("lr_techrep_token",
                       label    = "Biological-sample token(s) — select one or more to combine",
                       choices  = choices,
                       selected = sel,
                       inline   = TRUE)
  })
  
  # PCA and Stats group token UIs — read-only mirrors (editing happens in Group Preview only)
  # Rendering a second checkboxGroupInput with the same ID causes Shiny to duplicate values.
  # Instead show a text summary of the current selection.
  .token_mirror_ui <- function() {
    sel <- input$lr_group_token
    if (is.null(sel) || length(sel) == 0)
      return(helpText(tags$small(style = "color:#888;",
                                 "No token selected — go to the Group Preview tab to select group tokens.")))
    tagList(
      helpText(tags$small(
        tags$strong("Active group token(s): "),
        paste(sel, collapse = " + "),
        br(),
        tags$em("Edit in the Group Preview tab.")
      ))
    )
  }
  
  output$pca_group_token_ui   <- renderUI({ .token_mirror_ui() })
  output$stats_group_token_ui <- renderUI({ .token_mirror_ui() })
  output$path_group_token_ui  <- renderUI({ .token_mirror_ui() })
  output$net_group_token_ui   <- renderUI({ .token_mirror_ui() })
  
  # ── Central grouping reactive — single source of truth ───────────────────────
  # Returns a named list: list(group = vector, factorA = vector, factorB = vector)
  # for a given vector of sample names.
  # Priority: CSV > token > fallback (first token)
  get_active_grouping <- reactive({
    function(samples) {
      samples <- as.character(samples)
      n <- length(samples)
      
      # Guard: if no samples or no grouping inputs initialised yet, return sample-level grouping
      if (n == 0) return(list(group = character(0),
                              factorA = character(0),
                              factorB = character(0)))
      
      # 1. CSV override
      if (isTRUE(input$use_group_csv) && !is.null(input$group_csv)) {
        gm <- tryCatch(group_map(), error = function(e) NULL)
        if (!is.null(gm)) {
          grp <- prefer_csv_group_labels(samples, "Delimiter-based", "_", 1, "^([^_]+)", gm)
          fA  <- prefer_csv_factor_labels(samples, "A", "Delimiter-based", "_", 1, "^([^_]+)", gm)
          fB  <- prefer_csv_factor_labels(samples, "B", "Delimiter-based", "_", 2, "^([^_]+)", gm)
          return(list(group = grp, factorA = fA, factorB = fB))
        }
      }
      
      # Safe helper: get a single finite integer from an input, or NA
      safe_int <- function(x) {
        if (is.null(x) || length(x) == 0) return(NA_integer_)
        v <- suppressWarnings(as.integer(x[1]))
        if (is.na(v) || !is.finite(v)) NA_integer_ else v
      }
      
      # 2. Regex mode
      if (isTRUE(identical(input$lr_group_mode, "regex"))) {
        flags <- if (isTRUE(input$lr_regex_ignore_case)) "(?i)" else ""
        pat   <- paste0(flags, input$lr_group_regex %||% "^([^_]+)")
        m     <- tryCatch(stringr::str_match(samples, pat), error = function(e) NULL)
        if (!is.null(m) && ncol(m) >= 2) {
          grp <- dplyr::coalesce(m[, 2], NA_character_)
        } else {
          grp <- rep(NA_character_, length(samples))
        }
        grp[is.na(grp) | !nzchar(grp)] <- "Unassigned"
        return(list(group = grp,
                    factorA = rep(NA_character_, length(samples)),
                    factorB = rep(NA_character_, length(samples))))
      }
      
      # 3. Delimiter mode
      delim <- input$lr_delimiter %||% "_"
      if (!nzchar(delim)) delim <- "_"
      sep   <- input$lr_group_sep %||% "_"
      
      tok <- input$lr_group_token
      if (!is.null(tok) && length(tok) > 0) {
        tok_int <- suppressWarnings(as.integer(tok))
        tok_int <- tok_int[!is.na(tok_int) & is.finite(tok_int)]
      } else {
        tok_int <- integer(0)
      }
      
      if (length(tok_int) > 0) {
        grp <- extract_combined_tokens(samples,
                                       delimiter     = delim,
                                       token_indices = tok_int,
                                       sep           = sep)
      } else {
        grp <- extract_token(samples, delimiter = delim, token_index = 1L)
      }
      grp[is.na(grp) | !nzchar(grp)] <- "Unassigned"
      
      # Factor A / B — safe single-integer extraction
      tokA <- safe_int(input$lr_factorA_token)
      tokB <- safe_int(input$lr_factorB_token)
      
      fA <- if (!is.na(tokA))
        extract_token(samples, delimiter = delim, token_index = tokA)
      else rep(NA_character_, length(samples))
      
      fB <- if (!is.na(tokB))
        extract_token(samples, delimiter = delim, token_index = tokB)
      else rep(NA_character_, length(samples))
      
      list(group = grp, factorA = fA, factorB = fB)
    }
  })
  
  # ── Central technical-replicate parent-sample reactive — single source of truth ─
  # Returns a function(samples) -> character vector of "biological/parent sample"
  # labels, one per input sample. Priority: grouping-CSV 'bio_sample' column (if
  # selected & present) > sample-name token(s) selected in Group Preview > sample
  # itself unchanged (i.e. no collapsing).
  get_active_tech_rep_map <- reactive({
    function(samples) {
      samples <- as.character(samples)
      n <- length(samples)
      if (n == 0) return(character(0))
      
      # 1. CSV column, if that source is selected and the column exists
      if (identical(input$tech_rep_source, "csv") &&
          isTRUE(input$use_group_csv) && !is.null(input$group_csv)) {
        gm <- tryCatch(group_map(), error = function(e) NULL)
        if (!is.null(gm) && !is.null(gm$bio_sample)) {
          bio <- labels_from_map(samples, gm$bio_sample)
          missing <- is.na(bio) | !nzchar(bio)
          bio[missing] <- samples[missing]
          return(bio)
        }
        # No bio_sample column in the CSV — fall through to token method below
      }
      
      # 2. Sample-name token(s)
      delim <- input$lr_delimiter %||% "_"
      if (!nzchar(delim)) delim <- "_"
      sep <- input$lr_group_sep %||% "_"
      
      tok <- input$lr_techrep_token
      if (!is.null(tok) && length(tok) > 0) {
        tok_int <- suppressWarnings(as.integer(tok))
        tok_int <- tok_int[!is.na(tok_int) & is.finite(tok_int)]
      } else {
        tok_int <- integer(0)
      }
      
      if (length(tok_int) > 0) {
        subdelim <- input$lr_techrep_subdelim %||% ""
        if (nzchar(subdelim)) {
          # A single token encodes both bio & tech rep together (e.g. "1-1" = bio 1,
          # tech 1). Extract each selected token individually, keep only the part
          # before the sub-delimiter (the biological replicate index), then combine.
          bio_part_for <- function(tok_idx) {
            raw <- extract_token(samples, delimiter = delim, token_index = tok_idx)
            out <- stringr::str_split_fixed(raw, stringr::fixed(subdelim), 2)[, 1]
            out[is.na(raw)] <- NA_character_
            out
          }
          parts_mat <- vapply(tok_int, bio_part_for, character(length(samples)))
          if (is.null(dim(parts_mat))) parts_mat <- matrix(parts_mat, ncol = length(tok_int))
          bio <- apply(parts_mat, 1, function(r) paste(r[nzchar(r) & !is.na(r)], collapse = sep))
        } else {
          bio <- extract_combined_tokens(samples, delimiter = delim,
                                         token_indices = tok_int, sep = sep)
        }
      } else {
        # No token selected — each sample is its own biological sample (no averaging)
        bio <- samples
      }
      missing <- is.na(bio) | !nzchar(bio)
      bio[missing] <- samples[missing]
      bio
    }
  })
  
  # ── Update available groups whenever token selection or grouping mode changes ──
  .refresh_available_groups <- reactive({
    tbl <- tryCatch(lr_token_preview(), error = function(e) NULL)
    if (is.null(tbl) || nrow(tbl) == 0) return(character(0))
    
    # Exclude blank and iQC samples before deriving group labels
    samples <- tbl$sample
    is_blank <- grepl("^blank", samples, ignore.case = TRUE)
    is_iqc   <- grepl("iqc",   samples, ignore.case = TRUE)
    samples  <- samples[!is_blank & !is_iqc]
    if (length(samples) == 0) return(character(0))
    
    fn  <- get_active_grouping()
    res <- fn(samples)
    grps <- res$group
    
    # Drop NA, empty, Unassigned, and any group label that looks like Blank/iQC/QC/ISTD/ITSD
    grps <- grps[!is.na(grps) & nzchar(grps) & grps != "Unassigned"]
    grps <- grps[!grepl("^(blank|iqc|qc|istd|itsd)$", grps, ignore.case = TRUE)]
    sort(unique(grps))
  })
  
  observe({
    grps <- tryCatch(.refresh_available_groups(), error = function(e) character(0))
    pca_available_groups(grps)
    stats_available_groups(grps)
    met_available_groups(grps)
    class_all_available_groups(grps)
  })
  
  # ── Plot tab group reactiveVals ───────────────────────────────────────────────
  met_available_groups       <- reactiveVal(character(0))
  class_all_available_groups <- reactiveVal(character(0))
  net_available_groups       <- reactiveVal(character(0))
  
  # Network group selector — only groups that were included in statistics
  observe({
    all_grps  <- tryCatch(.refresh_available_groups(), error = function(e) character(0))
    sel_grps  <- input$stats_selected_groups
    # If stats groups are set, restrict to those; otherwise fall back to all
    net_grps  <- if (!is.null(sel_grps) && length(sel_grps) > 0)
      intersect(all_grps, sel_grps)
    else
      all_grps
    net_available_groups(net_grps)
  })
  
  # Plot single lipid — group selector
  output$met_group_select_ui <- renderUI({
    grps <- met_available_groups()
    if (length(grps) == 0)
      return(helpText(tags$small(style = "color:#888;",
                                 "Preview tokens in the Group Preview tab to populate groups.")))
    cur <- isolate(input$met_selected_groups)
    sel <- if (!is.null(cur) && length(cur) > 0 && all(cur %in% grps)) cur else grps
    tagList(
      tags$div(
        style = "display:flex; gap:6px; margin-bottom:4px;",
        actionButton("met_select_all",  "All",  class = "btn-xs btn-default",
                     style = "font-size:11px; padding:1px 7px;"),
        actionButton("met_select_none", "None", class = "btn-xs btn-default",
                     style = "font-size:11px; padding:1px 7px;")
      ),
      checkboxGroupInput("met_selected_groups", label = NULL, choices = grps, selected = sel)
    )
  })
  observeEvent(input$met_select_all,  {
    updateCheckboxGroupInput(session, "met_selected_groups",
                             choices = met_available_groups(), selected = met_available_groups())
  })
  observeEvent(input$met_select_none, {
    updateCheckboxGroupInput(session, "met_selected_groups",
                             choices = met_available_groups(), selected = character(0))
  })
  
  # Plot lipids by class — group selector (single group for the bar chart)
  output$class_all_group_select_ui <- renderUI({
    grps <- class_all_available_groups()
    if (length(grps) == 0)
      return(helpText(tags$small(style = "color:#888;",
                                 "Preview tokens in the Group Preview tab to populate groups.")))
    cur <- isolate(input$class_all_selected_group)
    sel <- if (!is.null(cur) && cur %in% grps) cur else grps[1]
    selectInput("class_all_selected_group", label = NULL, choices = grps, selected = sel)
  })
  
  # -------- Sample/Blank helpers --------
  normalize_sample_name <- function(x) {
    x <- if (is.null(x)) character(0) else x
    stringr::str_trim(stringr::str_to_lower(x))
  }
  
  filter_blanks <- function(df, exclude, sample_col = "sample", exact = TRUE) {
    if (!isTRUE(exclude) || is.null(df) || nrow(df) == 0) return(df)
    sn <- normalize_sample_name(df[[sample_col]])
    if (exact) {
      df[sn != "blank", , drop = FALSE]
    } else {
      df[!startsWith(sn, "blank"), , drop = FALSE]
    }
  }
  
  filter_blank_rownames <- function(mat, exclude, exact = TRUE) {
    if (!isTRUE(exclude) || is.null(mat) || nrow(mat) == 0) return(mat)
    rn <- normalize_sample_name(rownames(mat))
    if (exact) {
      mat[rn != "blank", , drop = FALSE]
    } else {
      mat[!startsWith(rn, "blank"), , drop = FALSE]
    }
  }
  
  #----ReactiveVal for stats
  
  stats_results_val <- reactiveVal(NULL)
  stats_results_all <- reactiveVal(NULL)  # always full unfiltered results for volcano
  
  observeEvent(input$run_stats, {
    # Full unfiltered results — always all tested features (used by volcano plot)
    all_stats <- .compute_stats_results(
      df = stats_input_long(),
      test_choice = input$stats_test,
      equal_var = isTRUE(input$ttest_equal_var),
      tw_filter = input$twoway_effect_filter,
      alpha = input$alpha %||% 0.05,
      padj_method = input$padj_method,
      use_adj = isTRUE(input$use_adj_threshold),
      show_all_rows_flag = TRUE  # always keep all rows
    )
    stats_results_all(all_stats)
    
    # Filtered results — respects show_all_rows toggle (used by bar plots, heatmap, etc.)
    stats <- .compute_stats_results(
      df = stats_input_long(),
      test_choice = input$stats_test,
      equal_var = isTRUE(input$ttest_equal_var),
      tw_filter = input$twoway_effect_filter,
      alpha = input$alpha %||% 0.05,
      padj_method = input$padj_method,
      use_adj = isTRUE(input$use_adj_threshold),
      show_all_rows_flag = isTRUE(input$show_all_rows)
    )
    stats_results_val(stats)
  })
  
  
  output$uploadedNames <- renderText({
    req(input$msdial_txts)
    paste("Selected files:", paste(input$msdial_txts$name, collapse = ", "))
  })
  
  output$loadedFiles <- renderText({
    req(uploaded_data())
    files <- attr(uploaded_data(), "source_files")
    paste("Loaded files:", paste(files, collapse = ", "))
  })
  
  
  pick_measure_col <- function(kind) if (identical(kind, "value_bs")) "value_bs" else "norm"
  axis_label <- function(kind, mode, unit_label = NULL) {
    if (identical(mode, "percent")) return("Percent of class total (%)")
    if (identical(kind, "value_bs")) return("Background-subtracted (a.u.)")
    # kind == "norm"
    ulab <- unit_label %||% "pmol"
    paste0("Lipid amount (", ulab, ")")
  }
  
  uploaded_data <- eventReactive(
    input$load_data,
    {
      validate(need(is.numeric(input$ISTD_vol) && !is.na(input$ISTD_vol) && input$ISTD_vol > 0,
                    "Please enter a positive ISTD amount (pmol) before loading."))
      validate(need(!is.null(input$msdial_txts) && nrow(input$msdial_txts) >= 1,
                    "Upload at least one MS-DIAL .txt file."))
      validate(need(nrow(input$msdial_txts) <= 2,
                    "Upload no more than two MS-DIAL .txt files."))
      
      withProgress(message = "Loading MS-DIAL files…", value = 0, {
        incProgress(0.20, detail = "Reading uploaded .txt file(s)")
        dat <- load_msdial_uploaded(input$msdial_txts)
        
        incProgress(0.80, detail = "Preparing tidy format")
        dat
      })
    },
    ignoreInit = TRUE
  )
  
  # ---- ISTD CSV reader & validator ----
  istd_map_df <- reactive({
    req(input$istd_map_csv)
    df <- readr::read_csv(input$istd_map_csv$datapath, show_col_types = FALSE)
    # Basic column validation
    required <- c("ISTD", "amount", "units")
    missing  <- setdiff(required, names(df))
    validate(need(length(missing) == 0,
                  paste0("ISTD CSV must have columns: ", paste(required, collapse=", "), 
                         ". Missing: ", paste(missing, collapse=", "))))
    # Coerce types and trim names
    df <- df %>%
      dplyr::mutate(
        ISTD   = stringr::str_trim(as.character(.data$ISTD)),
        amount = suppressWarnings(as.numeric(.data$amount)),
        units  = stringr::str_trim(as.character(.data$units))
      )
    validate(need(all(is.finite(df$amount)),
                  "ISTD CSV: 'amount' column must be numeric (all rows)."))
    # Remove obvious empty rows
    df <- df %>% dplyr::filter(ISTD != "", !is.na(amount), units != "")
    df
  })
  
  # ---- ISTD CSV summary against current dataset ----
  output$istd_map_summary <- renderText({
    if (is.null(input$istd_map_csv)) return("No ISTD CSV uploaded.")
    df <- tryCatch(istd_map_df(), error = function(e) NULL)
    if (is.null(df)) return("ISTD CSV error — check columns: ISTD, amount, units.")
    is_table <- tryCatch(is_table_reactive(), error = function(e) NULL)
    if (is.null(is_table) || nrow(is_table) == 0) {
      return(paste0(
        "ISTD CSV ok. Rows: ", nrow(df),
        " • Units present: ", paste(sort(unique(df$units)), collapse = ", "),
        "\n(Load MS-DIAL files first to see matching stats.)"
      ))
    }
    if (!"istd_name" %in% names(is_table)) {
      return(paste0(
        "ISTD CSV ok. Rows: ", nrow(df),
        " • Units present: ", paste(sort(unique(df$units)), collapse = ", ")
      ))
    }
    is_names  <- sort(unique(is_table$istd_name))
    matched   <- sum(is_names %in% df$ISTD)
    unmatched <- setdiff(is_names, df$ISTD)
    paste0(
      "ISTD CSV ok. Rows: ", nrow(df),
      " • Units present: ", paste(sort(unique(df$units)), collapse = ", "),
      "\nMatched [IS]: ", matched, " / ", length(is_names),
      if (length(unmatched) > 0)
        paste0("\nUnmatched (fallback to numeric): ",
               paste(unmatched, collapse = ", ")) else ""
    )
  })
  
  # ── Promoted reactive: avoids repeated attr(bg_norm_long(), "is_table") calls ──
  is_table_reactive <- reactive({
    req(bg_norm_long())
    attr(bg_norm_long(), "is_table")
  })
  
  # ---- Protein normalisation CSV reader ----
  protein_df <- reactive({
    req(input$protein_csv)
    df <- readr::read_csv(input$protein_csv$datapath, show_col_types = FALSE)
    required <- c("sample", "protein")
    missing  <- setdiff(required, names(df))
    validate(need(length(missing) == 0,
                  paste0("Protein CSV must have columns: sample, protein. Missing: ",
                         paste(missing, collapse = ", "))))
    df <- df %>%
      dplyr::mutate(
        sample  = stringr::str_trim(as.character(.data$sample)),
        protein = suppressWarnings(as.numeric(.data$protein))
      )
    validate(need(all(is.finite(df$protein)),
                  "Protein CSV: 'protein' column must be numeric (all rows)."))
    df
  })
  
  output$protein_csv_summary <- renderText({
    if (is.null(input$protein_csv)) return("No protein CSV uploaded.")
    df <- tryCatch(protein_df(), error = function(e) NULL)
    if (is.null(df)) return("Protein CSV error — check columns: sample, protein.")
    paste0(
      "Protein CSV ok. Samples: ", nrow(df),
      " • Range: ", round(min(df$protein, na.rm = TRUE), 3),
      " – ", round(max(df$protein, na.rm = TRUE), 3)
    )
  })
  
  # ---- Protein normalisation sample-matching diagnostic ----
  # Compares distinct imported MS-DIAL samples (sample_norm) against the
  # uploaded protein CSV, independent of whether "Apply protein normalisation"
  # is currently ticked, so mismatches can be spotted before enabling it.
  protein_match_table <- reactive({
    req(bg_norm_long())
    
    msdial_samples <- bg_norm_long() %>%
      dplyr::distinct(sample, sample_norm) %>%
      dplyr::filter(!is_protected_sample(sample)) %>%
      dplyr::arrange(sample)
    
    base_missing_csv <- msdial_samples %>%
      dplyr::transmute(
        Sample            = sample,
        `In MS-DIAL Data` = TRUE,
        `In Protein CSV`  = FALSE,
        `Protein Value`   = NA_real_,
        Status            = "No protein CSV uploaded"
      )
    
    if (is.null(input$protein_csv)) return(base_missing_csv)
    
    prot <- tryCatch(protein_df(), error = function(e) NULL)
    if (is.null(prot)) {
      return(dplyr::mutate(base_missing_csv, Status = "Protein CSV error — check columns: sample, protein."))
    }
    
    prot_norm <- prot %>%
      dplyr::filter(!is_protected_sample(sample)) %>%
      dplyr::mutate(sample_norm = normalize_sample_name(sample)) %>%
      dplyr::rename(protein_sample = sample)
    
    msdial_samples %>%
      dplyr::full_join(prot_norm, by = "sample_norm") %>%
      dplyr::mutate(
        Sample            = dplyr::coalesce(sample, protein_sample),
        `In MS-DIAL Data` = !is.na(sample),
        `In Protein CSV`  = !is.na(protein_sample),
        `Protein Value`   = protein,
        Status = dplyr::case_when(
          `In MS-DIAL Data` & `In Protein CSV`  ~ "Matched",
          `In MS-DIAL Data` & !`In Protein CSV` ~ "Missing from protein CSV",
          !`In MS-DIAL Data` & `In Protein CSV` ~ "Extra in protein CSV (no matching sample)",
          TRUE ~ "Unknown"
        )
      ) %>%
      dplyr::select(Sample, `In MS-DIAL Data`, `In Protein CSV`, `Protein Value`, Status) %>%
      dplyr::arrange(Status != "Matched", Sample)
  })
  
  output$protein_match_summary <- renderText({
    df <- protein_match_table()
    req(df)
    
    if (is.null(input$protein_csv)) return("Upload a protein CSV (Background & Normalisation sidebar) to check sample matching.")
    
    n_total   <- sum(df$`In MS-DIAL Data`)
    n_matched <- sum(df$Status == "Matched")
    n_missing <- sum(df$Status == "Missing from protein CSV")
    n_extra   <- sum(df$Status == "Extra in protein CSV (no matching sample)")
    
    paste0(
      n_matched, " of ", n_total, " imported samples matched to protein data.",
      if (n_missing > 0) paste0("\n\u26a0 ", n_missing, " sample(s) have no protein value — will NOT be protein-normalised.") else "",
      if (n_extra   > 0) paste0("\n\u26a0 ", n_extra,   " protein CSV row(s) do not match any imported sample — check for typos.") else "",
      if (n_missing == 0 && n_extra == 0) "\n\u2713 All samples matched." else ""
    )
  })
  
  output$proteinMatchTable <- DT::renderDT({
    df <- protein_match_table()
    
    DT::datatable(
      df,
      rownames = FALSE,
      options  = list(pageLength = 15, order = list()),
      class    = "stripe hover"
    ) %>%
      DT::formatStyle(
        "Status",
        target = "row",
        backgroundColor = DT::styleEqual(
          c("Matched",
            "Missing from protein CSV",
            "Extra in protein CSV (no matching sample)",
            "No protein CSV uploaded",
            "Protein CSV error — check columns: sample, protein."),
          c("#e6f4ea", "#fdecea", "#fdecea", "#f5f5f5", "#fdecea")
        )
      )
  })
  
  # ---- Outlier Detection tab outputs ----
  .build_outlier_sample_plot <- function(fsz = 13) {
    df <- sample_outlier_flags()
    validate(need(nrow(df) > 0,
                  "Not enough samples/features for PCA-based outlier detection (need \u2265 5 non-blank samples and \u2265 2 informative features)."))
    df <- df %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(sample = factor(sample, levels = sample))
    
    ggplot2::ggplot(df, ggplot2::aes(x = sample, y = score, fill = is_outlier)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = df$threshold[1], linetype = "dashed", colour = "red") +
      ggplot2::scale_fill_manual(values = c(`FALSE` = "#4c72b0", `TRUE` = "#c44e52"),
                                 labels = c("Normal", "Flagged"), name = NULL) +
      ggplot2::labs(x = NULL, y = "Hotelling's T\u00b2") +
      ggplot2::theme_minimal(base_size = fsz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  .build_outlier_iqc_plot <- function(fsz = 13) {
    df <- iqc_outlier_flags()
    validate(need(nrow(df) > 0, "Fewer than 3 iQC samples detected \u2014 iQC deviation check skipped."))
    df <- df %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(sample = factor(sample, levels = sample))
    
    ggplot2::ggplot(df, ggplot2::aes(x = sample, y = score, fill = is_outlier)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = df$threshold[1], linetype = "dashed", colour = "red") +
      ggplot2::scale_fill_manual(values = c(`FALSE` = "#4c72b0", `TRUE` = "#c44e52"),
                                 labels = c("Normal", "Flagged"), name = NULL) +
      ggplot2::labs(x = NULL, y = "Median relative deviation") +
      ggplot2::theme_minimal(base_size = fsz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  output$outlierSamplePlot <- renderPlot({
    .build_outlier_sample_plot(input$outlier_export_fontsize %||% 13)
  })
  
  output$outlierIqcPlot <- renderPlot({
    .build_outlier_iqc_plot(input$outlier_export_fontsize %||% 13)
  })
  
  # Export helpers
  .outlier_export_dims <- function() {
    px_w  <- input$outlier_export_width  %||% 1200
    px_h  <- input$outlier_export_height %||% 500
    dpi   <- input$outlier_export_dpi    %||% 300
    scale <- input$outlier_export_scale  %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  output$download_outlier_sample_png <- downloadHandler(
    filename = function() paste0("outlier_sample_hotelling_", Sys.Date(), ".png"),
    content = function(file) {
      dims <- isolate(.outlier_export_dims())
      fsz  <- isolate(input$outlier_export_fontsize %||% 13)
      p    <- isolate(.build_outlier_sample_plot(fsz))
      validate(need(!is.null(p), "No plot to export."))
      ggplot2::ggsave(file, plot = p,
                      width = dims$w, height = dims$h,
                      dpi = dims$dpi, device = "png")
    }
  )
  
  output$download_outlier_sample_svg <- downloadHandler(
    filename = function() paste0("outlier_sample_hotelling_", Sys.Date(), ".svg"),
    content = function(file) {
      dims <- isolate(.outlier_export_dims())
      fsz  <- isolate(input$outlier_export_fontsize %||% 13)
      p    <- isolate(.build_outlier_sample_plot(fsz))
      validate(need(!is.null(p), "No plot to export."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  output$download_outlier_iqc_png <- downloadHandler(
    filename = function() paste0("outlier_iqc_deviation_", Sys.Date(), ".png"),
    content = function(file) {
      dims <- isolate(.outlier_export_dims())
      fsz  <- isolate(input$outlier_export_fontsize %||% 13)
      p    <- isolate(.build_outlier_iqc_plot(fsz))
      validate(need(!is.null(p), "No plot to export."))
      ggplot2::ggsave(file, plot = p,
                      width = dims$w, height = dims$h,
                      dpi = dims$dpi, device = "png")
    }
  )
  
  output$download_outlier_iqc_svg <- downloadHandler(
    filename = function() paste0("outlier_iqc_deviation_", Sys.Date(), ".svg"),
    content = function(file) {
      dims <- isolate(.outlier_export_dims())
      fsz  <- isolate(input$outlier_export_fontsize %||% 13)
      p    <- isolate(.build_outlier_iqc_plot(fsz))
      validate(need(!is.null(p), "No plot to export."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  output$outlier_feature_summary <- renderText({
    df <- feature_outlier_flags()
    if (nrow(df) == 0) return("No feature-level outliers flagged with current settings.")
    paste0(
      nrow(df), " flagged data point(s) across ",
      dplyr::n_distinct(df$`Metabolite name`), " feature(s) and ",
      dplyr::n_distinct(df$group), " group(s)."
    )
  })
  
  output$outlierFeatureTable <- DT::renderDT({
    df <- feature_outlier_flags()
    validate(need(nrow(df) > 0, "No feature-level outliers flagged with current settings."))
    
    keys <- .feature_key(df$sample, df[["Metabolite name"]])
    
    df %>%
      dplyr::mutate(
        `Manual override` = dplyr::case_when(
          keys %in% manual_feature_keep()    ~ "Kept",
          keys %in% manual_feature_exclude() ~ "Excluded",
          TRUE ~ ""
        ),
        `Effective status` = dplyr::case_when(
          keys %in% manual_feature_keep()         ~ "Included",
          keys %in% manual_feature_exclude()       ~ "Excluded",
          isTRUE(input$outlier_exclude_features)   ~ "Excluded",
          TRUE ~ "Included"
        )
      ) %>%
      dplyr::transmute(
        Sample              = sample,
        `Metabolite name`   = `Metabolite name`,
        Class               = class,
        Group               = group,
        Value               = signif(value, 4),
        Score               = round(score, 2),
        Method              = method,
        `Manual override`,
        `Effective status`
      ) %>%
      DT::datatable(rownames = FALSE, selection = "multiple", options = list(pageLength = 15)) %>%
      DT::formatStyle(
        "Effective status", target = "row",
        backgroundColor = DT::styleEqual(c("Excluded", "Included"), c("#fdecea", "white"))
      )
  })
  
  # ---- Read uploaded CSV and build mapping ----
  group_csv_df <- reactive({
    req(input$group_csv)
    readr::read_csv(input$group_csv$datapath, show_col_types = FALSE)
  })
  
  group_map <- reactive({
    req(input$group_csv)
    build_group_map(group_csv_df())
  })
  
  
  output$group_csv_summary <- renderText({
    if (is.null(input$group_csv)) return("No CSV uploaded.")
    gm <- try(group_map(), silent = TRUE)
    if (inherits(gm, "try-error")) return(paste("CSV error:", gm))
    
    # Extract unique values for each factor (case insensitive)
    group_levels  <- if (!is.null(gm$group))   unique(tolower(na.omit(gm$group)))   else NULL
    factorA_levels <- if (!is.null(gm$factorA)) unique(tolower(na.omit(gm$factorA))) else NULL
    factorB_levels <- if (!is.null(gm$factorB)) unique(tolower(na.omit(gm$factorB))) else NULL
    
    
    # Compose summary
    summary <- paste0(
      "CSV columns → sample: '", gm$colnames$sample,
      "', group: '", gm$colnames$group,
      if (!is.na(gm$colnames$factorA)) paste0("', factorA: '", gm$colnames$factorA) else "",
      if (!is.na(gm$colnames$factorB)) paste0("', factorB: '", gm$colnames$factorB) else "",
      "'."
    )
    if (!is.null(group_levels)) {
      summary <- paste0(summary, "\nGroups found: ", paste(group_levels, collapse = ", "))
    }
    if (!is.null(factorA_levels)) {
      summary <- paste0(summary, "\nFactorA levels: ", paste(factorA_levels, collapse = ", "))
    }
    if (!is.null(factorB_levels)) {
      summary <- paste0(summary, "\nFactorB levels: ", paste(factorB_levels, collapse = ", "))
    }
    summary
  })
  
  
  # ---- Clean the data ----
  data_clean <- reactive({
    req(uploaded_data())
    df <- uploaded_data()
    
    keep_meta <- c("Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type")
    
    df_keep <- df %>%
      dplyr::select(
        dplyr::any_of(keep_meta),
        dplyr::matches("(?i)(_pos|_neg)(?:$|_)", perl = TRUE)
      ) %>%
      dplyr::mutate(
        class = dplyr::case_when(
          stringr::str_detect(`Metabolite name`, "^TG") ~ "TG",
          stringr::str_detect(`Metabolite name`, "^d5 TG") ~ "TG", # for d5 ISTD
          TRUE ~ NA_character_
        ),
        class = dplyr::coalesce(
          class,
          stringr::str_extract(`Metabolite name`, "^[A-Za-z]{2,3}(?=\\s*\\(?\\d)"),
          stringr::str_extract(`Metabolite name`, "^[A-Za-z]{2,3}")
        ),
        `Adduct type` = stringr::str_trim(`Adduct type`),
        ion.mode = dplyr::case_when(
          stringr::str_detect(`Adduct type`, "\\+$") ~ "positive",
          stringr::str_detect(`Adduct type`, "\\-$") ~ "negative",
          TRUE ~ NA_character_
        )
      )
    
    keep_meta2 <- c(keep_meta, "class", "ion.mode")
    
    df_long <- df_keep %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(keep_meta2),
        names_to  = "sample_with_mode",
        values_to = "value",
        values_drop_na = FALSE
      ) %>%
      dplyr::mutate(
        mode = dplyr::case_when(
          stringr::str_detect(sample_with_mode, "(?i)_pos(?:$|_)") ~ "pos",
          stringr::str_detect(sample_with_mode, "(?i)_neg(?:$|_)") ~ "neg",
          TRUE ~ NA_character_
        ),
        # strip mode and any trailing suffix (e.g., '.1')
        sample = stringr::str_replace(sample_with_mode, "(?i)(_pos|_neg).*", "")
      ) %>%
      dplyr::select(dplyr::all_of(keep_meta2), sample, mode, value)
    
    df_dedup <- df_long %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(keep_meta2)), sample) %>%
      dplyr::arrange(factor(mode, levels = c("pos", "neg")), .by_group = TRUE) %>%
      dplyr::summarise(
        value = {
          v <- value
          v <- v[!is.na(v) & v != ""]
          dplyr::first(v, default = NA_character_)
        },
        .groups = "drop"
      )
    
    df_wide <- df_dedup %>%
      tidyr::pivot_wider(
        names_from  = sample,
        values_from = value
      )
    
    df_wide
  })
  
  
  # ---------- Background subtraction + normalisation (IS from RAW intensity) ----------
  bg_norm_long <- reactive({
    req(data_clean())
    
    keep_meta  <- c("Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type", "class", "ion.mode")
    df_wide    <- data_clean()
    sample_cols <- setdiff(colnames(df_wide), keep_meta)
    
    validate(need(length(sample_cols) > 0, "No sample columns found after cleaning."))
    
    # Long with raw intensity 'value'
    df_long <- df_wide %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(sample_cols),
        names_to  = "sample",
        values_to = "value",
        values_drop_na = FALSE
      ) %>%
      dplyr::mutate(
        value = suppressWarnings(as.numeric(value)),
        sample_norm = normalize_sample_name(sample)
      )
    
    # Identify blanks flexibly (blank, blank_1, blank-...)
    is_blank  <- stringr::str_detect(df_long$sample_norm, "^blank(?:$|[\\-_])") | df_long$sample_norm == "blank"
    has_blank <- any(is_blank, na.rm = TRUE)
    
    # Compute blank per metabolite × ion.mode (mean of any blanks present)
    blanks <- df_long %>%
      dplyr::filter(is_blank) %>%
      dplyr::group_by(`Metabolite name`, ion.mode) %>%
      dplyr::summarise(blank_value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Background-subtracted analyte signal (NEVER used for IS)
    df_bs <- df_long %>%
      dplyr::left_join(blanks, by = c("Metabolite name", "ion.mode")) %>%
      dplyr::mutate(
        blank_value = dplyr::coalesce(blank_value, 0),
        value_bs    = pmax(value - blank_value, 0)
      )
    
    # IS table from RAW 'value' (NOT value_bs)
    # Matched by class + ion.mode + adduct type
    # IS rows are always kept from whichever mode detected — never deduplicated
    is_table <- df_long %>%
      dplyr::filter(stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(
        istd_name   = stringr::str_trim(`Metabolite name`),
        istd_adduct = stringr::str_trim(`Adduct type`)
      ) %>%
      dplyr::group_by(sample, sample_norm, class, ion.mode, istd_adduct) %>%
      dplyr::summarise(
        IS_value  = dplyr::first(stats::na.omit(value)),
        istd_name = dplyr::first(na.omit(istd_name)),
        .groups   = "drop"
      )
    
    # Normalise analyte: join IS by class + ion.mode + adduct type
    # If ISTD CSV provided, also join by adduct — produces one row per matched adduct
    df_norm <- df_bs %>%
      dplyr::left_join(
        is_table %>% dplyr::select(sample, sample_norm, class, ion.mode,
                                   istd_adduct, IS_value, istd_name),
        by       = c("sample", "sample_norm", "class", "ion.mode",
                     "Adduct type" = "istd_adduct")
      )
    
    # ── Force Shiny to track toggle dependencies unconditionally ─────────────
    # Reading these here ensures the reactive re-runs when either toggle changes,
    # even if the if() branch below would not have evaluated them.
    use_csv        <- isTRUE(input$use_istd_map_csv)
    istd_csv_file  <- input$istd_map_csv        # register file input as dependency
    istd_vol       <- input$ISTD_vol            # register ISTD vol as dependency
    df_csv_pre     <- if (use_csv) tryCatch(istd_map_df(), error = function(e) NULL) else NULL
    
    if (use_csv && !is.null(df_csv_pre) && nrow(df_csv_pre) > 0) {
      df_norm <- df_norm %>%
        dplyr::left_join(df_csv_pre, by = c("istd_name" = "ISTD")) %>%
        dplyr::mutate(
          istd_amount_use = dplyr::coalesce(.data$amount, as.numeric(istd_vol)),
          istd_units_use  = dplyr::coalesce(.data$units, "pmol"),
          norm       = dplyr::if_else(!is.na(IS_value) & IS_value > 0,
                                      (value_bs / IS_value) * istd_amount_use, NA_real_),
          norm_units = dplyr::if_else(!is.na(norm), istd_units_use, NA_character_)
        )
    } else {
      df_norm <- df_norm %>%
        dplyr::mutate(
          norm       = dplyr::if_else(!is.na(IS_value) & IS_value > 0,
                                      (value_bs / IS_value) * istd_vol, NA_real_),
          norm_units = dplyr::if_else(!is.na(norm), "pmol", NA_character_)
        )
    }
    
    # Preserve IS table on output (already done in your code)
    attr(df_norm, "is_table") <- is_table
    attr(df_norm, "has_blank") <- has_blank
    df_norm
  })
  
  # ---------- Ion mode resolution (generalised for all classes) ----------
  # Rules (user-configurable via Settings tab):
  #   - IS rows: always kept from whatever mode detected — never touched
  #   - Prefer-negative classes: keep neg if present; include pos only if neg absent,
  #     or if a pos species has no neg counterpart (different metabolite name)
  #   - Prefer-positive classes: mirror of above
  #   - All other classes: keep both modes as separate rows
  # Metabolite names are kept exactly as MS-DIAL exports them.
  resolve_ion_modes <- function(df, pref_neg, pref_pos) {
    req_cols <- c("class", "ion.mode", "Metabolite name", "sample")
    missing  <- setdiff(req_cols, names(df))
    if (length(missing)) {
      warning("resolve_ion_modes(): missing column(s): ",
              paste(missing, collapse = ", "), " — returning input unmodified.")
      return(df)
    }
    
    df <- df %>%
      dplyr::mutate(
        class  = tidyr::replace_na(class, "Unknown"),
        is_is  = stringr::str_detect(`Metabolite name`, "\\[IS\\]"),
        # base class strips -O suffix for lookup
        base_class = sub("-O$", "", class)
      )
    
    # IS rows are never touched — pass straight through
    df_is     <- df %>% dplyr::filter(is_is)
    df_analyte <- df %>% dplyr::filter(!is_is)
    
    # Determine mode preference per row
    df_analyte <- df_analyte %>%
      dplyr::mutate(
        mode_pref = dplyr::case_when(
          class %in% pref_neg | base_class %in% pref_neg ~ "neg",
          class %in% pref_pos | base_class %in% pref_pos ~ "pos",
          TRUE ~ "both"
        )
      )
    
    # Split by preference group
    df_both <- df_analyte %>% dplyr::filter(mode_pref == "both")
    df_pref  <- df_analyte %>% dplyr::filter(mode_pref %in% c("neg", "pos"))
    
    if (nrow(df_pref) > 0) {
      df_pref <- df_pref %>%
        dplyr::group_by(`Metabolite name`, sample, mode_pref) %>%
        dplyr::mutate(
          has_preferred = any(
            (mode_pref == "neg" & ion.mode == "negative") |
              (mode_pref == "pos" & ion.mode == "positive"),
            na.rm = TRUE
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::filter(
          # Keep row if: it IS the preferred mode, OR there is no preferred mode present
          (mode_pref == "neg" & ion.mode == "negative") |
            (mode_pref == "pos" & ion.mode == "positive") |
            !has_preferred
        ) %>%
        dplyr::select(-has_preferred)
    }
    
    dplyr::bind_rows(
      df_is,
      df_both,
      df_pref
    ) %>%
      dplyr::select(-is_is, -base_class,
                    -dplyr::any_of(c("mode_pref", "has_preferred")))
  }
  
  # ---------- Cross-mode identity reconciliation (neg identity -> pos quant) ----------
  # For designated classes (default PC, PE, PI, PS, CL) the positive-mode feature is
  # the quantitation source (it has an adduct-matched positive IS), while its
  # acyl-resolved identity is transferred from the negative-mode feature of the same
  # species, matched by sum-composition shorthand (text before "|") + retention time
  # within `rt_tol` minutes. Negative mode resolves individual fatty acids that
  # positive mode reports only as sum composition (or, for CL, as the two combined
  # halves), so identity is authoritative in negative and quantitation in positive.
  #
  #   - Positive feature WITH a negative match : positive sample values are kept
  #       (quant); Metabolite name is replaced by the negative acyl identity.
  #   - Positive feature WITHOUT a negative match : retained and positive-quantified,
  #       shown at the depth positive can determine — sum composition for
  #       glycerophospholipids, positive substructure kept for CL.
  #   - Negative feature WITHOUT a positive match : retained and quantified from its
  #       own (negative) mode IS.
  #   - Multiple negative isomers matching one positive : nearest-RT wins the identity;
  #       any losing isomer with no positive match of its own is retained as neg-only.
  #   - IS rows and all non-designated classes pass through untouched.
  #
  # The main long data frame is returned with an unchanged schema (only rows dropped
  # and Metabolite name rewritten). A per-decision audit tibble is attached as
  # attr(., "cross_mode_map") for display and reproducibility.
  # `pos_substructure_classes`: classes for which the POSITIVE annotation carries
  # genuine (if partial) structure worth keeping when unmatched — e.g. CL, whose
  # positive feature resolves the two combined halves. For all other cross-mode
  # classes (the glycerophospholipids) a positive acyl assignment is unreliable, so
  # an unmatched positive is collapsed to sum composition (the depth positive can
  # actually determine). A matched positive always takes the negative acyl identity.
  resolve_cross_mode <- function(df, cross_classes, rt_tol = 0.1,
                                 pos_substructure_classes = c("CL")) {
    empty_map <- tibble::tibble(
      class = character(0), shorthand = character(0),
      pos_name = character(0), pos_rt = numeric(0), pos_adduct = character(0),
      neg_identity = character(0), neg_rt = numeric(0),
      delta_rt = numeric(0), final_identity = character(0),
      quant_mode = character(0), status = character(0)
    )
    req_cols <- c("class", "ion.mode", "Metabolite name", "Adduct type",
                  "Average Rt(min)", "sample")
    if (length(cross_classes) == 0 || any(!req_cols %in% names(df))) {
      attr(df, "cross_mode_map") <- empty_map
      return(df)
    }
    if (is.na(rt_tol) || rt_tol < 0) rt_tol <- 0.1
    
    shorthand_of <- function(x) stringr::str_trim(sub("\\|.*$", "", x))
    feat_key     <- function(nm, ad, rt)
      paste(nm, ad, sprintf("%.4f", as.numeric(rt)), sep = "\u0001")
    
    df <- df %>%
      dplyr::mutate(.is_is = stringr::str_detect(`Metabolite name`, "\\[IS\\]"))
    in_scope <- (df$class %in% cross_classes) & !df$.is_is
    df_cm   <- df[in_scope, , drop = FALSE]
    df_rest <- df[!in_scope, , drop = FALSE]
    
    if (nrow(df_cm) == 0) {
      out <- df_rest %>% dplyr::select(-.is_is)
      attr(out, "cross_mode_map") <- empty_map
      return(out)
    }
    
    # Feature-level view: one row per physical feature (not per sample).
    # `Average Rt(min)` can arrive as character from the MS-DIAL export, so carry a
    # numeric copy (.rt) for the RT arithmetic and leave the original column untouched.
    feats <- df_cm %>%
      dplyr::distinct(class, ion.mode, `Metabolite name`, `Adduct type`, `Average Rt(min)`) %>%
      dplyr::mutate(shorthand = shorthand_of(`Metabolite name`),
                    .rt = suppressWarnings(as.numeric(`Average Rt(min)`)))
    pos_f <- feats %>% dplyr::filter(ion.mode == "positive")
    neg_f <- feats %>% dplyr::filter(ion.mode == "negative")
    
    matches <- pos_f %>%
      dplyr::inner_join(neg_f, by = c("class", "shorthand"),
                        suffix = c("_pos", "_neg"),
                        relationship = "many-to-many") %>%
      dplyr::filter(!is.na(.rt_pos), !is.na(.rt_neg)) %>%
      dplyr::mutate(delta_rt = abs(.rt_pos - .rt_neg)) %>%
      dplyr::filter(delta_rt <= rt_tol) %>%
      dplyr::group_by(class, `Metabolite name_pos`, `Adduct type_pos`, `Average Rt(min)_pos`) %>%
      dplyr::slice_min(delta_rt, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    # Rewrite matched positive features to their negative acyl identity
    rename_tbl <- if (nrow(matches) > 0) {
      matches %>%
        dplyr::transmute(
          .k = feat_key(`Metabolite name_pos`, `Adduct type_pos`, `Average Rt(min)_pos`),
          new_name = `Metabolite name_neg`
        ) %>%
        dplyr::distinct(.k, .keep_all = TRUE)
    } else {
      tibble::tibble(.k = character(0), new_name = character(0))
    }
    
    # Negative donor features (consumed as identity source -> drop their sample rows)
    neg_consumed <- if (nrow(matches) > 0) {
      feat_key(matches$`Metabolite name_neg`, matches$`Adduct type_neg`,
               matches$`Average Rt(min)_neg`)
    } else character(0)
    
    df_cm <- df_cm %>%
      dplyr::mutate(.k = feat_key(`Metabolite name`, `Adduct type`, `Average Rt(min)`)) %>%
      dplyr::left_join(rename_tbl, by = ".k") %>%
      dplyr::mutate(.shorthand = shorthand_of(`Metabolite name`))
    keep <- !(df_cm$ion.mode == "negative" & df_cm$.k %in% neg_consumed)
    df_cm <- df_cm[keep, , drop = FALSE] %>%
      dplyr::mutate(
        `Metabolite name` = dplyr::case_when(
          # matched positive -> negative acyl identity
          ion.mode == "positive" & !is.na(new_name)              ~ new_name,
          # unmatched positive, class keeps positive substructure (e.g. CL) -> as-is
          ion.mode == "positive" & class %in% pos_substructure_classes ~ `Metabolite name`,
          # unmatched positive glycerophospholipid -> sum composition only
          ion.mode == "positive"                                 ~ .shorthand,
          # negatives (matched donors already dropped) -> keep full acyl identity
          TRUE                                                   ~ `Metabolite name`
        )
      ) %>%
      dplyr::select(-.k, -new_name, -.shorthand)
    
    out <- dplyr::bind_rows(df_rest, df_cm) %>% dplyr::select(-.is_is)
    
    # ---- Audit map ----
    matched_map <- if (nrow(matches) > 0) {
      matches %>%
        dplyr::transmute(
          class, shorthand,
          pos_name = `Metabolite name_pos`, pos_rt = .rt_pos,
          pos_adduct = `Adduct type_pos`,
          neg_identity = `Metabolite name_neg`, neg_rt = .rt_neg,
          delta_rt,
          final_identity = `Metabolite name_neg`, quant_mode = "positive",
          status = "pos quant / neg id"
        )
    } else empty_map
    
    matched_pos_keys <- if (nrow(matches) > 0)
      feat_key(matches$`Metabolite name_pos`, matches$`Adduct type_pos`,
               matches$`Average Rt(min)_pos`) else character(0)
    pos_unmatched <- pos_f %>%
      dplyr::mutate(.k = feat_key(`Metabolite name`, `Adduct type`, `Average Rt(min)`)) %>%
      dplyr::filter(!.k %in% matched_pos_keys) %>%
      dplyr::transmute(
        class, shorthand, pos_name = `Metabolite name`, pos_rt = .rt,
        pos_adduct = `Adduct type`, neg_identity = NA_character_, neg_rt = NA_real_,
        delta_rt = NA_real_,
        final_identity = dplyr::if_else(class %in% pos_substructure_classes,
                                        `Metabolite name`, shorthand),
        quant_mode = "positive",
        status = dplyr::if_else(class %in% pos_substructure_classes,
                                "pos only (positive substructure)",
                                "pos only (sum composition)")
      )
    neg_unmatched <- neg_f %>%
      dplyr::mutate(.k = feat_key(`Metabolite name`, `Adduct type`, `Average Rt(min)`)) %>%
      dplyr::filter(!.k %in% neg_consumed) %>%
      dplyr::transmute(
        class, shorthand, pos_name = NA_character_, pos_rt = NA_real_,
        pos_adduct = NA_character_, neg_identity = `Metabolite name`,
        neg_rt = .rt, delta_rt = NA_real_,
        final_identity = `Metabolite name`, quant_mode = "negative",
        status = "neg only (neg-quantified)"
      )
    
    attr(out, "cross_mode_map") <-
      dplyr::bind_rows(matched_map, pos_unmatched, neg_unmatched) %>%
      dplyr::arrange(class, shorthand, dplyr::desc(status))
    out
  }
  
  append_duplicate_counters <- function(df) {
    feature_keys <- c("Average Rt(min)", "Average Mz", "Adduct type", "class", "ion.mode", "Metabolite name")
    feature_map <- df %>%
      dplyr::distinct(dplyr::across(all_of(feature_keys))) %>%
      dplyr::arrange(`Metabolite name`, `Average Rt(min)`, `Average Mz`, `Adduct type`, ion.mode, class) %>%
      dplyr::group_by(`Metabolite name`) %>%
      dplyr::mutate(
        dup_n = dplyr::n(),
        dup_idx = dplyr::row_number(),
        name_suffixed = dplyr::if_else(
          dup_n > 1 & !stringr::str_detect(`Metabolite name`, "\\[IS\\]"),
          paste0(`Metabolite name`, "_", dup_idx),
          `Metabolite name`
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(dplyr::all_of(feature_keys), name_suffixed)
    
    df %>%
      dplyr::left_join(feature_map, by = feature_keys) %>%
      dplyr::mutate(`Metabolite name` = dplyr::coalesce(name_suffixed, `Metabolite name`)) %>%
      dplyr::select(-name_suffixed)
  }
  
  add_plot_class_variant <- function(df) {
    df %>%
      dplyr::mutate(
        plot_class = dplyr::case_when(
          stringr::str_detect(`Metabolite name`, "^[A-Za-z]{2,3}\\s*-?O\\b") ~ paste0(class, "-O"),
          TRUE ~ class
        ),
        plot_class = tidyr::replace_na(plot_class, "Unknown")
      )
  }
  
  # ---- Manual outlier override state ----
  # Individual add/remove of specific samples or feature points, independent of
  # (and taking precedence over) the bulk "exclude all flagged" toggles.
  manual_sample_exclude  <- reactiveVal(character(0))
  manual_sample_keep     <- reactiveVal(character(0))
  manual_feature_exclude <- reactiveVal(character(0))
  manual_feature_keep    <- reactiveVal(character(0))
  
  .feature_key <- function(sample, metabolite) paste(sample, metabolite, sep = "\r__\r")
  
  observeEvent(input$sample_manual_exclude_btn, {
    sel <- input$sampleReviewTable_rows_selected
    req(length(sel) > 0)
    tbl <- isolate(sample_review_table())
    picked <- tbl$Sample[sel]
    manual_sample_exclude(union(manual_sample_exclude(), picked))
    manual_sample_keep(setdiff(manual_sample_keep(), picked))
  })
  
  observeEvent(input$sample_manual_keep_btn, {
    sel <- input$sampleReviewTable_rows_selected
    req(length(sel) > 0)
    tbl <- isolate(sample_review_table())
    picked <- tbl$Sample[sel]
    manual_sample_keep(union(manual_sample_keep(), picked))
    manual_sample_exclude(setdiff(manual_sample_exclude(), picked))
  })
  
  observeEvent(input$sample_manual_clear_btn, {
    manual_sample_exclude(character(0))
    manual_sample_keep(character(0))
  })
  
  observeEvent(input$feature_manual_exclude_btn, {
    sel <- input$outlierFeatureTable_rows_selected
    req(length(sel) > 0)
    tbl <- isolate(feature_outlier_flags())
    picked <- .feature_key(tbl$sample[sel], tbl[["Metabolite name"]][sel])
    manual_feature_exclude(union(manual_feature_exclude(), picked))
    manual_feature_keep(setdiff(manual_feature_keep(), picked))
  })
  
  observeEvent(input$feature_manual_keep_btn, {
    sel <- input$outlierFeatureTable_rows_selected
    req(length(sel) > 0)
    tbl <- isolate(feature_outlier_flags())
    picked <- .feature_key(tbl$sample[sel], tbl[["Metabolite name"]][sel])
    manual_feature_keep(union(manual_feature_keep(), picked))
    manual_feature_exclude(setdiff(manual_feature_exclude(), picked))
  })
  
  observeEvent(input$feature_manual_clear_btn, {
    manual_feature_exclude(character(0))
    manual_feature_keep(character(0))
  })
  
  bg_norm_pre_outlier <- reactive({
    req(bg_norm_long())
    
    # Cross-mode identity reconciliation (neg identity -> pos quant) runs first,
    # for its designated classes only. Those classes are then excluded from the
    # prefer-neg/prefer-pos lists so resolve_ion_modes() does not touch them again;
    # when the toggle is off they revert to their normal (prefer-negative) handling.
    cm_on      <- isTRUE(input$cross_mode_enable)
    cm_classes <- if (cm_on) cross_mode_classes() else character(0)
    cm_rt_tol  <- suppressWarnings(as.numeric(input$cross_mode_rt_tol %||% 0.1))
    if (length(cm_rt_tol) != 1 || is.na(cm_rt_tol) || cm_rt_tol < 0) cm_rt_tol <- 0.1
    
    df_src <- bg_norm_long()
    if (cm_on && length(cm_classes) > 0)
      df_src <- resolve_cross_mode(df_src, cm_classes, rt_tol = cm_rt_tol)
    
    df_res <- resolve_ion_modes(df_src,
                                pref_neg = setdiff(pref_neg_classes(), cm_classes),
                                pref_pos = setdiff(pref_pos_classes(), cm_classes))
    df_res <- append_duplicate_counters(df_res)
    df_res <- add_plot_class_variant(df_res)
    df_res <- df_res %>%
      dplyr::mutate(
        sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample)),
        is_iqc = is_iqc_sample(sample),
        # Preserve pre-protein-normalisation values so individual tabs (e.g. PCA)
        # can opt out of protein normalisation even when the global toggle is on —
        # useful since iQC samples are rarely present in the protein CSV and would
        # otherwise sit on a different scale to protein-normalised biological samples.
        norm_raw     = norm,
        value_bs_raw = value_bs
      )
    
    # ── Force Shiny to track protein norm toggle unconditionally ─────────────
    use_prot      <- isTRUE(input$use_protein_norm)
    prot_csv_file <- input$protein_csv            # register as dependency
    prot_units_in <- input$protein_units           # register as dependency
    prot_data     <- if (use_prot && !is.null(prot_csv_file))
      tryCatch(protein_df(), error = function(e) NULL)
    else NULL
    
    if (use_prot && !is.null(prot_data) && nrow(prot_data) > 0) {
      prot_join <- prot_data %>%
        dplyr::mutate(sample_norm = normalize_sample_name(sample)) %>%
        dplyr::select(sample_norm, protein)
      # Label appended to norm_units so downstream y-axis labels are accurate
      # (e.g. "pmol" -> "pmol/mg protein", "nmol" -> "nmol/ug protein").
      protein_unit_label <- switch(
        prot_units_in %||% "mg",
        "ug" = "ug protein",
        "mg" = "mg protein",
        "g"  = "g protein",
        "mg protein"
      )
      df_res <- df_res %>%
        dplyr::left_join(prot_join, by = "sample_norm") %>%
        dplyr::mutate(
          norm     = dplyr::if_else(!is.na(norm)     & !is.na(protein) & protein > 0,
                                    norm     / protein, norm),
          value_bs = dplyr::if_else(!is.na(value_bs) & !is.na(protein) & protein > 0,
                                    value_bs / protein, value_bs),
          protein_norm_applied = !is.na(protein),
          norm_units = dplyr::if_else(
            protein_norm_applied & !is.na(norm_units) & nzchar(norm_units),
            paste0(norm_units, "/", protein_unit_label),
            norm_units
          )
        ) %>%
        dplyr::select(-protein)
    }
    
    df_res
  })
  
  # ---- Sample-level outlier detection: PCA Hotelling's T² ----
  # Applied to non-blank samples (iQC included) on IS-normalised values, log-transformed
  # and unit-variance scaled, mirroring the PCA tab's default treatment but kept
  # self-contained so this tab doesn't depend on the PCA tab's own selections.
  sample_outlier_flags <- reactive({
    req(bg_norm_pre_outlier())
    
    empty_result <- tibble::tibble(
      sample = character(0), score = numeric(0), threshold = numeric(0),
      is_outlier = logical(0)
    )
    
    df0 <- bg_norm_pre_outlier() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"))
    df0 <- df0[!is_protected_sample(df0$sample), , drop = FALSE]
    if (nrow(df0) == 0) return(empty_result)
    
    mat <- df0 %>%
      dplyr::group_by(sample, `Metabolite name`) %>%
      dplyr::summarise(value = mean(norm, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = `Metabolite name`, values_from = value, values_fill = NA_real_) %>%
      as.data.frame()
    
    if (is.null(mat$sample) || nrow(mat) < 5) return(empty_result)
    
    rownames(mat) <- mat$sample
    mat$sample <- NULL
    mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
    
    keep <- vapply(mat, function(x) {
      sdv <- stats::sd(x, na.rm = TRUE)
      sum(!is.na(x)) >= 3 && is.finite(sdv) && sdv > 0
    }, logical(1))
    mat <- mat[, keep, drop = FALSE]
    if (ncol(mat) < 2 || nrow(mat) < 5) return(empty_result)
    
    mat <- log1p(pmax(mat, 0, na.rm = TRUE))
    for (j in seq_len(ncol(mat))) {
      v <- mat[[j]]
      if (anyNA(v)) v[is.na(v)] <- stats::median(v, na.rm = TRUE)
      mat[[j]] <- v
    }
    
    conf <- suppressWarnings(as.numeric(input$outlier_sample_conf %||% "0.975"))
    if (!is.finite(conf)) conf <- 0.975
    n <- nrow(mat)
    k <- min(5, ncol(mat) - 1, n - 2)
    if (k < 2) return(empty_result)
    
    pca <- tryCatch(stats::prcomp(mat, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (is.null(pca)) return(empty_result)
    
    scores <- pca$x[, seq_len(k), drop = FALSE]
    eig    <- (pca$sdev[seq_len(k)])^2
    eig[eig <= 0] <- .Machine$double.eps
    
    t2        <- rowSums(sweep(scores^2, 2, eig, "/"))
    f_crit    <- stats::qf(conf, k, n - k)
    t2_thresh <- (k * (n - 1) / (n - k)) * f_crit
    
    tibble::tibble(
      sample     = rownames(mat),
      score      = as.numeric(t2),
      threshold  = t2_thresh,
      is_outlier = t2 > t2_thresh
    )
  })
  
  # ---- Sample-level outlier detection: iQC replicate deviation ----
  # For each iQC sample, the median relative deviation from the cross-replicate
  # median across all features — a robust measure of "how far off is this run".
  iqc_outlier_flags <- reactive({
    req(bg_norm_pre_outlier())
    
    empty_result <- tibble::tibble(
      sample = character(0), score = numeric(0), threshold = numeric(0),
      is_outlier = logical(0)
    )
    
    df0 <- bg_norm_pre_outlier() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"), is_iqc_sample(sample))
    if (nrow(df0) == 0) return(empty_result)
    
    mat <- df0 %>%
      dplyr::group_by(sample, `Metabolite name`) %>%
      dplyr::summarise(value = mean(norm, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = `Metabolite name`, values_from = value, values_fill = NA_real_) %>%
      as.data.frame()
    
    if (is.null(mat$sample) || nrow(mat) < 3) return(empty_result)
    
    rownames(mat) <- mat$sample
    mat$sample <- NULL
    mat_m <- as.matrix(mat)
    mat_m[] <- suppressWarnings(as.numeric(mat_m))
    
    mult <- suppressWarnings(as.numeric(input$outlier_iqc_mad %||% "3"))
    if (!is.finite(mult)) mult <- 3
    
    med     <- apply(mat_m, 2, stats::median, na.rm = TRUE)
    rel_dev <- abs(sweep(mat_m, 2, med, "-"))
    rel_dev <- sweep(rel_dev, 2, abs(med), "/")
    rel_dev[!is.finite(rel_dev)] <- NA_real_
    
    score  <- apply(rel_dev, 1, stats::median, na.rm = TRUE)
    center <- stats::median(score, na.rm = TRUE)
    spread <- stats::mad(score, na.rm = TRUE)
    thresh <- center + mult * spread
    
    tibble::tibble(
      sample     = rownames(mat_m),
      score      = as.numeric(score),
      threshold  = thresh,
      is_outlier = score > thresh
    )
  })
  
  # ---- Combined sample review table for manual override ----
  # Lists every non-protected (ISTD/Blank/iQC excluded) sample with its PCA
  # Hotelling's T² score/flag and the current manual override state, so a
  # specific sample can be individually excluded or kept regardless of the
  # bulk "exclude all flagged" toggle.
  sample_review_table <- reactive({
    req(bg_norm_pre_outlier())
    
    base <- bg_norm_pre_outlier() %>%
      dplyr::distinct(sample) %>%
      dplyr::filter(!is_protected_sample(sample)) %>%
      dplyr::arrange(sample)
    
    validate(need(nrow(base) > 0, "No non-ISTD/Blank/iQC samples available."))
    
    flags <- sample_outlier_flags()
    
    base %>%
      dplyr::left_join(
        flags %>% dplyr::select(sample, t2_score = score, t2_flag = is_outlier),
        by = "sample"
      ) %>%
      dplyr::mutate(
        t2_flag          = dplyr::coalesce(t2_flag, FALSE),
        `Manual override` = dplyr::case_when(
          sample %in% manual_sample_keep()    ~ "Kept",
          sample %in% manual_sample_exclude() ~ "Excluded",
          TRUE ~ ""
        ),
        `Effective status` = dplyr::case_when(
          sample %in% manual_sample_keep()               ~ "Included",
          sample %in% manual_sample_exclude()             ~ "Excluded",
          isTRUE(input$outlier_exclude_samples) & t2_flag ~ "Excluded",
          TRUE ~ "Included"
        )
      ) %>%
      dplyr::transmute(
        Sample          = sample,
        `T² score`      = round(t2_score, 2),
        `Auto flagged`  = t2_flag,
        `Manual override`,
        `Effective status`
      ) %>%
      dplyr::arrange(`Effective status` == "Included", dplyr::desc(`Auto flagged`), Sample)
  })
  
  output$sampleReviewTable <- DT::renderDT({
    df <- sample_review_table()
    DT::datatable(
      df, rownames = FALSE, selection = "multiple",
      options = list(pageLength = 10), class = "stripe hover"
    ) %>%
      DT::formatStyle(
        "Effective status", target = "row",
        backgroundColor = DT::styleEqual(c("Excluded", "Included"), c("#fdecea", "white"))
      )
  })
  
  # ---- Feature-level outlier detection: per lipid × group ----
  # Flags individual sample values within a Metabolite name × group combination
  # (blanks and iQC excluded, groups from the active grouping scheme).
  feature_outlier_flags <- reactive({
    req(bg_norm_pre_outlier())
    
    empty_result <- tibble::tibble(
      sample = character(0), `Metabolite name` = character(0), class = character(0),
      group = character(0), value = numeric(0), score = numeric(0),
      method = character(0)
    )
    
    df0 <- bg_norm_pre_outlier() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"))
    df0 <- df0[!is_protected_sample(df0$sample), , drop = FALSE]
    if (nrow(df0) == 0) return(empty_result)
    
    fn        <- get_active_grouping()
    df0$group <- fn(df0$sample)$group
    df0       <- df0[!is.na(df0$group), , drop = FALSE]
    if (nrow(df0) == 0) return(empty_result)
    
    method <- input$outlier_feature_method %||% "mad"
    mult   <- suppressWarnings(as.numeric(input$outlier_feature_mult %||% "3"))
    if (!is.finite(mult)) mult <- 3
    
    if (identical(method, "iqr")) {
      res <- df0 %>%
        dplyr::group_by(`Metabolite name`, class, group) %>%
        dplyr::filter(sum(!is.na(norm)) >= 4) %>%
        dplyr::mutate(
          q1             = stats::quantile(norm, 0.25, na.rm = TRUE),
          q3             = stats::quantile(norm, 0.75, na.rm = TRUE),
          iqr_val        = q3 - q1,
          threshold_low  = q1 - mult * iqr_val,
          threshold_high = q3 + mult * iqr_val,
          is_outlier     = !is.na(norm) & (norm < threshold_low | norm > threshold_high),
          score          = dplyr::if_else(iqr_val > 0,
                                          pmax((q1 - norm) / iqr_val, (norm - q3) / iqr_val),
                                          NA_real_)
        ) %>%
        dplyr::ungroup()
    } else {
      res <- df0 %>%
        dplyr::group_by(`Metabolite name`, class, group) %>%
        dplyr::filter(sum(!is.na(norm)) >= 4) %>%
        dplyr::mutate(
          med_val    = stats::median(norm, na.rm = TRUE),
          mad_val    = stats::mad(norm, na.rm = TRUE),
          score      = dplyr::if_else(mad_val > 0, abs(norm - med_val) / mad_val, NA_real_),
          is_outlier = !is.na(score) & score > mult
        ) %>%
        dplyr::ungroup()
    }
    
    res %>%
      dplyr::filter(is_outlier) %>%
      dplyr::transmute(
        sample, `Metabolite name`, class, group,
        value  = norm,
        score,
        method = if (identical(method, "iqr")) "IQR rule" else "Modified Z-score (MAD)"
      ) %>%
      dplyr::arrange(dplyr::desc(score))
  })
  
  bg_norm_long_resolved <- reactive({
    df <- bg_norm_pre_outlier()
    
    # ── Sample-level exclusion ────────────────────────────────────────────────
    # Auto-flagged (PCA T²) samples if the bulk toggle is on, plus/minus any
    # individual manual overrides. iQC deviation is informational only and never
    # contributes to exclusion. ISTD/Blank/iQC samples are never excluded.
    auto_flagged_samples <- if (isTRUE(input$outlier_exclude_samples)) {
      sf <- sample_outlier_flags()
      sf$sample[sf$is_outlier]
    } else character(0)
    
    excluded_samples <- union(auto_flagged_samples, manual_sample_exclude())
    excluded_samples <- setdiff(excluded_samples, manual_sample_keep())
    excluded_samples <- excluded_samples[!is_protected_sample(excluded_samples)]
    
    if (length(excluded_samples) > 0) {
      df <- df %>% dplyr::filter(!sample %in% excluded_samples)
    }
    
    # ── Feature-level exclusion ───────────────────────────────────────────────
    auto_flagged_keys <- if (isTRUE(input$outlier_exclude_features)) {
      ff <- feature_outlier_flags()
      .feature_key(ff$sample, ff[["Metabolite name"]])
    } else character(0)
    
    excluded_keys <- union(auto_flagged_keys, manual_feature_exclude())
    excluded_keys <- setdiff(excluded_keys, manual_feature_keep())
    
    if (length(excluded_keys) > 0) {
      key_df <- .feature_key(df$sample, df[["Metabolite name"]])
      hit    <- key_df %in% excluded_keys
      if (any(hit)) {
        df$norm[hit]         <- NA_real_
        df$value_bs[hit]     <- NA_real_
        df$norm_raw[hit]     <- NA_real_
        df$value_bs_raw[hit] <- NA_real_
      }
    }
    
    df
  })
  
  # ---------- Technical replicate averaging ----------
  # Collapses multiple technical-injection rows down to one row per
  # (feature x biological sample) by averaging the numeric measurement
  # columns, when the "Average technical replicates" toggle is on. iQC,
  # ISTD, and Blank samples are always passed through untouched — each
  # injection is a distinct QC event and should never be averaged away.
  bg_norm_long_avg <- reactive({
    df <- bg_norm_long_resolved()
    
    if (!isTRUE(input$use_tech_rep_avg)) return(df)
    
    fn  <- get_active_tech_rep_map()
    bio <- fn(df$sample)
    
    protect <- is_protected_sample(df$sample)
    df$bio_sample <- dplyr::if_else(protect, df$sample, bio)
    
    # Nothing to collapse (no CSV column / no token selected) — skip the group-by
    if (all(df$bio_sample == df$sample)) {
      df$bio_sample <- NULL
      return(df)
    }
    
    id_cols <- intersect(
      c("Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type",
        "class", "ion.mode", "plot_class", "mode", "istd_name"),
      names(df)
    )
    avg_cols <- intersect(
      c("value", "blank_value", "value_bs", "IS_value", "norm",
        "norm_raw", "value_bs_raw"),
      names(df)
    )
    other_cols <- setdiff(names(df), c(id_cols, avg_cols, "sample", "sample_norm", "bio_sample"))
    
    safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    safe_first <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA else x[1] }
    
    df_avg <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(id_cols)), bio_sample) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(avg_cols), safe_mean),
        dplyr::across(dplyr::all_of(other_cols), safe_first),
        n_tech_reps = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::rename(sample = bio_sample) %>%
      dplyr::mutate(sample_norm = normalize_sample_name(sample))
    
    df_avg
  })
  
  # ---------- IS classes & plots ----------
  is_classes_available <- reactive({
    req(bg_norm_long())
    is_tbl <- attr(bg_norm_long(), "is_table")
    if (is.null(is_tbl) || nrow(is_tbl) == 0) return(character(0))
    is_tbl %>%
      dplyr::mutate(class = tidyr::replace_na(class, "Unknown")) %>%
      dplyr::distinct(class) %>%
      dplyr::arrange(class) %>%
      dplyr::pull(class)
  })
  
  observeEvent(bg_norm_long(), {
    classes <- is_classes_available()
    updateSelectInput(
      session, "is_class",
      choices  = if (length(classes) > 0) classes else "No IS found",
      selected = if (length(classes) > 0) classes[1] else "No IS found"
    )
  }, ignoreInit = FALSE)
  
  is_plot_data <- reactive({
    req(bg_norm_long(), input$is_class, input$is_mode)
    is_table <- attr(bg_norm_long(), "is_table")
    validate(need(!is.null(is_table) && nrow(is_table) > 0, "No internal standards data available."))
    
    df <- is_table %>%
      dplyr::mutate(
        class = tidyr::replace_na(class, "Unknown"),
        ion.mode = dplyr::coalesce(ion.mode, "unknown"),
        sample_norm = coalesce(sample_norm, normalize_sample_name(sample)),
        istd_name = dplyr::coalesce(istd_name, "unknown IS")
      ) %>%
      dplyr::filter(class == input$is_class)
    
    
    # Single, consistent blank filter
    df <- filter_blanks(df, input$exclude_blank_is, sample_col = "sample_norm", exact = FALSE)
    
    # Filter out iQC
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    
    validate(need(nrow(df) > 0, "No IS values available for the selected class."))
    
    
    df <- df %>%
      dplyr::mutate(group = resolve_group_labels(
        samples = sample, method = input$grouping_method_is,
        delim = input$group_delim_is, tokens = input$group_tokens_is,
        regex = input$group_regex_is, group_csv_input = input$group_csv,
        use_csv_flag = input$use_group_csv,
        group_map = if (!is.null(input$group_csv) && isTRUE(input$use_group_csv))
          group_map() else NULL
      ))
    
    mode_sel <- input$is_mode
    if (mode_sel == "neg") {
      df <- df %>% dplyr::filter(ion.mode == "negative")
    } else if (mode_sel == "pos") {
      df <- df %>% dplyr::filter(ion.mode == "positive")
    } else if (mode_sel == "auto") {
      if (any(df$class == "CL", na.rm = TRUE)) {
        df <- df %>%
          dplyr::group_by(group, sample) %>%
          dplyr::mutate(
            has_neg = any(ion.mode == "negative", na.rm = TRUE),
            has_pos = any(ion.mode == "positive", na.rm = TRUE)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::filter(!(has_neg & ion.mode == "positive")) %>%
          dplyr::select(-has_neg, -has_pos)
      }
    }
    
    if (mode_sel == "both") {
      df_sum <- df %>%
        dplyr::group_by(group, ion.mode) %>%
        dplyr::summarise(
          mean = mean(IS_value, na.rm = TRUE),
          sd   = sd(IS_value, na.rm = TRUE),
          n    = dplyr::n(),
          sem  = sd / sqrt(pmax(n, 1)),
          ci95 = 1.96 * sem,
          .groups = "drop"
        )
    } else {
      df_sum <- df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(
          mean = mean(IS_value, na.rm = TRUE),
          sd   = sd(IS_value, na.rm = TRUE),
          n    = dplyr::n(),
          sem  = sd / sqrt(pmax(n, 1)),
          ci95 = 1.96 * sem,
          .groups = "drop"
        )
    }
    
    list(points = df, summary = df_sum, mode = mode_sel)
  })
  
  # ---------- Classes for plots and export ----------
  met_classes_available <- reactive({
    req(bg_norm_long_avg())
    bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::distinct(plot_class) %>%
      dplyr::arrange(plot_class) %>%
      dplyr::pull(plot_class)
  })
  
  observeEvent(bg_norm_long_avg(), {
    updateSelectInput(
      session, "met_class",
      choices  = met_classes_available(),
      selected = if (length(met_classes_available()) > 0) met_classes_available()[1] else character(0)
    )
    classes <- met_classes_available()
    updateCheckboxGroupInput(
      session,
      "export_classes",
      choices  = classes,
      selected = classes
    )
  }, ignoreInit = FALSE)
  
  # Shared class + ion-mode slice for the single-lipid tab (drives both the adduct
  # and name selectors so they always agree with what will actually plot).
  met_class_df <- reactive({
    req(bg_norm_long_avg(), input$met_class)
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"),
                    plot_class == input$met_class)
    mode_sel <- input$met_ion_mode %||% "all"
    if (!identical(mode_sel, "all") && "ion.mode" %in% names(df))
      df <- df %>% dplyr::filter(ion.mode == mode_sel)
    df
  })
  
  met_names_available <- reactive({
    df <- met_class_df()
    adduct_sel <- input$met_adduct_filter %||% "all"
    if (!is.null(adduct_sel) && adduct_sel != "all" && nzchar(adduct_sel))
      df <- df %>% dplyr::filter(`Adduct type` == adduct_sel)
    df %>%
      dplyr::distinct(`Metabolite name`) %>%
      dplyr::arrange(`Metabolite name`) %>%
      dplyr::pull(`Metabolite name`)
  })
  
  # Name list cascades off class + ion mode + adduct; preserve the current
  # selection when it is still valid, otherwise fall back to the class total.
  observeEvent(list(input$met_class, input$met_ion_mode, input$met_adduct_filter), {
    nm <- tryCatch(met_names_available(), error = function(e) character(0))
    choices <- c("Total (class sum)", nm)
    cur <- isolate(input$met_name)
    sel <- if (!is.null(cur) && cur %in% choices) cur else "Total (class sum)"
    updateSelectInput(session, "met_name",
                      choices  = if (length(choices) > 0) choices else "No metabolites in this class",
                      selected = if (length(choices) > 0) sel else "No metabolites in this class")
  }, ignoreInit = FALSE)
  
  # ── Adduct filter reactives ────────────────────────────────────────────────────
  # Returns adducts available for the current class + ion mode
  met_adducts_available <- reactive({
    adducts <- met_class_df() %>%
      dplyr::distinct(`Adduct type`) %>%
      dplyr::pull(`Adduct type`)
    sort(unique(adducts[!is.na(adducts) & nzchar(adducts)]))
  })
  
  class_all_adducts_available <- reactive({
    req(bg_norm_long_avg(), input$class_all)
    adducts <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"),
                    plot_class == input$class_all) %>%
      dplyr::distinct(`Adduct type`) %>%
      dplyr::pull(`Adduct type`)
    sort(unique(adducts[!is.na(adducts) & nzchar(adducts)]))
  })
  
  # Update adduct selector when class or ion mode changes (reset to All)
  observeEvent(list(input$met_class, input$met_ion_mode), {
    adducts <- tryCatch(met_adducts_available(), error = function(e) character(0))
    choices  <- c("All" = "all", setNames(adducts, adducts))
    updateSelectInput(session, "met_adduct_filter",
                      choices = choices, selected = "all")
  }, ignoreInit = FALSE)
  
  observeEvent(input$class_all, {
    adducts <- tryCatch(class_all_adducts_available(), error = function(e) character(0))
    choices  <- c("All" = "all", setNames(adducts, adducts))
    updateSelectInput(session, "class_all_adduct_filter",
                      choices = choices, selected = "all")
  }, ignoreInit = FALSE)
  
  met_plot_data <- reactive({
    req(
      bg_norm_long_avg(),
      input$met_class, input$met_name,
      input$plot_value_type_met, input$display_mode_met
    )
    
    # --- Base DF: resolve blanks & select class, exclude IS
    df <- bg_norm_long_avg() %>%
      dplyr::mutate(
        sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample))
      )
    
    # Single, consistent blank filter
    df <- filter_blanks(
      df,
      exclude = input$exclude_blank_met,
      sample_col = "sample_norm",
      exact = FALSE
    )
    
    #filter out iQc
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    
    df <- df %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::filter(plot_class == input$met_class)
    
    # --- Ion-mode filter (cascades with the name/adduct selectors)
    mode_sel <- input$met_ion_mode %||% "all"
    if (!identical(mode_sel, "all") && "ion.mode" %in% names(df))
      df <- df %>% dplyr::filter(ion.mode == mode_sel)
    
    # --- Adduct filter (if not "All")
    adduct_sel <- input$met_adduct_filter
    if (!is.null(adduct_sel) && adduct_sel != "all" && nzchar(adduct_sel))
      df <- df %>% dplyr::filter(`Adduct type` == adduct_sel)
    measure_col <- pick_measure_col(input$plot_value_type_met)
    # Ensure norm_units exists (avoids "unknown column" warnings when a filtered
    # slice or a non-normalised value type has no units column)
    if (!"norm_units" %in% names(df)) df$norm_units <- NA_character_
    
    # --- Class totals per sample — used for percent mode
    class_totals <- df %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        class_total = sum(.data[[measure_col]], na.rm = TRUE),
        .groups = "drop"
      )
    
    # --- Build per-metabolite data (or class total)
    if (identical(input$met_name, "Total (class sum)")) {
      # Sum across all metabolites in the class per sample
      df_met <- df %>%
        dplyr::group_by(sample, plot_class) %>%
        dplyr::summarise(
          value = sum(.data[[measure_col]], na.rm = TRUE),
          # carry representative normalized units if applicable
          norm_units = dplyr::first(stats::na.omit(norm_units)),
          .groups = "drop"
        )
      
      if (identical(input$display_mode_met, "percent")) {
        df_met <- df_met %>%
          dplyr::left_join(class_totals, by = "sample") %>%
          dplyr::mutate(
            value = dplyr::if_else(class_total > 0, 100, NA_real_)
          ) %>%
          dplyr::select(-class_total)
      }
      
      df_met <- df_met %>%
        dplyr::mutate(`Metabolite name` = "Total (class sum)")
      
    } else {
      # Single metabolite: keep norm_units when plotting normalized values
      df_met <- df %>%
        dplyr::filter(`Metabolite name` == input$met_name) %>%
        dplyr::transmute(
          sample, plot_class, `Metabolite name`,
          value = .data[[measure_col]],
          norm_units = if (measure_col == "norm") norm_units else NA_character_
        ) %>%
        dplyr::left_join(class_totals, by = "sample")
      
      if (identical(input$display_mode_met, "percent")) {
        df_met <- df_met %>%
          dplyr::mutate(
            value = dplyr::if_else(class_total > 0, 100 * value / class_total, NA_real_)
          ) %>%
          dplyr::select(-class_total)
      } else {
        df_met <- df_met %>% dplyr::select(-class_total)
      }
    }
    
    # --- Validate we have point data
    validate(need(nrow(df_met) > 0, "No data available for the selected option."))
    
    # --- Assign group labels via central grouping reactive
    fn  <- get_active_grouping()
    res <- fn(df_met$sample)
    f_met <- df_met %>% dplyr::mutate(group = res$group)
    
    # --- Filter to selected groups
    sel_grps <- input$met_selected_groups
    if (!is.null(sel_grps) && length(sel_grps) > 0)
      f_met <- f_met %>% dplyr::filter(.data$group %in% sel_grps)
    
    # Extra sanity checks (surface issues early)
    validate(need("group" %in% names(f_met), "Grouping failed: no 'group' column on points data."))
    validate(need(nrow(f_met) > 0, "No point data after grouping."))
    
    # --- Group-level summary stats (MUST use f_met which contains 'group')
    df_sum <- f_met %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean = mean(value, na.rm = TRUE),
        sd   = stats::sd(value, na.rm = TRUE),
        n    = dplyr::n(),
        sem  = sd / sqrt(pmax(n, 1)),
        ci95 = 1.96 * sem,
        .groups = "drop"
      )
    
    # --- Unit label for y-axis (only relevant for normalized, non-percent)
    units_label <- NULL
    if (identical(input$plot_value_type_met, "norm") &&
        !identical(input$display_mode_met, "percent")) {
      u <- stats::na.omit(unique(f_met$norm_units))
      if (length(u) == 1) {
        units_label <- u
      } else if (length(u) > 1) {
        showNotification(
          "Metabolite plot: multiple ISTD units detected across data (mixed units).",
          type = "warning", duration = 6
        )
        units_label <- "(mixed units)"
      }
      # if length(u) == 0 → leave NULL to fall back to default in axis_label()
    }
    
    # --- Return points, summary, and computed y-axis label
    list(
      points  = f_met,  # ✅ grouped points data
      summary = df_sum, # ✅ summary by group
      ylab    = axis_label(
        kind       = input$plot_value_type_met,
        mode       = input$display_mode_met,
        unit_label = units_label
      )
    )
  })
  
  
  # ---------- Class (all metabolites) ----------
  observeEvent(bg_norm_long_avg(), {
    updateSelectInput(session, "class_all",
                      choices = met_classes_available(),
                      selected = if (length(met_classes_available()) > 0) met_classes_available()[1] else character(0))
  })
  
  all_groups_available <- reactive({
    tbl <- tryCatch(lr_token_preview(), error = function(e) NULL)
    if (is.null(tbl) || nrow(tbl) == 0) return(character(0))
    fn   <- get_active_grouping()
    res  <- fn(tbl$sample)
    sort(unique(na.omit(res$group[nzchar(res$group)])))
  })
  
  # Update class_all_selected_group when groups change
  observeEvent(all_groups_available(), {
    choices <- all_groups_available()
    updateSelectInput(session, "class_all_selected_group",
                      choices  = if (length(choices) > 0) choices else "No groups found",
                      selected = if (length(choices) > 0) choices[1] else "No groups found")
  }, ignoreInit = FALSE)
  
  class_all_plot_data <- reactive({
    req(bg_norm_long_avg(), input$class_all,
        input$plot_value_type_all, input$display_mode_all)
    
    withProgress(message = "Building class plot…", value = 0, {
      incProgress(0.15, detail = "Filtering data")
      df <- bg_norm_long_avg() %>%
        dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
        dplyr::filter(plot_class == input$class_all) %>%
        dplyr::mutate(sample_norm = coalesce(sample_norm, normalize_sample_name(sample)))
      
      # Adduct filter
      adduct_sel <- input$class_all_adduct_filter
      if (!is.null(adduct_sel) && adduct_sel != "all" && nzchar(adduct_sel))
        df <- df %>% dplyr::filter(`Adduct type` == adduct_sel)
      
      # consistent blank filter
      df <- filter_blanks(df, input$exclude_blank_all, sample_col = "sample_norm", exact = FALSE)
      
      #Filter out iQC
      df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
      
      incProgress(0.25, detail = "Grouping & totals")
      fn  <- get_active_grouping()
      res <- fn(df$sample)
      df  <- df %>% dplyr::mutate(group = res$group)
      
      # Filter to selected group (single group for this plot)
      sel_grp <- input$class_all_selected_group
      if (!is.null(sel_grp) && nzchar(sel_grp))
        df <- df %>% dplyr::filter(group == sel_grp)
      
      validate(need(nrow(df) > 0, "No data available for the selected group."))
      
      measure_col <- pick_measure_col(input$plot_value_type_all)
      class_totals <- df %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(class_total = sum(.data[[measure_col]], na.rm = TRUE), .groups = "drop")
      
      incProgress(0.30, detail = "Summarising means & errors")
      df_pts <- df %>%
        dplyr::transmute(sample, group, plot_class, `Metabolite name`, value = .data[[measure_col]]) %>%
        dplyr::left_join(class_totals, by = "sample")
      
      if (identical(input$display_mode_all, "percent")) {
        df_pts <- df_pts %>%
          dplyr::mutate(value = dplyr::if_else(class_total > 0, 100 * value / class_total, NA_real_))
      }
      df_pts <- df_pts %>% dplyr::select(-class_total)
      
      df_sum <- df_pts %>%
        dplyr::group_by(`Metabolite name`) %>%
        dplyr::summarise(
          mean = mean(value, na.rm = TRUE),
          sd   = sd(value, na.rm = TRUE),
          n    = dplyr::n(),
          sem  = sd / sqrt(pmax(n, 1)),
          ci95 = 1.96 * sem,
          .groups = "drop"
        )
      
      # After computing df_pts and df_sum ...
      
      units_label <- {
        if (identical(input$plot_value_type_all, "norm")) {
          u <- na.omit(unique(df_pts$norm_units))
          if (length(u) == 1) u else {
            if (length(u) > 1) {
              showNotification("Class-all plot: multiple ISTD units detected across data (mixed units).", type = "warning", duration = 6)
            }
            "(mixed units)"
          }
        } else NULL
      }
      
      list(
        points = df_pts,
        summary = df_sum,
        ylab = axis_label(input$plot_value_type_all, input$display_mode_all, unit_label = units_label)
      )
      
      order_mode <- input$class_all_order
      if (identical(order_mode, "alpha_asc")) {
        df_sum <- df_sum %>% dplyr::arrange(`Metabolite name`)
      } else {
        df_sum <- df_sum %>% dplyr::arrange(dplyr::desc(mean), `Metabolite name`)
      }
      
      incProgress(0.30, detail = "Finalising")
      list(points = df_pts, summary = df_sum,
           ylab = axis_label(input$plot_value_type_all, input$display_mode_all))
    })
  })
  
  # ---------- EXPORT ----------
  export_wide_data <- reactive({
    req(bg_norm_long_avg())
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(sample_norm = coalesce(sample_norm, normalize_sample_name(sample)))
    
    # Single, consistent blank filter
    df <- filter_blanks(df, input$exclude_blank_export, sample_col = "sample_norm", exact = FALSE)
    
    # Filter out iQC
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    
    classes <- input$export_classes
    if (!is.null(classes) && length(classes) > 0) {
      df <- df %>% dplyr::filter(plot_class %in% classes)
    }
    
    fn  <- get_active_grouping()
    res <- fn(df$sample)
    df  <- df %>% dplyr::mutate(group = res$group)
    
    measure_col <- pick_measure_col(input$export_value_type)
    
    if (identical(input$export_display_mode, "percent")) {
      totals <- df %>%
        dplyr::group_by(sample, plot_class) %>%
        dplyr::summarise(class_total = sum(.data[[measure_col]], na.rm = TRUE), .groups = "drop")
      df <- df %>%
        dplyr::left_join(totals, by = c("sample", "plot_class")) %>%
        dplyr::mutate(export_value = dplyr::if_else(class_total > 0, 100 * .data[[measure_col]] / class_total, NA_real_)) %>%
        dplyr::select(-class_total)
    } else {
      df <- df %>% dplyr::mutate(export_value = .data[[measure_col]])
    }
    
    meta_cols <- c("plot_class", "Metabolite name", "Average Rt(min)", "Average Mz", "Adduct type", "ion.mode")
    
    wide <- df %>%
      dplyr::select(dplyr::all_of(meta_cols), sample, export_value) %>%
      tidyr::pivot_wider(
        names_from = sample,
        values_from = export_value
      ) %>%
      dplyr::arrange(plot_class, `Metabolite name`)
    
    wide
  })
  
  observeEvent(bg_norm_long_avg(), {
    classes <- met_classes_available()
    updateCheckboxGroupInput(
      session,
      "export_classes",
      choices  = classes,
      selected = classes
    )
  })
  
  observeEvent(input$toggle_classes, {
    all_classes <- met_classes_available()
    current <- isolate(input$export_classes)
    if (is.null(current)) current <- character(0)
    new_selection <- if (length(current) < length(all_classes)) all_classes else character(0)
    updateCheckboxGroupInput(session, "export_classes", selected = new_selection)
  })
  
  # ---- PCA ----
  # iQC samples are handled as FactoMineR "supplementary individuals": the PCA
  # axes/loadings are computed from biological samples only, then iQC are
  # projected into that space afterward without influencing it. Since iQC are
  # rarely present in the protein CSV (so bg_norm_pre_outlier() leaves their
  # norm/value_bs un-normalised), they're scaled here by the *median* protein
  # content of the biological samples — an assumed, not measured, value —
  # purely so the projection lands somewhere visually meaningful.
  pca_data <- reactive({
    req(bg_norm_long_avg())
    measure   <- req(input$pca_measure)
    units     <- req(input$pca_units)
    excl      <- isTRUE(input$exclude_blank_pca)
    show_iqc  <- isTRUE(input$pca_show_iqc)
    show_istd <- isTRUE(input$pca_show_istd)
    
    df0 <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"))
    
    # ── Group selection filter ────────────────────────────────────────────────
    # iQC and ISTD samples are distinct: iQC is a pooled biological QC sample;
    # a sample named ISTD/ITSD is typically an internal-standard-only injection
    # with no biological matrix, so it must never be treated as equivalent to
    # iQC (mixing them was producing a spurious extra/dominant point). Both are
    # excluded from group filtering (neither has a real group label) and, if
    # shown at all, are added back as SUPPLEMENTARY points — never active.
    sel_grps   <- input$pca_selected_groups
    df0_iqc    <- df0 %>% dplyr::filter(is_iqc_sample(sample))
    df0_istd   <- df0 %>% dplyr::filter(is_istd_sample(sample))
    df0_noiqc  <- df0 %>% dplyr::filter(!is_qc_type_sample(sample))
    
    if (!is.null(sel_grps) && length(sel_grps) > 0) {
      fn      <- get_active_grouping()
      grp_vec <- fn(df0_noiqc$sample)$group
      df0_noiqc <- df0_noiqc[!is.na(grp_vec) & grp_vec %in% sel_grps, , drop = FALSE]
    }
    
    # Scale both iQC and ISTD by the median biological protein content — neither
    # is typically present in the protein CSV — purely so their projection lands
    # somewhere comparable; never affects the fit itself.
    if (isTRUE(input$use_protein_norm) && (nrow(df0_iqc) > 0 || nrow(df0_istd) > 0)) {
      prot <- tryCatch(protein_df(), error = function(e) NULL)
      if (!is.null(prot) && nrow(prot) > 0) {
        prot_median <- stats::median(prot$protein, na.rm = TRUE)
        if (is.finite(prot_median) && prot_median > 0) {
          .scale_supp <- function(d) {
            d %>% dplyr::mutate(
              norm     = dplyr::if_else(!is.na(norm),     norm     / prot_median, norm),
              value_bs = dplyr::if_else(!is.na(value_bs), value_bs / prot_median, value_bs)
            )
          }
          if (nrow(df0_iqc)  > 0) df0_iqc  <- .scale_supp(df0_iqc)
          if (nrow(df0_istd) > 0) df0_istd <- .scale_supp(df0_istd)
        }
      }
    }
    
    # Recombine: filtered biological samples always included; iQC/ISTD added as
    # supplementary rows only if their respective toggle is on
    df0 <- df0_noiqc
    if (show_iqc)  df0 <- dplyr::bind_rows(df0, df0_iqc)
    if (show_istd) df0 <- dplyr::bind_rows(df0, df0_istd)
    
    if (identical(units, "percent")) {
      totals <- df0 %>%
        dplyr::group_by(sample, plot_class) %>%
        dplyr::summarise(class_total = sum(.data[[measure]], na.rm = TRUE), .groups = "drop")
      
      df0 <- df0 %>%
        dplyr::left_join(totals, by = c("sample", "plot_class")) %>%
        dplyr::mutate(value_for_pca = dplyr::if_else(class_total > 0, 100 * .data[[measure]] / class_total, NA_real_)) %>%
        dplyr::select(sample, `Metabolite name`, value_for_pca)
    } else {
      df0 <- df0 %>% dplyr::transmute(sample, `Metabolite name`, value_for_pca = .data[[measure]])
    }
    
    df <- df0 %>%
      dplyr::group_by(sample, `Metabolite name`) %>%
      dplyr::summarise(value = mean(value_for_pca, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from  = `Metabolite name`,
        values_from = value,
        values_fill = NA_real_
      ) %>%
      as.data.frame()
    
    if (!is.null(df$sample)) {
      rownames(df) <- df$sample
      df$sample <- NULL
    }
    
    if (excl && nrow(df) > 0) {
      df <- filter_blank_rownames(df, exclude = TRUE, exact = FALSE)
    }
    
    validate(need(ncol(df) > 0, "No metabolite columns available for PCA."))
    
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
    if (ncol(df) > 0) {
      all_na_cols <- colSums(!is.na(df)) == 0
      if (any(all_na_cols)) df <- df[, !all_na_cols, drop = FALSE]
    }
    if (ncol(df) > 0) {
      zero_var_cols <- vapply(df, function(x) {
        s <- stats::sd(x, na.rm = TRUE)
        isTRUE(is.nan(s)) || isTRUE(s == 0)
      }, logical(1))
      if (any(zero_var_cols)) df <- df[, !zero_var_cols, drop = FALSE]
    }
    
    validate(
      need(sum(!is_qc_type_sample(rownames(df))) >= 3, "Not enough biological samples for PCA (need \u2265 3)."),
      need(ncol(df) >= 2, "Not enough variables for PCA (need \u2265 2 with variance).")
    )
    
    if (anyNA(df)) {
      for (j in seq_len(ncol(df))) {
        v <- df[[j]]
        if (anyNA(v) && !all(is.na(v))) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          df[[j]] <- v
        }
      }
    }
    
    # FactoMineR's ind.sup mechanism requires supplementary rows to be
    # identifiable by index — put biological (active) rows first, iQC/ISTD
    # (supplementary) rows last, and record the split as an attribute.
    is_supp     <- is_qc_type_sample(rownames(df))
    ord         <- order(is_supp)
    df          <- df[ord, , drop = FALSE]
    attr(df, "is_supp") <- is_supp[ord]
    
    df
  })
  
  pca_groups <- reactive({
    req(pca_data())
    samples <- rownames(pca_data())
    fn  <- get_active_grouping()
    res <- fn(samples)
    factor(res$group)
  })
  
  pca_result <- reactive({
    req(pca_data())
    withProgress(message = "Running PCA…", value = 0, {
      incProgress(0.20, detail = "Scaling features")
      
      mat     <- pca_data()
      is_supp <- attr(mat, "is_supp")
      scaling <- input$pca_scaling %||% "uv"
      
      # Impute any remaining NA/zero with half-min per feature
      mat_imp <- apply(mat, 2, function(x) {
        x[!is.finite(x) | x <= 0] <- NA_real_
        hm <- min(x, na.rm = TRUE) / 2
        x[is.na(x)] <- hm
        x
      })
      
      # Apply scaling — FactoMineR will mean-centre internally; we pre-scale here
      scaled_mat <- switch(scaling,
                           uv = {
                             # Standard unit-variance: scale.unit=TRUE in FactoMineR handles this
                             mat_imp
                           },
                           pareto = {
                             # Divide each feature by square root of its SD (after centering)
                             sds <- apply(mat_imp, 2, sd, na.rm = TRUE)
                             sds[sds == 0] <- 1
                             sweep(mat_imp, 2, sqrt(sds), "/")
                           },
                           log_uv = {
                             # Log transform (log1p to handle zeros), then unit variance
                             log1p(mat_imp)
                           },
                           none = {
                             mat_imp
                           },
                           mat_imp
      )
      
      # For UV and log+UV, tell FactoMineR to scale to unit variance
      do_scale <- scaling %in% c("uv", "log_uv")
      
      incProgress(0.60, detail = "FactoMineR::PCA")
      n_active <- sum(!is_supp)
      ind_sup  <- if (any(is_supp)) which(is_supp) else NULL
      ncp      <- max(2, min(10, ncol(scaled_mat), n_active - 1))
      res <- FactoMineR::PCA(scaled_mat, graph = FALSE,
                             scale.unit = do_scale, ncp = ncp,
                             ind.sup = ind_sup)
      
      incProgress(0.20, detail = "Preparing loadings")
      res
    })
  })
  
  
  
  # ── Shared plot reactives for bulk export ────────────────────────────────────
  # These build the ggplot/pheatmap objects reused by both renderPlot and download_stats_pdf
  
  # Shared palette helper — returns a 101-colour vector for pheatmap
  .heatmap_palette <- function(palette_id, n = 101) {
    switch(palette_id,
           "viridis"  = viridis::viridis(n,  option = "viridis"),
           "magma"    = viridis::viridis(n,  option = "magma"),
           "plasma"   = viridis::viridis(n,  option = "plasma"),
           "inferno"  = viridis::viridis(n,  option = "inferno"),
           "cividis"  = viridis::viridis(n,  option = "cividis"),
           "rocket"   = viridis::viridis(n,  option = "rocket"),
           "mako"     = viridis::viridis(n,  option = "mako"),
           "turbo"    = viridis::viridis(n,  option = "turbo"),
           "bwr"      = colorRampPalette(c("#2166ac", "#f7f7f7", "#d6604d"))(n),
           "gwp"      = colorRampPalette(c("#1b7837", "#f7f7f7", "#762a83"))(n),
           viridis::viridis(n, option = "viridis")  # default fallback
    )
  }
  
  heatmap_plot_obj <- reactive({
    hd <- tryCatch(sig_heatmap_data(), error = function(e) NULL)
    if (is.null(hd)) return(NULL)
    cols <- .heatmap_palette(input$hm_palette %||% "viridis")
    pheatmap::pheatmap(
      hd$mat,
      color            = cols,
      cluster_rows     = isTRUE(input$hm_cluster_rows),
      cluster_cols     = isTRUE(input$hm_cluster_cols),
      clustering_distance_rows = input$hm_dist %||% "correlation",
      clustering_distance_cols = input$hm_dist %||% "correlation",
      clustering_method = input$hm_linkage %||% "complete",
      annotation_col   = if (!is.null(hd$anno_col) && ncol(hd$anno_col) > 0) hd$anno_col else NULL,
      show_rownames    = isTRUE(input$hm_show_rownames),
      show_colnames    = isTRUE(input$hm_show_colnames),
      fontsize_row     = 8,
      fontsize_col     = 9,
      border_color     = NA,
      main             = paste0("Significant features heatmap — ", hd$class),
      silent           = TRUE
    )
  })
  
  volcano_plot_obj <- reactive({
    tryCatch({
      req(stats_results())
      res <- stats_results()
      if (is.null(res) || nrow(res) == 0) return(NULL)
      # Build simple volcano from stats_results
      df_v <- res %>%
        dplyr::filter(!is.na(p_adj) & !is.na(log2FC)) %>%
        dplyr::mutate(sig = p_adj < (input$alpha %||% 0.05))
      ggplot2::ggplot(df_v, ggplot2::aes(x = log2FC, y = -log10(p_adj),
                                         colour = sig, label = `Metabolite name`)) +
        ggplot2::geom_point(alpha = 0.7, size = 2) +
        ggplot2::scale_colour_manual(values = c("TRUE" = "#d7191c", "FALSE" = "#aaaaaa")) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(x = "log2 fold change", y = "-log10(adj p-value)",
                      title = paste0("Volcano — ", input$stats_class))
    }, error = function(e) NULL)
  })
  
  enrichment_plot_obj <- reactive({
    tryCatch({
      req(enrichment_results())
      res <- enrichment_results()
      if (is.null(res) || nrow(res) == 0) return(NULL)
      top <- res %>% dplyr::slice_min(padj, n = 20, with_ties = FALSE)
      ggplot2::ggplot(top,
                      ggplot2::aes(x = NES %||% ES, y = stats::reorder(pathway, NES %||% ES),
                                   fill = padj)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_viridis_c(direction = -1) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::labs(x = "Enrichment score", y = NULL, fill = "adj p",
                      title = "Lipid set enrichment (top 20)")
    }, error = function(e) NULL)
  })
  
  pathway_heatmap_obj <- reactive({
    tryCatch({
      req(path_scores_long())
      sc <- path_scores_long()
      if (is.null(sc) || nrow(sc) == 0) return(NULL)
      mat <- sc %>%
        dplyr::select(sample, score_name, score) %>%
        tidyr::pivot_wider(names_from = sample, values_from = score) %>%
        tibble::column_to_rownames("score_name") %>%
        as.matrix()
      mat <- mat[, !apply(mat, 2, function(x) all(is.na(x))), drop = FALSE]
      pheatmap::pheatmap(
        mat,
        cluster_rows = TRUE, cluster_cols = TRUE,
        show_colnames = isTRUE(input$hm_show_colnames), show_rownames = isTRUE(input$hm_show_rownames),
        fontsize_row = 8, fontsize_col = 8,
        border_color = NA,
        main = "Synthesis pathway scores",
        silent = TRUE
      )
    }, error = function(e) NULL)
  })
  
  is_plot_obj <- reactive({
    tryCatch({
      req(is_plot_data())
      pd       <- is_plot_data()
      df_pts   <- pd$points
      df_sum   <- pd$summary
      mode_sel <- pd$mode
      err_vec  <- df_sum$sem
      df_sum$ymin <- df_sum$mean - err_vec
      df_sum$ymax <- df_sum$mean + err_vec
      if (mode_sel == "both") {
        ggplot() +
          geom_col(data = df_sum, aes(x = group, y = mean, fill = ion.mode),
                   position = position_dodge(0.6), width = 0.6) +
          geom_errorbar(data = df_sum,
                        aes(x = group, ymin = ymin, ymax = ymax, colour = ion.mode),
                        position = position_dodge(0.6), width = 0.2) +
          theme_minimal(base_size = 12) +
          labs(title = paste0("IS Values — ", input$is_class), x = "Group", y = "IS Value")
      } else {
        ggplot() +
          geom_col(data = df_sum, aes(x = group, y = mean), fill = "#72B7B2", width = 0.7) +
          geom_errorbar(data = df_sum, aes(x = group, ymin = ymin, ymax = ymax),
                        width = 0.2, colour = "#333") +
          theme_minimal(base_size = 12) +
          labs(title = paste0("IS Values — ", input$is_class), x = "Group", y = "IS Value")
      }
    }, error = function(e) NULL)
  })
  
  # Shared PCA plot reactive (used by renderPlot and bulk export)
  # Manual ggplot2 build (rather than factoextra::fviz_pca_ind) so we have full
  # control over the legend — specifically, so "iQC" can appear as a real,
  # colour-coded legend entry even though it's a supplementary group — and so
  # sample-ID labels and font size can be toggled cleanly for both the on-screen
  # plot and PNG/SVG export.
  .build_pca_plot <- function(fsz = 14, show_labels = FALSE) {
    req(pca_result(), pca_data())
    
    res      <- pca_result()
    mat_pca  <- pca_data()
    is_supp  <- attr(mat_pca, "is_supp")
    if (is.null(is_supp)) is_supp <- rep(FALSE, nrow(mat_pca))
    sample_names <- rownames(mat_pca)
    
    orig_grp <- pca_groups()
    if (length(orig_grp) != length(sample_names)) orig_grp <- rep("Other", length(sample_names))
    
    active_coord <- as.data.frame(res$ind$coord[, 1:2, drop = FALSE])
    colnames(active_coord) <- c("Dim1", "Dim2")
    active_coord$Sample <- sample_names[!is_supp]
    active_coord$Group  <- as.character(orig_grp[!is_supp])
    active_coord$Type   <- "Active"
    
    has_supp <- any(is_supp) && !is.null(res$ind.sup)
    plot_df  <- active_coord
    
    if (has_supp) {
      supp_coord <- as.data.frame(res$ind.sup$coord[, 1:2, drop = FALSE])
      colnames(supp_coord) <- c("Dim1", "Dim2")
      supp_names <- sample_names[is_supp]
      supp_coord$Sample <- supp_names
      supp_coord$Group  <- dplyr::case_when(
        is_iqc_sample(supp_names)  ~ "iQC",
        is_istd_sample(supp_names) ~ "ISTD",
        TRUE ~ "Other (supp.)"
      )
      supp_coord$Type <- "Supplementary"
      plot_df <- dplyr::bind_rows(plot_df, supp_coord)
    }
    
    supp_levels <- if (has_supp) intersect(c("iQC", "ISTD", "Other (supp.)"), unique(plot_df$Group)) else character(0)
    grp_levels  <- c(sort(unique(active_coord$Group)), supp_levels)
    plot_df$Group <- factor(plot_df$Group, levels = grp_levels)
    
    supp_colours <- c(iQC = "#444444", ISTD = "#d9822b", `Other (supp.)` = "#999999")
    palette <- scales::hue_pal()(length(setdiff(grp_levels, names(supp_colours))))
    names(palette) <- setdiff(grp_levels, names(supp_colours))
    palette <- c(palette, supp_colours[intersect(names(supp_colours), grp_levels)])
    
    eig     <- res$eig
    pc1_pct <- round(eig[1, 2], 1)
    pc2_pct <- round(eig[2, 2], 1)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, colour = Group, shape = Type)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      ggplot2::geom_point(size = 4, alpha = 0.9) +
      ggplot2::scale_colour_manual(values = palette) +
      ggplot2::scale_shape_manual(values = c(Active = 16, Supplementary = 17), guide = "none") +
      ggplot2::labs(
        x      = paste0("Dim1 (", pc1_pct, "%)"),
        y      = paste0("Dim2 (", pc2_pct, "%)"),
        colour = "Group"
      ) +
      ggplot2::theme_minimal(base_size = fsz) +
      ggplot2::theme(
        legend.position  = "right",
        panel.grid.minor = ggplot2::element_blank()
      )
    
    if (has_supp) {
      p <- p + ggplot2::labs(caption = "Triangles = supplementary points (iQC/ISTD); projected after fitting, so they never influenced the PCA axes.")
    }
    
    if (isTRUE(show_labels)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = Sample),
        size = fsz / 4, show.legend = FALSE, max.overlaps = Inf
      )
    }
    
    p
  }
  
  output$pcaPlot <- renderPlot({
    .build_pca_plot(input$pca_export_fontsize %||% 14, isTRUE(input$pca_show_labels))
  })
  
  .pca_export_dims <- function() {
    px_w  <- input$pca_export_width  %||% 1200
    px_h  <- input$pca_export_height %||% 900
    dpi   <- input$pca_export_dpi    %||% 300
    scale <- input$pca_export_scale  %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  output$download_pca_png <- downloadHandler(
    filename = function() paste0("pca_plot_", Sys.Date(), ".png"),
    content = function(file) {
      dims <- isolate(.pca_export_dims())
      fsz  <- isolate(input$pca_export_fontsize %||% 14)
      lbl  <- isolate(isTRUE(input$pca_show_labels))
      p    <- isolate(.build_pca_plot(fsz, lbl))
      validate(need(!is.null(p), "No plot to export."))
      ggplot2::ggsave(file, plot = p,
                      width = dims$w, height = dims$h,
                      dpi = dims$dpi, device = "png")
    }
  )
  
  output$download_pca_svg <- downloadHandler(
    filename = function() paste0("pca_plot_", Sys.Date(), ".svg"),
    content = function(file) {
      dims <- isolate(.pca_export_dims())
      fsz  <- isolate(input$pca_export_fontsize %||% 14)
      lbl  <- isolate(isTRUE(input$pca_show_labels))
      p    <- isolate(.build_pca_plot(fsz, lbl))
      validate(need(!is.null(p), "No plot to export."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  
  output$pcaLoadingsTable <- renderDT({
    req(pca_result())
    res  <- pca_result()
    load <- as.data.frame(res$var$coord)
    cos2 <- as.data.frame(res$var$cos2)
    ctrb <- as.data.frame(res$var$contrib)
    
    k <- min(5, ncol(load))
    out <- cbind(
      Variable = rownames(load),
      load[, seq_len(k), drop = FALSE],
      cos2[, seq_len(min(2, ncol(cos2))), drop = FALSE],
      ctrb[, seq_len(min(2, ncol(ctrb))), drop = FALSE]
    )
    
    colnames(out)[2:(1 + k)] <- paste0("Loading_PC", seq_len(k))
    next_start <- 2 + k
    if (ncol(cos2) >= 1) {
      c2 <- min(2, ncol(cos2))
      colnames(out)[next_start:(next_start + c2 - 1)] <- paste0("cos2_PC", seq_len(c2))
      next_start <- next_start + c2
    }
    if (ncol(ctrb) >= 1) {
      cc <- min(2, ncol(ctrb))
      colnames(out)[next_start:(next_start + cc - 1)] <- paste0("contrib_PC", seq_len(cc))
    }
    
    num_cols <- vapply(out, is.numeric, logical(1))
    out[num_cols] <- lapply(out[num_cols], function(x) round(x, 4))
    datatable(out, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$exportWideDT <- renderDT({
    req(export_wide_data())
    datatable(export_wide_data(), options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$download_export_wide <- downloadHandler(
    filename = function() {
      mode <- if (identical(input$export_display_mode, "percent")) "percent" else "absolute"
      kind <- if (identical(input$export_value_type, "value_bs")) "bgsub" else "quant"
      paste0("export_wide_", kind, "_", mode, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      readr::write_csv(export_wide_data(), file)
    }
  )
  
  # ---------- UI renders ----------
  output$wideData <- renderDT({
    req(data_clean())
    datatable(data_clean(), options = list(scrollX = TRUE))
  })
  
  output$ISTDData <- renderDT({
    req(is_table_reactive())
    datatable(is_table_reactive(), options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$isPlot <- renderPlotly({
    req(is_plot_data())
    pd       <- is_plot_data()
    df_pts   <- pd$points
    df_sum   <- pd$summary
    mode_sel <- pd$mode
    
    err_vec <- switch(input$error_type_is,
                      "SEM"    = df_sum$sem,
                      "SD"     = df_sum$sd,
                      "95% CI" = df_sum$ci95,
                      df_sum$sem)
    df_sum$ymin <- df_sum$mean - err_vec
    df_sum$ymax <- df_sum$mean + err_vec
    
    if (mode_sel == "both") {
      p <- ggplot() +
        geom_col(data = df_sum,
                 aes(x = group, y = mean, fill = ion.mode,
                     text = paste0("Group: ", group, "<br>Ion mode: ", ion.mode,
                                   "<br>Mean: ", round(mean, 2))),
                 position = position_dodge(width = 0.6), width = 0.6) +
        geom_errorbar(data = df_sum,
                      aes(x = group, ymin = ymin, ymax = ymax, colour = ion.mode),
                      position = position_dodge(width = 0.6), width = 0.2) +
        geom_point(data = df_pts,
                   aes(x = group, y = IS_value, colour = ion.mode,
                       text = paste0("IS: ", istd_name,
                                     "<br>Sample: ", sample,
                                     "<br>Ion mode: ", ion.mode,
                                     "<br>IS value: ", round(IS_value, 2))),
                   position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6),
                   alpha = 0.7, size = 2) +
        labs(title = paste0("IS Values (Class: ", input$is_class, ") — Both ion modes"),
             x = "Group", y = "IS Value") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.grid.minor = element_blank())
    } else {
      p <- ggplot() +
        geom_col(data = df_sum,
                 aes(x = group, y = mean,
                     text = paste0("Group: ", group, "<br>Mean IS: ", round(mean, 2))),
                 fill = "#72B7B2", width = 0.7) +
        geom_errorbar(data = df_sum,
                      aes(x = group, ymin = ymin, ymax = ymax),
                      width = 0.2, colour = "#333") +
        geom_jitter(data = df_pts,
                    aes(x = group, y = IS_value,
                        text = paste0("IS: ", istd_name,
                                      "<br>Sample: ", sample,
                                      "<br>IS value: ", round(IS_value, 2))),
                    width = 0.12, height = 0, alpha = 0.6, size = 2, colour = "#E45756") +
        labs(title = paste0("IS Values (Class: ", input$is_class, ") — ",
                            switch(mode_sel,
                                   "auto" = "Auto (preferred mode)",
                                   "neg"  = "Negative only",
                                   "pos"  = "Positive only")),
             x = "Group", y = "IS Value") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.grid.minor = element_blank())
    }
    
    ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        legend = list(orientation = "v"),
        margin = list(b = 80)
      )
  })
  
  
  output$metPlot <- renderPlotly({
    req(met_plot_data())
    pd <- met_plot_data()
    df_pts <- pd$points
    df_sum <- pd$summary
    
    # Ensure group is a factor
    df_sum$group <- factor(df_sum$group)
    df_pts$group <- factor(df_pts$group, levels = levels(df_sum$group))
    group_levels <- levels(df_sum$group)
    
    # Numeric conversion
    df_sum$group_num <- as.numeric(df_sum$group)
    df_pts$group_num <- as.numeric(df_pts$group)
    
    # Remove NA rows
    df_sum <- df_sum[!is.na(df_sum$group_num) & is.finite(df_sum$mean), ]
    df_pts <- df_pts[!is.na(df_pts$group_num) & is.finite(df_pts$value), ]
    
    # Error bars
    err_vec <- switch(input$error_type_met,
                      "SEM" = df_sum$sem,
                      "SD" = df_sum$sd,
                      "95% CI" = df_sum$ci95,
                      df_sum$sem
    )
    err_vec <- suppressWarnings(as.numeric(err_vec))
    err_vec[!is.finite(err_vec)] <- 0
    ymin <- df_sum$mean - err_vec
    ymax <- df_sum$mean + err_vec
    
    # Hover text
    bar_hover_text <- switch(input$error_type_met,
                             "SEM" = paste0("Mean: ", signif(df_sum$mean, 5), "\nSEM: ", signif(df_sum$sem, 5)),
                             "SD" = paste0("Mean: ", signif(df_sum$mean, 5), "\nSD: ", signif(df_sum$sd, 5)),
                             "95% CI" = paste0("Mean: ", signif(df_sum$mean, 5), "\n95% CI: ", signif(df_sum$ci95, 5)),
                             paste0("Mean: ", signif(df_sum$mean, 5), "\nSEM: ", signif(df_sum$sem, 5))
    )
    df_pts$val_fmt <- signif(df_pts$value, 5)
    pt_hover_text <- paste0("Sample: ", as.character(df_pts$sample), "\nValue: ", df_pts$val_fmt)
    
    # Try ggplotly first
    p <- ggplot() +
      geom_col(data = df_sum, aes(x = group, y = mean, text = bar_hover_text), fill = "#54A24B", width = 0.7) +
      geom_errorbar(data = transform(df_sum, ymin = ymin, ymax = ymax), aes(x = group, ymin = ymin, ymax = ymax), width = 0.2, colour = "#333") +
      geom_jitter(data = df_pts, aes(x = group, y = value, text = pt_hover_text), width = 0.12, height = 0, alpha = 0.6, size = 2, colour = "#E45756") +
      labs(
        title = paste0("Metabolite: ", input$met_name, " (Class: ", input$met_class, ")"),
        x = "Group",
        y = pd$ylab
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), panel.grid.minor = element_blank())
    
    tryCatch({
      ggplotly(p, tooltip = "text") %>%
        layout(showlegend = FALSE)
    }, error = function(e) {
      plt <- plotly::plot_ly()
      
      # Bar plot for group means
      plt <- plt %>%
        add_bars(
          data = df_sum,
          x = ~group_num,
          y = ~mean,
          text = ~bar_hover_text,
          hoverinfo = "text+y",
          marker = list(color = "#54A24B"),
          showlegend = FALSE,
          textposition = "none"
        )
      
      # Capped error bars
      cap_width <- 0.3
      for (i in seq_len(nrow(df_sum))) {
        # Vertical line
        plt <- plt %>%
          add_segments(
            x = df_sum$group_num[i], xend = df_sum$group_num[i],
            y = ymin[i], yend = ymax[i],
            line = list(color = "#333", width = 2),
            showlegend = FALSE,
            hoverinfo = "skip"
          )
        # Lower cap
        plt <- plt %>%
          add_segments(
            x = df_sum$group_num[i] - cap_width/2, xend = df_sum$group_num[i] + cap_width/2,
            y = ymin[i], yend = ymin[i],
            line = list(color = "#333", width = 2),
            showlegend = FALSE,
            hoverinfo = "skip"
          )
        # Upper cap
        plt <- plt %>%
          add_segments(
            x = df_sum$group_num[i] - cap_width/2, xend = df_sum$group_num[i] + cap_width/2,
            y = ymax[i], yend = ymax[i],
            line = list(color = "#333", width = 2),
            showlegend = FALSE,
            hoverinfo = "skip"
          )
      }
      
      # Overlay points for individual samples
      plt <- plt %>%
        add_markers(
          data = df_pts,
          x = ~group_num,
          y = ~value,
          text = ~pt_hover_text,
          hoverinfo = "text+y",
          marker = list(color = "#E45756", size = 8, opacity = 0.7),
          showlegend = FALSE
        )
      
      # Set x-axis labels
      if (!is.null(group_levels) && length(group_levels) > 0) {
        plt <- plt %>%
          layout(
            title = paste0("Metabolite: ", input$met_name, " (Class: ", input$met_class, ")"),
            xaxis = list(
              title = "Group",
              tickvals = seq_along(group_levels),
              ticktext = group_levels
            ),
            yaxis = list(title = pd$ylab),
            barmode = "overlay",
            showlegend = FALSE
          )
      } else {
        plt <- plt %>%
          layout(
            title = paste0("Metabolite: ", input$met_name, " (Class: ", input$met_class, ")"),
            xaxis = list(title = "Group"),
            yaxis = list(title = pd$ylab),
            barmode = "overlay",
            showlegend = FALSE
          )
      }
      plt
    })
  })
  
  
  output$metPlotSummary <- renderDT({
    req(met_plot_data())
    df <- met_plot_data()$summary %>% dplyr::mutate(across(where(is.numeric), ~round(., 4)))
    datatable(df, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$classAllPlot <- renderPlotly({
    req(class_all_plot_data())
    pd     <- class_all_plot_data()
    df_pts <- pd$points
    df_sum <- pd$summary
    
    # --- Error vector and whiskers (unchanged)
    err_vec <- switch(input$error_type_all,
                      "SEM"    = df_sum$sem,
                      "SD"     = df_sum$sd,
                      "95% CI" = df_sum$ci95,
                      df_sum$sem)
    ymin <- df_sum$mean - err_vec
    ymax <- df_sum$mean + err_vec
    
    # --- Keep your ordering logic (already done in class_all_plot_data)
    # Factor the metabolite names for stable ordering in coord_flip
    df_sum <- df_sum %>%
      dplyr::mutate(`Metabolite name` = factor(`Metabolite name`,
                                               levels = rev(`Metabolite name`)))
    df_pts <- df_pts %>%
      dplyr::mutate(`Metabolite name` = factor(`Metabolite name`,
                                               levels = levels(df_sum$`Metabolite name`)))
    
    # --- Build hover text content
    # Format numbers for cleaner display
    df_sum$mean_fmt <- signif(df_sum$mean, 5)
    df_sum$sem_fmt  <- signif(df_sum$sem, 5)
    df_sum$sd_fmt   <- signif(df_sum$sd, 5)
    df_sum$ci_fmt   <- signif(df_sum$ci95, 5)
    
    err_lab <- switch(input$error_type_all,
                      "SEM"    = "SEM",
                      "SD"     = "SD",
                      "95% CI" = "95% CI",
                      "SEM")
    
    # Bar hover: show mean plus chosen error metric
    bar_hover_text <- switch(input$error_type_all,
                             "SEM"    = paste0("Mean: ", df_sum$mean_fmt, "\n", err_lab, ": ", df_sum$sem_fmt),
                             "SD"     = paste0("Mean: ", df_sum$mean_fmt, "\n", err_lab, ": ", df_sum$sd_fmt),
                             "95% CI" = paste0("Mean: ", df_sum$mean_fmt, "\n", err_lab, ": ", df_sum$ci_fmt),
                             paste0("Mean: ", df_sum$mean_fmt, "\n", err_lab, ": ", df_sum$sem_fmt)
    )
    
    # Points hover: sample name + metabolite name + value
    df_pts$val_fmt <- signif(df_pts$value, 5)
    pt_hover_text  <- paste0("Sample: ", as.character(df_pts$sample),
                             "\nMetabolite: ", as.character(df_pts$`Metabolite name`),
                             "\nValue: ", df_pts$val_fmt)
    
    # --- Plot
    p <- ggplot() +
      # Bars with our mean-based hover text
      geom_col(
        data = df_sum,
        aes(x = `Metabolite name`, y = mean, text = bar_hover_text),
        fill = "#4C78A8", width = 0.7
      ) +
      # Error bars (hover off—bars already expose the values)
      geom_errorbar(
        data = transform(df_sum, ymin = ymin, ymax = ymax),
        aes(x = `Metabolite name`, ymin = ymin, ymax = ymax),
        width = 0.2, colour = "#333"
      ) +
      # Points (individual sample values) with per-point hover
      geom_jitter(
        data = df_pts,
        aes(x = `Metabolite name`, y = value, text = pt_hover_text),
        width = 0.15, height = 0, alpha = 0.5, size = 1.8, colour = "#F58518"
      ) +
      coord_flip() +
      labs(
        title = paste0(
          "All Metabolites in Class: ", input$class_all,
          " (Group: ", input$class_all_selected_group, ", Order: ",
          if (identical(input$class_all_order, "alpha_asc")) "Alphabetical" else "Abundance",
          ")"
        ),
        x = "Metabolite",
        y = pd$ylab
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
    
    plot_height <- max(400, length(unique(df_sum$`Metabolite name`)) * 30)
    
    ggplotly(p, tooltip = "text", height = plot_height)
  })
  
  
  output$classAllSummary <- renderDT({
    req(class_all_plot_data())
    df <- class_all_plot_data()$summary %>% dplyr::mutate(across(where(is.numeric), ~round(., 4)))
    datatable(df, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$bgNormTable <- renderDT({
    req(bg_norm_long())
    df <- bg_norm_long() %>%
      dplyr::select(class, ion.mode, `Metabolite name`, `Adduct type`,
                    sample, value, value_bs, IS_value, norm)
    datatable(df, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  # =========================
  # ===== STATISTICS TAB ====
  # =========================
  
  .safe_sd <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s)) NA_real_ else s
  }
  cohens_d_unpaired <- function(x, g) {
    if (length(levels(g)) != 2) return(NA_real_)
    g <- droplevels(g)
    x1 <- x[g == levels(g)[1]]; x2 <- x[g == levels(g)[2]]
    n1 <- sum(!is.na(x1)); n2 <- sum(!is.na(x2))
    if (n1 < 2 || n2 < 2) return(NA_real_)
    m1 <- mean(x1, na.rm = TRUE); m2 <- mean(x2, na.rm = TRUE)
    sd1 <- .safe_sd(x1); sd2 <- .safe_sd(x2)
    if (!is.finite(sd1) || !is.finite(sd2)) return(NA_real_)
    sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    (m2 - m1) / sp
  }
  eta2_from_aov_table <- function(aov_tbl, effect_row_name) {
    if (is.null(aov_tbl) || nrow(aov_tbl) == 0) return(NA_real_)
    if (!("Sum Sq" %in% colnames(aov_tbl))) return(NA_real_)
    ss_total <- sum(aov_tbl[,"Sum Sq"], na.rm = TRUE)
    if (!effect_row_name %in% rownames(aov_tbl)) return(NA_real_)
    ss_eff <- aov_tbl[effect_row_name, "Sum Sq"]
    if (!is.finite(ss_total) || ss_total == 0) return(NA_real_)
    as.numeric(ss_eff / ss_total)
  }
  
  
  observeEvent(bg_norm_long_avg(), {
    classes <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::distinct(plot_class) %>%
      dplyr::arrange(plot_class) %>%
      dplyr::pull(plot_class)
    
    classes <- c("All classes", classes)  # ⬅️ add the “All” option
    updateSelectInput(session, "stats_class",
                      choices = classes,
                      selected = classes[1]
    )
  })
  
  
  
  observeEvent(input$stats_class, {
    req(bg_norm_long_avg(), input$stats_class)
    
    base <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]"))
    
    mets <- if (identical(input$stats_class, "All classes")) {
      base %>%
        dplyr::distinct(`Metabolite name`) %>%
        dplyr::arrange(`Metabolite name`) %>%
        dplyr::pull(`Metabolite name`)
    } else {
      base %>%
        dplyr::filter(plot_class == input$stats_class) %>%
        dplyr::distinct(`Metabolite name`) %>%
        dplyr::arrange(`Metabolite name`) %>%
        dplyr::pull(`Metabolite name`)
    }
    
    updateSelectInput(session, "stats_met_name",
                      choices  = if (length(mets) > 0) mets else "No metabolites in this class",
                      selected = if (length(mets) > 0) mets[1] else "No metabolites in this class")
  })
  
  
  # Build dataset to test
  stats_input_long <- reactive({
    req(bg_norm_long_avg(), input$stats_class, input$stats_value_type, input$stats_units)
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      { if (!identical(input$stats_class, "All classes")) dplyr::filter(., plot_class == input$stats_class) else . } %>%
      dplyr::mutate(sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample)))
    
    # Single, consistent blank & iQC filter
    df <- filter_blanks(df, input$exclude_blank_stats, sample_col = "sample_norm", exact = FALSE)
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    measure_col <- if (identical(input$stats_value_type, "value_bs")) "value_bs" else "norm"
    if (identical(input$stats_units, "percent")) {
      totals <- df %>%
        dplyr::group_by(sample, plot_class) %>%
        dplyr::summarise(class_total = sum(.data[[measure_col]], na.rm = TRUE), .groups = "drop")
      df <- df %>%
        dplyr::left_join(totals, by = c("sample", "plot_class")) %>%
        dplyr::mutate(value = dplyr::if_else(class_total > 0, 100 * .data[[measure_col]] / class_total, NA_real_)) %>%
        dplyr::select(sample, plot_class, `Metabolite name`, value)
    } else {
      df <- df %>%
        dplyr::transmute(sample, plot_class, `Metabolite name`, value = .data[[measure_col]])
    }
    if (identical(input$stats_scope, "single")) {
      df <- df %>% dplyr::filter(`Metabolite name` == input$stats_met_name)
    }
    if (identical(input$stats_scope, "total")) {
      df <- df %>%
        dplyr::group_by(sample, plot_class) %>%
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(`Metabolite name` = "Total (class sum)")
    }
    
    # ── Apply central grouping ────────────────────────────────────────────────
    fn  <- get_active_grouping()
    res <- fn(df$sample)
    df$group   <- res$group
    df$factorA <- res$factorA
    df$factorB <- res$factorB
    
    # ── Group selection filter (from Stats tab checkbox) ──────────────────────
    sel_grps <- input$stats_selected_groups
    if (!is.null(sel_grps) && length(sel_grps) > 0 && "group" %in% names(df)) {
      df <- df %>% dplyr::filter(.data$group %in% sel_grps)
    }
    
    df
  })
  .compute_stats_results <- function(df, test_choice, equal_var, tw_filter, alpha, padj_method, use_adj, show_all_rows_flag) {
    spl <- split(df, df$`Metabolite name`)
    out <- vector("list", length(spl))
    i <- 0L
    
    withProgress(message = "Running statistics...", value = 0, {
      for (nm in names(spl)) {
        i <- i + 1L
        incProgress(1 / max(1, length(spl)), detail = nm)
        dd <- spl[[nm]] %>% dplyr::filter(!is.na(value), !is.na(group))
        if (nrow(dd) < 3) { out[[i]] <- NULL; next }
        
        if (dplyr::n_distinct(dd$group) == 2) {
          lv <- sort(unique(dd$group))
          dd$group <- factor(dd$group, levels = lv)
        } else {
          dd$group <- factor(dd$group)
        }
        k <- nlevels(dd$group)
        
        used <- test_choice
        if (identical(used, "auto")) used <- if (k == 2) "ttest" else if (k >= 3) "oneway" else NA_character_
        if (is.na(used)) { out[[i]] <- NULL; next }
        
        res_row <- NULL
        
        if (identical(used, "ttest")) {
          if (k != 2) { out[[i]] <- NULL; next }
          tt <- tryCatch(stats::t.test(value ~ group, data = dd, var.equal = isTRUE(equal_var)), error = function(e) NULL)
          if (is.null(tt)) { out[[i]] <- NULL; next }
          pval <- as.numeric(tt$p.value)
          
          means <- dd %>% dplyr::group_by(group) %>%
            dplyr::summarise(mean = mean(value, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
          
          d_eff <- tryCatch(cohens_d_unpaired(dd$value, dd$group), error = function(e) NA_real_)
          
          g1 <- levels(dd$group)[1]; g2 <- levels(dd$group)[2]
          m1 <- means$mean[means$group == g1]; m2 <- means$mean[means$group == g2]
          pseudo <- 1e-9
          log2fc <- log2((m2 + pseudo) / (m1 + pseudo))
          
          res_row <- tibble::tibble(
            plot_class = dd$plot_class[1],
            metabolite = nm,
            test       = "ttest",
            effect     = "group",
            p          = pval,
            cohen_d    = d_eff,
            log2FC     = log2fc,
            grp_ref    = g1,
            grp_comp   = g2,
            groups     = paste0(paste0(means$group, " (n=", means$n, ", mean=", round(means$mean, 4), ")"), collapse = " | ")
          )
          
        } else if (identical(used, "oneway")) {
          if (k < 3) { out[[i]] <- NULL; next }
          fit <- tryCatch(stats::aov(value ~ group, data = dd), error = function(e) NULL)
          if (is.null(fit)) { out[[i]] <- NULL; next }
          sm <- tryCatch(summary(fit)[[1]], error = function(e) NULL)
          if (is.null(sm)) { out[[i]] <- NULL; next }
          pval <- tryCatch(sm[["Pr(>F)"]][1], error = function(e) NA_real_)
          eta2 <- tryCatch(eta2_from_aov_table(sm, "group"), error = function(e) NA_real_)
          
          means <- dd %>% dplyr::group_by(group) %>%
            dplyr::summarise(mean = mean(value, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
          
          res_row <- tibble::tibble(
            plot_class = dd$plot_class[1],
            metabolite = nm,
            test       = "oneway",
            effect     = "group",
            p          = as.numeric(pval),
            eta2       = eta2,
            groups     = paste0(paste0(means$group, " (n=", means$n, ", mean=", round(means$mean, 4), ")"), collapse = " | ")
          )
          
        } else if (identical(used, "twoway")) {
          dd2 <- dd %>% dplyr::filter(!is.na(factorA), !is.na(factorB))
          if (dplyr::n_distinct(dd2$factorA) != 2 || dplyr::n_distinct(dd2$factorB) != 2) { out[[i]] <- NULL; next }
          dd2$factorA <- factor(dd2$factorA); dd2$factorB <- factor(dd2$factorB)
          fit <- tryCatch(stats::aov(value ~ factorA * factorB, data = dd2), error = function(e) NULL)
          if (is.null(fit)) { out[[i]] <- NULL; next }
          sm <- tryCatch(summary(fit)[[1]], error = function(e) NULL)
          if (is.null(sm)) { out[[i]] <- NULL; next }
          
          pA  <- suppressWarnings(sm["factorA", "Pr(>F)"]);  eA  <- tryCatch(eta2_from_aov_table(sm, "factorA"), error = function(e) NA_real_)
          pB  <- suppressWarnings(sm["factorB", "Pr(>F)"]);  eB  <- tryCatch(eta2_from_aov_table(sm, "factorB"), error = function(e) NA_real_)
          pAB <- suppressWarnings(sm["factorA:factorB", "Pr(>F)"]); eAB <- tryCatch(eta2_from_aov_table(sm, "factorA:factorB"), error = function(e) NA_real_)
          
          sel <- tw_filter
          chosen_p  <- NA_real_; chosen_eff <- NA_character_; chosen_eta2 <- NA_real_
          if (sel == "A")  { chosen_p <- pA;  chosen_eff <- "factorA";         chosen_eta2 <- eA  }
          if (sel == "B")  { chosen_p <- pB;  chosen_eff <- "factorB";         chosen_eta2 <- eB  }
          if (sel == "AxB"){ chosen_p <- pAB; chosen_eff <- "factorA:factorB"; chosen_eta2 <- eAB }
          if (sel == "any" || !is.finite(chosen_p)) {
            allp <- c(A = pA, B = pB, `A:B` = pAB)
            ok <- is.finite(allp)
            if (any(ok)) {
              idx <- which.min(allp[ok])
              chosen_p <- allp[ok][idx]
              chosen_eff <- names(allp[ok])[idx]
              chosen_eta2 <- c(A = eA, B = eB, `A:B` = eAB)[chosen_eff]
            }
          }
          if (!is.finite(chosen_p)) { out[[i]] <- NULL; next }
          
          dd2 <- dd2 %>% dplyr::mutate(AB = interaction(factorA, factorB, drop = TRUE))
          means <- dd2 %>% dplyr::group_by(AB) %>%
            dplyr::summarise(mean = mean(value, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
          
          res_row <- tibble::tibble(
            plot_class = dd2$plot_class[1],
            metabolite = nm,
            test       = "twoway",
            effect     = chosen_eff,
            p          = as.numeric(chosen_p),
            eta2       = chosen_eta2,
            groups     = paste0(paste0(means$AB, " (n=", means$n, ", mean=", round(means$mean, 4), ")"), collapse = " | ")
          )
        }
        
        out[[i]] <- res_row
      }
    })
    
    res <- dplyr::bind_rows(out)
    if (nrow(res) == 0) {
      return(tibble::tibble(plot_class=character(0), metabolite=character(0),
                            test=character(0), effect=character(0), p=numeric(0),
                            p_adj=numeric(0), groups=character(0)))
    }
    
    res <- res %>%
      dplyr::mutate(p_adj = stats::p.adjust(p, method = if (identical(padj_method, "none")) "none" else padj_method))
    
    sig_mask <- if (!isTRUE(use_adj)) res$p < alpha else res$p_adj < alpha
    if (!isTRUE(show_all_rows_flag)) {
      res <- res[sig_mask, , drop = FALSE]
    }
    
    res %>% dplyr::arrange(dplyr::coalesce(p_adj, p), p)
  }
  
  
  output$stats_header <- renderText({
    res <- stats_results_val()
    if (is.null(res) || nrow(res) == 0) {
      "No results yet — adjust options to see results."
    } else {
      paste0("Results: ", nrow(res), " feature(s) returned",
             if (!isTRUE(input$show_all_rows)) " (filtered by significance)" else "")
    }
  })
  
  output$statsResults <- renderDT({
    res <- stats_results_val()
    validate(
      need(!is.null(res), 
           "Statistics have not been run yet. Please click 'Run statistics' to begin. You can let 'Auto' choose the test, or pick a test and then run."),
      need(nrow(res) > 0, 
           "No results to display. Try adjusting your test or thresholds.")
    )
    numeric_candidates <- c("p","p_adj","cohen_d","log2FC","eta2")
    present_num <- intersect(numeric_candidates, names(res))
    if (length(present_num) > 0) {
      res[present_num] <- lapply(res[present_num], function(x) round(x, 6))
    }
    
    tests_present <- unique(res$test)
    base_cols <- c("plot_class","metabolite","test","effect","p","p_adj")
    cols <- base_cols
    if ("ttest" %in% tests_present) {
      cols <- c(cols, intersect(c("cohen_d","log2FC","grp_ref","grp_comp"), names(res)))
    }
    if (any(c("oneway","twoway") %in% tests_present)) {
      cols <- c(cols, intersect("eta2", names(res)))
    }
    cols <- c(cols, intersect("groups", names(res)))
    
    res_disp <- res[, cols, drop = FALSE]
    datatable(res_disp, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$download_stats <- downloadHandler(
    filename = function() {
      paste0("stats_", input$stats_class, "_", input$stats_value_type, "_", input$stats_units, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- stats_results_val()
      readr::write_csv(res, file)
    }
  )
  
  # ---- Post-hoc results table (single metabolite focus) ----
  posthoc_results <- reactive({
    req(stats_input_long(), input$posthoc)
    if (length(input$posthoc) == 0L) return(NULL)
    if (!identical(input$stats_scope, "single")) return(NULL)
    
    dd <- stats_input_long()
    if (nrow(dd) == 0) return(NULL)
    met <- unique(dd$`Metabolite name`)
    if (length(met) != 1) return(NULL)
    
    out <- list()
    has_rstatix <- requireNamespace("rstatix", quietly = TRUE)
    k <- dplyr::n_distinct(dd$group)
    
    if ((identical(input$stats_test, "oneway") || identical(input$stats_test, "auto")) && k >= 3) {
      fit <- tryCatch(stats::aov(value ~ group, data = dd), error = function(e) NULL)
      if (!is.null(fit)) {
        if ("tukey" %in% input$posthoc) {
          if (has_rstatix) {
            tk <- tryCatch(rstatix::tukey_hsd(fit), error = function(e) NULL)
            if (!is.null(tk) && nrow(tk) > 0) {
              tk <- tk %>% dplyr::mutate(method = "TukeyHSD (one-way)") %>%
                dplyr::rename(p.adj = adj.p.value)
              out <- c(out, list(tk))
            }
          } else {
            tk <- tryCatch(stats::TukeyHSD(fit, "group"), error = function(e) NULL)
            if (!is.null(tk)) {
              otk <- as.data.frame(tk$group)
              otk$contrast <- rownames(otk); rownames(otk) <- NULL
              colnames(otk) <- sub("^p\\.adj$", "p.adj", colnames(otk))
              otk$method <- "TukeyHSD (one-way, base R)"
              out <- c(out, list(otk))
            }
          }
        }
        if ("pairwise" %in% input$posthoc) {
          if (has_rstatix) {
            pw <- tryCatch(rstatix::pairwise_t_test(value ~ group, dd, p.adjust.method = "holm"), error = function(e) NULL)
            if (!is.null(pw) && nrow(pw) > 0) {
              pw$method <- "Pairwise t-tests (Holm) — rstatix"
              out <- c(out, list(pw))
            }
          } else {
            pw <- tryCatch(stats::pairwise.t.test(dd$value, dd$group, p.adjust.method = "holm"), error = function(e) NULL)
            if (!is.null(pw) && !is.null(pw$p.value)) {
              mat <- pw$p.value
              mm <- tibble::tibble()
              for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
                if (!is.na(mat[i, j])) {
                  mm <- dplyr::bind_rows(mm, tibble::tibble(
                    contrast = paste0(rownames(mat)[i], " - ", colnames(mat)[j]),
                    p.adj = mat[i, j],
                    method = "Pairwise t-tests (Holm) — base R"
                  ))
                }
              }
              if (nrow(mm) > 0) out <- c(out, list(mm))
            }
          }
        }
      }
    }
    
    if (identical(input$stats_test, "twoway")) {
      dd2 <- dd %>% dplyr::filter(!is.na(factorA), !is.na(factorB))
      if (dplyr::n_distinct(dd2$factorA) == 2 && dplyr::n_distinct(dd2$factorB) == 2) {
        dd2 <- dd2 %>% dplyr::mutate(AB = interaction(factorA, factorB, drop = TRUE))
        fit2 <- tryCatch(stats::aov(value ~ factorA * factorB, data = dd2), error = function(e) NULL)
        
        if (!is.null(fit2) && "tukey" %in% input$posthoc) {
          if (has_rstatix) {
            fit_ab <- tryCatch(stats::aov(value ~ AB, data = dd2), error = function(e) NULL)
            tk2 <- if (!is.null(fit_ab)) tryCatch(rstatix::tukey_hsd(fit_ab), error = function(e) NULL) else NULL
            if (!is.null(tk2) && nrow(tk2) > 0) {
              tk2 <- tk2 %>% dplyr::mutate(method = "TukeyHSD on A:B — rstatix") %>%
                dplyr::rename(p.adj = adj.p.value)
              out <- c(out, list(tk2))
            }
          } else {
            tk2 <- tryCatch(stats::TukeyHSD(fit2, "factorA:factorB"), error = function(e) NULL)
            if (!is.null(tk2)) {
              otk2 <- as.data.frame(tk2$`factorA:factorB`)
              otk2$contrast <- rownames(otk2); rownames(otk2) <- NULL
              colnames(otk2) <- sub("^p\\.adj$", "p.adj", colnames(otk2))
              otk2$method <- "TukeyHSD on A:B — base R"
              out <- c(out, list(otk2))
            }
          }
        }
        
        if ("pairwise" %in% input$posthoc) {
          if (has_rstatix) {
            pw2 <- tryCatch(rstatix::pairwise_t_test(value ~ AB, dd2, p.adjust.method = "holm"), error = function(e) NULL)
            if (!is.null(pw2) && nrow(pw2) > 0) {
              pw2$method <- "Pairwise t-tests (Holm) on A:B — rstatix"
              out <- c(out, list(pw2))
            }
          } else {
            pw2 <- tryCatch(stats::pairwise.t.test(dd2$value, dd2$AB, p.adjust.method = "holm"), error = function(e) NULL)
            if (!is.null(pw2) && !is.null(pw2$p.value)) {
              mat <- pw2$p.value
              mm2 <- tibble::tibble()
              for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
                if (!is.na(mat[i, j])) {
                  mm2 <- dplyr::bind_rows(mm2, tibble::tibble(
                    contrast = paste0(rownames(mat)[i], " - ", colnames(mat)[j]),
                    p.adj = mat[i, j],
                    method = "Pairwise t-tests (Holm) on A:B — base R"
                  ))
                }
              }
              if (nrow(mm2) > 0) out <- c(out, list(mm2))
            }
          }
        }
      }
    }
    
    if (length(out) == 0) return(NULL)
    res <- dplyr::bind_rows(out, .id = "block")
    res$metabolite <- met
    num_cols <- vapply(res, is.numeric, logical(1))
    res[num_cols] <- lapply(res[num_cols], function(x) round(x, 6))
    res
  })
  
  output$posthocResults <- renderDT({
    ph <- posthoc_results()
    validate(need(!is.null(ph) && nrow(ph) > 0, "No post-hoc results to display."))
    datatable(ph, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  # ---- Per-class detection & significance summary ----
  class_stats_summary <- reactive({
    # Detected: all species in bg_norm_long_resolved (IS excluded)
    df_all <- tryCatch(bg_norm_long_avg(), error = function(e) NULL)
    if (is.null(df_all)) return(NULL)
    detected <- df_all %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::distinct(`Metabolite name`, plot_class) %>%
      dplyr::count(plot_class, name = "Detected")
    
    # Tested & significant: from stats_results_all / stats_results_val
    tested_df  <- stats_results_all()
    sig_df     <- stats_results_val()
    alpha      <- input$alpha %||% 0.05
    use_adj    <- isTRUE(input$use_adj_threshold) &&
      !identical(input$padj_method %||% "BH", "none")
    
    tested <- if (!is.null(tested_df) && nrow(tested_df) > 0)
      tested_df %>% dplyr::count(plot_class, name = "Tested")
    else
      tibble::tibble(plot_class = character(0), Tested = integer(0))
    
    sig <- if (!is.null(sig_df) && nrow(sig_df) > 0) {
      p_col <- if (use_adj && "p_adj" %in% names(sig_df)) "p_adj" else "p"
      sig_df %>%
        dplyr::filter(!is.na(.data[[p_col]]), .data[[p_col]] < alpha) %>%
        dplyr::count(plot_class, name = "Significant")
    } else {
      tibble::tibble(plot_class = character(0), Significant = integer(0))
    }
    
    out <- detected %>%
      dplyr::left_join(tested,  by = "plot_class") %>%
      dplyr::left_join(sig,     by = "plot_class") %>%
      dplyr::mutate(
        Tested      = tidyr::replace_na(Tested,      0L),
        Significant = tidyr::replace_na(Significant, 0L),
        `% significant` = dplyr::if_else(
          Tested > 0,
          paste0(round(Significant / Tested * 100, 1), "%"),
          "—"
        )
      ) %>%
      dplyr::arrange(plot_class) %>%
      dplyr::rename(Class = plot_class)
    
    out
  })
  
  output$classStatsSummary <- DT::renderDT({
    df <- class_stats_summary()
    validate(need(!is.null(df) && nrow(df) > 0,
                  "Load data to see the per-class summary."))
    DT::datatable(
      df,
      rownames  = FALSE,
      selection = "none",
      options   = list(
        dom        = "t",
        pageLength = nrow(df),
        scrollX    = TRUE,
        columnDefs = list(list(className = "dt-center",
                               targets   = 1:4))
      ),
      class = "stripe hover compact"
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleInterval(0, c("white", "#d4edda"))
      )
  })
  
  output$download_class_summary <- downloadHandler(
    filename = function() paste0("LipiRich_class_summary_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- class_stats_summary()
      validate(need(!is.null(df), "No data to export."))
      readr::write_csv(df, file)
    }
  )
  
  # ---- Bar plot data for significant features ----
  stats_sig_plot_data <- reactive({
    res <- stats_results_val()
    req(res, nrow(res) > 0)
    df_all <- stats_input_long()
    req(df_all)
    
    sort_mode <- input$stats_plot_sort
    if (sort_mode == "log2fc" && !"log2FC" %in% names(res)) {
      sort_mode <- "p_adj"
    }
    
    res_ranked <- res
    if (sort_mode == "p_adj") {
      res_ranked <- res_ranked %>% dplyr::arrange(dplyr::coalesce(p_adj, p), p)
    } else if (sort_mode == "p") {
      res_ranked <- res_ranked %>% dplyr::arrange(p)
    } else if (sort_mode == "log2fc" && "log2FC" %in% names(res_ranked)) {
      res_ranked <- res_ranked %>% dplyr::arrange(dplyr::desc(abs(log2FC)), dplyr::coalesce(p_adj, p))
    }
    
    topn <- input$stats_plot_topn %||% 20
    feats <- head(res_ranked$metabolite, topn)
    if (length(feats) == 0) return(NULL)
    
    df <- df_all %>% dplyr::filter(`Metabolite name` %in% feats)
    
    # --- NA-safe, scalar checks for two-way factors ---
    show_two_way <- any(res_ranked$test == "twoway")
    
    # Do we have at least one non-NA level in both A and B across the current df?
    hasA <- any(!is.na(df$factorA))
    hasB <- any(!is.na(df$factorB))
    has_two_way_factors <- hasA && hasB
    
    if (isTRUE(show_two_way) && isTRUE(has_two_way_factors)) {
      bar_group <- df %>%
        dplyr::mutate(bar_group = interaction(factorA, factorB, drop = TRUE))
    } else {
      bar_group <- df %>%
        dplyr::mutate(bar_group = factor(group))
    }
    
    # (optional) defensively drop any rows where bar_group is NA, just in case
    bar_group <- bar_group %>% dplyr::filter(!is.na(bar_group))
    
    
    df <- bar_group %>%
      dplyr::group_by(`Metabolite name`, bar_group) %>%
      dplyr::summarise(
        mean = mean(value, na.rm = TRUE),
        sd   = stats::sd(value, na.rm = TRUE),
        n    = dplyr::n(),
        sem  = sd / sqrt(pmax(n, 1)),
        ci95 = 1.96 * sem,
        .groups = "drop"
      )
    
    df$`Metabolite name` <- factor(df$`Metabolite name`, levels = feats)
    
    err <- switch(input$stats_error_type, "SEM" = df$sem, "SD" = df$sd, "95% CI" = df$ci95, df$sem)
    df$ymin <- df$mean - err
    df$ymax <- df$mean + err
    
    df
  })
  
  # ---- Post-hoc annotations for the bar plot ----
  # Shared helper: compute pairwise comparisons for a set of features.
  # df_long must have columns: Metabolite name, bar_group (or group), value.
  # Returns a data frame with Metabolite name, group1, group2, p, p.adj,
  # p.adj.signif — or NULL if nothing significant.
  .compute_pairwise <- function(df_long, group_col = "bar_group",
                                padj_method = "BH", alpha = 0.05,
                                use_adj = TRUE, filter_alpha = alpha) {
    normalize_pcols <- function(x) {
      nm <- names(x)
      if ("adj.p.value" %in% nm && !"p.adj" %in% nm) x <- dplyr::rename(x, p.adj = `adj.p.value`)
      if ("p_adj"       %in% nm && !"p.adj" %in% nm) x <- dplyr::rename(x, p.adj = p_adj)
      if ("p.signif"    %in% nm && !"p.adj.signif" %in% nm) x <- dplyr::rename(x, p.adj.signif = p.signif)
      x
    }
    use_rstatix <- requireNamespace("rstatix", quietly = TRUE)
    feats       <- unique(df_long$`Metabolite name`)
    ann_list    <- list()
    
    for (met in feats) {
      dmet <- df_long %>%
        dplyr::filter(`Metabolite name` == met, !is.na(.data[[group_col]]),
                      !is.na(value))
      if (!identical(group_col, "bar_group"))
        dmet <- dmet %>% dplyr::rename(bar_group = dplyr::all_of(group_col))
      k <- dplyr::n_distinct(dmet$bar_group)
      if (k < 2) next
      comps <- NULL
      
      if (k == 2) {
        if (use_rstatix) {
          comps <- tryCatch({
            normalize_pcols(
              rstatix::t_test(value ~ bar_group, dmet) |>
                rstatix::adjust_pvalue(method = if (identical(padj_method, "none")) "none" else padj_method) |>
                rstatix::add_significance()
            )
          }, error = function(e) NULL)
          if (!is.null(comps) && nrow(comps) > 0) {
            if (!all(c("group1","group2") %in% names(comps))) {
              lv <- levels(dmet$bar_group); comps$group1 <- lv[1]; comps$group2 <- lv[2]
            }
            comps <- dplyr::transmute(comps,
                                      group1 = as.character(group1), group2 = as.character(group2),
                                      p      = dplyr::coalesce(as.numeric(p), as.numeric(p.adj)),
                                      p.adj  = dplyr::coalesce(as.numeric(p.adj), as.numeric(p)),
                                      p.adj.signif = as.character(p.adj.signif))
          } else comps <- NULL
        }
        if (is.null(comps)) {
          tt <- tryCatch(stats::t.test(value ~ bar_group, dmet), error = function(e) NULL)
          if (!is.null(tt)) {
            lv <- sort(unique(as.character(dmet$bar_group)))
            pval <- tt$p.value
            comps <- tibble::tibble(
              group1 = lv[1], group2 = lv[2],
              p = pval, p.adj = pval,
              p.adj.signif = dplyr::case_when(
                pval < 0.001 ~ "***", pval < 0.01 ~ "**",
                pval < 0.05  ~ "*",  TRUE ~ "ns"))
          }
        }
      } else {
        # ≥3 groups: Tukey HSD
        if (use_rstatix) {
          comps <- tryCatch({
            fit <- stats::aov(value ~ bar_group, data = dmet)
            tk  <- rstatix::tukey_hsd(fit)
            normalize_pcols(tk) |>
              dplyr::transmute(
                group1 = as.character(group1), group2 = as.character(group2),
                p      = dplyr::coalesce(as.numeric(p), as.numeric(p.adj)),
                p.adj  = dplyr::coalesce(as.numeric(p.adj), as.numeric(p)),
                p.adj.signif = as.character(p.adj.signif))
          }, error = function(e) NULL)
        }
        if (is.null(comps)) {
          comps <- tryCatch({
            fit <- stats::aov(value ~ bar_group, data = dmet)
            tk  <- as.data.frame(stats::TukeyHSD(fit, "bar_group")$bar_group)
            tk$contrast <- rownames(tk); rownames(tk) <- NULL
            parts <- strsplit(tk$contrast, "-")
            tibble::tibble(
              group1 = vapply(parts, `[`, "", 1L),
              group2 = vapply(parts, `[`, "", 2L),
              p      = tk[["p adj"]], p.adj = tk[["p adj"]],
              p.adj.signif = dplyr::case_when(
                p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                p.adj < 0.05  ~ "*",  TRUE ~ "ns"))
          }, error = function(e) NULL)
        }
      }
      
      if (is.null(comps) || nrow(comps) == 0) next
      p_col <- if (use_adj && "p.adj" %in% names(comps)) "p.adj" else "p"
      comps <- comps %>%
        dplyr::filter(.data[[p_col]] < filter_alpha) %>%
        dplyr::mutate(`Metabolite name` = met)
      if (nrow(comps) > 0) ann_list[[length(ann_list) + 1]] <- comps
    }
    
    if (length(ann_list) == 0) return(NULL)
    dplyr::bind_rows(ann_list)
  }
  
  posthoc_for_plot <- reactive({
    req(stats_sig_plot_data())
    req(stats_input_long())
    
    sdf    <- stats_sig_plot_data()
    df_all <- stats_input_long()
    
    show_two_way <- any(stats_results_val()$test == "twoway")
    hasA <- any(!is.na(df_all$factorA))
    hasB <- any(!is.na(df_all$factorB))
    df_all <- df_all %>%
      dplyr::mutate(bar_group = if (isTRUE(show_two_way) && hasA && hasB)
        interaction(factorA, factorB, drop = TRUE)
        else factor(group)) %>%
      dplyr::filter(!is.na(bar_group))
    
    feats <- unique(c(levels(sdf$`Metabolite name`), sdf$`Metabolite name`))
    df_feats <- df_all %>%
      dplyr::filter(`Metabolite name` %in% feats) %>%
      dplyr::select(`Metabolite name`, bar_group, value)
    
    .compute_pairwise(
      df_long     = df_feats,
      group_col   = "bar_group",
      padj_method = input$padj_method %||% "BH",
      alpha       = input$alpha %||% 0.05,
      use_adj     = isTRUE(input$use_adj_threshold) &&
        !identical(input$padj_method %||% "BH", "none")
    )
  })
  
  current_plot_index <- reactiveVal(1)
  observeEvent(input$next_plot, {
    n <- nrow(stats_results_val())
    if (n > 0) {
      new_index <- current_plot_index() + 1
      if (new_index > n) new_index <- 1  # wrap around
      current_plot_index(new_index)
    }
  })
  
  # Reset the plot index when key inputs change
  observeEvent(
    list(input$stats_scope,
         input$stats_class,
         input$stats_units,
         input$stats_value_type,
         input$padj_method,
         input$stats_plot_sort),
    {
      current_plot_index(1)
    },
    ignoreInit = TRUE
  )
  
  
  observeEvent(input$prev_plot, {
    n <- nrow(stats_results_val())
    if (n > 0) {
      new_index <- current_plot_index() - 1
      if (new_index < 1) new_index <- n  # wrap around
      current_plot_index(new_index)
    }
  })
  
  output$plotIndexInfo <- renderText({
    res <- stats_results_val()
    if (is.null(res) || nrow(res) == 0) return("")
    paste("Showing", current_plot_index(), "of", nrow(res), "features")
  })
  
  output$statsPlotUI <- renderUI({
    plotOutput("singleStatsPlot", height = "600px")
  })
  
  # ── Shared helper: build one stats bar plot for a single metabolite ──────────
  # Used by both singleStatsPlot (screen) and download_stats_pdf (PDF export).
  # Arguments mirror the current UI inputs so both callers pass the same values.
  .build_stats_barplot <- function(met, idx, n_total, df_long, res_all,
                                   err_type, show_pts, show_two_way,
                                   alpha, padj_method, use_adj, posthoc_choices,
                                   base_size = 14) {
    df_all <- df_long %>% dplyr::filter(`Metabolite name` == met)
    if (nrow(df_all) == 0) return(NULL)
    
    hasA <- any(!is.na(df_all$factorA))
    hasB <- any(!is.na(df_all$factorB))
    if (isTRUE(show_two_way) && hasA && hasB) {
      df_all <- df_all %>% dplyr::mutate(bar_group = interaction(factorA, factorB, drop = TRUE))
    } else {
      df_all <- df_all %>% dplyr::mutate(bar_group = factor(group))
    }
    df_all <- df_all %>% dplyr::filter(!is.na(bar_group))
    if (nrow(df_all) == 0) return(NULL)
    
    df_sum <- df_all %>%
      dplyr::group_by(bar_group) %>%
      dplyr::summarise(
        mean = mean(value, na.rm = TRUE),
        sd   = stats::sd(value, na.rm = TRUE),
        n    = dplyr::n(),
        sem  = sd / sqrt(pmax(n, 1)),
        ci95 = 1.96 * sem,
        .groups = "drop"
      )
    
    err_vec <- switch(err_type,
                      "SEM"    = df_sum$sem,
                      "SD"     = df_sum$sd,
                      "95% CI" = df_sum$ci95,
                      df_sum$sem)
    df_sum$ymin <- df_sum$mean - err_vec
    df_sum$ymax <- df_sum$mean + err_vec
    
    p <- ggplot2::ggplot(df_sum, ggplot2::aes(x = bar_group, y = mean)) +
      ggplot2::geom_col(fill = "#4C78A8", width = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                             width = 0.2, colour = "#333") +
      ggplot2::labs(
        title = if (!is.null(idx)) paste0("(", idx, "/", n_total, ")  ", met) else met,
        x = NULL,
        y = paste0("Mean \u00b1 ", err_type)
      ) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    if (isTRUE(show_pts)) {
      p <- p + ggplot2::geom_jitter(
        data = df_all,
        ggplot2::aes(x = bar_group, y = value),
        width = 0.15, height = 0, alpha = 0.6, size = 2, colour = "#F58518"
      )
    }
    
    # Significance brackets
    if (requireNamespace("ggpubr", quietly = TRUE)) {
      tryCatch({
        k <- dplyr::n_distinct(df_all$bar_group)
        comps <- NULL
        
        if (k == 2) {
          lv <- levels(df_all$bar_group)
          if (length(lv) == 2 && requireNamespace("rstatix", quietly = TRUE)) {
            tt <- try(
              rstatix::t_test(value ~ bar_group, data = df_all, var.equal = FALSE) |>
                rstatix::adjust_pvalue(method = if (identical(padj_method, "none")) "none" else padj_method) |>
                rstatix::add_significance(),
              silent = TRUE)
            if (!inherits(tt, "try-error") && !is.null(tt) && nrow(tt) > 0) {
              comps <- tt |>
                dplyr::mutate(
                  group1  = as.character(group1),
                  group2  = as.character(group2),
                  p       = dplyr::coalesce(as.numeric(.data$p),     as.numeric(.data$p.adj)),
                  p.adj   = dplyr::coalesce(as.numeric(.data$p.adj), as.numeric(.data$p)),
                  p_label = dplyr::coalesce(.data$p.adj.signif, sprintf("p = %.3g", p.adj))
                ) |>
                dplyr::select(group1, group2, p, p.adj, p_label)
            } else comps <- NULL
          }
        } else if (k >= 3) {
          want_tukey    <- "tukey"    %in% posthoc_choices
          want_pairwise <- "pairwise" %in% posthoc_choices
          if (want_tukey && requireNamespace("rstatix", quietly = TRUE)) {
            fit <- try(stats::aov(value ~ bar_group, data = df_all), silent = TRUE)
            tk  <- if (!inherits(fit, "try-error"))
              try(rstatix::tukey_hsd(fit), silent = TRUE) else NULL
            if (!inherits(tk, "try-error") && !is.null(tk) && nrow(tk) > 0) {
              nm_tk <- names(tk)
              if ("adj.p.value" %in% nm_tk && !"p.adj" %in% nm_tk)
                tk <- dplyr::rename(tk, p.adj = adj.p.value)
              if (!"p" %in% names(tk)) tk <- dplyr::mutate(tk, p = p.adj)
              if (!"p.adj.signif" %in% names(tk))
                tk <- dplyr::mutate(tk, p.adj.signif = dplyr::case_when(
                  p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                  p.adj < 0.05  ~ "*",   TRUE ~ "ns"))
              comps <- tk |>
                dplyr::mutate(p_label = dplyr::coalesce(.data$p.adj.signif,
                                                        sprintf("p=%.3g", p.adj))) |>
                dplyr::select(group1, group2, p, p.adj, p_label)
            }
          } else if (want_pairwise && requireNamespace("rstatix", quietly = TRUE)) {
            pw <- try(rstatix::pairwise_t_test(value ~ bar_group, df_all,
                                               p.adjust.method = "holm"), silent = TRUE)
            if (!inherits(pw, "try-error") && !is.null(pw) && nrow(pw) > 0) {
              comps <- pw |>
                dplyr::mutate(
                  p     = dplyr::coalesce(.data$p, .data$p.adj),
                  p.adj = dplyr::coalesce(.data$p.adj, .data$p),
                  p_label = dplyr::coalesce(.data$p.adj.signif, sprintf("p=%.3g", p.adj))
                ) |>
                dplyr::select(group1, group2, p, p.adj, p_label)
            }
          }
        }
        
        if (!is.null(comps) && nrow(comps) > 0) {
          use_adj_f <- isTRUE(use_adj) && !identical(padj_method, "none")
          comps <- comps %>%
            dplyr::mutate(p_keep = dplyr::coalesce(.data$p.adj, .data$p)) %>%
            { if (use_adj_f) dplyr::filter(., p_keep < alpha)
              else dplyr::filter(., p < alpha | p_keep < alpha) }
          if (nrow(comps) > 0) {
            top  <- max(df_sum$ymax, na.rm = TRUE)
            if (!is.finite(top)) top <- 1
            step <- 0.07 * top
            ann  <- comps %>%
              dplyr::filter(group1 %in% levels(df_sum$bar_group),
                            group2 %in% levels(df_sum$bar_group)) %>%
              dplyr::mutate(
                `Metabolite name` = met,
                y.position = top + step * dplyr::row_number()
              )
            if (nrow(ann) > 0) {
              p <- p +
                ggpubr::stat_pvalue_manual(
                  ann, label = "p_label",
                  xmin = "group1", xmax = "group2",
                  y.position = "y.position",
                  tip.length = 0.01, hide.ns = TRUE
                ) +
                ggplot2::expand_limits(y = max(ann$y.position, na.rm = TRUE) * 1.05)
            }
          }
        }
      }, error = function(e) NULL)
    }
    
    p
  }
  
  output$singleStatsPlot <- renderPlot({
    res <- stats_results_val()
    validate(
      need(!is.null(res),
           "Statistics have not been run yet. Please click 'Run statistics' to begin."),
      need(nrow(res) > 0,
           "No significant features to plot. Try adjusting your test or thresholds.")
    )
    idx <- max(1, min(current_plot_index(), nrow(res)))
    met <- res$metabolite[idx]
    p <- .build_stats_barplot(
      met           = met,
      idx           = idx,
      n_total       = nrow(res),
      df_long       = stats_input_long(),
      res_all       = res,
      err_type      = input$stats_error_type %||% "SEM",
      show_pts      = isTRUE(input$stats_plot_points),
      show_two_way  = any(res$test == "twoway"),
      alpha         = input$alpha %||% 0.05,
      padj_method   = input$padj_method %||% "BH",
      use_adj       = isTRUE(input$use_adj_threshold),
      posthoc_choices = input$posthoc %||% character(0),
      base_size     = 14
    )
    validate(need(!is.null(p), paste0("No data for metabolite '", met, "'.")))
    p
  })
  
  # ---------- Heatmap (significant features) ----------
  sig_heatmap_data <- reactive({
    # 1) Significant features
    res <- stats_results_val()
    validate(need(!is.null(res) && nrow(res) > 0,
                  "No results available. Compute statistics first."))
    
    alpha   <- input$alpha %||% 0.05
    use_adj <- isTRUE(input$use_adj_threshold) && !identical(input$padj_method, "none")
    
    is_sig <- if (use_adj) (!is.na(res$p_adj) & (res$p_adj < alpha)) else (!is.na(res$p) & (res$p < alpha))
    res_sig <- res[is_sig, , drop = FALSE]
    validate(need(nrow(res_sig) > 0,
                  "No statistically significant features under current thresholds."))
    
    res_sig <- res_sig %>% dplyr::arrange(dplyr::coalesce(p_adj, p), p)
    feats   <- head(res_sig$metabolite, input$hm_topn %||% 50)
    validate(need(length(feats) >= 2, "Need at least 2 significant features for a heatmap."))
    
    # 2) Tested values aligned to Statistics settings
    df_all <- stats_input_long() %>%
      dplyr::filter(`Metabolite name` %in% feats)
    
    # 3) Re-assert grouping via central reactive
    fn  <- get_active_grouping()
    res <- fn(df_all$sample)
    df_all <- df_all %>%
      dplyr::mutate(group = res$group) %>%
      dplyr::filter(!is.na(group))
    
    # --- New: decide if two-way ANOVA is actually active (among significant features)
    twoway_active <- any(res_sig$test == "twoway")
    
    # 4) Factor A/B from central grouping
    if (isTRUE(twoway_active)) {
      df_all$factorA <- res$factorA[match(df_all$sample, df_all$sample)]
      df_all$factorB <- res$factorB[match(df_all$sample, df_all$sample)]
    } else {
      df_all$factorA <- NA_character_
      df_all$factorB <- NA_character_
    }
    
    # helper: replicate suffix detection for samples
    looks_like_replicate_suffix <- function(x) {
      grepl("([_\\-])(rep|plate|r)[_\\-]?[0-9]+$", x, ignore.case = TRUE)
    }
    
    # 5) Build matrix: group means vs sample-level
    if (isTRUE(input$hm_use_group_means)) {
      # Decide A×B strictly: only if two-way is active AND both factors have >=2 levels
      nA <- dplyr::n_distinct(stats::na.omit(df_all$factorA))
      nB <- dplyr::n_distinct(stats::na.omit(df_all$factorB))
      rep_suffix_present <- any(looks_like_replicate_suffix(df_all$sample))
      use_axb <- isTRUE(twoway_active) &&
        isTRUE(input$hm_use_axb) &&
        (nA >= 2) && (nB >= 2) &&
        !rep_suffix_present
      
      if (use_axb) {
        df_all <- df_all %>% dplyr::mutate(col_group = interaction(factorA, factorB, drop = TRUE))
        col_lab <- "A:B"
      } else {
        df_all <- df_all %>% dplyr::mutate(col_group = factor(group))
        col_lab <- "Group"
      }
      
      mat_df <- df_all %>%
        dplyr::group_by(`Metabolite name`, col_group) %>%
        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = col_group, values_from = value) %>%
        as.data.frame()
      
      anno_col <- data.frame(mat_df[, -1, drop = FALSE] %>% names(), check.names = FALSE)
      colnames(anno_col) <- col_lab
      rownames(anno_col) <- anno_col[[col_lab]]
      
    } else {
      # --- Sample-level heatmap: NEVER collapse by A×B, only per-sample columns
      df_ann <- df_all %>%
        dplyr::distinct(sample, group, factorA, factorB) %>%
        dplyr::arrange(sample)
      
      mat_df <- df_all %>%
        dplyr::group_by(`Metabolite name`, sample) %>%
        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = sample, values_from = value) %>%
        as.data.frame()
      
      all_cols <- colnames(mat_df)[-1]
      df_ann   <- df_ann[match(all_cols, df_ann$sample), , drop = FALSE]
      
      # Only annotate Group; include FactorA/B columns ONLY if they exist and are non-NA (after two-way guard)
      grp <- df_ann$group
      anno_col <- data.frame(Group = grp, check.names = FALSE)
      
      if (isTRUE(twoway_active) && any(!is.na(df_ann$factorA))) anno_col$FactorA <- df_ann$factorA
      if (isTRUE(twoway_active) && any(!is.na(df_ann$factorB))) anno_col$FactorB <- df_ann$factorB
      
      rownames(anno_col) <- all_cols
      
      # Drop empty annotation columns
      if ("FactorA" %in% names(anno_col) && all(is.na(anno_col$FactorA))) anno_col$FactorA <- NULL
      if ("FactorB" %in% names(anno_col) && all(is.na(anno_col$FactorB))) anno_col$FactorB <- NULL
    }
    
    # 6) Matrix + impute missing by row median
    validate(need(ncol(mat_df) >= 2, "No matrix columns found."))
    rn  <- mat_df$`Metabolite name`
    mat <- as.matrix(mat_df[, -1, drop = FALSE])
    rownames(mat) <- rn
    
    if (anyNA(mat)) mat <- impute_row_median(mat)
    
    # 7) Row scaling
    scale_mode <- input$hm_scale %||% "z"
    if (identical(scale_mode, "z")) {
      mat <- t(scale(t(mat), center = TRUE, scale = TRUE))
      mat[!is.finite(mat)] <- 0
    } else if (identical(scale_mode, "center")) {
      mat <- t(scale(t(mat), center = TRUE, scale = FALSE))
      mat[!is.finite(mat)] <- 0
    }
    
    # 8) Return payload
    list(
      mat      = mat,
      anno_col = anno_col,
      n_feats  = nrow(mat),
      n_cols   = ncol(mat),
      class    = input$stats_class
    )
  })
  
  
  observeEvent(sig_heatmap_data(), {
    hd <- sig_heatmap_data()
    # Rough scaling: 0.25in per feature, capped
    w_default <- 8
    h_default <- max(6, min(25, 0.25 * hd$n_feats))
    updateNumericInput(session, "hm_export_width",  value = round(w_default * 300))
    updateNumericInput(session, "hm_export_height", value = round(h_default * 300))
  })
  
  
  
  
  output$sigHeatmapInfo <- renderText({
    res <- stats_results_val()
    if (is.null(res) || nrow(res) == 0) return("")
    tests <- unique(res$test)
    if (all(tests == "ttest")) {
      "Note: column clustering is disabled for t-test results (only 2 columns)."
    } else ""
  })
  
  output$sigHeatmap <- renderPlot({
    hd  <- sig_heatmap_data()
    mat <- hd$mat
    anno_col   <- hd$anno_col
    cols       <- .heatmap_palette(input$hm_palette %||% "viridis")
    dist_opt   <- input$hm_dist    %||% "correlation"
    method_opt <- input$hm_linkage %||% "complete"
    pheatmap::pheatmap(
      mat,
      color            = cols,
      cluster_rows     = isTRUE(input$hm_cluster_rows),
      cluster_cols     = isTRUE(input$hm_cluster_cols) && ncol(mat) >= 3,
      clustering_distance_rows = dist_opt,
      clustering_distance_cols = dist_opt,
      clustering_method = method_opt,
      annotation_col   = if (ncol(anno_col) > 0) anno_col else NULL,
      show_rownames    = isTRUE(input$hm_show_rownames),
      show_colnames    = isTRUE(input$hm_show_colnames),
      fontsize_row     = 8,
      fontsize_col     = 9,
      border_color     = NA,
      main             = paste0("Significant features heatmap — ", hd$class)
    )
  }, height = function() {
    hd <- tryCatch(sig_heatmap_data(), error = function(e) NULL)
    if (is.null(hd)) return(400)
    pmin(1600, pmax(400, 18 * (hd$n_feats %||% 10)))
  })
  
  
  
  output$download_sig_hm_matrix <- downloadHandler(
    filename = function() {
      paste0("heatmap_significant_", input$stats_class, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      hd <- sig_heatmap_data()
      df <- cbind(Feature = rownames(hd$mat), as.data.frame(hd$mat, check.names = FALSE))
      readr::write_csv(df, file)
    }
  )
  
  
  
  .hm_export_dims <- function() {
    px_w  <- input$hm_export_width    %||% 1100
    px_h  <- input$hm_export_height   %||% 900
    dpi   <- input$hm_export_dpi      %||% 300
    scale <- input$hm_export_scale    %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  .build_sig_heatmap <- function(fsz = 10) {
    hd <- sig_heatmap_data()
    if (is.null(hd)) return(NULL)
    pheatmap::pheatmap(
      hd$mat,
      color = .heatmap_palette(input$hm_palette %||% "viridis"),
      cluster_rows = isTRUE(input$hm_cluster_rows),
      cluster_cols = isTRUE(input$hm_cluster_cols) && ncol(hd$mat) >= 3,
      clustering_distance_rows = input$hm_dist    %||% "correlation",
      clustering_distance_cols = input$hm_dist    %||% "correlation",
      clustering_method        = input$hm_linkage %||% "complete",
      annotation_col = if (ncol(hd$anno_col) > 0) hd$anno_col else NULL,
      show_rownames = isTRUE(input$hm_show_rownames),
      show_colnames = isTRUE(input$hm_show_colnames),
      fontsize_row = fsz, fontsize_col = fsz + 1,
      border_color = NA,
      main = paste0("Significant features heatmap — ", hd$class),
      silent = TRUE
    )
  }
  
  output$download_hm_png <- downloadHandler(
    filename = function() paste0("heatmap_significant_", isolate(input$stats_class), "_", Sys.Date(), ".png"),
    content = function(file) {
      dims <- isolate(.hm_export_dims())
      fsz  <- isolate(input$hm_export_fontsize %||% 10)
      hm   <- isolate(.build_sig_heatmap(fsz))
      validate(need(!is.null(hm), "No heatmap to export. Run statistics first."))
      grDevices::png(file,
                     width  = round(dims$w * dims$dpi),
                     height = round(dims$h * dims$dpi),
                     res    = dims$dpi)
      on.exit(grDevices::dev.off(), add = TRUE)
      grid::grid.newpage()
      grid::grid.draw(hm$gtable)
    }
  )
  
  output$download_hm_svg <- downloadHandler(
    filename = function() paste0("heatmap_significant_", isolate(input$stats_class), "_", Sys.Date(), ".svg"),
    content = function(file) {
      dims <- isolate(.hm_export_dims())
      fsz  <- isolate(input$hm_export_fontsize %||% 10)
      hm   <- isolate(.build_sig_heatmap(fsz))
      validate(need(!is.null(hm), "No heatmap to export. Run statistics first."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      grid::grid.newpage()
      grid::grid.draw(hm$gtable)
    }
  )
  
  
  
  
  # ================================
  # Lipid Set Enrichment (LSEA) — server block
  # ================================
  
  # ---- Parsing utilities (sum vs resolved FA) ----
  split_sum_resolved <- function(name) {
    # Accept either "PC 34:0 | PC 16:0_18:1" or "PC 34:0\nPC 16:0_18:1"
    parts <- unlist(strsplit(name, "\\||\n", perl = TRUE))
    parts <- stringr::str_trim(parts)
    parts <- parts[nchar(parts) > 0]
    sum_part <- if (length(parts) >= 1) parts[1] else NA_character_
    resolved_part <- if (length(parts) >= 2) parts[2] else NA_character_
    list(sum = sum_part, resolved = resolved_part)
  }
  
  # Extract FA tokens "16:0", "18:1" etc. from the resolved part (after '|').
  # Extract FA tokens from the resolved part (after '|'), including ether "O-XX:Y" tokens.
  # Accept "_" "/" "," delimiters; strip class token if present (e.g., "PC ").
  extract_fa_list <- function(resolved_part) {
    if (is.na(resolved_part) || resolved_part == "") return(character(0))
    rhs <- stringr::str_remove(resolved_part, "^[A-Za-z0-9]+\\s+")  # drop class token like "PC "
    rhs <- gsub("[()\\[\\]\\s]", "", rhs)                           # clean brackets/spaces
    toks <- unlist(strsplit(rhs, "[_/ ,]"))
    toks <- toks[nchar(toks) > 0]
    # Allow plain "XX:Y" and ether "O-XX:Y" (optionally with trailing ";keto/oxo counts")
    keep <- grepl("^(?:O-)?\\d+:(\\d+)(?:;\\d+)*$", toks, perl = TRUE)
    toks[keep]
  }
  
  
  parse_total_c_db <- function(sum_name) {
    m <- stringr::str_match(sum_name, "(\\d+):(\\d+)")
    if (!any(is.na(m))) return(list(total_c = as.numeric(m[2]), total_db = as.numeric(m[3])))
    list(total_c = NA_real_, total_db = NA_real_)
  }
  
  sat_bucket <- function(total_db) {
    dplyr::case_when(
      is.na(total_db) ~ NA_character_,
      total_db == 0   ~ "SFA",
      total_db == 1   ~ "MUFA",
      total_db >= 2   ~ "PUFA",
      TRUE            ~ NA_character_
    )
  }
  
  # ---- Base annotations used by both parsing modes ----
  lipid_annotations <- reactive({
    req(bg_norm_long_avg())
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::distinct(`Metabolite name`, plot_class)
    
    # Split SUM | RESOLVED
    split_sum_resolved <- function(name) {
      parts <- unlist(strsplit(name, "\\||\n", perl = TRUE))
      parts <- stringr::str_trim(parts)
      parts <- parts[nchar(parts) > 0]
      sum_part <- if (length(parts) >= 1) parts[1] else NA_character_
      resolved_part <- if (length(parts) >= 2) parts[2] else NA_character_
      list(sum = sum_part, resolved = resolved_part)
    }
    
    sp <- lapply(df$`Metabolite name`, split_sum_resolved)
    sum_names <- vapply(sp, `[[`, "", "sum")
    res_parts <- vapply(sp, `[[`, "", "resolved")
    
    # Parse total carbons / total double bonds from SUM side
    parse_total_c_db <- function(sum_name) {
      m <- stringr::str_match(sum_name, "(\\d+):(\\d+)")
      if (!any(is.na(m))) {
        list(total_c = as.numeric(m[2]), total_db = as.numeric(m[3]))
      } else {
        list(total_c = NA_real_, total_db = NA_real_)
      }
    }
    sums <- lapply(sum_names, parse_total_c_db)
    df$total_c  <- vapply(sums, function(x) x$total_c,  numeric(1))
    df$total_db <- vapply(sums, function(x) x$total_db, numeric(1))
    
    # Resolved FA tokens (includes O-XX:Y)
    df$fa_list <- lapply(res_parts, extract_fa_list)
    
    # Ether class tagging from plot_class
    df$is_ether_class <- stringr::str_detect(df$plot_class, "-O$")
    df$base_class <- sub("-O$", "", df$plot_class)
    
    df
  })
  
  # ---- Ether subclass mapping (plasmanyl / plasmenyl) from resolved FA ----
  ether_subclass_lookup <- reactive({
    ann <- lipid_annotations()
    validate(need(nrow(ann) > 0, "No lipid annotations available."))
    
    # Pull the DB count from the first O-XX:YY token in the resolved FA list (if present)
    get_o_db <- function(fa_vec) {
      if (length(fa_vec) == 0) return(NA_real_)
      o_tokens <- fa_vec[grepl("^O-\\d+:\\d+", fa_vec)]
      if (length(o_tokens) == 0) return(NA_real_)
      m <- stringr::str_match(o_tokens[1], "^O-(\\d+):(\\d+)")
      suppressWarnings(as.numeric(m[, 3]))
    }
    
    o_db <- vapply(ann$fa_list, get_o_db, numeric(1))
    
    subclass <- ifelse(
      ann$plot_class %in% c("PE-O", "PC-O") & !is.na(o_db),
      ifelse(o_db == 0, "plasmanyl",
             ifelse(o_db >= 1, "plasmenyl", NA_character_)),
      NA_character_
    )
    
    # IMPORTANT: keep the space in "Metabolite name"
    data.frame(
      `Metabolite name` = ann$`Metabolite name`,
      plot_class        = ann$plot_class,
      ether_subclass    = subclass,
      stringsAsFactors  = FALSE,
      check.names       = FALSE   # <— preserve "Metabolite name" exactly
    )
  })
  
  
  
  # ---- Sets from sum composition (adds class membership) ----
  # --- REPLACE: build_sets_sum() ---
  build_sets_sum <- reactive({
    ann <- lipid_annotations()
    validate(need(nrow(ann) > 0, "No lipid annotations available."))
    
    sets <- list()
    
    # 1) Global: class membership
    for (cl in unique(ann$plot_class)) {
      sets[[paste0("All: class ", cl)]] <- ann$`Metabolite name`[ann$plot_class == cl]
    }
    
    # 2) Global: exact total carbons (C)
    for (tc in sort(unique(stats::na.omit(ann$total_c)))) {
      idx <- which(ann$total_c == tc)
      if (length(idx)) sets[[paste0("All: totalC ", tc)]] <- ann$`Metabolite name`[idx]
    }
    
    # 3) Global: exact total double bonds (DB)
    for (db in sort(unique(stats::na.omit(ann$total_db)))) {
      idx <- which(ann$total_db == db)
      if (length(idx)) sets[[paste0("All: totalDB ", db)]] <- ann$`Metabolite name`[idx]
    }
    
    # 4) Global: combined exact totals C:DB
    key <- paste0("C", ann$total_c, ":DB", ann$total_db)
    for (lev in sort(unique(key[is.finite(ann$total_c) & is.finite(ann$total_db)]))) {
      idx <- which(key == lev)
      if (length(idx)) sets[[paste0("All: ", lev)]] <- ann$`Metabolite name`[idx]
    }
    
    # 5) Per-class variants for totalC, totalDB, and C:DB
    for (cl in unique(ann$plot_class)) {
      sub <- ann[ann$plot_class == cl, , drop = FALSE]
      
      for (tc in sort(unique(stats::na.omit(sub$total_c)))) {
        idx <- which(sub$total_c == tc)
        if (length(idx)) sets[[paste0(cl, ": totalC ", tc)]] <- sub$`Metabolite name`[idx]
      }
      
      for (db in sort(unique(stats::na.omit(sub$total_db)))) {
        idx <- which(sub$total_db == db)
        if (length(idx)) sets[[paste0(cl, ": totalDB ", db)]] <- sub$`Metabolite name`[idx]
      }
      
      sub_key <- paste0("C", sub$total_c, ":DB", sub$total_db)
      keep_keys <- sub_key[is.finite(sub$total_c) & is.finite(sub$total_db)]
      for (lev in sort(unique(keep_keys))) {
        idx <- which(sub_key == lev)
        if (length(idx)) sets[[paste0(cl, ": ", lev)]] <- sub$`Metabolite name`[idx]
      }
    }
    
    # Clean, de-dup, enforce min set size
    sets <- lapply(sets, function(v) sort(unique(stats::na.omit(v))))
    minsz <- input$lsea_min_set; if (is.null(minsz)) minsz <- 3
    sets[vapply(sets, length, integer(1)) >= minsz]
  })
  
  
  # ---- Sets from resolved FA (FA identity + FA categories; adds class membership)
  #      NOW WITH length classes: LCFA (13–21), VLCFA (22–25), ULCFA (>25) ----
  # --- REPLACE: build_sets_resolved() ---
  build_sets_resolved <- reactive({
    ann <- lipid_annotations()
    validate(need(nrow(ann) > 0, "No lipid annotations available."))
    
    # Keep species with resolved FA lists or (optionally) append SUM sets later if fallback is ON
    has_res <- vapply(ann$fa_list, function(x) length(x) > 0, logical(1))
    keep_idx <- has_res | isTRUE(input$fa_resolved_fallback_sum)
    ann <- ann[keep_idx, , drop = FALSE]
    
    sets <- list()
    
    # ------------------------------------------------------------
    # A) CLASS-LEVEL ENRICHMENT, including explicit Ether vs Diacyl
    # ------------------------------------------------------------
    
    # Global class membership (unchanged)
    for (cl in unique(ann$plot_class)) {
      sets[[paste0("All: class ", cl)]] <- ann$`Metabolite name`[ann$plot_class == cl]
    }
    
    # NEW: Global Ether vs Diacyl sets using class naming
    if (any(ann$is_ether_class, na.rm = TRUE)) {
      sets[["All: class Ether (−O)"]]  <- ann$`Metabolite name`[ann$is_ether_class %in% TRUE]
    }
    if (any(!ann$is_ether_class, na.rm = TRUE)) {
      sets[["All: class Diacyl (non‑O)"]] <- ann$`Metabolite name`[ann$is_ether_class %in% FALSE]
    }
    
    # NEW: Per-base-class Ether vs Diacyl sets (e.g., PC: Ether vs PC: Diacyl)
    for (bc in unique(ann$base_class)) {
      sub <- ann[ann$base_class == bc, , drop = FALSE]
      if (nrow(sub) == 0) next
      
      ether_idx  <- which(sub$is_ether_class %in% TRUE)
      diacyl_idx <- which(sub$is_ether_class %in% FALSE)
      
      if (length(ether_idx))
        sets[[paste0(bc, ": class Ether (−O)")]]  <- sub$`Metabolite name`[ether_idx]
      if (length(diacyl_idx))
        sets[[paste0(bc, ": class Diacyl (non‑O)")]] <- sub$`Metabolite name`[diacyl_idx]
    }
    
    # ------------------------------------------------------------
    # B) FA-RESOLVED ENRICHMENT: Identity + Categories (including Ether FA subtypes)
    # ------------------------------------------------------------
    if (any(has_res)) {
      # Expand to long (lipid, class, base_class, FA)
      all_pairs <- do.call(
        rbind,
        lapply(seq_len(nrow(ann)), function(i) {
          if (length(ann$fa_list[[i]]) == 0) return(NULL)
          data.frame(
            lipid      = ann$`Metabolite name`[i],
            plot_class = ann$plot_class[i],
            base_class = ann$base_class[i],
            fa         = ann$fa_list[[i]],
            stringsAsFactors = FALSE
          )
        })
      )
      
      if (!is.null(all_pairs) && nrow(all_pairs) > 0) {
        # (1) FA identity sets
        for (fa in unique(all_pairs$fa)) {
          sets[[paste0("All: FA ", fa)]] <- unique(all_pairs$lipid[all_pairs$fa == fa])
        }
        for (cl in unique(all_pairs$plot_class)) {
          sub <- all_pairs[all_pairs$plot_class == cl, , drop = FALSE]
          for (fa in unique(sub$fa)) {
            sets[[paste0(cl, ": FA ", fa)]] <- unique(sub$lipid[sub$fa == fa])
          }
        }
        
        # Build FA meta (len, db, saturation, bins), supporting O- tokens
        fa_info <- unique(all_pairs["fa"])
        m <- stringr::str_match(fa_info$fa, "^(?:O-)?(\\d+):(\\d+)")
        fa_info$len <- suppressWarnings(as.numeric(m[, 2]))   # carbon number
        fa_info$db  <- suppressWarnings(as.numeric(m[, 3]))   # double bonds
        
        # SFA/MUFA/PUFA
        sat_bucket <- function(total_db) {
          dplyr::case_when(
            is.na(total_db) ~ NA_character_,
            total_db == 0   ~ "SFA",
            total_db == 1   ~ "MUFA",
            total_db >= 2   ~ "PUFA",
            TRUE            ~ NA_character_
          )
        }
        fa_info$sat <- sat_bucket(fa_info$db)
        
        # FA length bins (fine)
        fa_info$len_bin <- cut(
          fa_info$len,
          breaks = c(-Inf, 16, 18, 20, 22, 24, Inf),
          labels = c("≤16", "18", "20", "22", "24", "≥26")
        )
        
        # Length classes (LCFA/VLCFA/ULCFA)
        fa_info$len_class <- dplyr::case_when(
          is.na(fa_info$len)                    ~ NA_character_,
          fa_info$len >= 13 & fa_info$len <= 21 ~ "LCFA (13–21)",
          fa_info$len >= 22 & fa_info$len <= 25 ~ "VLCFA (22–25)",
          fa_info$len > 25                      ~ "ULCFA (>25)",
          TRUE                                  ~ "Other (≤12)"
        )
        
        # (NEW) Ether FA subtype tagging from FA tokens (O-XX:0 vs O-XX:1)
        fa_info$is_ether_fa <- grepl("^O-\\d+:(\\d+)", fa_info$fa)
        fa_info$ether_type <- NA_character_
        fa_info$ether_type[fa_info$is_ether_fa & fa_info$db == 0] <- "plasmanyl"
        fa_info$ether_type[fa_info$is_ether_fa & fa_info$db >= 1] <- "plasmenyl"
        # Note: Only O-XX:0 or O-XX:1 are mapped as requested; others remain NA.
        
        # Map FA token -> lipids
        fa2lip <- split(all_pairs$lipid, all_pairs$fa)
        
        # (2) Global FA categories: SFA/MUFA/PUFA (resolved mode only)
        for (k in na.omit(unique(fa_info$sat))) {
          famem <- fa_info$fa[fa_info$sat == k]
          if (length(famem)) {
            sets[[paste0("All: FAcat ", k)]] <- unique(unlist(fa2lip[famem], use.names = FALSE))
          }
        }
        
        # (3) Global FA length bins
        for (lb in na.omit(levels(fa_info$len_bin))) {
          famem <- fa_info$fa[as.character(fa_info$len_bin) == lb]
          if (length(famem)) {
            sets[[paste0("All: FAlen ", lb)]] <- unique(unlist(fa2lip[famem], use.names = FALSE))
          }
        }
        
        # (4) Global length classes
        for (lc in na.omit(unique(fa_info$len_class))) {
          famem <- fa_info$fa[fa_info$len_class == lc]
          if (length(famem)) {
            sets[[paste0("All: ", lc)]] <- unique(unlist(fa2lip[famem], use.names = FALSE))
          }
        }
        
        # (NEW) (5) Global Ether FA subtypes (plasmanyl / plasmenyl)
        for (etype in c("plasmanyl", "plasmenyl")) {
          famem <- fa_info$fa[fa_info$ether_type == etype]
          if (length(famem)) {
            sets[[paste0("All: Ether FA ", etype)]] <- unique(unlist(fa2lip[famem], use.names = FALSE))
          }
        }
        
        # Per-class (by base_class) variants for FAcat / FAlen / length classes / Ether FA types
        for (bc in unique(all_pairs$base_class)) {
          sub_pairs <- all_pairs[all_pairs$base_class == bc, , drop = FALSE]
          if (nrow(sub_pairs) == 0) next
          sub_map <- split(sub_pairs$lipid, sub_pairs$fa)
          
          # SFA/MUFA/PUFA per base class
          for (k in na.omit(unique(fa_info$sat))) {
            famem <- fa_info$fa[fa_info$sat == k]
            lipids <- unique(unlist(sub_map[intersect(names(sub_map), famem)], use.names = FALSE))
            if (length(lipids)) sets[[paste0(bc, ": FAcat ", k)]] <- lipids
          }
          
          # FAlen bins per base class
          for (lb in na.omit(levels(fa_info$len_bin))) {
            famem <- fa_info$fa[as.character(fa_info$len_bin) == lb]
            lipids <- unique(unlist(sub_map[intersect(names(sub_map), famem)], use.names = FALSE))
            if (length(lipids)) sets[[paste0(bc, ": FAlen ", lb)]] <- lipids
          }
          
          # Length classes per base class
          for (lc in na.omit(unique(fa_info$len_class))) {
            famem <- fa_info$fa[fa_info$len_class == lc]
            lipids <- unique(unlist(sub_map[intersect(names(sub_map), famem)], use.names = FALSE))
            if (length(lipids)) sets[[paste0(bc, ": ", lc)]] <- lipids
          }
          
          # (NEW) Ether FA subtypes per base class
          for (etype in c("plasmanyl", "plasmenyl")) {
            famem <- fa_info$fa[fa_info$ether_type == etype]
            lipids <- unique(unlist(sub_map[intersect(names(sub_map), famem)], use.names = FALSE))
            if (length(lipids)) sets[[paste0(bc, ": Ether FA ", etype)]] <- lipids
          }
        }
      }
    }
    
    # Optional fallback: append SUM-mode sets (totalC/totalDB/C:DB) if user enables it
    if (isTRUE(input$fa_resolved_fallback_sum)) {
      sets <- c(sets, build_sets_sum())
    }
    
    # Clean & size filter
    sets <- lapply(sets, function(v) sort(unique(stats::na.omit(v))))
    minsz <- input$lsea_min_set; if (is.null(minsz)) minsz <- 3
    sets[vapply(sets, length, integer(1)) >= minsz]
  })
  
  
  
  
  # ---- Switch sets by mode and include/exclude per-class ----
  lsea_sets <- reactive({
    s <- if (identical(input$fa_parsing_mode, "resolved")) build_sets_resolved() else build_sets_sum()
    if (!isTRUE(input$lsea_per_class)) {
      s <- s[startsWith(names(s), "All:")]  # keep only global sets
    }
    validate(need(length(s) > 0, "No lipid sets available under current options."))
    s
  })
  
  # ---- Base long DF aligned with your Statistics settings (for per-group) ----
  enrich_base_df <- reactive({
    req(bg_norm_long_avg())
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample)))
    
    # Blank handling
    df <- filter_blanks(df, input$exclude_blank_stats, sample_col = "sample_norm", exact = FALSE)
    # IQC handling
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    
    # Measure/units (unchanged) ...
    measure_col <- if (identical(input$stats_value_type, "value_bs")) "value_bs" else "norm"
    if (identical(input$stats_units, "percent")) {
      totals <- df %>%
        dplyr::group_by(sample, plot_class) %>%
        dplyr::summarise(class_total = sum(.data[[measure_col]], na.rm = TRUE), .groups = "drop")
      df <- df %>%
        dplyr::left_join(totals, by = c("sample","plot_class")) %>%
        dplyr::mutate(value = dplyr::if_else(class_total > 0, 100 * .data[[measure_col]] / class_total, NA_real_)) %>%
        dplyr::select(sample, plot_class, `Metabolite name`, value)
    } else {
      df <- df %>% dplyr::transmute(sample, plot_class, `Metabolite name`, value = .data[[measure_col]])
    }
    
    # Scope: current class vs all
    if (identical(input$lsea_univ_scope, "current")) {
      req(input$stats_class)
      df <- df %>% dplyr::filter(plot_class == input$stats_class)
      validate(need(nrow(df) > 0,
                    sprintf("No rows found for current class '%s' under the present filters.",
                            input$stats_class)))
    }
    
    fn  <- get_active_grouping()
    res <- fn(df$sample)
    df$group <- res$group
    
    # Respect the group selection from the Statistics tab
    sel_grps <- input$stats_selected_groups
    if (!is.null(sel_grps) && length(sel_grps) > 0) {
      df <- df %>% dplyr::filter(group %in% sel_grps)
    }
    
    # Ensure at least 2 groups for “one-vs-rest”
    if (identical(input$lsea_scope, "group")) {
      validate(need(dplyr::n_distinct(df$group) >= 2,
                    "Per-group mode requires ≥2 groups in the selected universe."))
    }
    df
  })
  
  
  # ---- Populate the group selector (server-side) ----
  available_groups_for_enrich <- reactive({
    req(enrich_base_df())
    unique(enrich_base_df()$group)
  })
  observeEvent(enrich_base_df(), {
    choices <- available_groups_for_enrich()
    sel <- if (length(choices) > 0) head(choices, 1) else NULL
    updateSelectizeInput(
      session, "lsea_group_select",
      choices = choices, selected = sel, server = TRUE
    )
  }, ignoreInit = FALSE)
  
  
  # ---- Universe (features considered) ----
  lsea_universe <- reactive({
    if (identical(input$lsea_scope, "overall")) {
      df <- stats_input_long();  validate(need(nrow(df) > 0, "No tested features found."))
      unique(df$`Metabolite name`)
    } else {
      df <- enrich_base_df();    validate(need(nrow(df) > 0, "No tested features found."))
      unique(df$`Metabolite name`)
    }
  })
  
  # ---- Hits for ORA ----
  lsea_hits_overall <- reactive({
    res <- stats_results_val(); validate(need(nrow(res) > 0, "No significant results."))
    unique(res$metabolite)
  })
  
  # One-vs-rest test per group → p, log2FC per lipid
  compute_ovr_stats <- function(df, groups) {
    feats <- unique(df$`Metabolite name`)
    out <- lapply(groups, function(g){
      dd <- df %>% dplyr::mutate(is_g = group == g)
      rows <- lapply(feats, function(m){
        x <- dd %>% dplyr::filter(`Metabolite name` == m)
        if (length(unique(x$is_g)) < 2) return(NULL)
        tt <- try(stats::t.test(value ~ is_g, data = x), silent = TRUE)
        if (inherits(tt, "try-error") || is.null(tt)) return(NULL)
        m1 <- mean(x$value[x$is_g], na.rm = TRUE); m0 <- mean(x$value[!x$is_g], na.rm = TRUE)
        pseudo <- 1e-9
        data.frame(
          metabolite = m,
          for_group  = g,
          p          = as.numeric(tt$p.value),
          log2FC     = log2((m1 + pseudo) / (m0 + pseudo)),
          stringsAsFactors = FALSE
        )
      })
      res <- dplyr::bind_rows(rows)
      if (!is.null(res) && nrow(res) > 0) res$p_adj <- p.adjust(res$p, "BH")
      res
    })
    dplyr::bind_rows(out)
  }
  
  lsea_hits_per_group <- reactive({
    df <- enrich_base_df()
    groups <- input$lsea_group_select
    validate(need(length(groups) > 0, "Select at least one group."))
    
    # one-vs-rest t-tests (p, log2FC per lipid per group)
    ovr <- compute_ovr_stats(df, groups)
    validate(need(!is.null(ovr) && nrow(ovr) > 0, "No one-vs-rest results."))
    
    alpha   <- input$alpha %||% 0.05
    use_adj <- isTRUE(input$use_adj_threshold) && !identical(input$padj_method, "none")
    
    hits_list <- lapply(groups, function(g) {
      gg <- ovr %>% dplyr::filter(for_group == g)
      if (nrow(gg) == 0) return(NULL)
      
      gg <- if (use_adj) dplyr::filter(gg, p_adj < alpha) else dplyr::filter(gg, p < alpha)
      if (isTRUE(input$ora_up_only)) gg <- dplyr::filter(gg, log2FC > 0)
      
      if (nrow(gg) == 0) return(NULL)  # <-- IMPORTANT: nothing to contribute for this group
      
      gg %>% dplyr::transmute(for_group = g, metabolite)
    })
    
    res <- dplyr::bind_rows(hits_list)
    validate(need(!is.null(res) && nrow(res) > 0,
                  "No hits for selected group(s) under current thresholds."))
    res
  })
  
  
  # ---- Ranks for FGSEA ----
  lsea_rank_overall <- reactive({
    res <- stats_results_val(); validate(need(nrow(res) > 0, "No results to build a ranking."))
    pcol <- if ("p_adj" %in% names(res)) "p_adj" else if ("p.adj" %in% names(res)) "p.adj" else "p"
    p_use <- res[[pcol]];  p_use[!is.finite(p_use)] <- NA_real_;  p_use <- pmin(p_use, 1, na.rm = TRUE)
    
    if ("log2FC" %in% names(res) && any(is.finite(res$log2FC))) {
      w <- 0.05
      r <- sign(res$log2FC) * (abs(res$log2FC) + w * -log10(p_use + 1e-12))
    } else {
      r <- -log10(p_use + 1e-12)
    }
    names(r) <- res$metabolite
    if (any(duplicated(r))) { set.seed(123); r <- r + rnorm(length(r), sd = 1e-9) }
    if (sum(is.finite(r)) < 30) {
      showNotification("FGSEA ranking is short. Consider turning ON 'Show all rows' in Statistics.", type = "warning", duration = 7)
    }
    r
  })
  
  lsea_rank_per_group <- reactive({
    df <- enrich_base_df()
    groups <- input$lsea_group_select
    validate(need(length(groups) > 0, "Select at least one group."))
    ovr <- compute_ovr_stats(df, groups)
    validate(need(nrow(ovr) > 0, "No one-vs-rest results."))
    
    ranks <- lapply(groups, function(g) {
      gg <- ovr %>% dplyr::filter(for_group == g)
      p_use <- gg$p; p_use[!is.finite(p_use)] <- NA_real_; p_use <- pmin(p_use, 1, na.rm = TRUE)
      w <- 0.05
      r <- sign(gg$log2FC) * (abs(gg$log2FC) + w * -log10(p_use + 1e-12))
      names(r) <- gg$metabolite
      if (any(duplicated(r))) { set.seed(123); r <- r + rnorm(length(r), sd = 1e-9) }
      r
    })
    names(ranks) <- groups
    ranks
  })
  
  # ---- ORA & FGSEA engines ----
  # ---- OVER-REPRESENTATION ANALYSIS (ORA) ----
  run_ora <- function(hit_list, universe, sets, minSetSize = 1) {
    results <- lapply(names(sets), function(sn) {
      S <- intersect(sets[[sn]], universe)
      if (length(S) < minSetSize) return(NULL)  # allow singletons if minSetSize = 1
      k <- sum(hit_list %in% S)
      m <- length(S)
      n <- length(universe) - m
      K <- length(hit_list)
      # Upper-tail hypergeometric (enrichment)
      p <- phyper(q = k - 1, m = m, n = n, k = K, lower.tail = FALSE)
      data.frame(set = sn, k = k, K = K, m = m, n = n, p = p, stringsAsFactors = FALSE)
    })
    out <- dplyr::bind_rows(results)
    if (is.null(out) || nrow(out) == 0) return(NULL)
    out$p_adj <- p.adjust(out$p, method = "BH")
    dplyr::arrange(out, p_adj, p)
  }
  
  
  # ---- Main results (Overall or Per-group scopes; ORA or FGSEA) ----
  
  
  
  run_fgsea_safe <- function(ranks, sets, minSize = 3, maxSize = 5000, eps = 1e-6) {
    # Preconditions
    if (!requireNamespace("fgsea", quietly = TRUE)) {
      stop("Package 'fgsea' is not installed/loaded.")
    }
    if (is.null(names(ranks)) || !is.numeric(ranks)) {
      stop("FGSEA: 'ranks' must be a named numeric vector.")
    }
    
    # Keep only finite ranks, preserve names
    ranks <- ranks[is.finite(ranks)]
    if (length(ranks) < minSize) return(NULL)
    
    # Intersect each set with available rank names and drop small sets
    clean_sets <- lapply(sets, function(s) {
      intersect(unique(as.character(stats::na.omit(s))), names(ranks))
    })
    keep <- vapply(clean_sets, length, integer(1)) >= minSize
    clean_sets <- clean_sets[keep]
    if (length(clean_sets) == 0) return(NULL)
    
    # Use fgseaMultilevel (recommended) — no nperm, no parallel
    res <- fgsea::fgseaMultilevel(
      pathways = clean_sets,
      stats    = ranks,
      minSize  = minSize,
      maxSize  = maxSize,
      eps      = eps
    )
    
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    # Normalize output
    out <- as.data.frame(res)
    # Sort by adjusted p-value, then raw p-value
    out <- dplyr::arrange(out, padj, pval)
    out
  }
  
  lsea_results <- reactive({
    # 1) Get sets and universe; fail fast if empty
    sets <- lsea_sets()
    validate(need(length(sets) > 0, "No lipid sets available under current options."))
    
    univ <- lsea_universe()
    validate(need(length(univ) > 0, "No tested features found for the selected universe."))
    
    # 2) Branch by scope: Overall vs Per-group (one-vs-rest)
    if (identical(input$lsea_scope, "overall")) {
      
      # ---- Overall scope ----
      if (identical(input$lsea_method, "ora")) {
        # ORA on overall hits
        res <- try(run_ora(lsea_hits_overall(), univ, sets), silent = TRUE)
        validate(need(!inherits(res, "try-error") && !is.null(res) && nrow(res) > 0,
                      "No ORA enrichment detected under current settings."))
        
        # Optional: only sets with ≥1 overlapping hit
        if (isTRUE(input$ora_only_overlap)) {
          res <- subset(res, k > 0)
          validate(need(nrow(res) > 0, "No sets with ≥1 hit under current settings."))
        }
        
        # Standardize columns for downstream plot/table
        res$size  <- res$m
        res$label <- res$set
        res$type  <- "ORA"
        res$for_group <- NA_character_
        return(res)
        
      } else {
        
        # Overall scope FGSEA branch
        r    <- lsea_rank_overall()
        sets <- lsea_sets()
        univ <- lsea_universe()
        
        
        res <- try(run_fgsea_safe(r, sets, minSize = input$lsea_min_set %||% 3), silent = TRUE)
        validate(need(!inherits(res, "try-error") && !is.null(res) && nrow(res) > 0,
                      "No FGSEA enrichment detected (insufficient rank coverage or sets)."))
        res$label     <- res$pathway
        res$type      <- "FGSEA"
        res$for_group <- NA_character_
        return(res)
      }
      
    } else {
      
      # ---- Per-group (one-vs-rest) scope ----
      groups <- input$lsea_group_select
      validate(need(length(groups) > 0, "Select at least one group."))
      
      if (identical(input$lsea_method, "ora")) {
        # Build hits per group with robust empties handling
        hits <- lsea_hits_per_group()
        validate(need(!is.null(hits) && nrow(hits) > 0,
                      "No hits for selected group(s) under current thresholds."))
        
        out <- lapply(groups, function(g) {
          # Subset this group's hits; skip if none
          hg <- hits[hits$for_group == g, , drop = FALSE]
          if (is.null(hg) || nrow(hg) == 0) return(NULL)
          
          # Run ORA for this group; skip on error/empty
          resg <- try(run_ora(hg$metabolite, univ, sets), silent = TRUE)
          if (inherits(resg, "try-error") || is.null(resg) || nrow(resg) == 0) return(NULL)
          
          # Optional: show only sets with ≥1 overlapping hit; skip if filtered to empty
          if (isTRUE(input$ora_only_overlap)) {
            resg <- subset(resg, k > 0)
            if (nrow(resg) == 0) return(NULL)
          }
          
          # Standardize columns
          resg$size  <- resg$m
          resg$label <- resg$set
          resg$type  <- "ORA"
          resg$for_group <- g
          resg
        })
        
        res <- dplyr::bind_rows(out)
        validate(need(!is.null(res) && nrow(res) > 0,
                      "No ORA enrichment detected for selected group(s)."))
        return(res)
        
      } else {
        
        
        ranks <- lsea_rank_per_group()
        sets  <- lsea_sets()
        validate(need(length(ranks) > 0, "No ranks available for selected group(s)."))
        
        out <- lapply(names(ranks), function(g) {
          rg <- ranks[[g]]
          if (is.null(rg) || length(rg) == 0) return(NULL)
          
          resg <- try(run_fgsea_safe(rg, sets, minSize = input$lsea_min_set %||% 3), silent = TRUE)
          if (inherits(resg, "try-error") || is.null(resg) || nrow(resg) == 0) return(NULL)
          
          resg$label     <- resg$pathway
          resg$type      <- "FGSEA"
          resg$for_group <- g
          resg
        })
        
        res <- dplyr::bind_rows(out)
        validate(need(!is.null(res) && nrow(res) > 0, "No FGSEA enrichment detected for selected group(s)."))
        return(res)
      }
    }
  })
  
  
  # ---- Outputs: context, plot, table, download, references ----
  
  output$lseaContextInfo <- renderText({
    show_loading_modal("Preparing enrichment…")
    on.exit(hide_loading_modal(), add = TRUE)
    
    univ_n <- length(lsea_universe())
    sets_n <- length(lsea_sets())
    scope  <- if (identical(input$lsea_scope, "group")) "Per-group (one-vs-rest)" else "Overall"
    meth   <- if (identical(input$lsea_method, "ora")) "ORA" else "FGSEA"
    paste0("Scope: ", scope, " • Universe: ", univ_n, " lipid species • Sets: ", sets_n, " • Method: ", meth)
  })
  
  
  lsea_plot_obj <- reactive({
    res  <- lsea_results()
    topn <- input$lsea_topn; if (is.null(topn)) topn <- 30
    
    if (identical(res$type[1], "FGSEA")) {
      top <- res %>% dplyr::group_by(for_group) %>% dplyr::slice_head(n = topn) %>% dplyr::ungroup()
      p <- ggplot(top, aes(x = NES, y = reorder(label, NES), size = size, color = -log10(padj))) +
        geom_point() +
        labs(x = "Normalized Enrichment Score (NES)", y = "Set",
             color = "-log10(adj p)", size = "Set size") +
        theme_minimal(base_size = 12)
    } else {
      top <- res %>% dplyr::group_by(for_group) %>% dplyr::slice_head(n = topn) %>% dplyr::ungroup()
      xval <- if (identical(input$ora_x_axis, "p")) -log10(pmax(top$p,     1e-300))
      else                                 -log10(pmax(top$p_adj, 1e-300))
      top$xval <- xval
      p <- ggplot(top, aes(x = xval, y = reorder(label, xval), size = size, color = -log10(p))) +
        geom_point() +
        labs(x = if (identical(input$ora_x_axis, "p")) "-log10(p)" else "-log10(adj p)",
             y = "Set", size = "Set size", color = "-log10(p)") +
        theme_minimal(base_size = 12)
    }
    
    if (identical(input$lsea_scope, "group") && isTRUE(input$lsea_facet_groups) &&
        any(!is.na(res$for_group)) && length(unique(na.omit(res$for_group))) > 1) {
      p <- p + facet_wrap(~ for_group, scales = "free_y", ncol = 1)
    }
    p
  })
  
  output$lseaDotPlot <- renderPlot({
    show_loading_modal("Running enrichment…", "Building enrichment plot")
    on.exit(hide_loading_modal(), add = TRUE)
    lsea_plot_obj()
  })
  
  output$lseaTable <- renderDT({
    # Show loading while FGSEA/ORA table prepares
    show_loading_modal("Running enrichment…", "Preparing enrichment table")
    on.exit(hide_loading_modal(), add = TRUE)
    
    res <- lsea_results()  # heavy reactive
    if (identical(res$type[1], "FGSEA")) {
      keep <- res[, c("for_group","label","NES","pval","padj","size"), drop = FALSE]
      colnames(keep) <- c("Group","Set","NES","p","p_adj","Size")
    } else {
      res$neglog10_p     <- -log10(pmax(res$p, 1e-300))
      res$neglog10_p_adj <- -log10(pmax(res$p_adj, 1e-300))
      keep <- res[, c("for_group","label","k","K","m","n","p","p_adj","neglog10_p","neglog10_p_adj","size"), drop = FALSE]
      colnames(keep) <- c("Group","Set","k (hits in set)","K (hits total)","m (set size)","n (bg rest)",
                          "p","p_adj","-log10(p)","-log10(p_adj)","Size")
    }
    
    DT::datatable(keep, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  
  output$download_lsea_table <- downloadHandler(
    filename = function() paste0("lsea_", input$lsea_scope, "_", input$lsea_method, "_", Sys.Date(), ".csv"),
    content = function(file) { readr::write_csv(lsea_results(), file) }
  )
  
  
  
  .lsea_export_dims <- function() {
    px_w  <- input$lsea_export_width    %||% 1200
    px_h  <- input$lsea_export_height   %||% 700
    dpi   <- input$lsea_export_dpi      %||% 300
    scale <- input$lsea_export_scale    %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  output$download_lsea_png <- downloadHandler(
    filename = function() paste0("enrichment_", isolate(input$lsea_method), "_", Sys.Date(), ".png"),
    content = function(file) {
      dims <- isolate(.lsea_export_dims())
      p    <- isolate(lsea_plot_obj())
      validate(need(!is.null(p), "No enrichment plot to export. Run enrichment first."))
      ggplot2::ggsave(file, plot = p,
                      width = dims$w, height = dims$h,
                      dpi = dims$dpi, device = "png")
    }
  )
  
  output$download_lsea_svg <- downloadHandler(
    filename = function() paste0("enrichment_", isolate(input$lsea_method), "_", Sys.Date(), ".svg"),
    content = function(file) {
      dims <- isolate(.lsea_export_dims())
      p    <- isolate(lsea_plot_obj())
      validate(need(!is.null(p), "No enrichment plot to export. Run enrichment first."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  
  observeEvent(input$show_refs_lsea, {
    showModal(modalDialog(
      title = "References (Lipid Enrichment)",
      easyClose = TRUE, footer = modalButton("Close"),
      HTML(paste(
        "<p><strong>lipidr</strong> — LSEA on class, total chain length and unsaturation (Bioconductor). ",
        "Docs: <a href='https://www.bioconductor.org/packages/lipidr.html' target='_blank'>link</a>; ",
        "Vignette: <a href='https://bioconductor.org/packages/release/bioc/vignettes/lipidr/inst/doc/workflow.html' target='_blank'>workflow</a>.</p>",
        "<p><strong>LION/web</strong> — Lipid ontology & biophysical property enrichment. ",
        "GigaScience 2019: <a href='https://academic.oup.com/gigascience/article/8/6/giz061/5505544' target='_blank'>paper</a>; ",
        "Preprint: <a href='https://www.biorxiv.org/content/10.1101/398040v2.full.pdf' target='_blank'>bioRxiv</a>.</p>",
        "<p><strong>FGSEA</strong> — Fast GSEA implementation for ranked lists (R). ",
        "CRAN/Bioconductor.</p>"
      ))
    ))
  })
  
  # ---- Pathway class sets for ORA (curated; editable in one place) ----
  pathway_sets_static <- reactive({
    # Each set is a vector of class names as they appear in 'plot_class'
    # (including "-O" where relevant and present in your data).
    list(
      # Kennedy pathway branches (product-side class membership)
      "Kennedy-PC (DAG→PC)"                 = c("PC", "PC-O"),
      "Kennedy-PE (DAG→PE)"                 = c("PE", "PE-O"),
      
      # NEW: Specific glycerophospholipid synthesis steps (product classes)
      "PG synthesis (DAG→PG)"               = c("PG"),
      "CL synthesis (PG→CL)"                = c("CL"),
      
      # Storage & CDP-DAG branch
      "DGAT (DAG→TG; storage)"              = c("TG"),
      "CDP-DAG branch (PI/PG/CL)"           = c("PI", "PG", "CL"),
      "PA node (upstream glycerolipid)"     = c("PA"),
      
      # Sphingolipid branch points
      "Sphingomyelin synthase (Cer→SM)"     = c("SM"),
      "Glycosphingolipid branch (Cer→HexCer)" = c("HexCer"),
      "Ceramide node (upstream sphingolipid)" = c("Cer"),
      
      # Ether pathways (product-side class membership)
      "Ether pathway to PC-O"               = c("PC-O"),
      "Ether pathway to PE-O"               = c("PE-O"),
      
      # Lands cycle lyso classes
      "Lands cycle lyso-PC"                 = c("LPC"),
      "Lands cycle lyso-PE"                 = c("LPE"),
      "Lands cycle lyso-PI"                 = c("LPI"),
      "Lands cycle lyso-PS"                 = c("LPS"),
      "Lands cycle lyso-PG"                 = c("LPG")
    )
  })
  
  
  
  
  # ---- Pathway scores (ratios/fractions) ----
  # type = "ratio": (sum(num) / sum(den)); type = "fraction": num / (num + others in den)
  pathway_scores_defs <- reactive({
    list(
      
      # --- NEW: Ether lipid synthesis scores (FA-resolved where available) ---
      "TG-O synthesis (TG-O/DG-O)" = list(
        type = "ratio",
        num  = c("TG-O"),
        den  = c("DG-O")
      ),
      "plasmanyl-PE synthesis (PE-O XX:0/DG-O)" = list(
        type = "ratio",
        num  = c("PE-O plasmanyl"),
        den  = c("DG-O")
      ),
      "plasmanyl-PC synthesis (PC-O XX:0/PE-O XX:0)" = list(
        type = "ratio",
        num  = c("PC-O plasmanyl"),
        den  = c("PE-O plasmanyl")
      ),
      "plasmanyl-PC synthesis, direct (PC-O XX:0/DG-O)" = list(
        type = "ratio",
        num  = c("PC-O plasmanyl"),
        den  = c("DG-O")
      ),
      "plasmenyl-PE synthesis (PE-O XX:>=1/PE-O XX:0)" = list(
        type = "ratio",
        num  = c("PE-O plasmenyl"),
        den  = c("PE-O plasmanyl")
      ),
      "plasmenyl-PC synthesis (PC-O XX:>=1/PE-O XX:>=1)" = list(
        type = "ratio",
        num  = c("PC-O plasmenyl"),
        den  = c("PE-O plasmenyl")
      ),
      
      # --- NEW: Additional glycerophospholipid pathways ---
      "PG synthesis (PG/DAG)" = list(
        type = "ratio",
        num  = c("PG"),
        den  = c("DAG")
      ),
      "CL synthesis (CL/PG)" = list(
        type = "ratio",
        num  = c("CL"),
        den  = c("PG")
      ),
      
      # --- EXISTING SCORES (unchanged) ---
      "TAG storage (TG/DAG)" = list(
        type = "ratio",
        num  = c("TG"),
        den  = c("DAG")
      ),
      "DGAT activity (TG/(TG + DAG))" = list(
        type = "fraction",
        num  = c("TG"),
        den  = c("TG", "DAG")
      ),
      "PC synthesis (PC/DAG)" = list(
        type = "ratio",
        num  = c("PC", "PC-O"),
        den  = c("DAG")
      ),
      "PE synthesis (PE/DAG)" = list(
        type = "ratio",
        num  = c("PE", "PE-O"),
        den  = c("DAG")
      ),
      "PEMT (PC/PE)" = list(
        type = "ratio",
        num  = c("PC", "PC-O"),
        den  = c("PE", "PE-O")
      ),
      "SM synthase (SM/Cer)" = list(
        type = "ratio",
        num  = c("SM"),
        den  = c("Cer")
      ),
      "Glycosphingolipids (HexCer/Cer)" = list(
        type = "ratio",
        num  = c("HexCer"),
        den  = c("Cer")
      ),
      "PI via CDP-DAG (PI/PA)" = list(
        type = "ratio",
        num  = c("PI"),
        den  = c("PA")
      ),
      "PG+CL via CDP-DAG ((PG+CL)/PA)" = list(
        type = "ratio",
        num  = c("PG", "CL"),
        den  = c("PA")
      ),
      "Ether index PC (PC-O/(PC-O + PC))" = list(
        type = "fraction",
        num  = c("PC-O"),
        den  = c("PC-O", "PC")
      ),
      "Ether index PE (PE-O/(PE-O + PE))" = list(
        type = "fraction",
        num  = c("PE-O"),
        den  = c("PE-O", "PE")
      ),
      "Phospholipid to storage ((PC+PE+PI+PS+PG)/TG)" = list(
        type = "ratio",
        num  = c("PC", "PC-O", "PE", "PE-O", "PI", "PS", "PG"),
        den  = c("TG")
      )
    )
  })
  
  # Sum per sample × class, with unit option "absolute" or "% of total lipids"
  path_class_totals <- reactive({
    req(bg_norm_long_avg(), input$path_value_type, input$path_units)
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample))) %>%
      filter_blanks(exclude = input$exclude_blank_path, sample_col = "sample_norm", exact = FALSE) 
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    
    # Respect the group selection from the Statistics tab
    fn <- tryCatch(get_active_grouping(), error = function(e) NULL)
    if (!is.null(fn)) {
      grp_res <- tryCatch(fn(df$sample), error = function(e) NULL)
      if (!is.null(grp_res)) {
        df$group <- grp_res$group
        sel_grps <- input$stats_selected_groups
        if (!is.null(sel_grps) && length(sel_grps) > 0) {
          df <- df %>% dplyr::filter(group %in% sel_grps)
        }
        df$group <- NULL  # remove temp column before summarise
      }
    }
    
    measure_col <- pick_measure_col(input$path_value_type)
    
    class_tot <- df %>%
      dplyr::group_by(sample, plot_class) %>%
      dplyr::summarise(value = sum(.data[[measure_col]], na.rm = TRUE), .groups = "drop")
    
    if (identical(input$path_units, "percent_total")) {
      totals <- class_tot %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(tot = sum(value, na.rm = TRUE), .groups = "drop")
      class_tot <- class_tot %>%
        dplyr::left_join(totals, by = "sample") %>%
        dplyr::mutate(value = dplyr::if_else(tot > 0, 100 * value / tot, NA_real_)) %>%
        dplyr::select(-tot)
    }
    class_tot
  })
  
  # ---- Class totals + FA-filtered ether sub-class totals for scoring ----
  path_totals_extended <- reactive({
    req(bg_norm_long_avg(), input$path_value_type, input$path_units)
    
    # Base long DF with chosen measure & blank handling (same as your pipeline)
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample))) %>%
      filter_blanks(exclude = input$exclude_blank_path, sample_col = "sample_norm", exact = FALSE) %>% 
      filter_iqc(include_iqc = FALSE, sample_col = "sample")
    
    measure_col <- pick_measure_col(input$path_value_type)
    
    # Respect the group selection from the Statistics tab
    fn_grp <- tryCatch(get_active_grouping(), error = function(e) NULL)
    if (!is.null(fn_grp)) {
      grp_res <- tryCatch(fn_grp(df$sample), error = function(e) NULL)
      if (!is.null(grp_res)) {
        df$group <- grp_res$group
        sel_grps <- input$stats_selected_groups
        if (!is.null(sel_grps) && length(sel_grps) > 0) {
          df <- df %>% dplyr::filter(group %in% sel_grps)
        }
        df$group <- NULL
      }
    }
    
    df_vals <- df %>%
      dplyr::transmute(sample, plot_class, `Metabolite name`, value = .data[[measure_col]])
    
    # Total lipid per sample (used when Units = % of total lipids)
    tot_by_sample <- df_vals %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(tot = sum(value, na.rm = TRUE), .groups = "drop")
    
    # Base class totals (unchanged)
    class_tot <- df_vals %>%
      dplyr::group_by(sample, plot_class) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
    
    # Ether sub-class totals (requires resolved FA)
    ann <- try(ether_subclass_lookup(), silent = TRUE)
    ether_tot <- NULL
    if (!inherits(ann, "try-error") && !is.null(ann) && nrow(ann) > 0) {
      df_join <- df_vals %>%
        dplyr::left_join(ann, by = c("Metabolite name" = "Metabolite name", "plot_class" = "plot_class"))
      
      # PE-O/PC-O plasmanyl & plasmenyl
      ether_tot <- df_join %>%
        dplyr::filter(plot_class %in% c("PE-O", "PC-O"), !is.na(ether_subclass)) %>%
        dplyr::mutate(label = paste0(plot_class, " ", ether_subclass)) %>%
        dplyr::group_by(sample, label) %>%
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
      
      # Convenience total: PE-O >=0 (all PE-O) to match your ">=0" denominator
      peo_total <- df_join %>%
        dplyr::filter(plot_class == "PE-O") %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(label = "PE-O total")
      
      ether_tot <- dplyr::bind_rows(ether_tot, peo_total %>% dplyr::select(sample, label, value))
    }
    
    # Combine: base classes use label = plot_class
    base_long <- class_tot %>% dplyr::rename(label = plot_class)
    ext_long  <- if (!is.null(ether_tot)) dplyr::bind_rows(base_long, ether_tot) else base_long
    
    # Apply unit conversion if requested (% of total lipids)
    if (identical(input$path_units, "percent_total")) {
      ext_long <- ext_long %>%
        dplyr::left_join(tot_by_sample, by = "sample") %>%
        dplyr::mutate(value = dplyr::if_else(tot > 0, 100 * value / tot, NA_real_)) %>%
        dplyr::select(-tot)
    }
    
    ext_long
  })
  
  
  observeEvent(path_class_totals(), {
    df <- path_class_totals()
    # Build sample->group mapping safely (token UI may not exist yet)
    samples <- unique(df$sample)
    fn  <- tryCatch(get_active_grouping(), error = function(e) NULL)
    grps_raw <- if (!is.null(fn))
      tryCatch(fn(samples)$group, error = function(e) rep("Unassigned", length(samples)))
    else
      rep("Unassigned", length(samples))
    sample_map <- data.frame(sample = samples, group = grps_raw, stringsAsFactors = FALSE)
    # Score choices
    sc_names <- names(pathway_scores_defs())
    updateSelectizeInput(session, "path_scores_select",
                         choices = sc_names,
                         selected = sc_names)
  })
  
  path_scores_long <- reactive({
    req(path_totals_extended())
    class_tot <- path_totals_extended()
    
    # Wide: sample × label (label = class OR ether-subclass)
    wide <- tidyr::pivot_wider(
      class_tot,
      names_from  = label,
      values_from = value,
      values_fill = NA_real_  # keep NA so we can detect missing FA-resolved info
    )
    
    stopifnot("sample" %in% names(wide))
    
    defs <- pathway_scores_defs()
    eps  <- 1e-12
    
    # helper to add missing columns; here we add NA (not zeros) to reflect "not observed"
    add_missing_cols <- function(df, cols) {
      miss <- setdiff(cols, names(df))
      if (length(miss)) for (c in miss) df[[c]] <- NA_real_
      df
    }
    
    out <- list()
    for (nm in names(defs)) {
      d <- defs[[nm]]
      need_cols <- unique(c(d$num, d$den))
      tmp <- add_missing_cols(wide, need_cols)
      
      num_val <- rowSums(tmp[, intersect(names(tmp), d$num), drop = FALSE], na.rm = TRUE)
      den_val <- rowSums(tmp[, intersect(names(tmp), d$den), drop = FALSE], na.rm = TRUE)
      
      val <- if (identical(d$type, "ratio")) {
        # NA when both numerator and denominator are zero/absent
        ifelse(num_val == 0 & den_val == 0, NA_real_, (num_val + eps) / (den_val + eps))
      } else {
        # fraction
        ifelse(den_val == 0 & num_val == 0, NA_real_, num_val / pmax(den_val, eps))
      }
      
      out[[nm]] <- as.numeric(val)
    }
    
    scores <- cbind(sample = wide$sample, as.data.frame(out, check.names = FALSE))
    
    # Attach group labels (unchanged)
    fn <- get_active_grouping()
    scores$group <- fn(scores$sample)$group
    
    
    
    scores %>% tidyr::pivot_longer(cols = -c(sample, group), names_to = "score", values_to = "value")
  })
  
  path_heatmap_data <- reactive({
    req(path_scores_long())
    df <- path_scores_long()
    sel <- input$path_scores_select
    if (!is.null(sel) && length(sel) > 0)
      df <- df %>% dplyr::filter(score %in% sel)
    
    if (isTRUE(input$path_use_group_means)) {
      mat_df <- df %>%
        dplyr::group_by(score, group) %>%
        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = group, values_from = value) %>%
        as.data.frame()
    } else {
      mat_df <- df %>%
        tidyr::pivot_wider(names_from = sample, values_from = value) %>%
        as.data.frame()
    }
    
    rn <- mat_df$score
    mat <- as.matrix(mat_df[, -1, drop = FALSE])
    rownames(mat) <- rn
    
    # z-score by row for visual comparability
    mat_sc <- t(scale(t(mat)))
    mat_sc[!is.finite(mat_sc)] <- 0
    list(mat = mat_sc)
  })
  
  output$pathHeatmap <- renderPlot({
    hd <- path_heatmap_data()
    cols <- .heatmap_palette(input$path_hm_palette %||% "viridis")
    pheatmap::pheatmap(
      hd$mat, color = cols,
      cluster_rows = TRUE, cluster_cols = TRUE,
      border_color = NA,
      show_rownames = TRUE, show_colnames = TRUE,
      fontsize_row = 9, fontsize_col = 10,
      main = "Pathway scores (z-scored by row)"
    )
  }, height = function() {
    hd <- path_heatmap_data()
    pmin(1600, pmax(400, 18 * (nrow(hd$mat) %||% 10)))
  })
  
  output$pathScoresTable <- DT::renderDT({
    sc <- path_scores_long()
    sel <- input$path_scores_select
    if (!is.null(sel) && length(sel) > 0) sc <- sc %>% dplyr::filter(score %in% sel)
    DT::datatable(sc, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$download_path_scores <- downloadHandler(
    filename = function() paste0("pathway_scores_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(path_scores_long(), file)
  )
  
  # ---- Synthesis Pathways plot export ----
  .path_export_dims <- function() {
    px_w  <- input$path_export_width    %||% 1400
    px_h  <- input$path_export_height   %||% 900
    dpi   <- input$path_export_dpi      %||% 300
    scale <- input$path_export_scale    %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  .build_path_plot <- function(which, fsz = 11) {
    if (which == "scorestats") which <- "scorestats"  # explicit alias
    if (which == "heatmap") {
      hd <- path_heatmap_data()
      if (is.null(hd)) return(NULL)
      cols <- .heatmap_palette(input$path_hm_palette %||% "viridis")
      pheatmap::pheatmap(
        hd$mat, color = cols,
        cluster_rows = TRUE, cluster_cols = TRUE,
        border_color = NA,
        show_rownames = TRUE, show_colnames = TRUE,
        fontsize_row = fsz, fontsize_col = fsz + 1,
        main = "Pathway scores (z-scored by row)",
        silent = TRUE
      )
    } else {
      # Score statistics lollipop — rebuild from path_score_stats
      res <- path_score_stats()
      if (is.null(res) || nrow(res) == 0) return(NULL)
      alpha    <- input$path_alpha %||% 0.05
      use_padj <- !identical(input$path_padj_method %||% "BH", "none")
      p_col    <- if (use_padj) "p_adj" else "p"
      x_label  <- if (use_padj) expression(-log[10](p[adj])) else expression(-log[10](p))
      res <- res %>%
        dplyr::mutate(
          neg_log10_p = -log10(pmax(.data[[p_col]], 1e-300)),
          sig         = .data[[p_col]] < alpha,
          score_fct   = forcats::fct_reorder(score, neg_log10_p)
        )
      ggplot2::ggplot(res, ggplot2::aes(
        x = neg_log10_p, y = score_fct, colour = direction, fill = direction
      )) +
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = neg_log10_p, yend = score_fct),
          linewidth = 0.7, alpha = 0.5
        ) +
        ggplot2::geom_point(ggplot2::aes(shape = sig), size = 3.5) +
        ggplot2::scale_shape_manual(
          values = c("TRUE" = 19, "FALSE" = 1),
          labels = c("TRUE" = paste0("p < ", alpha), "FALSE" = "NS"),
          name   = "Significance"
        ) +
        ggplot2::geom_vline(xintercept = -log10(alpha), linetype = "dashed", colour = "grey50") +
        ggplot2::labs(x = x_label, y = NULL, colour = "Higher in", fill = "Higher in",
                      title = "Pathway score group comparisons") +
        ggplot2::theme_minimal(base_size = fsz) +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "right")
    }
  }
  
  output$download_path_png <- downloadHandler(
    filename = function() {
      sel <- isolate(input$path_export_plot_sel %||% "heatmap")
      paste0("pathway_", sel, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      dims <- isolate(.path_export_dims())
      fsz  <- isolate(input$path_export_fontsize %||% 11)
      sel  <- isolate(input$path_export_plot_sel %||% "heatmap")
      obj  <- isolate(.build_path_plot(sel, fsz))
      validate(need(!is.null(obj), "No plot to export."))
      if (inherits(obj, "pheatmap") || inherits(obj, "list")) {
        grDevices::png(file,
                       width  = round(dims$w * dims$dpi),
                       height = round(dims$h * dims$dpi),
                       res    = dims$dpi)
        on.exit(grDevices::dev.off(), add = TRUE)
        grid::grid.newpage()
        grid::grid.draw(obj$gtable)
      } else {
        ggplot2::ggsave(file, plot = obj,
                        width = dims$w, height = dims$h,
                        dpi = dims$dpi, device = "png")
      }
    }
  )
  
  output$download_path_svg <- downloadHandler(
    filename = function() {
      sel <- isolate(input$path_export_plot_sel %||% "heatmap")
      paste0("pathway_", sel, "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      dims <- isolate(.path_export_dims())
      fsz  <- isolate(input$path_export_fontsize %||% 11)
      sel  <- isolate(input$path_export_plot_sel %||% "heatmap")
      obj  <- isolate(.build_path_plot(sel, fsz))
      validate(need(!is.null(obj), "No plot to export."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (inherits(obj, "pheatmap") || inherits(obj, "list")) {
        grid::grid.newpage()
        grid::grid.draw(obj$gtable)
      } else {
        print(obj)
      }
    }
  )
  
  # One-vs-rest t-tests on class totals to define "hits" per group (classes as features)
  # ========== SYNTHESIS PATHWAYS — SCORE STATISTICS (Option D) ==========
  
  # Reactive: run t-test / one-way ANOVA on each score's sample values
  path_score_stats <- reactive({
    req(path_scores_long())
    sc   <- path_scores_long()
    sel  <- input$path_scores_select
    if (!is.null(sel) && length(sel) > 0) sc <- sc %>% dplyr::filter(score %in% sel)
    
    # Drop samples with no group assignment
    sc <- sc %>% dplyr::filter(!is.na(group) & nzchar(group) & group != "Unassigned")
    
    grps <- unique(sc$group)
    if (length(grps) < 2) {
      showNotification(
        "Score statistics: need ≥2 groups. Check Group Preview tab.",
        type = "warning", duration = 6
      )
      return(NULL)
    }
    
    padj_method  <- input$path_padj_method  %||% "BH"
    equal_var    <- isTRUE(input$path_stats_equal_var)
    n_grps       <- length(grps)
    
    rows <- lapply(unique(sc$score), function(s) {
      df_s <- sc %>%
        dplyr::filter(score == s) %>%
        dplyr::filter(!is.na(value))
      
      # Need ≥2 groups with ≥1 observation each
      grp_counts <- table(df_s$group)
      valid_grps <- names(grp_counts[grp_counts >= 1])
      if (length(valid_grps) < 2) return(NULL)
      
      df_s <- df_s %>% dplyr::filter(group %in% valid_grps)
      groups_present <- sort(unique(df_s$group))
      
      # Group means for direction / fold-change
      means <- df_s %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")
      
      if (length(groups_present) == 2) {
        # t-test
        g1 <- df_s$value[df_s$group == groups_present[1]]
        g2 <- df_s$value[df_s$group == groups_present[2]]
        tt <- tryCatch(
          stats::t.test(g1, g2, var.equal = equal_var),
          error = function(e) NULL
        )
        if (is.null(tt)) return(NULL)
        m1 <- mean(g1, na.rm = TRUE)
        m2 <- mean(g2, na.rm = TRUE)
        log2fc <- log2((m1 + 1e-12) / (m2 + 1e-12))
        data.frame(
          score      = s,
          test       = "t-test",
          groups     = paste(groups_present, collapse = " vs "),
          statistic  = round(as.numeric(tt$statistic), 3),
          log2FC     = round(log2fc, 3),
          direction  = if (log2fc > 0) groups_present[1] else groups_present[2],
          p          = as.numeric(tt$p.value),
          stringsAsFactors = FALSE
        )
      } else {
        # one-way ANOVA
        fit <- tryCatch(
          stats::aov(value ~ group, data = df_s),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        sm  <- summary(fit)[[1]]
        Fval <- sm[["F value"]][1]
        pval <- sm[["Pr(>F)"]][1]
        if (is.na(pval)) return(NULL)
        # Direction = group with highest mean
        top_grp <- means$group[which.max(means$mean)]
        data.frame(
          score      = s,
          test       = "one-way ANOVA",
          groups     = paste(groups_present, collapse = " / "),
          statistic  = round(Fval, 3),
          log2FC     = NA_real_,
          direction  = top_grp,
          p          = as.numeric(pval),
          stringsAsFactors = FALSE
        )
      }
    })
    
    out <- dplyr::bind_rows(rows)
    if (is.null(out) || nrow(out) == 0) return(NULL)
    
    # BH correction across all scores
    out$p_adj <- p.adjust(out$p, method = if (identical(padj_method, "none")) "none" else padj_method)
    out <- out %>% dplyr::arrange(p_adj, p)
    out
  })
  
  # Plot: horizontal lollipop of -log10(p_adj), coloured by direction
  output$pathStatsPlot <- renderPlot({
    res <- path_score_stats()
    validate(need(!is.null(res) && nrow(res) > 0,
                  "No score statistics yet. Check grouping in Group Preview tab."))
    
    alpha    <- input$path_alpha %||% 0.05
    use_padj <- !identical(input$path_padj_method %||% "BH", "none")
    p_col    <- if (use_padj) "p_adj" else "p"
    x_label  <- if (use_padj) expression(-log[10](p[adj])) else expression(-log[10](p))
    
    res <- res %>%
      dplyr::mutate(
        neg_log10_p = -log10(pmax(.data[[p_col]], 1e-300)),
        sig         = .data[[p_col]] < alpha,
        score_fct   = forcats::fct_reorder(score, neg_log10_p)
      )
    
    direction_colours <- setNames(
      scales::hue_pal()(length(unique(res$direction))),
      unique(res$direction)
    )
    
    ggplot2::ggplot(res, ggplot2::aes(
      x     = neg_log10_p,
      y     = score_fct,
      colour = direction,
      fill   = direction
    )) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = neg_log10_p, yend = score_fct),
        linewidth = 0.7, alpha = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(shape = sig),
        size = 3.5
      ) +
      ggplot2::scale_shape_manual(
        values = c("TRUE" = 19, "FALSE" = 1),
        labels = c("TRUE" = paste0("p < ", alpha), "FALSE" = "NS"),
        name   = "Significance"
      ) +
      ggplot2::geom_vline(
        xintercept = -log10(alpha),
        linetype   = "dashed",
        colour     = "grey50"
      ) +
      ggplot2::labs(
        x      = x_label,
        y      = NULL,
        colour = "Higher in",
        fill   = "Higher in",
        title  = "Pathway score group comparisons"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position    = "right"
      )
  })
  
  # Table: clean formatted results
  output$pathStatsTable <- DT::renderDT({
    res <- path_score_stats()
    validate(need(!is.null(res) && nrow(res) > 0, "No results yet."))
    alpha    <- input$path_alpha %||% 0.05
    use_padj <- !identical(input$path_padj_method %||% "BH", "none")
    display  <- res %>%
      dplyr::mutate(
        p      = signif(p,     3),
        p_adj  = signif(p_adj, 3),
        log2FC = ifelse(is.na(log2FC), "—", as.character(round(log2FC, 2)))
      ) %>%
      dplyr::select(score, test, groups, statistic, log2FC, direction, p, p_adj)
    colnames(display) <- c("Score", "Test", "Groups", "Statistic", "log2FC", "Higher in", "p", "p_adj")
    DT::datatable(
      display,
      rownames = FALSE,
      options  = list(scrollX = TRUE, pageLength = 20),
      class    = "stripe hover compact"
    ) %>%
      DT::formatStyle(
        "p_adj",
        backgroundColor = DT::styleInterval(
          alpha,
          c("#d4edda", "white")
        )
      )
  })
  
  # ========== CLASS BAR PLOTS TAB ==========
  
  # Class selector — sources classes from bg_norm_long_resolved so it is
  # populated as soon as data is loaded, not only after stats are run
  output$cbp_class_ui <- renderUI({
    df <- tryCatch(bg_norm_long_avg(), error = function(e) NULL)
    classes <- if (!is.null(df) && "plot_class" %in% names(df))
      sort(unique(df$plot_class[!stringr::str_detect(df$`Metabolite name`, "\\[IS\\]")]))
    else character(0)
    selectInput("cbp_class", "Lipid class",
                choices  = classes,
                selected = if (length(classes)) classes[1] else NULL)
  })
  
  # Core reactive: filtered species for the selected class
  cbp_data <- reactive({
    req(bg_norm_long_avg(), stats_results_all())
    # Use stats_results_all so all tested species are available regardless of
    # the significance filter; sig_only toggle is applied within this reactive
    res <- stats_results_all()
    validate(need(nrow(res) > 0, "No features for this class. Run statistics first."))
    
    # Build long df from bg_norm_long_resolved — same filtering as stats_input_long
    df_long <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample)))
    df_long <- filter_blanks(df_long, input$exclude_blank_stats,
                             sample_col = "sample_norm", exact = FALSE)
    df_long <- filter_iqc(df_long, include_iqc = FALSE, sample_col = "sample")
    
    # Apply grouping and restrict to stats-selected groups
    fn  <- get_active_grouping()
    grp_res <- fn(df_long$sample)
    df_long$group <- grp_res$group
    sel_grps <- input$stats_selected_groups
    if (!is.null(sel_grps) && length(sel_grps) > 0) {
      df_long <- df_long %>% dplyr::filter(group %in% sel_grps)
    }
    
    # Value type — column names are "norm" and "value_bs" in bg_norm_long_avg
    measure_col <- if (identical(input$cbp_value_type %||% "norm", "value_bs")) "value_bs" else "norm"
    if (!measure_col %in% names(df_long)) measure_col <- names(df_long)[grep("^norm$|^value_bs$", names(df_long))[1]]
    df_long <- df_long %>%
      dplyr::mutate(value_plot = .data[[measure_col]])
    
    # Auto-detect the correct y-axis units from norm_units (set from the ISTD CSV,
    # optionally suffixed with the protein unit when protein normalisation is on)
    # rather than letting the user pick an arbitrary, possibly-mismatched label.
    cbp_units_label <- NULL
    if (identical(measure_col, "norm")) {
      u <- stats::na.omit(unique(df_long$norm_units))
      if (length(u) == 1) {
        cbp_units_label <- u
      } else if (length(u) > 1) {
        showNotification(
          "Class bar plot: multiple ISTD units detected across data (mixed units).",
          type = "warning", duration = 6
        )
        cbp_units_label <- "(mixed units)"
      }
    }
    cbp_ylab <- axis_label(kind = measure_col, mode = "absolute", unit_label = cbp_units_label)
    
    # Restrict to features present in stats results
    df_long <- df_long %>% dplyr::filter(`Metabolite name` %in% res$metabolite)
    
    # Class filter — always filter to the selected class
    cls <- input$cbp_class
    req(nzchar(cls %||% ""))
    res     <- res     %>% dplyr::filter(plot_class == cls)
    df_long <- df_long %>% dplyr::filter(plot_class == cls)
    validate(need(nrow(res) > 0 && nrow(df_long) > 0,
                  "No features for this class. Run statistics first."))
    
    # Significant-only filter
    alpha   <- input$alpha %||% 0.05
    use_adj <- isTRUE(input$use_adj_threshold) && !identical(input$padj_method, "none")
    sig_mask <- if (use_adj) (!is.na(res$p_adj) & res$p_adj < alpha) else (!is.na(res$p) & res$p < alpha)
    sig_feats <- res$metabolite[sig_mask]
    
    if (isTRUE(input$cbp_sig_only)) {
      res     <- res[sig_mask, , drop = FALSE]
      df_long <- df_long %>% dplyr::filter(`Metabolite name` %in% sig_feats)
      validate(need(nrow(res) > 0, "No significant features under current thresholds."))
    }
    
    # Abundance range filter — based on mean value_plot across all samples
    ab_range <- input$cbp_abundance_range %||% c(0, 100)
    if (ab_range[1] > 0 || ab_range[2] < 100) {
      means_per_feat <- df_long %>%
        dplyr::group_by(`Metabolite name`) %>%
        dplyr::summarise(mean_val = mean(value_plot, na.rm = TRUE), .groups = "drop")
      max_val <- max(means_per_feat$mean_val, na.rm = TRUE)
      if (is.finite(max_val) && max_val > 0) {
        pct  <- means_per_feat$mean_val / max_val * 100
        keep <- means_per_feat$`Metabolite name`[pct >= ab_range[1] & pct <= ab_range[2]]
        df_long <- df_long %>% dplyr::filter(`Metabolite name` %in% keep)
        res     <- res     %>% dplyr::filter(metabolite %in% keep)
        sig_feats <- sig_feats[sig_feats %in% keep]
      }
      validate(need(nrow(df_long) > 0,
                    "No features remain after abundance filter. Adjust the range."))
    }
    
    df_long <- df_long %>%
      dplyr::mutate(significant = `Metabolite name` %in% sig_feats)
    
    # Use .compute_pairwise directly on cbp_data's own filtered long data.
    # Pass a clean minimal frame to avoid column name conflicts.
    ph_class <- tryCatch(
      .compute_pairwise(
        df_long = df_long %>%
          dplyr::transmute(
            `Metabolite name` = `Metabolite name`,
            bar_group         = group,
            value             = value_plot
          ),
        group_col   = "bar_group",
        padj_method = input$padj_method %||% "BH",
        alpha       = input$alpha %||% 0.05,
        use_adj     = isTRUE(input$use_adj_threshold) &&
          !identical(input$padj_method %||% "BH", "none")
      ),
      error = function(e) NULL
    )
    
    list(df = df_long, res = res, sig_feats = sig_feats,
         measure_col = measure_col, posthoc = ph_class, ylab = cbp_ylab)
  })
  
  # Info text
  output$cbpInfoUI <- renderUI({
    d <- tryCatch(cbp_data(), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    n_total <- dplyr::n_distinct(d$df$`Metabolite name`)
    n_sig   <- length(d$sig_feats)
    tags$p(tags$small(
      paste0(n_total, " species shown (", n_sig, " significant).")
    ))
  })
  
  # Build the class bar plot — single grouped bar plot, species on x-axis
  .build_cbp <- function(fsz = 11) {
    d        <- cbp_data()
    df       <- d$df
    err_type <- input$cbp_error_type  %||% "SEM"
    show_pts <- isTRUE(input$cbp_show_points)
    
    # Compute summary stats per feature × group
    summ <- df %>%
      dplyr::group_by(`Metabolite name`, group) %>%
      dplyr::summarise(
        mean_val = mean(value_plot, na.rm = TRUE),
        sd_val   = sd(value_plot, na.rm = TRUE),
        n_val    = sum(!is.na(value_plot)),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        se  = sd_val / sqrt(pmax(n_val, 1)),
        ci  = qt(0.975, df = pmax(n_val - 1, 1)) * se,
        err = switch(err_type, SEM = se, SD = sd_val, CI95 = ci, se)
      )
    
    # Species order — alphabetical or by total abundance
    sort_order <- input$cbp_sort_order %||% "alpha"
    feat_order <- if (identical(sort_order, "abundance")) {
      summ %>%
        dplyr::group_by(`Metabolite name`) %>%
        dplyr::summarise(total = sum(mean_val, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(total)) %>%
        dplyr::pull(`Metabolite name`)
    } else {
      sort(unique(as.character(summ$`Metabolite name`)))
    }
    summ <- summ %>%
      dplyr::mutate(
        species   = factor(`Metabolite name`, levels = feat_order),
        is_sig    = `Metabolite name` %in% d$sig_feats
      )
    
    # Group colour palette
    n_groups     <- dplyr::n_distinct(summ$group)
    group_levels <- sort(unique(summ$group))
    palette_id   <- input$cbp_palette %||% "OkabeIto"
    bar_alpha    <- input$cbp_bar_alpha  %||% 0.85
    pt_size      <- input$cbp_point_size %||% 2.2
    
    # Colour resolution
    okabe_ito  <- c("#E69F00","#56B4E9","#009E73","#F0E442",
                    "#0072B2","#D55E00","#CC79A7","#000000")
    tableau10  <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
    viridis_opts <- c("viridis","plasma","inferno","magma","cividis")
    brewer_qual  <- c("Set1","Set2","Dark2","Paired","Accent","Pastel1",
                      "RdBu","PuOr","BrBG","Blues","Greens","Purples")
    
    pal <- if (palette_id == "OkabeIto") {
      okabe_ito[seq_len(min(n_groups, length(okabe_ito)))]
    } else if (palette_id == "Tableau10") {
      tableau10[seq_len(min(n_groups, length(tableau10)))]
    } else if (palette_id %in% viridis_opts) {
      viridis::viridis(n_groups, option = palette_id, end = 0.85)
    } else if (palette_id %in% brewer_qual) {
      max_cols <- RColorBrewer::brewer.pal.info[palette_id, "maxcolors"]
      if (n_groups <= max_cols) {
        RColorBrewer::brewer.pal(max(3, n_groups), palette_id)[seq_len(n_groups)]
      } else {
        colorRampPalette(RColorBrewer::brewer.pal(max_cols, palette_id))(n_groups)
      }
    } else {
      # Fallback
      scales::hue_pal()(n_groups)
    }
    pal <- setNames(pal, group_levels)
    
    # Slightly darker border per bar — 55% luminance of fill
    darken <- function(col, fac = 0.55) {
      m <- col2rgb(col) / 255
      grDevices::rgb(m[1]*fac, m[2]*fac, m[3]*fac)
    }
    border_pal <- setNames(sapply(pal, darken), group_levels)
    
    dodge_w <- 0.8
    
    p <- ggplot2::ggplot(
      summ,
      ggplot2::aes(x = species, y = mean_val, fill = group, colour = group)
    ) +
      ggplot2::geom_col(
        position  = ggplot2::position_dodge(width = dodge_w),
        width     = 0.72,
        linewidth = 0.4,
        alpha     = bar_alpha
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(0, mean_val - err), ymax = mean_val + err),
        position  = ggplot2::position_dodge(width = dodge_w),
        width     = 0.18,
        linewidth = 0.5,
        colour    = "grey20"
      ) +
      ggplot2::scale_fill_manual(values = pal,        name = NULL) +
      ggplot2::scale_colour_manual(values = border_pal, guide = "none") +
      ggplot2::labs(
        x       = NULL,
        y       = d$ylab,
        caption = paste0("Error bars: ", err_type)
      ) +
      ggplot2::theme_minimal(base_size = fsz) +
      ggplot2::theme(
        axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
        axis.line.x        = ggplot2::element_line(colour = "grey70", linewidth = 0.4),
        axis.line.y        = ggplot2::element_line(colour = "grey70", linewidth = 0.4),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.35),
        legend.position    = "right",
        legend.key.size    = ggplot2::unit(0.8, "lines"),
        plot.margin        = ggplot2::margin(t = 8, r = 16, b = 8, l = 48, unit = "pt")
      ) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 1.2))
    
    # Individual points — fill and outline controlled by UI
    if (show_pts && nrow(df) > 0) {
      pt_fill_mode    <- input$cbp_point_fill    %||% "palette"
      pt_outline_mode <- input$cbp_point_outline %||% "white"
      
      df_pts <- df %>%
        dplyr::mutate(
          species = factor(`Metabolite name`, levels = feat_order),
          pt_fill = if (pt_fill_mode == "black") "#222222"
                    else if (pt_fill_mode == "white") "#ffffff"
                    else pal[as.character(group)]
        )
      
      # Shape 21 always — outline=none achieved by matching colour to fill
      pt_stroke <- if (identical(pt_outline_mode, "none")) 0 else 0.65
      
      p <- p + ggplot2::geom_point(
        data     = df_pts,
        ggplot2::aes(x = species, y = value_plot, group = group),
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.07,
          dodge.width  = dodge_w,
          seed         = 42
        ),
        shape       = 21,
        size        = pt_size,
        stroke      = pt_stroke,
        fill        = df_pts$pt_fill,
        colour      = if (pt_outline_mode == "none") df_pts$pt_fill
                      else if (pt_outline_mode == "black") "#222222"
                      else "white",
        alpha       = 0.92,
        inherit.aes = FALSE
      )
    }
    
    # Significance brackets — drawn from posthoc_for_plot data if available,
    # otherwise from overall stats p-value (for t-test / ANOVA without posthoc)
    if (length(d$sig_feats) > 0) {
      tryCatch({
        ph       <- d$posthoc   # may be NULL if no posthoc was run
        alpha_br <- input$alpha %||% 0.05
        use_adj_f <- isTRUE(input$use_adj_threshold) &&
          !identical(input$padj_method %||% "BH", "none")
        
        y_global_top <- max(summ$mean_val + summ$err, na.rm = TRUE)
        if (!is.finite(y_global_top)) y_global_top <- 1
        tip <- 0.025 * y_global_top
        bracket_tops <- numeric(0)
        
        # Helper: draw one bracket between two x positions
        draw_bracket <- function(p, x1, x2, y_br, y_cap, sym) {
          p +
            ggplot2::annotate("segment",
                              x = x1, xend = x1, y = y_br - tip, yend = y_cap,
                              colour = "grey30", linewidth = 0.4) +
            ggplot2::annotate("segment",
                              x = x2, xend = x2, y = y_br - tip, yend = y_cap,
                              colour = "grey30", linewidth = 0.4) +
            ggplot2::annotate("segment",
                              x = x1, xend = x2, y = y_cap, yend = y_cap,
                              colour = "grey30", linewidth = 0.4) +
            ggplot2::annotate("text",
                              x = (x1 + x2) / 2, y = y_cap + 0.01 * y_global_top,
                              label = sym, vjust = 0, size = fsz * 0.32, colour = "grey20")
        }
        
        # x offset of each bar within a dodged species group
        # With position_dodge(width=dodge_w) and n groups evenly spaced:
        # bar_centres[i] = sp_idx + dodge_w * ((i - 1) / (n-1) - 0.5) * (n-1)/n
        # Simplify: centres are sp_idx + seq(-0.5, 0.5, len=n)*dodge_w*(n-1)/n
        bar_offsets <- function(n) {
          if (n == 1) return(0)
          dodge_w * seq(-(n-1), (n-1), by = 2) / (2*n)
        }
        grp_offsets <- setNames(bar_offsets(length(group_levels)), group_levels)
        
        for (sp in d$sig_feats) {
          sp_idx  <- which(feat_order == sp)
          if (length(sp_idx) == 0) next
          sp_summ <- summ %>% dplyr::filter(`Metabolite name` == sp)
          y_top   <- max(sp_summ$mean_val + sp_summ$err, na.rm = TRUE)
          if (!is.finite(y_top)) next
          
          # Posthoc comparisons for this species
          sp_ph <- if (!is.null(ph))
            ph %>% dplyr::filter(`Metabolite name` == sp)
          else NULL
          
          if (!is.null(sp_ph) && nrow(sp_ph) > 0) {
            # Draw one bracket per significant pairwise comparison, stacking upward
            sp_ph <- sp_ph %>%
              dplyr::mutate(
                p_use = if (use_adj_f) dplyr::coalesce(p.adj, p) else p,
                sym   = dplyr::case_when(
                  p_use < 0.001 ~ "***",
                  p_use < 0.01  ~ "**",
                  p_use < 0.05  ~ "*",
                  TRUE          ~ "ns"
                )
              ) %>%
              dplyr::filter(sym != "ns",
                            group1 %in% group_levels,
                            group2 %in% group_levels) %>%
              dplyr::arrange(p_use)
            
            y_step <- 0.10 * y_global_top
            for (bi in seq_len(nrow(sp_ph))) {
              row   <- sp_ph[bi, ]
              x1    <- sp_idx + grp_offsets[as.character(row$group1)]
              x2    <- sp_idx + grp_offsets[as.character(row$group2)]
              y_br  <- y_top + 0.05 * y_global_top + (bi - 1) * y_step
              y_cap <- y_br + 0.04 * y_global_top
              bracket_tops <- c(bracket_tops, y_cap)
              p <- draw_bracket(p, x1, x2, y_br, y_cap, row$sym)
            }
            
          } else {
            # No posthoc — single bracket spanning all groups (or star for >2)
            sp_res  <- d$res %>% dplyr::filter(metabolite == sp)
            pval    <- if (nrow(sp_res) > 0) {
              if (use_adj_f && !is.na(sp_res$p_adj[1])) sp_res$p_adj[1]
              else sp_res$p[1]
            } else NA_real_
            sym <- dplyr::case_when(
              is.na(pval)  ~ "",
              pval < 0.001 ~ "***",
              pval < 0.01  ~ "**",
              pval < 0.05  ~ "*",
              TRUE         ~ "ns"
            )
            if (!nzchar(sym) || sym == "ns") next
            
            if (length(group_levels) == 2) {
              x1   <- sp_idx + grp_offsets[group_levels[1]]
              x2   <- sp_idx + grp_offsets[group_levels[2]]
              y_br  <- y_top + 0.05 * y_global_top
              y_cap <- y_br + 0.04 * y_global_top
              bracket_tops <- c(bracket_tops, y_cap)
              p <- draw_bracket(p, x1, x2, y_br, y_cap, sym)
            } else {
              # ≥3 groups, no posthoc: single star above species
              y_cap <- y_top + 0.08 * y_global_top
              bracket_tops <- c(bracket_tops, y_cap)
              p <- p + ggplot2::annotate("text",
                                         x = sp_idx, y = y_cap,
                                         label = sym, vjust = 0, size = fsz * 0.38, colour = "grey20")
            }
          }
        }
        
        if (length(bracket_tops) > 0) {
          p <- p + ggplot2::expand_limits(
            y = max(bracket_tops, na.rm = TRUE) * 1.1
          )
        }
      }, error = function(e) NULL)
    }
    
    p
  }
  
  # Fixed height for single plot
  cbp_height <- reactive({ 550 })
  
  output$classBarPlot <- renderPlot({
    validate(need(!is.null(stats_results_val()) && nrow(stats_results_val()) > 0,
                  "Run statistics first."))
    fsz <- input$cbp_export_fontsize %||% 11
    .build_cbp(fsz)
  }, height = function() cbp_height())
  
  # Export helpers
  .cbp_export_dims <- function() {
    px_w  <- input$cbp_export_width    %||% 1600
    px_h  <- input$cbp_export_height   %||% 1200
    dpi   <- input$cbp_export_dpi      %||% 300
    scale <- input$cbp_export_scale    %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  output$download_cbp_png <- downloadHandler(
    filename = function() {
      cls <- isolate(input$cbp_class %||% "all")
      paste0("class_barplots_", cls, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      dims <- isolate(.cbp_export_dims())
      fsz  <- isolate(input$cbp_export_fontsize %||% 11)
      p    <- isolate(.build_cbp(fsz))
      validate(need(!is.null(p), "No plot to export."))
      ggplot2::ggsave(file, plot = p,
                      width = dims$w, height = dims$h,
                      dpi = dims$dpi, device = "png")
    }
  )
  
  output$download_cbp_svg <- downloadHandler(
    filename = function() {
      cls <- isolate(input$cbp_class %||% "all")
      paste0("class_barplots_", cls, "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      dims <- isolate(.cbp_export_dims())
      fsz  <- isolate(input$cbp_export_fontsize %||% 11)
      p    <- isolate(.build_cbp(fsz))
      validate(need(!is.null(p), "No plot to export."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  # ========== VOLCANO PLOT ==========
  # Shared ggplot reactive — used by both renderPlotly (screen) and download handlers
  # Populate comparison selector — only shown for ANOVA results
  output$volcano_comparison_ui <- renderUI({
    res <- stats_results_all()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    is_anova <- any(res$test %in% c("oneway", "twoway"), na.rm = TRUE)
    if (!is_anova) return(NULL)
    
    ph <- tryCatch(posthoc_for_plot(), error = function(e) NULL)
    if (is.null(ph) || nrow(ph) == 0) {
      return(helpText(tags$small(
        tags$b("ANOVA results detected."),
        " Run post-hoc tests on the Statistics tab to enable pairwise volcano plots."
      )))
    }
    
    pairs <- ph %>%
      dplyr::mutate(pair = paste(group1, "vs", group2)) %>%
      dplyr::distinct(pair) %>%
      dplyr::pull(pair)
    
    tagList(
      selectInput("volcano_comparison",
                  "Pairwise comparison (ANOVA):",
                  choices  = c("Overall (ANOVA p-value)" = "overall", pairs),
                  selected = pairs[1])
    )
  })
  
  # Helper: build plot data for a specific pairwise comparison
  # Uses all tested features (filter_alpha = Inf) so the volcano shows
  # non-significant points too, not just those that passed posthoc filtering.
  .volcano_pairwise_data <- function(pair_str) {
    df_long <- tryCatch(stats_input_long(), error = function(e) NULL)
    if (is.null(df_long)) return(NULL)
    
    parts <- strsplit(pair_str, " vs ")[[1]]
    g1    <- trimws(parts[1])
    g2    <- trimws(parts[2])
    
    df_pair <- df_long %>%
      dplyr::filter(group %in% c(g1, g2))
    if (nrow(df_pair) == 0) return(NULL)
    
    # Compute group means and log2FC for every feature
    means <- df_pair %>%
      dplyr::group_by(`Metabolite name`, plot_class, group) %>%
      dplyr::summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = group, values_from = mean_val)
    
    if (!all(c(g1, g2) %in% names(means))) return(NULL)
    
    means <- means %>%
      dplyr::mutate(
        log2FC = log2(((.data[[g1]] + 1e-12) / (.data[[g2]] + 1e-12)))
      ) %>%
      dplyr::rename(metabolite = `Metabolite name`)
    
    # Pairwise stats for ALL features — no significance filter
    df_feats <- df_pair %>%
      dplyr::rename(bar_group = group) %>%
      dplyr::select(`Metabolite name`, bar_group, value)
    
    ph_all <- tryCatch(
      .compute_pairwise(
        df_long     = df_feats,
        group_col   = "bar_group",
        padj_method = input$padj_method %||% "BH",
        alpha       = input$alpha       %||% 0.05,
        use_adj     = isTRUE(input$use_adj_threshold) &&
          !identical(input$padj_method %||% "BH", "none"),
        filter_alpha = Inf   # return all features, not just significant
      ),
      error = function(e) NULL
    )
    if (is.null(ph_all) || nrow(ph_all) == 0) return(NULL)
    
    ph_pair <- ph_all %>%
      dplyr::filter(
        (as.character(group1) == g1 & as.character(group2) == g2) |
          (as.character(group1) == g2 & as.character(group2) == g1)
      ) %>%
      dplyr::mutate(metabolite = as.character(`Metabolite name`)) %>%
      dplyr::select(metabolite, p, p.adj, p.adj.signif)
    
    if (nrow(ph_pair) == 0) return(NULL)
    
    means %>%
      dplyr::inner_join(ph_pair, by = "metabolite") %>%
      dplyr::select(metabolite, plot_class, log2FC, p, p.adj)
  }
  
  # Core volcano plot builder — called by both the reactive and export handlers.
  # add_repel_labels = TRUE  → ggrepel labels (for PNG/SVG export)
  # add_repel_labels = FALSE → no label layer (for plotly; annotations added separately)
  .build_volcano_ggplot <- function(add_repel_labels = TRUE) {
    res    <- stats_results_all()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    alpha  <- input$volcano_alpha  %||% 0.05
    fc_thr <- input$volcano_log2fc %||% 1
    use_adj <- isTRUE(input$volcano_use_padj)
    
    # Determine if we should use pairwise mode
    is_anova  <- any(res$test %in% c("oneway", "twoway"), na.rm = TRUE)
    pair_sel  <- input$volcano_comparison %||% "overall"
    use_pair  <- is_anova && !identical(pair_sel, "overall") && nzchar(pair_sel)
    
    if (use_pair) {
      plot_df <- .volcano_pairwise_data(pair_sel)
      validate(need(!is.null(plot_df) && nrow(plot_df) > 0,
                    "No pairwise data for this comparison. Run post-hoc tests on the Statistics tab."))
      p_col    <- if (use_adj && "p.adj" %in% names(plot_df)) "p.adj" else "p"
      subtitle <- paste0("Comparison: ", pair_sel)
    } else {
      if (!"log2FC" %in% names(res)) {
        validate(need(FALSE,
                      "Volcano plot requires log2FC — only available for 2-group (t-test) results. ",
                      "Select a pairwise comparison above for ANOVA data."))
      }
      plot_df  <- res
      p_col    <- if (use_adj && "p_adj" %in% names(plot_df)) "p_adj" else "p"
      subtitle <- NULL
    }
    
    plot_df <- plot_df %>%
      dplyr::mutate(
        neg_log10_p = -log10(pmax(.data[[p_col]], 1e-300)),
        direction   = dplyr::case_when(
          .data[[p_col]] < alpha & log2FC >  fc_thr ~ "Up",
          .data[[p_col]] < alpha & log2FC < -fc_thr ~ "Down",
          TRUE ~ "NS"
        ),
        hover_text = paste0(
          "Feature: ", metabolite,
          "\nClass: ", plot_class,
          "\nlog2FC: ", round(log2FC, 3),
          "\n", p_col, ": ", signif(.data[[p_col]], 3)
        )
      )
    
    # --- Axis capping ---
    # When enabled, features beyond the cap are plotted AT the cap with a
    # triangle marker so the compressed region remains readable while no
    # data are hidden. The original values are preserved in hover_text.
    cap_axes <- isTRUE(input$volcano_cap_axes)
    cap_y    <- if (cap_axes) (input$volcano_cap_y %||% 10)  else Inf
    cap_x    <- if (cap_axes) (input$volcano_cap_x %||% 5)   else Inf
    
    plot_df <- plot_df %>%
      dplyr::mutate(
        capped = cap_axes & (neg_log10_p > cap_y | abs(log2FC) > cap_x),
        # Append capping note to hover text for capped points
        hover_text = dplyr::if_else(
          capped,
          paste0(hover_text, "\n[axis-capped for display]"),
          hover_text
        ),
        neg_log10_p = pmin(neg_log10_p, cap_y),
        log2FC      = pmax(pmin(log2FC,  cap_x), -cap_x),
        pt_shape    = dplyr::if_else(capped, "capped", "normal")
      )
    
    shape_vals <- c(normal = 16, capped = 17)  # circle / filled triangle
    
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = log2FC, y = neg_log10_p,
                   colour = direction, shape = pt_shape, text = hover_text)
    ) +
      ggplot2::geom_point(alpha = 0.7, size = 1.8) +
      ggplot2::scale_shape_manual(
        values = shape_vals,
        guide  = if (cap_axes) ggplot2::guide_legend(title = NULL) else "none"
      ) +
      ggplot2::geom_vline(xintercept = c(-fc_thr, fc_thr),
                          linetype = "dashed", colour = "#888888") +
      ggplot2::geom_hline(yintercept = -log10(alpha),
                          linetype = "dashed", colour = "#888888") +
      ggplot2::scale_colour_manual(
        values = c(Up = "#E45756", Down = "#4C78A8", NS = "#AAAAAA")
      ) +
      ggplot2::labs(
        x        = "log2 Fold Change",
        y        = paste0("-log10(", p_col, ")"),
        colour   = "Direction",
        subtitle = subtitle
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::coord_cartesian(clip = "off")
    
    if (isTRUE(input$volcano_label_sig) && add_repel_labels) {
      top_n  <- input$volcano_topn_labels %||% 20
      lab_df <- plot_df %>%
        dplyr::filter(direction != "NS") %>%
        dplyr::arrange(.data[[p_col]]) %>%
        dplyr::slice_head(n = top_n)
      if (nrow(lab_df) > 0) {
        # Expand plot limits to give ggrepel room to place labels without
        # clipping. When axis capping is active, use the cap as the boundary
        # so labels don't try to position themselves beyond the cap line.
        x_range <- range(plot_df$log2FC,      na.rm = TRUE)
        y_range <- range(plot_df$neg_log10_p, na.rm = TRUE)
        x_max   <- if (cap_axes) cap_x  else x_range[2]
        x_min   <- if (cap_axes) -cap_x else x_range[1]
        y_max   <- if (cap_axes) cap_y  else y_range[2]
        x_pad   <- diff(c(x_min, x_max)) * 0.15
        y_pad   <- diff(c(y_range[1], y_max)) * 0.15
        
        p <- p + ggrepel::geom_text_repel(
          data           = lab_df,
          ggplot2::aes(label = metabolite),
          size           = 2.8,
          max.overlaps   = Inf,       # never silently drop labels
          box.padding    = 0.4,
          point.padding  = 0.3,
          force          = 2,
          force_pull     = 0.5,
          direction      = "both",
          segment.colour = "grey60",
          segment.size   = 0.3,
          segment.curvature = 0.1,
          show.legend    = FALSE,
          seed           = 42,
          xlim           = c(x_min - x_pad, x_max + x_pad),
          ylim           = c(y_range[1],    y_max + y_pad)
        )
      }
    }
    
    # Attach plot_df and p_col as attributes so renderPlotly can add annotations
    attr(p, "plot_df") <- plot_df
    attr(p, "p_col")   <- p_col
    p
  }
  
  volcano_ggplot <- reactive({
    .build_volcano_ggplot(add_repel_labels = TRUE)
  })
  
  output$volcanoPlot <- renderPlotly({
    # Build base plot without ggrepel (not supported by plotly)
    p <- .build_volcano_ggplot(add_repel_labels = FALSE)
    validate(need(!is.null(p), "Run statistics first (click Run statistics button)."))
    
    pl <- plotly::ggplotly(p, tooltip = "text")
    
    # Add plotly-native text annotations for top-N significant features
    if (isTRUE(input$volcano_label_sig)) {
      plot_df <- attr(p, "plot_df")
      p_col   <- attr(p, "p_col")
      top_n   <- input$volcano_topn_labels %||% 20
      
      if (!is.null(plot_df) && !is.null(p_col)) {
        lab_df <- plot_df %>%
          dplyr::filter(direction != "NS") %>%
          dplyr::arrange(.data[[p_col]]) %>%
          dplyr::slice_head(n = top_n)
        
        if (nrow(lab_df) > 0) {
          # Fan arrow offsets so labels spread out rather than all pointing the
          # same direction — prevents off-screen clipping near plot edges.
          n_lab   <- nrow(lab_df)
          angles  <- seq(0, 2 * pi, length.out = n_lab + 1)[seq_len(n_lab)]
          ax_vals <- round(cos(angles) * 40)
          ay_vals <- round(sin(angles) * -40)  # negative = up in plotly coords
          
          pl <- pl %>%
            plotly::add_annotations(
              x          = lab_df$log2FC,
              y          = lab_df$neg_log10_p,
              text       = lab_df$metabolite,
              xref       = "x", yref = "y",
              showarrow  = TRUE,
              arrowhead  = 2,
              arrowsize  = 0.5,
              arrowcolor = "grey60",
              ax         = ax_vals,
              ay         = ay_vals,
              font       = list(size = 9),
              bgcolor    = "rgba(255,255,255,0.7)",
              borderpad  = 2
            )
        }
      }
    }
    pl
  })
  
  
  .volc_export_dims <- function() {
    px_w  <- input$volc_export_width  %||% 1200
    px_h  <- input$volc_export_height %||% 700
    dpi   <- input$volc_export_dpi    %||% 300
    scale <- input$volc_export_scale  %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  output$download_volcano_png <- downloadHandler(
    filename = function() paste0("volcano_", isolate(input$stats_class), "_", Sys.Date(), ".png"),
    content = function(file) {
      dims <- isolate(.volc_export_dims())
      p    <- isolate(volcano_ggplot())
      validate(need(!is.null(p), "No plot to export. Run statistics first."))
      ggplot2::ggsave(file, plot = p,
                      width = dims$w, height = dims$h,
                      dpi = dims$dpi, device = "png")
    }
  )
  
  output$download_volcano_svg <- downloadHandler(
    filename = function() paste0("volcano_", isolate(input$stats_class), "_", Sys.Date(), ".svg"),
    content = function(file) {
      dims <- isolate(.volc_export_dims())
      p    <- isolate(volcano_ggplot())
      validate(need(!is.null(p), "No plot to export. Run statistics first."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  # ---- Statistics plot export helpers ----
  # Shared helper: build dimensions from the new MetaboRich-style inputs
  .stats_export_dims <- function() {
    px_w  <- input$stats_export_width    %||% 1200
    px_h  <- input$stats_export_height   %||% 700
    dpi   <- input$stats_export_dpi      %||% 300
    scale <- input$stats_export_scale    %||% 1.0
    list(w = (px_w / dpi) * scale,
         h = (px_h / dpi) * scale,
         dpi = dpi)
  }
  
  # Helper: rebuild the currently displayed single-feature stats plot
  .current_stats_plot <- function() {
    res <- stats_results_val()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    idx <- max(1, min(isolate(current_plot_index()), nrow(res)))
    met <- res$metabolite[idx]
    tryCatch(
      .build_stats_barplot(
        met             = met,
        idx             = idx,
        n_total         = nrow(res),
        df_long         = stats_input_long(),
        res_all         = res,
        err_type        = input$stats_error_type %||% "SEM",
        show_pts        = isTRUE(input$stats_plot_points),
        show_two_way    = any(res$test == "twoway"),
        alpha           = input$alpha %||% 0.05,
        padj_method     = input$padj_method %||% "BH",
        use_adj         = isTRUE(input$use_adj_threshold),
        posthoc_choices = input$posthoc %||% character(0),
        base_size       = isolate(input$stats_export_fontsize %||% 14)
      ),
      error = function(e) NULL
    )
  }
  
  # PNG export of the currently displayed single-feature plot
  output$download_stats_png <- downloadHandler(
    filename = function() {
      paste0("LipiRich_stats_", isolate(input$stats_class), "_", Sys.Date(), ".png")
    },
    content = function(file) {
      dims <- isolate(.stats_export_dims())
      p    <- isolate(.current_stats_plot())
      validate(need(!is.null(p), "No plot to export. Run statistics first."))
      ggplot2::ggsave(file, plot = p,
                      width  = dims$w, height = dims$h,
                      dpi    = dims$dpi, device = "png")
    }
  )
  
  # SVG export of the currently displayed single-feature plot
  output$download_stats_svg <- downloadHandler(
    filename = function() {
      paste0("LipiRich_stats_", isolate(input$stats_class), "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      dims <- isolate(.stats_export_dims())
      p    <- isolate(.current_stats_plot())
      validate(need(!is.null(p), "No plot to export. Run statistics first."))
      svglite::svglite(file, width = dims$w, height = dims$h)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(p)
    }
  )
  
  # ---- Statistics PDF export — one significant plot per page ----
  output$download_stats_pdf <- downloadHandler(
    filename    = function() {
      paste0("LipiRich_significant_", isolate(input$stats_class), "_", Sys.Date(), ".pdf")
    },
    contentType = "application/pdf",
    content = function(file) {
      res <- isolate(stats_results_val())
      validate(need(!is.null(res) && nrow(res) > 0,
                    "No significant results to export. Run statistics first."))
      
      px_w      <- isolate(input$stats_export_width    %||% 1200)
      px_h      <- isolate(input$stats_export_height   %||% 700)
      dpi       <- isolate(input$stats_export_dpi      %||% 300)
      scale     <- isolate(input$stats_export_scale    %||% 1.0)
      base_sz   <- isolate(input$stats_export_fontsize %||% 14)
      layout    <- isolate(input$stats_pdf_scale %||% "single")
      per_page  <- switch(layout, single = 1L, two = 2L, four = 4L, 1L)
      # Convert px → inches for grDevices::pdf(), apply scale fraction
      pw        <- (px_w / dpi) * scale
      ph        <- (px_h / dpi) * scale
      err_type  <- isolate(input$stats_error_type %||% "SEM")
      show_pts  <- isTRUE(isolate(input$stats_plot_points))
      show_two_way <- any(res$test == "twoway")
      alpha     <- isolate(input$alpha %||% 0.05)
      padj_meth <- isolate(input$padj_method %||% "BH")
      use_adj   <- isTRUE(isolate(input$use_adj_threshold))
      phoc      <- isolate(input$posthoc %||% character(0))
      df_long   <- isolate(stats_input_long())
      n_total   <- nrow(res)
      
      # Page dimensions: for multi-plot layouts, expand the page
      page_w <- if (per_page >= 2) pw * 2 else pw
      page_h <- if (per_page == 4) ph * 2 else ph
      
      tmp_pdf <- tempfile(fileext = ".pdf")
      on.exit({
        tryCatch(grDevices::dev.off(), error = function(e) NULL)
        if (file.exists(tmp_pdf)) file.copy(tmp_pdf, file, overwrite = TRUE)
      }, add = TRUE)
      
      grDevices::pdf(tmp_pdf, width = page_w, height = page_h)
      
      withProgress(message = "Building PDF…", value = 0, {
        plots_buffer <- list()
        
        flush_page <- function(buf) {
          if (length(buf) == 0) return(invisible(NULL))
          if (length(buf) == 1) {
            tryCatch(print(buf[[1]]), error = function(e) NULL)
          } else {
            tryCatch(
              print(gridExtra::grid.arrange(grobs = buf,
                                            ncol  = if (per_page == 2L) 2L else 2L,
                                            nrow  = if (per_page == 4L) 2L else 1L)),
              error = function(e) NULL
            )
          }
        }
        
        for (i in seq_len(n_total)) {
          incProgress(i / n_total, detail = paste("Plot", i, "of", n_total))
          met <- res$metabolite[i]
          p   <- tryCatch(
            .build_stats_barplot(
              met           = met,
              idx           = i,
              n_total       = n_total,
              df_long       = df_long,
              res_all       = res,
              err_type      = err_type,
              show_pts      = show_pts,
              show_two_way  = show_two_way,
              alpha         = alpha,
              padj_method   = padj_meth,
              use_adj       = use_adj,
              posthoc_choices = phoc,
              base_size     = if (per_page == 1L) base_sz else max(6L, base_sz - 3L)
            ),
            error = function(e) NULL
          )
          if (!is.null(p)) {
            plots_buffer <- c(plots_buffer, list(p))
            if (length(plots_buffer) == per_page) {
              flush_page(plots_buffer)
              plots_buffer <- list()
            }
          }
        }
        # Flush any remaining plots
        if (length(plots_buffer) > 0) flush_page(plots_buffer)
      })
    }
  )
  
  # ========== LIPID NETWORK TAB ==========
  
  # Populate network group selector from the same central reactiveVal as PCA/Stats
  observe({
    grps <- net_available_groups()
    cur  <- isolate(input$net_group)
    updateSelectInput(session, "net_group",
                      choices  = grps,
                      selected = if (!is.null(cur) && cur %in% grps) cur
                      else if (length(grps) > 0) grps[1]
                      else NULL)
  })
  
  # Store edge list for download
  net_edge_data <- reactiveVal(NULL)
  
  observeEvent(input$net_run, {
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
      showNotification(
        "visNetwork is required. Install with: install.packages('visNetwork')",
        type = "error"); return()
    }
    
    # ── Require statistics have been run ──────────────────────────────────────
    sm_sig <- stats_results_val()
    if (is.null(sm_sig) || nrow(sm_sig) == 0) {
      showNotification("Run statistics first (Statistics tab).", type = "warning")
      return()
    }
    
    grp      <- input$net_group
    alpha    <- input$net_alpha    %||% 0.05
    r_thresh <- input$net_r_thresh %||% 0.7
    max_n    <- input$net_max_nodes %||% 50
    
    if (is.null(grp) || !nzchar(grp)) {
      showNotification("Select a group.", type = "warning"); return()
    }
    
    # ── Pick significant features ─────────────────────────────────────────────
    p_col <- if ("p_adj" %in% names(sm_sig)) "p_adj" else "p"
    sm_sig <- sm_sig %>%
      dplyr::filter(!is.na(.data[[p_col]]), .data[[p_col]] < alpha) %>%
      dplyr::arrange(.data[[p_col]]) %>%
      dplyr::slice_head(n = max_n)
    
    if (nrow(sm_sig) == 0) {
      showNotification("No significant features under current threshold.", type = "warning")
      return()
    }
    feats <- unique(sm_sig$metabolite)
    
    # ── Build wide matrix for selected group ──────────────────────────────────
    measure_col <- pick_measure_col(input$net_value_type)
    
    df <- bg_norm_long_avg() %>%
      dplyr::filter(!stringr::str_detect(`Metabolite name`, "\\[IS\\]")) %>%
      dplyr::mutate(sample_norm = dplyr::coalesce(sample_norm, normalize_sample_name(sample)))
    fn  <- get_active_grouping()
    res <- fn(df$sample)
    df  <- df %>% dplyr::mutate(group = res$group)
    df <- filter_blanks(df, isTRUE(input$net_exclude_blank),
                        sample_col = "sample_norm", exact = FALSE)
    df <- filter_iqc(df, include_iqc = FALSE, sample_col = "sample")
    
    wide <- df %>%
      dplyr::filter(`Metabolite name` %in% feats, group == grp) %>%
      dplyr::select(sample, `Metabolite name`, dplyr::all_of(measure_col)) %>%
      tidyr::pivot_wider(names_from  = `Metabolite name`,
                         values_from = dplyr::all_of(measure_col)) %>%
      tibble::column_to_rownames("sample")
    
    if (nrow(wide) < 3) {
      showNotification(
        sprintf("Need >=3 samples in group '%s' (got %d).", grp, nrow(wide)),
        type = "warning"); return()
    }
    
    M <- as.matrix(wide)
    M[!is.finite(M)] <- NA
    M <- M[, colSums(is.finite(M)) >= 3, drop = FALSE]
    # Drop zero-variance (constant) features so cor() does not warn or return NA
    if (ncol(M) > 0) {
      keep_var <- vapply(seq_len(ncol(M)),
                         function(j) { v <- stats::var(M[, j], na.rm = TRUE); is.finite(v) && v > 0 },
                         logical(1))
      M <- M[, keep_var, drop = FALSE]
    }
    
    if (ncol(M) < 2) {
      showNotification("Too few complete features for correlation.", type = "warning")
      return()
    }
    
    cr <- tryCatch(
      stats::cor(M, use = "pairwise.complete.obs"),
      error = function(e) NULL
    )
    if (is.null(cr)) {
      showNotification("Correlation matrix failed.", type = "error"); return()
    }
    
    # ── Build edge list ───────────────────────────────────────────────────────
    fn <- colnames(cr)
    edge_list <- do.call(rbind, Filter(Negate(is.null),
                                       lapply(seq_len(nrow(cr) - 1), function(i) {
                                         do.call(rbind, Filter(Negate(is.null),
                                                               lapply((i + 1):ncol(cr), function(j) {
                                                                 r <- cr[i, j]
                                                                 if (!is.na(r) && abs(r) >= r_thresh)
                                                                   data.frame(from = fn[i], to = fn[j],
                                                                              r = r, abs_r = abs(r),
                                                                              stringsAsFactors = FALSE)
                                                                 else NULL
                                                               })
                                         ))
                                       })
    ))
    
    if (is.null(edge_list) || nrow(edge_list) == 0) {
      showNotification(
        sprintf("No edges above |r| = %.2f. Try lowering the threshold.", r_thresh),
        type = "warning"); return()
    }
    net_edge_data(edge_list)
    
    # ── Node properties ───────────────────────────────────────────────────────
    node_feats <- unique(c(edge_list$from, edge_list$to))
    
    # Direction: determine up/down in selected group using log2FC if available,
    # else compare group mean to grand mean
    node_dirs <- vapply(node_feats, function(f) {
      row <- sm_sig[sm_sig$metabolite == f, ]
      if (nrow(row) == 0) return(NA_character_)
      if ("log2FC" %in% names(row) && is.finite(row$log2FC[1])) {
        # log2FC: comp/ref where ref = grp_ref, comp = grp_comp from stats
        fc <- row$log2FC[1]
        # Direction relative to selected group: if selected grp == grp_comp, fc>0 means up
        if (!is.na(fc)) return(if (fc > 0) "up" else "down")
      }
      # Fallback: compare group mean to grand mean
      gm <- df %>%
        dplyr::filter(`Metabolite name` == f) %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(m = mean(.data[[measure_col]], na.rm = TRUE), .groups = "drop")
      if (nrow(gm) == 0) return(NA_character_)
      grand <- mean(gm$m, na.rm = TRUE)
      grp_m <- gm$m[gm$group == grp]
      if (length(grp_m) == 0 || !is.finite(grp_m)) return(NA_character_)
      if (grp_m > grand) "up" else "down"
    }, character(1))
    
    node_padj <- sm_sig[[p_col]][match(node_feats, sm_sig$metabolite)]
    
    # Lipid class for each node
    node_class <- vapply(node_feats, function(f) {
      pc <- df$plot_class[df$`Metabolite name` == f]
      pc <- pc[!is.na(pc)]
      if (length(pc) > 0) pc[1] else "Unknown"
    }, character(1))
    
    node_short <- ifelse(
      nchar(node_feats) > 25,
      paste0(substr(node_feats, 1, 23), "..."),
      node_feats
    )
    
    node_col <- dplyr::case_when(
      node_dirs == "up"   ~ "#D85A30",
      node_dirs == "down" ~ "#1D9E75",
      TRUE                ~ "#B4B2A9"
    )
    
    # Node size: larger = more significant
    node_size <- pmax(8, 22 - pmin(12,
                                   -log10(dplyr::coalesce(node_padj, 1) + 1e-10) * 2
    ))
    
    nodes <- data.frame(
      id    = node_feats,
      label = node_short,
      title = paste0(
        "<b>", node_feats, "</b><br>",
        "Class: ", node_class, "<br>",
        "adj p: ", signif(dplyr::coalesce(node_padj, NA_real_), 3), "<br>",
        "Direction (", grp, "): ",
        dplyr::coalesce(node_dirs, "n/a")
      ),
      color = node_col,
      size  = node_size,
      stringsAsFactors = FALSE
    )
    
    vis_edges <- data.frame(
      from  = edge_list$from,
      to    = edge_list$to,
      value = edge_list$abs_r,
      color = ifelse(edge_list$r > 0, "#1D9E75", "#D85A30"),
      title = paste0("r = ", round(edge_list$r, 3)),
      stringsAsFactors = FALSE
    )
    
    # ── Render visNetwork ─────────────────────────────────────────────────────
    output$net_plot <- visNetwork::renderVisNetwork({
      visNetwork::visNetwork(nodes, vis_edges,
                             width = "100%", height = "600px") %>%
        visNetwork::visOptions(
          highlightNearest = list(enabled = TRUE, degree = 1),
          nodesIdSelection = FALSE
        ) %>%
        visNetwork::visPhysics(
          solver = "forceAtlas2Based",
          forceAtlas2Based = list(gravitationalConstant = -60)
        ) %>%
        visNetwork::visLayout(randomSeed = 42) %>%
        visNetwork::visEdges(
          smooth = FALSE,
          scaling = list(min = 1, max = 6)
        ) %>%
        visNetwork::visNodes(
          borderWidth = 0.5,
          font = list(size = 11)
        ) %>%
        visNetwork::visInteraction(
          dragNodes   = TRUE,
          zoomView    = TRUE,
          tooltipDelay = 150,
          navigationButtons = TRUE
        ) %>%
        visNetwork::visLegend(
          addNodes = data.frame(
            label = c(paste0("Up in ", grp), paste0("Down in ", grp), "No direction"),
            color = c("#D85A30", "#1D9E75", "#B4B2A9"),
            stringsAsFactors = FALSE
          ),
          addEdges = data.frame(
            label = c("Positive corr.", "Negative corr."),
            color = c("#1D9E75", "#D85A30"),
            stringsAsFactors = FALSE
          ),
          useGroups = FALSE
        )
    })
    
    showNotification(
      sprintf("Network: %d nodes, %d edges (|r| >= %.2f, group = %s).",
              nrow(nodes), nrow(vis_edges), r_thresh, grp),
      type = "message", duration = 6
    )
  })
  
  # Download UI (only shown after network is built)
  output$net_download_ui <- renderUI({
    if (!is.null(net_edge_data()))
      downloadButton("net_download_csv", "Download edge list (CSV)",
                     class = "btn-sm btn-default")
  })
  
  output$net_download_csv <- downloadHandler(
    filename = function() paste0("lipid_network_edges_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(net_edge_data(), file)
  )
  
  output$net_table <- DT::renderDT({
    req(net_edge_data())
    df <- net_edge_data() %>%
      dplyr::mutate(r = round(r, 4), abs_r = round(abs_r, 4)) %>%
      dplyr::arrange(dplyr::desc(abs_r)) %>%
      dplyr::rename(
        `Feature A` = from, `Feature B` = to,
        `Pearson r`  = r,   `|r|` = abs_r
      )
    DT::datatable(df, options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
}

shinyApp(ui = ui, server = server)