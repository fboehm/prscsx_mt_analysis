#!/usr/bin/env Rscript
# aou/04_stack_evaluate.R
# -----------------------
# Inside the All of Us Researcher Workbench: build every method from the
# candidate scores, fit ridge weights in a tuning set, and compare the methods
# in a held-out validation set, separately for each ancestry x outcome.
#
# Methods (same definitions as sim_sweep/compare_prsxtra.py):
#   prscsx           PRS-CSx, ancestry/trait-matched score, no tuning
#   mtag_prscsx      MTAG -> PRS-CSx, ancestry/trait-matched score, no tuning
#   prscsx_mt        PRS-CSx-MT, ancestry/trait-matched score, no tuning
#   prsxa            ridge over the PRS-CSx scores of the target trait (all ancestries)
#   prsxtra          ridge over all MTAG -> PRS-CSx scores (all traits x ancestries)
#   prscsx_mt_ridge  ridge over all PRS-CSx-MT scores (all traits x ancestries)
#
# Ridge: glmnet, alpha = 0, 10-fold CV, lambda.min (He et al.); covariates
# enter unpenalized (penalty.factor = 0) and the method's score is the PRS part
# of the linear predictor. All ridge methods share the tuning set and the CV
# folds. The untuned methods need an ancestry-matched score, so they are skipped
# for ancestries without a GWAS (e.g. MID, SAS) unless --pop_map maps them.
#
# Validation metrics:
#   continuous  incremental R2 = R2(y ~ cov + score) - R2(y ~ cov); paired
#               bootstrap over validation participants for contrasts
#   binary      AUC of the score (DeLong CI; paired DeLong test for contrasts),
#               OR per SD from y ~ score + cov
#
# Output (all aggregate; cells with any count < 20 are flagged for the All of Us
# dissemination policy and their counts masked):
#   <out_dir>/metrics.tsv  contrasts.tsv  ridge_coefs.tsv  cells.tsv  report.md
#
# Usage:
#   Rscript aou/04_stack_evaluate.R --scores scores/lipids.scores.tsv.gz \
#     --pheno pheno_lipids.tsv --manifest candidates/manifest.json --out_dir eval_lipids
#   options: --n_splits N --n_boot B --seed S --pop_map MID=EUR,SAS=EUR --threads T

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
  library(jsonlite)
})

# ── arguments ────────────────────────────────────────────────────────────────
parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    out[[key]] <- argv[i + 1]
    i <- i + 2
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
for (req in c("scores", "pheno", "manifest", "out_dir")) {
  if (is.null(args[[req]])) stop("missing --", req)
}
man <- fromJSON(args$manifest, simplifyVector = TRUE)
aou <- man$aou
opt <- function(name, default) {
  if (!is.null(args[[name]])) return(type.convert(args[[name]], as.is = TRUE))
  if (!is.null(aou[[name]])) return(aou[[name]])
  default
}
n_splits    <- opt("n_splits", 1)
n_boot      <- opt("n_boot", 1000)
seed        <- opt("seed", 2026)
tune_frac   <- opt("tune_frac", 0.7)
min_n       <- opt("min_group_n", 200)
min_cases   <- opt("min_group_cases", 50)
covariates  <- aou$covariates
outcomes    <- aou$outcomes
pops        <- man$populations
traits      <- man$traits
pop_map <- list()
if (!is.null(args$pop_map)) {
  for (kv in strsplit(strsplit(args$pop_map, ",")[[1]], "=")) pop_map[[kv[1]]] <- kv[2]
}
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
SUPPRESS <- 20

# ── data ─────────────────────────────────────────────────────────────────────
read_tsv <- function(path) {
  # .gz via gzip so data.table does not need R.utils
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)),
                                    colClasses = list(character = "person_id"))
  else fread(path, colClasses = list(character = "person_id"))
}
S <- read_tsv(args$scores)
P <- read_tsv(args$pheno)
D <- merge(P, S, by = "person_id")
cat(sprintf("%d participants with phenotypes and scores\n", nrow(D)))
missing_cov <- setdiff(covariates, names(D))
if (length(missing_cov)) stop("covariates missing from pheno: ", paste(missing_cov, collapse = ", "))

col_of <- function(fam, trait, pop) sprintf("%s__%s__%s", fam, trait, pop)
all_of <- function(fam) grep(sprintf("^%s__", fam), names(S), value = TRUE)

# ── helpers ──────────────────────────────────────────────────────────────────
zs <- function(x, mu = mean(x), s = sd(x)) if (s > 0) (x - mu) / s else x * 0

fit_ridge <- function(Xp, C, y, family, foldid) {
  # Xp: tuning PRS matrix (standardized); C: covariate matrix (standardized)
  X <- cbind(Xp, C)
  pf <- c(rep(1, ncol(Xp)), rep(0, ncol(C)))
  cv <- cv.glmnet(X, y, alpha = 0, family = family, foldid = foldid,
                  penalty.factor = pf, standardize = FALSE,
                  type.measure = if (family == "binomial") "auc" else "mse")
  b <- as.matrix(coef(cv, s = "lambda.min"))[, 1]
  list(beta = b[colnames(Xp)], lambda = cv$lambda.min)
}

r2 <- function(y, X) {
  f <- .lm.fit(cbind(1, X), y)
  1 - sum(f$residuals^2) / sum((y - mean(y))^2)
}

inc_r2 <- function(y, C, s) r2(y, cbind(C, s)) - r2(y, C)

# ── one ancestry x outcome x split ───────────────────────────────────────────
run_cell <- function(dat, outcome, spec, anc, split_id) {
  type <- spec$type
  family <- if (type == "binary") "binomial" else "gaussian"
  y <- dat[[outcome]]
  if (!is.null(spec$statin_adjust) && "statin" %in% names(dat)) {
    y <- ifelse(dat$statin == 1, y / spec$statin_adjust, y)
  }
  if (identical(spec$transform, "log")) y <- log(y)

  set.seed(seed + 1000 * split_id + match(anc, sort(unique(D$ancestry))))
  strat <- if (type == "binary") y else rep(0, length(y))
  tune <- unlist(lapply(split(seq_along(y), strat), function(ix)
    ix[sample.int(length(ix), round(tune_frac * length(ix)))]))
  is_tune <- seq_along(y) %in% tune
  val <- which(!is_tune)
  foldid <- sample(rep(1:10, length.out = sum(is_tune)))

  Cm <- as.matrix(dat[, ..covariates])
  cm <- colMeans(Cm[is_tune, , drop = FALSE]); cs <- apply(Cm[is_tune, , drop = FALSE], 2, sd)
  cs[cs == 0] <- 1
  Cz <- sweep(sweep(Cm, 2, cm), 2, cs, "/")

  tpop <- if (anc %in% pops) anc else pop_map[[anc]]
  mtrait <- spec$matched_trait
  stacks <- list(prsxa = col_of("prscsx", mtrait, pops),
                 prsxtra = all_of("mtag_prscsx"),
                 prscsx_mt_ridge = all_of("prscsx_mt"))
  untuned <- if (is.null(tpop)) list() else
    list(prscsx = col_of("prscsx", mtrait, tpop),
         mtag_prscsx = col_of("mtag_prscsx", mtrait, tpop),
         prscsx_mt = col_of("prscsx_mt", mtrait, tpop))

  scores <- list(); coefs <- list()
  for (m in names(untuned)) scores[[m]] <- dat[[untuned[[m]]]]
  for (m in names(stacks)) {
    X <- as.matrix(dat[, stacks[[m]], with = FALSE])
    mu <- colMeans(X[is_tune, , drop = FALSE]); sdv <- apply(X[is_tune, , drop = FALSE], 2, sd)
    sdv[sdv == 0] <- 1
    Xz <- sweep(sweep(X, 2, mu), 2, sdv, "/")
    fr <- fit_ridge(Xz[is_tune, , drop = FALSE], Cz[is_tune, , drop = FALSE], y[is_tune],
                    family, foldid)
    scores[[m]] <- drop(Xz %*% fr$beta)
    coefs[[m]] <- data.table(method = m, candidate = names(fr$beta), coef = fr$beta,
                             lambda = fr$lambda)
  }

  yv <- y[val]; Cv <- Cz[val, , drop = FALSE]
  sv <- lapply(scores, function(s) zs(s[val]))
  n_val <- length(val)
  n_cases_val <- if (type == "binary") sum(yv) else NA
  metric_rows <- list(); rocs <- list()
  for (m in names(sv)) {
    if (type == "binary") {
      rocs[[m]] <- roc(yv, sv[[m]], quiet = TRUE, direction = "<", levels = c(0, 1))
      ci <- as.numeric(ci.auc(rocs[[m]], method = "delong"))
      g <- glm(yv ~ sv[[m]] + Cv, family = binomial)
      b <- coef(summary(g))[2, ]
      metric_rows[[m]] <- data.table(method = m, metric = "AUC", estimate = ci[2],
                                     ci_lo = ci[1], ci_hi = ci[3],
                                     or_per_sd = exp(b[1]),
                                     or_lo = exp(b[1] - 1.96 * b[2]), or_hi = exp(b[1] + 1.96 * b[2]))
    } else {
      metric_rows[[m]] <- data.table(method = m, metric = "incR2",
                                     estimate = inc_r2(yv, Cv, sv[[m]]),
                                     ci_lo = NA_real_, ci_hi = NA_real_,
                                     or_per_sd = NA_real_, or_lo = NA_real_, or_hi = NA_real_)
    }
  }
  metrics <- rbindlist(metric_rows)

  # paired contrasts
  pairs <- list(c("prscsx_mt_ridge", "prsxtra"), c("prscsx_mt", "prsxtra"),
                c("prsxtra", "prsxa"), c("prscsx_mt", "prscsx"),
                c("prscsx_mt_ridge", "prscsx_mt"), c("mtag_prscsx", "prscsx"))
  pairs <- Filter(function(p) all(p %in% names(sv)), pairs)
  boot_idx <- if (type != "binary") replicate(n_boot, sample.int(n_val, replace = TRUE),
                                              simplify = FALSE) else NULL
  if (type != "binary") {
    # bootstrap incremental R2 once per method on shared resamples
    boot_m <- sapply(names(sv), function(m) vapply(boot_idx, function(ix)
      inc_r2(yv[ix], Cv[ix, , drop = FALSE], sv[[m]][ix]), numeric(1)))
    for (m in names(sv)) {
      q <- quantile(boot_m[, m], c(0.025, 0.975))
      metrics[method == m, `:=`(ci_lo = q[1], ci_hi = q[2])]
    }
  }
  con_rows <- lapply(pairs, function(p) {
    if (type == "binary") {
      t <- roc.test(rocs[[p[1]]], rocs[[p[2]]], method = "delong", paired = TRUE)
      d <- as.numeric(auc(rocs[[p[1]]]) - auc(rocs[[p[2]]]))
      ci <- if (!is.null(t$conf.int)) as.numeric(t$conf.int) else c(NA, NA)
      data.table(contrast = paste(p, collapse = " - "), metric = "AUC", delta = d,
                 ci_lo = ci[1], ci_hi = ci[2], p = t$p.value, test = "DeLong (paired)")
    } else {
      d <- metrics[method == p[1], estimate] - metrics[method == p[2], estimate]
      bd <- boot_m[, p[1]] - boot_m[, p[2]]
      q <- quantile(bd, c(0.025, 0.975))
      pv <- min(1, 2 * min(mean(bd <= 0), mean(bd >= 0)))
      data.table(contrast = paste(p, collapse = " - "), metric = "incR2", delta = d,
                 ci_lo = q[1], ci_hi = q[2], p = max(pv, 1 / n_boot), test = "paired bootstrap")
    }
  })
  contrasts <- rbindlist(con_rows)

  cell <- data.table(ancestry = anc, outcome = outcome, split = split_id,
                     score_pop = if (is.null(tpop)) NA_character_ else tpop,
                     n_tune = sum(is_tune), n_val = n_val,
                     n_cases_tune = if (type == "binary") sum(y[is_tune]) else NA,
                     n_cases_val = n_cases_val)
  tag <- function(x) cbind(cell[, .(ancestry, outcome, split)], x)
  list(metrics = tag(metrics), contrasts = tag(contrasts),
       coefs = tag(rbindlist(coefs)), cell = cell)
}

# ── loop ─────────────────────────────────────────────────────────────────────
res <- list(); skipped <- list()
for (outcome in names(outcomes)) {
  spec <- outcomes[[outcome]]
  if (!outcome %in% names(D)) { warning("outcome column missing: ", outcome); next }
  for (anc in sort(unique(D$ancestry))) {
    dat <- D[ancestry == anc & !is.na(get(outcome))]
    dat <- dat[complete.cases(dat[, ..covariates])]
    n <- nrow(dat); nc <- if (spec$type == "binary") sum(dat[[outcome]]) else NA
    if (n < min_n || (spec$type == "binary" && (nc < min_cases || n - nc < min_cases))) {
      skipped[[length(skipped) + 1]] <- data.table(ancestry = anc, outcome = outcome, n = n,
                                                   n_cases = nc, reason = "below minimum size")
      next
    }
    for (s in seq_len(n_splits)) {
      cat(sprintf("%s / %s / split %d (n = %d)\n", outcome, anc, s, n))
      res[[length(res) + 1]] <- run_cell(dat, outcome, spec, anc, s)
    }
  }
}
if (!length(res)) stop("no ancestry x outcome cell met the minimum size")

metrics   <- rbindlist(lapply(res, `[[`, "metrics"), fill = TRUE)
contrasts <- rbindlist(lapply(res, `[[`, "contrasts"), fill = TRUE)
coefs     <- rbindlist(lapply(res, `[[`, "coefs"), fill = TRUE)
cells     <- rbindlist(lapply(res, `[[`, "cell"), fill = TRUE)
cells[, small_count := n_tune < SUPPRESS | n_val < SUPPRESS |
        (!is.na(n_cases_val) & (n_cases_val < SUPPRESS | n_val - n_cases_val < SUPPRESS |
                                n_cases_tune < SUPPRESS))]
for (v in c("n_tune", "n_val", "n_cases_tune", "n_cases_val"))
  cells[small_count == TRUE, (v) := NA]
if (length(skipped)) {
  sk <- rbindlist(skipped)
  sk[, `:=`(n = ifelse(n < SUPPRESS, NA, n), n_cases = ifelse(!is.na(n_cases) & n_cases < SUPPRESS, NA, n_cases))]
  fwrite(sk, file.path(args$out_dir, "skipped_cells.tsv"), sep = "\t")
}

fwrite(metrics,   file.path(args$out_dir, "metrics.tsv"), sep = "\t")
fwrite(contrasts, file.path(args$out_dir, "contrasts.tsv"), sep = "\t")
fwrite(coefs,     file.path(args$out_dir, "ridge_coefs.tsv"), sep = "\t")
fwrite(cells,     file.path(args$out_dir, "cells.tsv"), sep = "\t")

# ── short markdown report ────────────────────────────────────────────────────
fmt <- function(x, d = 4) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
lines <- c(sprintf("# PRS-CSx-MT vs PRSxtra in All of Us — %s", man$analysis_name), "",
           sprintf("Splits: %d (tuning fraction %.2f); bootstrap resamples: %d. Metric: incremental R² over covariates (continuous) or AUC of the score (binary).",
                   n_splits, tune_frac, n_boot), "")
agg_m <- metrics[, .(estimate = mean(estimate), ci_lo = mean(ci_lo), ci_hi = mean(ci_hi)),
                 by = .(outcome, ancestry, method, metric)]
agg_c <- contrasts[, .(delta = mean(delta), ci_lo = mean(ci_lo), ci_hi = mean(ci_hi), p = median(p)),
                   by = .(outcome, ancestry, contrast, metric)]
for (o in unique(agg_m$outcome)) {
  lines <- c(lines, sprintf("## %s", o), "", "| ancestry | method | metric | estimate | 95% CI |", "|---|---|---|---|---|")
  for (r in seq_len(nrow(agg_m[outcome == o]))) {
    x <- agg_m[outcome == o][r]
    lines <- c(lines, sprintf("| %s | %s | %s | %s | [%s, %s] |", x$ancestry, x$method, x$metric,
                              fmt(x$estimate), fmt(x$ci_lo), fmt(x$ci_hi)))
  }
  lines <- c(lines, "", "| ancestry | contrast | Δ | 95% CI | p |", "|---|---|---|---|---|")
  for (r in seq_len(nrow(agg_c[outcome == o]))) {
    x <- agg_c[outcome == o][r]
    lines <- c(lines, sprintf("| %s | %s | %s | [%s, %s] | %s |", x$ancestry, x$contrast,
                              fmt(x$delta), fmt(x$ci_lo), fmt(x$ci_hi), formatC(x$p, format = "g", digits = 3)))
  }
  lines <- c(lines, "")
}
if (n_splits > 1) lines <- c(lines, "With several splits, estimates and CI bounds are averaged over splits and p is the median; splits are not independent replicates.", "")
writeLines(lines, file.path(args$out_dir, "report.md"))
cat("Wrote metrics.tsv, contrasts.tsv, ridge_coefs.tsv, cells.tsv, report.md ->", args$out_dir, "\n")
