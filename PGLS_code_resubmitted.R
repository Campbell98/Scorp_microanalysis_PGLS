
# Integrated PGLS analysis + residual diagnostics + Excel export 


suppressPackageStartupMessages({
  library(phytools)
  library(ape)
  library(nlme)
  library(geiger)
  library(ggplot2)
  library(devtools)
  library(reshape)
  library(calibrate)
  library(qpcR)
  library(dplyr)
  library(openxlsx)
  library(lmtest)
})

# -------------------------
# Input data and plot onto trimmed phylo tree
# -------------------------
tree <- read.tree("santibanez2022_mcmc_independent_timetree.tre")
plot(tree)

metal.ions_check <- read.csv("Test_1.csv", row.names = 1)
head(metal.ions_check)
metal.ions <- read.csv("Test_1_8.csv", row.names = 1)
head(metal.ions)

obj <- name.check(tree, metal.ions_check)
obj

setdiff(tree$tip.label, row.names(metal.ions))

trimmed.tree <- drop.tip(tree, obj$tree_not_data)
plot(trimmed.tree)

print(trimmed.tree$tip.label)
trimmed.tree$tip.label <- c("lvar" = "Lychas mucronatus","Htrill"="Hottentotta hottentotta","Bisra"="Buthus occitanus","Aaustralis"="Androctonus australis",
                            "LquiE"="Leiurus quinquestriatus","Parabuthus"="Parabuthus granulatus","Uvit" = "Uroplectes vittatus","Tarc"="Tityus trinitatis",
                            "Ccar"="Centruroides suffusus","Bburmeiste"="Bothriurus bonariensis","Brotheas"= "Brotheas sp","Hadr"= "Hadrurus arizonesis", "Pbae"="Paruroctonus boreus",
                            "Pspinig" ="Paravaejovis spinigerus", "Ddiablo"="Diplocentrus lindo","Pandinus"= "Pandinus imperator","Hpaucidens"= "Hadogenes phylladus","Opisasper"= "Opisthcanthus lepturus")
print(trimmed.tree$tip.label)
plot(trimmed.tree)

name.check(trimmed.tree, metal.ions)

metal.ions <- metal.ions[trimmed.tree$tip.label,]
head(metal.ions)
colnames(metal.ions)

# -------------------------
# Check data & ancestral reconstruction
# -------------------------

P_zinc <- metal.ions[,1]; names(P_zinc) <- row.names(metal.ions)
P_iron <- metal.ions[,2]; names(P_iron) <- row.names(metal.ions)
P_calcium <- metal.ions[,3]; names(P_calcium) <- row.names(metal.ions)
T_zinc <- metal.ions[,4]; names(T_zinc) <- row.names(metal.ions)
T_Mn <- metal.ions[,5]; names(T_Mn) <- row.names(metal.ions)
T_calcium <- metal.ions[,6]; names(T_calcium) <- row.names(metal.ions)
P_CAR <- metal.ions[,7]; names(P_CAR) <- row.names(metal.ions)

P_zinc.anc <- fastAnc(trimmed.tree, P_zinc, CI = TRUE)
P_iron.anc <- fastAnc(trimmed.tree, P_iron, CI = TRUE)
P_calcium.anc <- fastAnc(trimmed.tree, P_calcium, CI = TRUE)
T_zinc.anc <- fastAnc(trimmed.tree, T_zinc, CI = TRUE)
T_Mn.anc <- fastAnc(trimmed.tree, T_Mn, CI = TRUE)
T_calcium.anc <- fastAnc(trimmed.tree, T_calcium, CI = TRUE)
P_CAR.anc <- fastAnc(trimmed.tree, P_CAR, CI = TRUE)

parsettings <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), cex = 0.4)
Chela_Zn <- contMap(trimmed.tree, P_zinc, res = 30, fsize = c(2, 0.5))
Chela_Fe <- contMap(trimmed.tree, P_iron, res = 30, fsize = c(2, 0.5))
Chela_Ca <- contMap(trimmed.tree, P_calcium, res = 30, fsize = c(2, 0.5))
Chela_CAR <- contMap(trimmed.tree, P_CAR, res = 30, fsize = c(2, 0.5))
Telson_Zn <- contMap(trimmed.tree, T_zinc, res = 30, fsize = c(2, 0.5))
Telson_Mn <- contMap(trimmed.tree, T_Mn, res = 30, fsize = c(2, 0.5))
Telson_Ca <- contMap(trimmed.tree, T_calcium, res = 30, fsize = c(2, 0.5))
par(parsettings)

# Chela
pdf("Chela_Zn.pdf", width = 8, height = 10)
plot(setMap(Chela_Zn, invert = TRUE))
dev.off()

pdf("Chela_Fe.pdf", width = 8, height = 10)
plot(setMap(Chela_Fe, invert = TRUE))
dev.off()

pdf("Chela_Ca.pdf", width = 8, height = 10)
plot(setMap(Chela_Ca, invert = TRUE))
dev.off()

pdf("Chela_CAR.pdf", width = 8, height = 10)
plot(setMap(Chela_CAR, invert = TRUE))
dev.off()

# Telson
pdf("Telson_Zn.pdf", width = 8, height = 10)
plot(setMap(Telson_Zn, invert = TRUE))
dev.off()

pdf("Telson_Mn.pdf", width = 8, height = 10)
plot(setMap(Telson_Mn, invert = TRUE))
dev.off()

pdf("Telson_Ca.pdf", width = 8, height = 10)
plot(setMap(Telson_Ca, invert = TRUE))
dev.off()






# -------------------------
# run_gls_matrix and diagnostics
# -------------------------
run_gls_matrix <- function(response, predictor, data, phy, save_models = FALSE) {
  library(nlme)
  df <- data[, c(response, predictor)]
  df <- df[complete.cases(df), ]
  if (nrow(df) < 3) stop("Not enough complete cases to fit models.")
  safe_gls <- function(correlation_type) {
    tryCatch({
      switch(correlation_type,
             "Brownian" = gls(as.formula(paste(response, "~", predictor)), data = df,
                              correlation = corBrownian(value = 1, phy = phy), method = "ML"),
             "Blomberg" = gls(as.formula(paste(response, "~", predictor)), data = df,
                              correlation = corBlomberg(value = 1, phy = phy, fixed = TRUE), method = "ML"),
             "Martins" = gls(as.formula(paste(response, "~", predictor)), data = df,
                             correlation = corMartins(value = 1, phy = phy, fixed = FALSE), method = "ML"),
             "Pagel" = gls(as.formula(paste(response, "~", predictor)), data = df,
                           correlation = corPagel(value = 1, phy = phy, fixed = FALSE), method = "ML"))
    }, error = function(e) {
      message(paste("⚠️ Model failed for", correlation_type, ":", e$message))
      return(NULL)
    })
  }
  
  models <- list(
    "Brownian Motion"    = safe_gls("Brownian"),
    "Early Burst"        = safe_gls("Blomberg"),
    "Ornstein-Uhlenbeck" = safe_gls("Martins"),
    "Lambda"             = safe_gls("Pagel")
  )
  
  result_matrix <- matrix(NA, nrow = 4, ncol = 8, dimnames = list(
    names(models),
    c("log likelihood", "AIC", "Delta AIC", "AIC Weights", "Slope", "p-value", "Lambda", "Alpha")
  ))
  
  for (i in seq_along(models)) {
    m <- models[[i]]
    if (!is.null(m)) {
      sm <- summary(m)
      corStruct <- m$modelStruct$corStruct
      result_matrix[i, "log likelihood"] <- as.numeric(logLik(m))
      result_matrix[i, "AIC"] <- AIC(m)
      coefvals <- tryCatch(coef(m), error = function(e) NA)
      if (!is.null(coefvals) && length(coefvals) >= 2) result_matrix[i, "Slope"] <- coefvals[2]
      if (!is.null(sm$tTable) && nrow(sm$tTable) >= 2) result_matrix[i, "p-value"] <- sm$tTable[2, 4]
      if (!is.null(corStruct)) {
        if (inherits(corStruct, "corPagel")) result_matrix[i, "Lambda"] <- coef(corStruct, unconstrained = FALSE)
        if (inherits(corStruct, "corMartins")) result_matrix[i, "Alpha"] <- coef(corStruct, unconstrained = FALSE)
      }
    }
  }
  
  valid_aics <- !is.na(result_matrix[, "AIC"])
  if (any(valid_aics)) {
    aic_vals <- as.numeric(result_matrix[valid_aics, "AIC"])
    delta_aic <- aic_vals - min(aic_vals, na.rm = TRUE)
    weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
    result_matrix[valid_aics, "Delta AIC"] <- delta_aic
    result_matrix[valid_aics, "AIC Weights"] <- weights
  }
  
  result_matrix <- round(result_matrix, 4)
  sorted <- result_matrix
  if (any(valid_aics)) sorted <- result_matrix[order(result_matrix[, "AIC"]), , drop = FALSE]
  
  if (save_models) return(list(results = sorted, models = models, df = df)) else return(sorted)
}

check_gls_diagnostics <- function(mod, data = NULL, plot = TRUE, verbose = TRUE) {
  if (is.null(mod)) stop("mod is NULL")
  res <- tryCatch({
    if (inherits(mod, "gls") || inherits(mod, "lme")) residuals(mod, type = "normalized") else residuals(mod)
  }, error = function(e) residuals(mod))
  fit <- tryCatch(fitted(mod), error = function(e) fitted.values(mod))
  n <- length(res)
  if (verbose) {
    cat("Model class:", paste(class(mod), collapse = "/"), "\n")
    cat("Number of observations:", n, "\n")
  }
  res_mean <- mean(res, na.rm = TRUE); res_sd <- sd(res, na.rm = TRUE)
  if (verbose) cat(sprintf("Residual mean = %.4g, sd = %.4g\n", res_mean, res_sd))
  shapiro_res <- NULL
  if (n >= 3 && n <= 5000) {
    shapiro_res <- tryCatch(stats::shapiro.test(res), error = function(e) { warning(e$message); NULL })
    if (verbose && !is.null(shapiro_res)) print(shapiro_res)
  } else if (verbose) {
    cat("Skipping Shapiro-Wilk (n outside 3..5000)\n")
  }
  
  spearman_res <- tryCatch(cor.test(abs(res), fit, method = "spearman"), error = function(e) { warning(e$message); NULL })
  if (verbose && !is.null(spearman_res)) print(spearman_res)
  
  aux_lm <- tryCatch(lm(I(res^2) ~ fit), error = function(e) { warning(e$message); NULL })
  aux_pvalue <- NA
  bp_chisq <- NULL; bp_pval <- NULL; bp_df <- 1
  if (!is.null(aux_lm)) {
    saux <- summary(aux_lm)
    if (nrow(saux$coefficients) >= 2) aux_pvalue <- saux$coefficients[2, 4]
    if (verbose) {
      cat("Auxiliary regression (res^2 ~ fitted) summary:\n")
      print(saux)
    }
    r2 <- summary(aux_lm)$r.squared
    bp_chisq <- n * r2
    bp_pval <- pchisq(bp_chisq, df = bp_df, lower.tail = FALSE)
    if (verbose) cat(sprintf("Breusch-Pagan style statistic: chi2 = %.4g, df = %d, p = %.4g\n", bp_chisq, bp_df, bp_pval))
  } else {
    if (verbose) cat("Auxiliary regression failed; cannot compute BP-style test.\n")
  }
  
  lm_bptest <- NULL
  if (!is.null(data) && inherits(mod, "gls")) {
    fm <- tryCatch(formula(mod), error = function(e) NULL)
    if (!is.null(fm) && all(all.vars(fm) %in% colnames(data))) {
      lm_for_bp <- tryCatch(lm(fm, data = data), error = function(e) NULL)
      if (!is.null(lm_for_bp)) {
        lm_bptest <- tryCatch(lmtest::bptest(lm_for_bp), error = function(e) { warning(e$message); NULL })
        if (verbose && !is.null(lm_bptest)) {
          cat("Breusch-Pagan test on an equivalent lm (ignoring phylogeny):\n")
          print(lm_bptest)
        }
      }
    } else {
      if (verbose) cat("Cannot run bptest: formula variables not all present in 'data'\n")
    }
  }
  
  if (plot) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(2, 2))
    hist(res, main = "Histogram of residuals", xlab = "Residuals", col = "lightgray")
    qqnorm(res, main = "QQ-plot of residuals"); qqline(res, col = "red", lwd = 2)
    plot(fit, res, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
    abline(h = 0, col = "red", lwd = 2); lines(lowess(fit, res), col = "blue", lwd = 2)
    sqrt_abs_res <- sqrt(abs(res))
    plot(fit, sqrt_abs_res, xlab = "Fitted values", ylab = expression(sqrt("|residuals|")), main = "Scale-Location")
    lines(lowess(fit, sqrt_abs_res), col = "blue", lwd = 2)
  }
  
  out <- list(
    residuals = res,
    fitted = fit,
    n = n,
    residual_mean = res_mean,
    residual_sd = res_sd,
    shapiro = shapiro_res,
    spearman_absres_fitted = spearman_res,
    auxiliary_lm = aux_lm,
    auxiliary_pvalue = aux_pvalue,
    bp_style = if (!is.null(aux_lm)) list(chi2 = bp_chisq, df = bp_df, p.value = bp_pval) else NULL,
    bptest_on_lm = lm_bptest
  )
  invisible(out)
}

# -------------------------
# Comparisons to run
# -------------------------
comparisons <- list(
  list(name = "ChelaTelsonZn", response = "Chela_R_Zn", predictor = "Telson_R_Zn"),
  list(name = "TelsonZnMn", response = "Telson_R_Zn", predictor = "Telson_R_Mn"),
  list(name = "TelsonChelaCa", response = "Telson_R_Ca", predictor = "Chela_R_Ca"),
  list(name = "TelsonZnCa", response = "Telson_R_Zn", predictor = "Chela_R_Ca"),
  list(name = "TelsonMnCa", response = "Telson_R_Mn", predictor = "Chela_R_Ca"),
  list(name = "ChelaZnFe", response = "Chela_R_Zn", predictor = "Chela_R_Fe"),
  list(name = "ChelaFeCa", response = "Chela_R_Fe", predictor = "Chela_R_Ca"),
  list(name = "ChelaZnCa", response = "Chela_R_Zn", predictor = "Chela_R_Ca"),
  list(name = "ChelaZnCAR", response = "Chela_R_Zn", predictor = "CAR_ratio"),
  list(name = "ChelaFeCAR", response = "Chela_R_Fe", predictor = "CAR_ratio"),
  list(name = "ChelaCaCAR", response = "Chela_R_Ca", predictor = "CAR_ratio")
)

# Containers
all_results_list <- list()
all_models_list <- list()
all_diagnostics <- list()

# Workbooks
wb_aic <- createWorkbook()
wb_diag <- createWorkbook()  # NEW: workbook for diagnostics summary

for (cmp in comparisons) {
  nm <- cmp$name; resp <- cmp$response; pred <- cmp$predictor
  message("Running comparison: ", nm, " (", resp, " ~ ", pred, ")")
  out <- tryCatch(run_gls_matrix(resp, pred, metal.ions, trimmed.tree, save_models = TRUE),
                  error = function(e) { message("Comparison failed: ", e$message); return(NULL) })
  if (is.null(out)) next
  
  # Save AIC table
  addWorksheet(wb_aic, nm)
  writeData(wb_aic, nm, as.data.frame(out$results), rowNames = TRUE)
  all_results_list[[nm]] <- out$results
  all_models_list[[nm]] <- out$models
  
  # Run diagnostics for each model and build a summary table
  diag_cmp <- list()
  diag_rows <- list()
  for (model_name in names(out$models)) {
    m <- out$models[[model_name]]
    if (is.null(m)) {
      diag_rows[[model_name]] <- data.frame(
        Model = model_name,
        n = NA_integer_,
        ResidualMean = NA_real_,
        ResidualSD = NA_real_,
        Shapiro_W = NA_real_,
        Shapiro_p = NA_real_,
        Spearman_rho = NA_real_,
        Spearman_p = NA_real_,
        Aux_pvalue = NA_real_,
        BP_chi2 = NA_real_,
        BP_df = NA_real_,
        BP_pvalue = NA_real_,
        bptest_statistic = NA_real_,
        bptest_pvalue = NA_real_,
        stringsAsFactors = FALSE
      )
      diag_cmp[[model_name]] <- NULL
      next
    }
    
    message("  Diagnostics for model: ", model_name)
    diag_res <- tryCatch(check_gls_diagnostics(m, data = out$df, plot = FALSE, verbose = FALSE),
                         error = function(e) { message("    Diagnostics failed: ", e$message); NULL })
    
    # Extract summary fields robustly
    n <- if (!is.null(diag_res)) diag_res$n else NA_integer_
    resid_mean <- if (!is.null(diag_res)) diag_res$residual_mean else NA_real_
    resid_sd <- if (!is.null(diag_res)) diag_res$residual_sd else NA_real_
    
    # Shapiro
    sh_W <- NA_real_; sh_p <- NA_real_
    if (!is.null(diag_res) && !is.null(diag_res$shapiro)) {
      sh_W <- as.numeric(diag_res$shapiro$statistic)
      sh_p <- as.numeric(diag_res$shapiro$p.value)
    }
    
    # Spearman
    sp_rho <- NA_real_; sp_p <- NA_real_
    if (!is.null(diag_res) && !is.null(diag_res$spearman_absres_fitted)) {
      est <- diag_res$spearman_absres_fitted$estimate
      # estimate may be named; take numeric
      sp_rho <- as.numeric(ifelse(length(est) > 0, est[[1]], NA_real_))
      sp_p <- as.numeric(diag_res$spearman_absres_fitted$p.value)
    }
    
    # Auxiliary p-value and BP-style
    aux_p <- NA_real_; bp_chi2 <- NA_real_; bp_df <- NA_integer_; bp_p <- NA_real_
    if (!is.null(diag_res)) {
      aux_p <- if (!is.null(diag_res$auxiliary_pvalue)) diag_res$auxiliary_pvalue else NA_real_
      if (!is.null(diag_res$bp_style)) {
        bp_chi2 <- as.numeric(diag_res$bp_style$chi2)
        bp_df <- as.integer(diag_res$bp_style$df)
        bp_p <- as.numeric(diag_res$bp_style$p.value)
      }
    }
    
    # bptest on equivalent lm (ignoring phylogeny)
    bp_stat <- NA_real_; bp_p_lm <- NA_real_
    if (!is.null(diag_res) && !is.null(diag_res$bptest_on_lm)) {
      bp_stat <- as.numeric(diag_res$bptest_on_lm$statistic)
      bp_p_lm <- as.numeric(diag_res$bptest_on_lm$p.value)
    }
    
    diag_rows[[model_name]] <- data.frame(
      Model = model_name,
      n = n,
      ResidualMean = resid_mean,
      ResidualSD = resid_sd,
      Shapiro_W = sh_W,
      Shapiro_p = sh_p,
      Spearman_rho = sp_rho,
      Spearman_p = sp_p,
      Aux_pvalue = aux_p,
      BP_chi2 = bp_chi2,
      BP_df = bp_df,
      BP_pvalue = bp_p,
      bptest_statistic = bp_stat,
      bptest_pvalue = bp_p_lm,
      stringsAsFactors = FALSE
    )
    
    diag_cmp[[model_name]] <- diag_res
  } # end per-model loop
  
  # Combine per-model rows into a data.frame and write to diagnostics workbook
  diag_df <- do.call(rbind, diag_rows)
  # Ensure sheet name length is <= 31 (Excel sheet limit). Truncate if needed.
  sheet_name <- substr(nm, 1, 31)
  addWorksheet(wb_diag, sheet_name)
  writeData(wb_diag, sheet_name, diag_df)
  
  all_diagnostics[[nm]] <- diag_cmp
}

# Save AIC workbook
saveWorkbook(wb_aic, "GLS_Model_Results_AllComparisons_log.xlsx", overwrite = TRUE)
saveWorkbook(wb_diag, "GLS_Model_Diagnostics_Summary_log.xlsx", overwrite = TRUE)
saveRDS(list(aic_tables = all_results_list, models = all_models_list, diagnostics = all_diagnostics),
        file = "GLS_models_and_diagnostics_log.rds")


#note, ChatGPT was used to automate GLS-residuals analysis. 