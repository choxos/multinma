#!/usr/bin/env Rscript
# Export R datasets and reference values for Julia comparison tests

library(multinma)

basedir <- "/Users/choxos/Documents/GitHub/multinma/multinma.jl"
outdir <- file.path(basedir, "data")
if (!dir.exists(outdir)) dir.create(outdir)

# Export datasets to CSV
datasets <- c("smoking", "blocker", "thrombolytics", "parkinsons", "diabetes",
              "statins", "transfusion", "atrial_fibrillation", "dietary_fat")

for (ds in datasets) {
  d <- get(ds)
  write.csv(d, file.path(outdir, paste0(ds, ".csv")), row.names = FALSE)
  cat("Exported:", ds, "- dims:", nrow(d), "x", ncol(d), "\n")
}

# =============================================================================
# Extract reference values for each vignette example
# =============================================================================

refdir <- file.path(basedir, "test", "reference_values")
if (!dir.exists(refdir)) dir.create(refdir)

# --- Smoking ---
smknet <- set_agd_arm(smoking,
                      study = studyn,
                      trt = trtc,
                      r = r,
                      n = n,
                      trt_ref = "No intervention")

cat("\n=== Smoking ===\n")
cat("Treatments:", paste(levels(smknet$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(smknet$treatments), "\n")
cat("N studies:", nlevels(smknet$studies), "\n")
cat("Reference:", levels(smknet$treatments)[1], "\n")
cat("Has AgD arm:", !is.null(smknet$agd_arm), "\n")
cat("Has AgD contrast:", !is.null(smknet$agd_contrast), "\n")
cat("Has IPD:", !is.null(smknet$ipd), "\n")

# Node-splits
smk_ns <- get_nodesplits(smknet)
cat("N nodesplits:", length(smk_ns), "\n")

sink(file.path(refdir, "smoking.txt"))
cat("treatments:", paste(levels(smknet$treatments), collapse = "|"), "\n")
cat("n_treatments:", nlevels(smknet$treatments), "\n")
cat("n_studies:", nlevels(smknet$studies), "\n")
cat("reference:", levels(smknet$treatments)[1], "\n")
cat("n_nodesplits:", length(smk_ns), "\n")
cat("nodesplits:\n")
for (ns in smk_ns) {
  cat("  ", ns[1], "vs", ns[2], "\n")
}
cat("likelihood:", "binomial_1par", "\n")
cat("link:", "logit", "\n")
sink()

# --- Blocker ---
blocker_net <- set_agd_arm(blocker,
                           study = studyn,
                           trt = trtc,
                           r = r,
                           n = n,
                           trt_ref = "Control")

cat("\n=== Blocker ===\n")
cat("Treatments:", paste(levels(blocker_net$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(blocker_net$treatments), "\n")
cat("N studies:", nlevels(blocker_net$studies), "\n")
cat("Reference:", levels(blocker_net$treatments)[1], "\n")

sink(file.path(refdir, "blocker.txt"))
cat("treatments:", paste(levels(blocker_net$treatments), collapse = "|"), "\n")
cat("n_treatments:", nlevels(blocker_net$treatments), "\n")
cat("n_studies:", nlevels(blocker_net$studies), "\n")
cat("reference:", levels(blocker_net$treatments)[1], "\n")
cat("likelihood:", "binomial_1par", "\n")
cat("link:", "logit", "\n")
sink()

# --- Thrombolytics ---
thrombo_net <- set_agd_arm(thrombolytics,
                           study = studyn,
                           trt = trtc,
                           r = r,
                           n = n)

cat("\n=== Thrombolytics ===\n")
cat("Treatments:", paste(levels(thrombo_net$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(thrombo_net$treatments), "\n")
cat("N studies:", nlevels(thrombo_net$studies), "\n")
cat("Reference:", levels(thrombo_net$treatments)[1], "\n")

thrombo_ns <- get_nodesplits(thrombo_net)
cat("N nodesplits:", length(thrombo_ns), "\n")

sink(file.path(refdir, "thrombolytics.txt"))
cat("treatments:", paste(levels(thrombo_net$treatments), collapse = "|"), "\n")
cat("n_treatments:", nlevels(thrombo_net$treatments), "\n")
cat("n_studies:", nlevels(thrombo_net$studies), "\n")
cat("reference:", levels(thrombo_net$treatments)[1], "\n")
cat("n_nodesplits:", length(thrombo_ns), "\n")
cat("nodesplits:\n")
for (ns in thrombo_ns) {
  cat("  ", ns[1], "vs", ns[2], "\n")
}
cat("likelihood:", "binomial_1par", "\n")
cat("link:", "logit", "\n")
sink()

# --- Parkinsons (arm-based) ---
arm_net <- set_agd_arm(parkinsons,
                       study = studyn,
                       trt = trtn,
                       y = y,
                       se = se,
                       sample_size = n)

cat("\n=== Parkinsons (arm) ===\n")
cat("Treatments:", paste(levels(arm_net$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(arm_net$treatments), "\n")
cat("N studies:", nlevels(arm_net$studies), "\n")
cat("Reference:", levels(arm_net$treatments)[1], "\n")

# Parkinsons (contrast-based)
contr_net <- set_agd_contrast(parkinsons,
                              study = studyn,
                              trt = trtn,
                              y = diff,
                              se = se_diff,
                              sample_size = n)

cat("\n=== Parkinsons (contrast) ===\n")
cat("Treatments:", paste(levels(contr_net$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(contr_net$treatments), "\n")
cat("N studies:", nlevels(contr_net$studies), "\n")

# Parkinsons (mixed)
parkinsons_arm <- parkinsons[parkinsons$studyn %in% 1:3, ]
parkinsons_contr <- parkinsons[parkinsons$studyn %in% 4:7, ]

mix_arm_net <- set_agd_arm(parkinsons_arm,
                           study = studyn,
                           trt = trtn,
                           y = y,
                           se = se,
                           sample_size = n)

mix_contr_net <- set_agd_contrast(parkinsons_contr,
                                  study = studyn,
                                  trt = trtn,
                                  y = diff,
                                  se = se_diff,
                                  sample_size = n)

mix_net <- combine_network(mix_arm_net, mix_contr_net)

cat("\n=== Parkinsons (mixed) ===\n")
cat("Treatments:", paste(levels(mix_net$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(mix_net$treatments), "\n")
cat("N studies:", nlevels(mix_net$studies), "\n")

sink(file.path(refdir, "parkinsons.txt"))
cat("arm_treatments:", paste(levels(arm_net$treatments), collapse = "|"), "\n")
cat("arm_n_treatments:", nlevels(arm_net$treatments), "\n")
cat("arm_n_studies:", nlevels(arm_net$studies), "\n")
cat("arm_reference:", levels(arm_net$treatments)[1], "\n")
cat("contr_treatments:", paste(levels(contr_net$treatments), collapse = "|"), "\n")
cat("contr_n_treatments:", nlevels(contr_net$treatments), "\n")
cat("contr_n_studies:", nlevels(contr_net$studies), "\n")
cat("mix_treatments:", paste(levels(mix_net$treatments), collapse = "|"), "\n")
cat("mix_n_treatments:", nlevels(mix_net$treatments), "\n")
cat("mix_n_studies:", nlevels(mix_net$studies), "\n")
cat("arm_likelihood:", "normal", "\n")
cat("arm_link:", "identity", "\n")
sink()

# --- Diabetes ---
db_net <- set_agd_arm(diabetes,
                      study = studyc,
                      trt = trtc,
                      r = r,
                      n = n)

cat("\n=== Diabetes ===\n")
cat("Treatments:", paste(levels(db_net$treatments), collapse = ", "), "\n")
cat("N treatments:", nlevels(db_net$treatments), "\n")
cat("N studies:", nlevels(db_net$studies), "\n")
cat("Reference:", levels(db_net$treatments)[1], "\n")

sink(file.path(refdir, "diabetes.txt"))
cat("treatments:", paste(levels(db_net$treatments), collapse = "|"), "\n")
cat("n_treatments:", nlevels(db_net$treatments), "\n")
cat("n_studies:", nlevels(db_net$studies), "\n")
cat("reference:", levels(db_net$treatments)[1], "\n")
cat("likelihood:", "binomial_1par", "\n")
cat("link:", "logit", "\n")
sink()

cat("\n=== Done exporting all reference values ===\n")
