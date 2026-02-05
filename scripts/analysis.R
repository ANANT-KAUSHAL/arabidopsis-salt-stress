# -----------------------------------------------------------
# PROJECT: Arabidopsis Salt Stress Transcriptomics
# AUTHOR: Anant Kaushal
# DATASET: GSE16765
# -----------------------------------------------------------

# === SAFETY BLOCK: FOLDER CREATION ===
# This fixes the "Unable to open file" error automatically.
# It creates a 'figures' folder right where you are running this script.
if (!dir.exists("figures")) {
  dir.create("figures")
  print("Created 'figures' folder.")
} else {
  print("'figures' folder already exists. Saving there.")
}

# 1. LOAD LIBRARIES
library(GEOquery)
library(limma)

# 2. LOAD DATA
print("Loading Salt Stress Data (GSE16765)...")
gse_salt <- getGEO("GSE16765", GSEMatrix = TRUE)
salt_data <- gse_salt[[1]]

# 3. SAMPLE SELECTION (Col-0 Genotype Only)
col_samples <- salt_data[, c("GSM420236", "GSM420237", "GSM420238", "GSM420239")]

# 4. EXPERIMENTAL DESIGN
groups <- factor(c("Control", "Control", "Salt", "Salt"))
design <- model.matrix(~0 + groups)
colnames(design) <- c("Control", "Salt")

# 5. DIFFERENTIAL EXPRESSION (Linear Modeling)
fit <- lmFit(col_samples, design)
contrast_matrix <- makeContrasts(Diff = Salt - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
all_results <- topTable(fit2, number = Inf)

# -----------------------------------------------------------
# OUTPUT 1: VOLCANO PLOT
# -----------------------------------------------------------
print("Generating Volcano Plot...")

all_results$Color <- "grey"
all_results$Color[all_results$adj.P.Val < 0.05 & abs(all_results$logFC) > 1] <- "red"

# FIXED: Save specifically to the folder we created
png("figures/volcano_plot_salt.png", width=800, height=600)

plot(all_results$logFC, -log10(all_results$adj.P.Val),
     col = all_results$Color,
     pch = 16, cex = 0.6,
     main = "Volcano Plot: Salt Stress in Col-0",
     xlab = "Log Fold Change",
     ylab = "-Log10 Adjusted P-value")

abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)
top_gene <- all_results["247925_at", ]
text(top_gene$logFC, -log10(top_gene$adj.P.Val), labels="TCH4", pos=2, font=2)

dev.off() # Closes the image file

# -----------------------------------------------------------
# OUTPUT 2: HEATMAP
# -----------------------------------------------------------
print("Generating Heatmap...")

wrky_rows <- grep("WRKY", all_results$Gene.Symbol)
syndicate_data <- col_samples[rownames(all_results)[wrky_rows], ]

# FIXED: Save to the figures folder
png("figures/wrky_heatmap.png", width=1000, height=600)

coolmap(syndicate_data, 
        main = "The WRKY Syndicate: Salt vs Control",
        labRow = "", 
        labCol = c("Control 1", "Control 2", "Salt 1", "Salt 2"))

dev.off() # Closes the image file

# -----------------------------------------------------------
# OUTPUT 3: CSV DATA
# -----------------------------------------------------------
write.csv(all_results, "figures/Salt_Stress_Full_Results.csv")

print("SUCCESS! Go check the 'figures' folder. The images are there.")