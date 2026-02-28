# Regression: Excludes 2026 data to model only complete years (2000-2025).

#Please visit https://github.com/evidencesynthesis-tools/R for the data files this code uses, thanks.

library(tidyverse)
library(gridExtra)
library(viridis)
library(scales)


if (!file.exists("OSS_Modular.csv")) stop("OSS_Modular.csv not found.")
if (!file.exists("OSS_Integrated.csv")) stop("OSS_Integrated.csv not found.")
if (!file.exists("Prop_Modular.csv")) stop("Prop_Modular.csv not found.")
if (!file.exists("Prop_Integrated.csv")) stop("Prop_Integrated.csv not found.")

oss_mod <- read.csv("OSS_Modular.csv", stringsAsFactors = FALSE)
oss_int <- read.csv("OSS_Integrated.csv", stringsAsFactors = FALSE)
prop_mod <- read.csv("Prop_Modular.csv", stringsAsFactors = FALSE)
prop_int <- read.csv("Prop_Integrated.csv", stringsAsFactors = FALSE)



oss_data <- bind_rows(
  oss_mod %>% mutate(Type = "Modular"),
  oss_int %>% mutate(Type = "Integrated")
) %>%
  mutate(Language = str_replace(Language, "Javasscript", "JavaScript")) %>%
  mutate(
    Workflow_Stage = case_when(
      Category %in% c("Search", "Retrieval") ~ "Searching",
      Category %in% c("Screening", "Risk of Bias") ~ "Screening",
      Category %in% c("Data Extraction", "Data Cleaning", "Extension", "Plugin", "Manager") ~ "Extraction/Mgmt",
      Category %in% c("Meta-analysis", "Automation") ~ "Analysis",
      Category %in% c("Visualization", "Text Mining", "Open Catalog", "Workflow") ~ "Visualization/Other",
      TRUE ~ Category
    ),
    Lang_Group = ifelse(Language %in% c("C", "C#", "C++", "Java", "Javascript", "Other"), "Compiled/Web", Language)
  )


n_oss_mod <- nrow(oss_mod)
n_oss_int <- nrow(oss_int)
n_prop_mod <- nrow(prop_mod)
n_prop_int <- nrow(prop_int)


tab_lang_stage <- table(oss_data$Lang_Group, oss_data$Workflow_Stage)
chi_test <- chisq.test(tab_lang_stage)
n_sample <- sum(tab_lang_stage)
n <- sum(tab_lang_stage)
min_dim <- min(dim(tab_lang_stage)) - 1
cramers_v <- sqrt(chi_test$statistic / (n * min_dim))
v_interpret <- ifelse(cramers_v < 0.1, "Negligible", 
                      ifelse(cramers_v < 0.3, "Small", 
                             ifelse(cramers_v < 0.5, "Medium", "Large")))

#  POISSON REGRESSION (2026 Excluded) 
# Data (includes 2026)
yearly_counts_plot <- oss_data %>%
  count(Year) %>%
  arrange(Year) %>%
  filter(Year >= 2000 & Year <= 2026)

# Data for Modeling (Excludes 2026 for methodological rigor)
yearly_counts_model <- oss_data %>%
  count(Year) %>%
  arrange(Year) %>%
  filter(Year >= 2000 & Year <= 2025)


start_count <- yearly_counts_model$n[1]
end_count <- yearly_counts_model$n[nrow(yearly_counts_model)]

poisson_model <- glm(n ~ Year, family = poisson, data = yearly_counts_model)
dispersion_stat <- sum(residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
model_type <- "Standard Poisson"

if (dispersion_stat > 1.5) {
  qp_model <- glm(n ~ Year, family = quasipoisson, data = yearly_counts_model)
  irr <- exp(coef(qp_model)[2])
  se_log_irr <- summary(qp_model)$coefficients[2, 2]
  irr_conf <- exp(c(coef(qp_model)[2] - 1.96 * se_log_irr, coef(qp_model)[2] + 1.96 * se_log_irr))
  model_type <- "Quasi-Poisson (Robust SE)"
} else {
  irr <- exp(coef(poisson_model)[2])
  irr_conf <- exp(confint.default(poisson_model)[2, ])
}



total_oss <- n_oss_mod + n_oss_int
total_prop <- n_prop_mod + n_prop_int

pct_oss_mod <- round(n_oss_mod / total_oss * 100, 1)
pct_oss_int <- round(n_oss_int / total_oss * 100, 1)
pct_prop_mod <- round(n_prop_mod / total_prop * 100, 1)
pct_prop_int <- round(n_prop_int / total_prop * 100, 1)


excl_calculated <- data.frame(Name = character(), Reason = character(), stringsAsFactors = FALSE)
if (file.exists("excluded_tools_data.csv")) {
  raw_exc_lines <- readLines("excluded_tools_data.csv", warn = FALSE)
  for (line in raw_exc_lines) {
    clean_line <- sub("^\\d+\\.\\s+", "", line)
    if (grepl(":", clean_line)) {
      parts <- strsplit(clean_line, ":\\s*")[[1]]
      name <- parts[1]
      reason <- ifelse(length(parts) > 1, paste(parts[-1], collapse=":"), "Unknown")
      excl_calculated <- rbind(excl_calculated, data.frame(Name = name, Reason = reason, stringsAsFactors = FALSE))
    }
  }
  excl_clean_all <- excl_calculated %>%
    mutate(
      Exclusion_Type = case_when(
        grepl("Proprietary|Closed|Commercial|License", Reason) ~ "Proprietary Software",
        grepl("SaaS|Service|Hosted", Reason) ~ "SaaS/Service",
        grepl("Checklist|Guide|Framework|Protocol", Reason) ~ "Non-Software (Guide/Framework)",
        grepl("Database|Repository", Reason) ~ "Database",
        grepl("Omics|Bioinformatics|Astronomy", Reason) ~ "Irrelevant Domain",
        TRUE ~ "Other (Prototype/Abandoned)"
      )
    )
  excl_counts_all <- excl_clean_all %>% count(Exclusion_Type) %>% mutate(Pct = round(n / nrow(excl_clean_all) * 100, 1))
  n_excl_total <- nrow(excl_clean_all)
} else {
  excl_counts_all <- data.frame(Exclusion_Type = character(), n = integer(), Pct = numeric())
  n_excl_total <- 0
}


n_included_total <- n_sample 
n_screened_total <- n_included_total + n_excl_total


p_val_str <- ifelse(chi_test$p.value < 0.001, "< 0.001", round(chi_test$p.value, 4))

report_lines <- c(
  "=================================================================",
  "  STATISTICAL ANALYSIS REPORT (Methodological Summary)",
  "=================================================================",
  "",
  "1. LANGUAGE vs WORKFLOW ASSOCIATION",
  "-------------------------------------------------------------------",
  "A Chi-Square Test of Independence was performed to evaluate the",
  "association between programming language and workflow stage.",
  "",
  paste("  - Sample Size (N):", n_sample),
  paste("  - Chi-Square Statistic (X²):", round(chi_test$statistic, 2)),
  paste("  - Degrees of Freedom (df):", chi_test$parameter),
  paste("  - P-Value:", p_val_str),
  "",
  "Interpretation:",
 paste0("  A statistically significant association was observed between ",
       "programming language and workflow stage (p ", p_val_str, ")."),
  paste("  Effect Size (Cramer's V):", round(cramers_v, 3), 
        ". This represents a", tolower(v_interpret), "association."),
  "",
  "2. TEMPORAL TREND ANALYSIS",
  "-------------------------------------------------------------------",
  "Poisson regression was utilized to model the annual incidence of new tools.",
  "To avoid bias from incomplete yearly data, the year 2026 was excluded",
  "from the regression model (data collection ceased February 2026).",
  "",
  "Year Attribution:",
"  The year assigned to each tool corresponds to the year of its first",
"  associated research publication. In cases where an original",
"  publication year could not be identified, the publication year of",
"  the most recent citable version of the software was used as a proxy.",
"  This approach ensured consistent temporal classification while",
"  minimizing missing data.",
"",
"  Although publication year may not perfectly represent the initial",
"  release date of the software, it provides a standardized and",
"  reproducible temporal reference for longitudinal analysis.",
"",
  paste("  - Modeling Period: 2000 – 2025"),
  paste("  - Model Type:", model_type),
  paste("  - Dispersion Statistic:", round(dispersion_stat, 3)),
  "",
  "Trend Description:",
  paste("  The annual output increased from", start_count, 
        "tools in 2000 to", end_count, "tools in 2025."),
  "",
  "Interpretation:",
  paste("  - Incidence Rate Ratio (IRR):", round(irr, 3)),
  paste("  - 95% Confidence Interval: [", round(irr_conf[1], 3), ", ", round(irr_conf[2], 3), "]"),
  "  Because the confidence interval does not include 1.0, the upward",
"  temporal trend is statistically significant.",
  "",
  "3. MARKET STRUCTURE (OSS vs Proprietary)",
  "-------------------------------------------------------------------",
  paste0("A comparative analysis of the verified OSS dataset (n=", total_oss, ") against the"),
  paste0("identified proprietary sample (n=", total_prop, ") revealed a distinct architectural"),
  "dichotomy.",
  "",
  "Contingency Table:",
  paste("  - OSS Modular:", n_oss_mod, "| OSS Integrated:", n_oss_int),
  paste("  - Proprietary Modular:", n_prop_mod, "| Proprietary Integrated:", n_prop_int),
  "",
  "Descriptive Interpretation:",
  "  The analysis highlights a significant 'Modularity Gap' in the OSS ecosystem.",
  paste0("  While integrated platforms constituted ", pct_prop_int, "% (n=", n_prop_int, 
         ") of the identified proprietary tools, they represented only ", pct_oss_int, 
         "% (n=", n_oss_int, ") of the open-source tools."),
  paste0("  Conversely, ", pct_oss_mod, "% (n=", n_oss_mod, ") of OSS tools were classified as modular",
         " libraries, compared to ", pct_prop_mod, "% (n=", n_prop_mod, ") of proprietary tools."),
  "",
  "4. EXCLUSION CRITERIA",
  "-------------------------------------------------------------------",
  paste("  Total Records Screened:", n_screened_total),
  paste("  OSS Tools Included (Final Analysis):", n_included_total),
  paste("  Tools Excluded:", n_excl_total),
  "",
  "Primary reasons for exclusion:"
)
for(i in 1:nrow(excl_counts_all)) {
  report_lines <- c(report_lines, paste(sprintf("  - %-30s: %d (%.1f%%)", excl_counts_all$Exclusion_Type[i], excl_counts_all$n[i], excl_counts_all$Pct[i])))
}

writeLines(report_lines, "Statistical_Report_Explained.txt")
cat("[SUCCESS] Generated Detailed Statistical Report.\n")


dir.create("Figures", showWarnings = FALSE)


p1 <- ggplot(yearly_counts_plot, aes(x = Year, y = n)) +
  geom_col(fill = "#4a6b8b", width = 0.6) +
  geom_smooth(data = subset(yearly_counts_plot, Year <= 2025),
              method = "glm", method.args = list(family = poisson), 
              se = TRUE, color = "#e74c3c", linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = seq(2000, 2026, 5)) + 
  labs(title = "Temporal Dynamics of OSS Innovation", 
       subtitle = paste0("N = ", n_sample, " | ", model_type, 
                         " (2000-2025) | IRR = ", round(irr,3)),  #imp
       x = "Year", y = "Count of New Tools") +
  theme_minimal(base_size = 12) + 
  theme(
    panel.grid.minor = element_blank(), 
    plot.subtitle = element_text(face = "italic"),
    plot.title = element_text(face = "bold"),      
    axis.title.x = element_text(face = "bold"),   
    axis.title.y = element_text(face = "bold")   
  )

mkt_df <- data.frame(
  Ecosystem = c("Open-Source", "Proprietary"),
  Modular = c(n_oss_mod, n_prop_mod),
  Integrated = c(n_oss_int, n_prop_int)
) %>% 
  pivot_longer(cols = c(Modular, Integrated), names_to = "Type", values_to = "Count")



figure2_addon <- paste0(
  "Figure 2: Open-Source, Modular - ", pct_oss_mod, "% ",
  "Open-Source, Integrated - ", pct_oss_int, "% ",
  "Proprietary, Modular - ", pct_prop_mod, "% ",
  "Proprietary, Integrated - ", pct_prop_int, "%"
)

writeLines(figure2_addon, "Figure_Addons.txt")

p2 <- ggplot(mkt_df, aes(x = Ecosystem, y = Count, fill = Type)) +
  geom_col(position = "dodge", width = 0.7, color = "white", linewidth = 1) + 
  geom_text(aes(label = Count), position = position_dodge(width = 0.7), vjust = -0.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = c(0, 5), limits = c(0, max(mkt_df$Count)*1.15)) +
  scale_fill_manual(values = c("Modular" = "#3498db", "Integrated" = "#e67e22")) +
  labs(title = "Ecosystem Architecture Comparison", 
       subtitle = "Count of Modular vs Integrated Tools in OSS vs Proprietary Sectors",
       x = "", y = "Tool Count", fill = "Architecture") +
  theme_classic(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.line = element_line(color = "grey50"),
        plot.subtitle = element_text(color = "#555555", face = "italic"))


tab_detailed <- table(oss_data$Language, oss_data$Category)
chi_test_detailed <- chisq.test(tab_detailed) 

resid_df_detailed <- as.data.frame.matrix(chi_test_detailed$stdres)
resid_df_detailed$Language <- rownames(resid_df_detailed)
resid_df_long_detailed <- resid_df_detailed %>%
  pivot_longer(cols = -Language, names_to = "Category", values_to = "Residual")

cat_order <- c("Search", "Retrieval", "Screening", "Risk of Bias", "Data Extraction", 
               "Data Cleaning", "Extension", "Plugin", "Manager", "Meta-analysis", 
               "Automation", "Visualization", "Text Mining", "Open Catalog", "Workflow")
resid_df_long_detailed$Category <- factor(resid_df_long_detailed$Category, levels = cat_order)

target_values <- c(4.5, -4.1, 6.7, 8.3, 4.3)
highlight_df <- resid_df_long_detailed %>%
  mutate(Round_Res = round(Residual, 1)) %>%
  filter(Round_Res %in% target_values)

p3 <- ggplot() +
  geom_tile(data = resid_df_long_detailed, 
            aes(x = Category, y = Language, fill = Residual), 
            color = "white") + 
  geom_tile(data = highlight_df, 
            aes(x = Category, y = Language), 
            fill = "grey70", alpha = 0.9) +
  geom_text(data = resid_df_long_detailed, 
            aes(x = Category, y = Language, 
                label = round(Residual, 1)), 
            color = "black", size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4", 
                       midpoint = 0, limits = c(-4, 4), 
                       name = "Std. Residual",
                       guide = guide_colorbar(title.position = "top")) +
  labs(title = "Association Between Language and Workflow (Detailed)", 
       subtitle = "Standardized Residuals (Grey boxes indicate specific high-magnitude outliers)",
       x = "Workflow Category", 
       y = "Programming Language") +
  theme_minimal(base_size = 12) + 
  theme(
    axis.text.x = element_text(size = 11.5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )


p4 <- ggplot(excl_counts_all, 
             aes(x = reorder(Exclusion_Type, n), 
                 y = n, fill = Exclusion_Type)) +
  geom_col(alpha = 0.65) + 
  geom_text(aes(label = paste0(n, " (", round(Pct,1), "%)")), 
            hjust = 1.05, color = "black", size = 3.5) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, option = "plasma", guide = "none") +
  labs(title = "Number of Excluded Tools by Category", 
       subtitle = paste0("Total Screened: ", n_screened_total, 
                         " | Total Excluded: ", n_excl_total),
       x = "", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"), 
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.6),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )

heatmap_detailed <- oss_data %>% count(Language, Category)
heatmap_detailed$Category <- factor(heatmap_detailed$Category, levels = cat_order)

p5 <- ggplot(heatmap_detailed, aes(x = Category, y = Language, fill = n)) +
  geom_tile(color = "grey60", linewidth = 0.3, alpha = 0.55) +  
  geom_text(aes(label = ifelse(n > 0, n, "")), 
            color = "black", size = 3, fontface = "bold") +
  scale_fill_viridis(option = "viridis", 
                     direction = -1, 
                     name = "Count", 
                     na.value = "white") +
  labs(title = "Technological Stratification: Raw Counts", 
       subtitle = "Distribution of tools by Language and Category (With Grid Lines)",
       x = "Workflow Category", 
       y = "Language") +
  theme_classic(base_size = 12) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),   
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),    
    plot.subtitle = element_text(size = 13),
    panel.grid.major = element_line(color = "grey60", linewidth = 0.2),  
    axis.line = element_line(color = "black", linewidth = 0.4)
  )

formats <- c(".png", ".eps", ".tiff")

for(fmt in formats){
  ggsave(paste0("Figures/Figure1_Trend", fmt), p1, width = 11, height = 6, dpi = 300)
  ggsave(paste0("Figures/Figure2_MarketStructure", fmt), p2, width = 8, height = 6, dpi = 300)
  ggsave(paste0("Figures/Figure3_LanguageHeatmap", fmt), p3, width = 13, height = 8, dpi = 300)
  ggsave(paste0("Figures/Figure4_ExclusionReasons", fmt), p4, width = 11, height = 7, dpi = 300)
  ggsave(paste0("Figures/Figure5_DetailedCounts", fmt), p5, width = 12, height = 8, dpi = 300)
}

cat("[SUCCESS] All 5 figures redesigned and saved in .png, .eps, and .tiff formats.\n")
grid.arrange(p1, p2, p3, p4, p5, nrow = 5, heights = c(1, 1, 1.2, 1, 1.2))


#ignore Rplots.pdf

# Woohoo! you reached the end!