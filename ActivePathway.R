
library(ActivePathways)



scores <- read.table(
  system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways'), 
  header = TRUE, sep = '\t', row.names = 'Gene')
scores <- as.matrix(scores)
scores


scores[is.na(scores)] <- 1


gmt_file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
result <- ActivePathways(scores, gmt_file)






library(readxl)
data <- read_excel("top_150_drivers_Pathway.xlsx", sheet = "Sheet1")

data <- data[,c("...1", "p_dndscv", "p_dndsloc", "p_oncodrivefml", "p_chasmplusl", "p_oncodriveclustl")]

data <- as.matrix(data)
data[is.na(data)] <- 1
rownames(data) <- data[,"...1"]
data <- data[,-1]

rownames_data <- rownames(data)
data <- apply(data, 2, as.numeric)
rownames(data) <- rownames_data 


fname_GMT <- read.GMT(filename = "c2.cp.reactome.v2024.1.Hs.symbols.gmt")


enriched_pathways <- ActivePathways(data, fname_GMT, significant = 0.05, cytoscape_file_tag = "enrichmentMap__")

enriched_pathways$term_id <- gsub("^REACTOME_", "", enriched_pathways$term_id)


export_as_CSV(enriched_pathways, "ReactomeDB_2024.1_ActivePathways.csv")



