# Import Packages
library(phyloseq)

# Load Data
gs <- data(GlobalPatterns)

# Data help
?GlobalPatterns

# Run provided examples
example("GlobalPatterns", ask=FALSE)
# GlblPt> data(GlobalPatterns)
# GlblPt> plot_richness(GlobalPatterns, x="SampleType", measures=c("Observed", "Chao1", "Shannon"))

# Verify that there is phylogenetic tree 
is.null(phy_tree(GlobalPatterns))



plot_heatmap(gs, taxa.label =  'Genre')

diversity_data <- estimate_richness(gs, measures = "Shannon") |>
  rowid_to_column()

colnames(diversity_data) <- c('sample', 'measure')

diversity_data 

ggplot(diversity_data, aes(x = sample, y = measure)) +
  geom_point() +
  theme_minimal()

plot_richness(gs, color="SampleType", measures=c("Chao1", "Shannon"))
p + geom_point(size=5, alpha=0.7)








ps1 <- ps_filter(enterotype, eval(parse(text = "SeqTech")) %in% c("Sanger"))


prune_taxa(taxa_sums(physeq()) > input$min_count)

load_object <- readRDS('phyloseq-2023-07-15.rds')
