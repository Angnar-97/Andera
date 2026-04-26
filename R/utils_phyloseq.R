# Utilidades para trabajar con objetos phyloseq y helpers de descarga.

has_tree <- function(ps) {
  tryCatch(
    !is.null(phyloseq::phy_tree(ps)),
    error = function(e) FALSE
  )
}

# UniFrac requiere un árbol enraizado. Si llega no enraizado, lo enraizamos
# con midpoint rooting (ape::midpoint.root vía phangorn) — silencioso para no
# alarmar al usuario, pero es necesario para que UniFrac dé valores correctos.
ensure_rooted_tree <- function(ps) {
  if (!has_tree(ps)) return(ps)
  tree <- phyloseq::phy_tree(ps)
  if (ape::is.rooted(tree)) return(ps)
  rooted <- tryCatch(
    phangorn::midpoint(tree),
    error = function(e) ape::root(tree, outgroup = 1L, resolve.root = TRUE)
  )
  phyloseq::phy_tree(ps) <- rooted
  ps
}

load_example_dataset <- function(name) {
  env <- new.env()
  switch(
    name,
    "Global Patterns" = {
      utils::data("GlobalPatterns", package = "phyloseq", envir = env)
      env$GlobalPatterns
    },
    "Enterotype" = {
      utils::data("enterotype", package = "phyloseq", envir = env)
      env$enterotype
    },
    "Soil" = {
      utils::data("soilrep", package = "phyloseq", envir = env)
      env$soilrep
    },
    NULL
  )
}

# ----- Helpers de descarga -----

# downloadHandler para un ggplot que vive en un reactive/función sin argumentos.
# Si el reactive devuelve NULL, la descarga se detiene silenciosamente.
download_plot <- function(plot_reactive, prefix,
                          width = 8, height = 6, dpi = 150) {
  shiny::downloadHandler(
    filename = function() paste0(prefix, "-", Sys.Date(), ".png"),
    content  = function(file) {
      p <- plot_reactive()
      shiny::req(p)
      ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = dpi)
    }
  )
}

# downloadHandler para un data.frame reactive → CSV.
download_csv <- function(df_reactive, prefix) {
  shiny::downloadHandler(
    filename = function() paste0(prefix, "-", Sys.Date(), ".csv"),
    content  = function(file) {
      d <- df_reactive()
      shiny::req(d)
      utils::write.csv(as.data.frame(d), file, row.names = TRUE)
    }
  )
}

# downloadHandler para cualquier objeto R → RDS.
download_rds <- function(obj_reactive, prefix) {
  shiny::downloadHandler(
    filename = function() paste0(prefix, "-", Sys.Date(), ".rds"),
    content  = function(file) {
      obj <- obj_reactive()
      shiny::req(obj)
      saveRDS(obj, file)
    }
  )
}

# ---------------------------------------------------------------------------
# Distancias soportadas
# ---------------------------------------------------------------------------
#
# `ANDERA_DISTANCES` es el catálogo único usado por los selectInputs de los
# módulos de PCoA, PERMANOVA y Grafos. La etiqueta es la mostrada al usuario;
# el valor es el identificador interno aceptado por `compute_distance()`.
#
# - bray, jaccard       — `phyloseq::distance` (vegdist)
# - unifrac, wunifrac   — `phyloseq::distance` (requieren phy_tree)
# - aitchison           — manual: CLR (`microbiome::transform`) + Euclidean.
#                          Composicionalmente correcta; no requiere árbol.

ANDERA_DISTANCES <- c(
  "Bray–Curtis"                          = "bray",
  "Jaccard"                              = "jaccard",
  "UniFrac (no ponderada)"               = "unifrac",
  "UniFrac ponderada"                    = "wunifrac",
  "Generalized UniFrac (α = 0.5)"   = "gunifrac",
  "Aitchison (CLR + Euclidean)"          = "aitchison"
)

distance_requires_tree <- function(method) {
  method %in% c("unifrac", "wunifrac", "gunifrac")
}

distance_is_compositional <- function(method) {
  method %in% c("aitchison")
}

# Cómputo unificado de distancia entre muestras. Devuelve un objeto `dist`
# aceptable directamente por `phyloseq::ordinate`, `vegan::adonis2` y
# `phyloseq::plot_net` (todos aceptan `dist` además de strings).
#
# - aitchison: CLR (`microbiome::transform`) + Euclidean. Composicionalmente
#   correcta; no requiere árbol. La imputación de ceros la hace `microbiome`
#   con la mitad del mínimo no-cero (suficiente para exploración; para
#   inferencia fina, usar `zCompositions::cmultRepl` previamente).
# - gunifrac: Generalized UniFrac (Chen 2012) con alpha = 0.5, balance entre
#   weighted (alpha=1) y unweighted (alpha=0). Mejor potencia general.
#   Implementado vía `GUniFrac::GUniFrac`.
# - bray, jaccard, unifrac, wunifrac: delegado a `phyloseq::distance`.
compute_distance <- function(ps, method) {
  if (method == "aitchison") {
    # CLR + Euclidean = Aitchison. Validar que partimos de conteos
    # (proporciones harían el CLR matemáticamente diferente al "real").
    if (!otu_is_integer_counts(ps)) {
      warning("Aitchison sobre proporciones es matemáticamente válido pero no equivalente al CLR sobre conteos. Considera usar conteos enteros.",
              call. = FALSE)
    }
    ps_clr <- microbiome::transform(ps, "clr")
    otu    <- as(phyloseq::otu_table(ps_clr), "matrix")
    if (phyloseq::taxa_are_rows(ps_clr)) otu <- t(otu)
    return(stats::dist(otu, method = "euclidean"))
  }

  if (method == "gunifrac") {
    if (!has_tree(ps)) stop("Generalized UniFrac requiere phy_tree.")
    ps_r <- ensure_rooted_tree(ps)   # GUniFrac falla silencioso si unrooted
    otu  <- as(phyloseq::otu_table(ps_r), "matrix")
    if (phyloseq::taxa_are_rows(ps_r)) otu <- t(otu)
    res <- GUniFrac::GUniFrac(otu, phyloseq::phy_tree(ps_r),
                                alpha = c(0.5))
    # `unifracs` es un array [n_samples, n_samples, length(alpha)]. Para
    # alpha=c(0.5) algunas versiones devuelven dimnames "d_0.5", otras solo
    # índices. Tomamos el primer (y único) slice de forma defensiva.
    d_mat <- res$unifracs[, , 1L]
    return(stats::as.dist(d_mat))
  }

  if (method %in% c("unifrac", "wunifrac")) {
    ps <- ensure_rooted_tree(ps)
  }

  phyloseq::distance(ps, method = method)
}
