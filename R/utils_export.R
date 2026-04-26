# Generación de scripts R reproducibles desde el estado de los módulos.
#
# Cada función `r_code_<module>(meta)` toma los metadatos guardados por el
# módulo tras una ejecución y devuelve un vector de líneas de código R que
# reproduce el análisis. Las líneas se concatenan con `paste(..., collapse=
# "\n")` antes de mostrarse en pantalla y se vuelcan tal cual al archivo `.R`
# descargado.
#
# Convención: el código generado asume que existe un objeto `ps` (phyloseq)
# cargado por el usuario; no intentamos templar la carga porque depende del
# origen de los datos (archivo local vs dataset de ejemplo).

# ---- Header común ------------------------------------------------------

r_code_header <- function(title) {
  c(
    sprintf("# Andera · %s", title),
    sprintf("# Generado: %s", format(Sys.Date(), "%Y-%m-%d")),
    "# Reproduce el análisis ejecutado en la app sobre el objeto `ps`",
    "# (asume que `ps` ya está cargado: phyloseq con otu_table, tax_table,",
    "# sample_data y, opcionalmente, phy_tree).",
    ""
  )
}

# Pretty-print de la distancia: para Aitchison añadimos el wrapper CLR
# explícito; para gunifrac, la llamada GUniFrac.
r_code_distance <- function(method) {
  if (method == "aitchison") {
    c(
      "# Aitchison = CLR + Euclidean",
      "ps_clr <- microbiome::transform(ps, 'clr')",
      "otu_clr <- as(phyloseq::otu_table(ps_clr), 'matrix')",
      "if (phyloseq::taxa_are_rows(ps_clr)) otu_clr <- t(otu_clr)",
      "d <- stats::dist(otu_clr, method = 'euclidean')"
    )
  } else if (method == "gunifrac") {
    c(
      "# Generalized UniFrac (alpha = 0.5)",
      "otu <- as(phyloseq::otu_table(ps), 'matrix')",
      "if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)",
      "gu <- GUniFrac::GUniFrac(otu, phyloseq::phy_tree(ps), alpha = 0.5)",
      "d <- stats::as.dist(gu$unifracs[, , 'd_0.5'])"
    )
  } else {
    c(sprintf("d <- phyloseq::distance(ps, method = '%s')", method))
  }
}

# ---- mod_ordination ----------------------------------------------------

r_code_ordination <- function(meta) {
  if (is.null(meta)) return("")
  out <- c(
    r_code_header("Ordenación"),
    "library(phyloseq)",
    "library(vegan)",
    "library(ggplot2)",
    "",
    r_code_distance(meta$distance),
    ""
  )
  out <- c(out, switch(
    meta$method,
    "PCoA" = "ord <- phyloseq::ordinate(ps, method = 'PCoA', distance = d)",
    "NMDS" = "ord <- phyloseq::ordinate(ps, method = 'NMDS', distance = d)",
    "CAP"  = sprintf(
      "ord <- phyloseq::ordinate(ps, method = 'CAP', distance = d,\n                          formula = ~ %s)",
      meta$constraint_variable
    )
  ), "")

  color_arg <- if (is.null(meta$color)) "NULL" else sprintf("'%s'", meta$color)
  out <- c(out,
    sprintf("phyloseq::plot_ordination(ps, ord, type = 'samples', color = %s) +",
            color_arg),
    "  ggplot2::geom_point(size = 3, alpha = 0.85) +",
    "  ggplot2::stat_ellipse(type = 'norm', linetype = 'dashed', alpha = 0.4) +",
    "  ggplot2::theme_minimal()"
  )

  if (meta$method == "NMDS") {
    out <- c(out, "",
      "# Diagnóstico — stress (Kruskal): < 0.10 excelente · 0.10–0.20 aceptable · > 0.20 deficiente",
      "cat('NMDS stress:', ord$stress, '\\n')")
  } else if (meta$method == "CAP") {
    out <- c(out, "",
      "# Diagnóstico — fracción de varianza constrained y test de permutación",
      "cat('R^2 constrained:', ord$CCA$tot.chi / ord$tot.chi, '\\n')",
      "vegan::anova.cca(ord, permutations = 999)")
  }
  out
}

# ---- mod_permanova ------------------------------------------------------

r_code_permanova <- function(meta) {
  if (is.null(meta)) return("")
  out <- c(
    r_code_header("PERMANOVA + betadisper"),
    "library(phyloseq)",
    "library(vegan)",
    "",
    r_code_distance(meta$distance),
    "",
    "sdat <- as(phyloseq::sample_data(ps), 'data.frame')",
    "",
    "# PERMANOVA global (vegan::adonis2)",
    sprintf("res_adonis <- vegan::adonis2(d ~ %s, data = sdat,", meta$variable),
    sprintf("                              permutations = %d)", meta$permutations),
    "print(res_adonis)",
    "",
    "# Diagnóstico de homogeneidad de dispersiones (Anderson 2001).",
    "# Si betadisper también es significativo, no se puede atribuir el",
    "# efecto solo a centroides — puede ser dispersión heterogénea.",
    sprintf("groups <- factor(sdat[['%s']])", meta$variable),
    "bd     <- vegan::betadisper(d, groups)",
    sprintf("bd_test <- vegan::permutest(bd, permutations = %d)", meta$permutations),
    "print(bd_test)"
  )

  if (isTRUE(meta$has_pairwise)) {
    out <- c(out, "",
      "# Pairwise PERMANOVA post-hoc (BH-FDR)",
      "# Necesario cuando el factor tiene >2 niveles para identificar qué",
      "# pares concretos difieren.",
      "groups_lvls <- levels(droplevels(groups))",
      "pairs <- utils::combn(groups_lvls, 2, simplify = FALSE)",
      "pw <- do.call(rbind, lapply(pairs, function(p) {",
      sprintf("  keep <- sdat[['%s']] %%in%% p", meta$variable),
      "  d_p  <- as.dist(as.matrix(d)[keep, keep])",
      "  sd_p <- sdat[keep, , drop = FALSE]",
      sprintf("  sd_p[['%s']] <- factor(sd_p[['%s']], levels = p)", meta$variable, meta$variable),
      sprintf("  a <- vegan::adonis2(d_p ~ %s, data = sd_p,", meta$variable),
      sprintf("                       permutations = %d)", meta$permutations),
      "  data.frame(group_a = p[1], group_b = p[2],",
      "             R2 = a$R2[1], F = a$F[1],",
      "             p  = a$`Pr(>F)`[1])",
      "}))",
      "pw$p_adj_BH <- stats::p.adjust(pw$p, method = 'BH')",
      "pw <- pw[order(pw$p_adj_BH), ]",
      "print(pw)"
    )
  }
  out
}

# ---- mod_diffabund ------------------------------------------------------

r_code_diffabund <- function(meta) {
  if (is.null(meta)) return("")

  out <- c(
    r_code_header("Abundancia diferencial"),
    "library(phyloseq)",
    if ("DESeq2" %in% meta$methods) "library(DESeq2)" else NULL,
    if ("ALDEx2" %in% meta$methods) "library(ALDEx2)" else NULL,
    ""
  )

  out <- c(out,
    sprintf("# Comparación: %s = %s vs %s", meta$variable, meta$group_a, meta$group_b),
    "sdat <- as(phyloseq::sample_data(ps), 'data.frame')",
    sprintf("ps_2 <- phyloseq::prune_samples(sdat[['%s']] %%in%% c('%s', '%s'), ps)",
            meta$variable, meta$group_a, meta$group_b),
    ""
  )

  if (meta$glom_rank != "__none__") {
    out <- c(out,
      sprintf("ps_2 <- phyloseq::tax_glom(ps_2, taxrank = '%s', NArm = TRUE)",
              meta$glom_rank),
      ""
    )
  }

  if (meta$min_prevalence > 0) {
    out <- c(out,
      "# Filtrado por prevalencia",
      "otu  <- as(phyloseq::otu_table(ps_2), 'matrix')",
      "if (!phyloseq::taxa_are_rows(ps_2)) otu <- t(otu)",
      "prev <- rowMeans(otu > 0)",
      sprintf("ps_2 <- phyloseq::prune_taxa(prev >= %.2f, ps_2)", meta$min_prevalence),
      ""
    )
  }

  if ("DESeq2" %in% meta$methods) {
    out <- c(out,
      "# DESeq2 (NB GLM con shrinkage)",
      sprintf("sdat2 <- as(phyloseq::sample_data(ps_2), 'data.frame')"),
      sprintf("sdat2[['%s']] <- factor(sdat2[['%s']], levels = c('%s', '%s'))",
              meta$variable, meta$variable, meta$group_a, meta$group_b),
      "phyloseq::sample_data(ps_2) <- phyloseq::sample_data(sdat2)",
      sprintf("dds <- phyloseq::phyloseq_to_deseq2(ps_2, ~ %s)", meta$variable),
      "# Trick para sparse data (Anders 2010): geometric mean robusta a ceros",
      "gm <- function(x) exp(sum(log(x[x > 0])) / length(x))",
      "geo <- apply(DESeq2::counts(dds), 1, gm)",
      "geo[!is.finite(geo)] <- 0",
      "dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geo)",
      "dds <- DESeq2::DESeq(dds, fitType = 'local', quiet = TRUE)",
      sprintf("res_deseq2 <- DESeq2::results(dds, contrast = c('%s', '%s', '%s'))",
              meta$variable, meta$group_b, meta$group_a),
      "print(head(as.data.frame(res_deseq2)[order(res_deseq2$padj), ], 20))",
      ""
    )
  }

  if ("ALDEx2" %in% meta$methods) {
    out <- c(out,
      "# ALDEx2 (Monte Carlo Dirichlet + CLR)",
      "otu2 <- as(phyloseq::otu_table(ps_2), 'matrix')",
      "if (!phyloseq::taxa_are_rows(ps_2)) otu2 <- t(otu2)",
      "storage.mode(otu2) <- 'integer'",
      "conds <- as.character(",
      sprintf("  as(phyloseq::sample_data(ps_2), 'data.frame')[['%s']]", meta$variable),
      ")",
      "res_aldex2 <- ALDEx2::aldex(reads = otu2, conditions = conds,",
      "                              test  = 't', effect = TRUE,",
      "                              denom = 'all', verbose = FALSE)",
      "print(head(res_aldex2[order(res_aldex2$we.eBH), ], 20))",
      ""
    )
  }

  out <- c(out,
    sprintf("# Umbrales aplicados: p_adj <= %.3f%s",
            meta$padj_threshold,
            if (meta$log2fc_threshold > 0)
              sprintf(", |log2FC| >= %.1f (DESeq2)", meta$log2fc_threshold)
            else "")
  )
  out
}

# ---- mod_composition ----------------------------------------------------

r_code_composition <- function(meta) {
  if (is.null(meta)) return("")
  out <- c(
    r_code_header("Composición taxonómica"),
    "library(phyloseq)",
    "library(ggplot2)",
    "",
    sprintf("ps_glom <- phyloseq::tax_glom(ps, taxrank = '%s', NArm = TRUE)", meta$rank),
    "ps_rel  <- phyloseq::transform_sample_counts(ps_glom,",
    "             function(x) if (sum(x) == 0) x else x / sum(x))",
    "df <- phyloseq::psmelt(ps_rel)",
    "",
    sprintf("# Top %d taxones por abundancia media · resto = 'Other'", meta$top_n),
    sprintf("means <- aggregate(df$Abundance, by = list(t = df[['%s']]), mean)",
            meta$rank),
    sprintf("top_n <- head(means[order(-means$x), 't'], %d)", meta$top_n),
    sprintf("df$.fill <- ifelse(df[['%s']] %%in%% top_n,", meta$rank),
    sprintf("                    as.character(df[['%s']]), 'Other')", meta$rank),
    "df$.fill <- factor(df$.fill, levels = c(top_n, 'Other'))",
    ""
  )

  if (meta$aggregation == "group") {
    out <- c(out,
      sprintf("agg <- aggregate(df$Abundance, by = list(group = df[['%s']],",
              meta$group_variable),
      "                                     taxon = df$.fill), mean)",
      "ggplot2::ggplot(agg, ggplot2::aes(group, x, fill = taxon)) +",
      "  ggplot2::geom_col(position = 'fill') +",
      sprintf("  ggplot2::labs(x = '%s', y = 'Abundancia relativa media',", meta$group_variable),
      sprintf("                 fill = '%s') +", meta$rank),
      "  ggplot2::theme_minimal()"
    )
  } else {
    out <- c(out,
      "ggplot2::ggplot(df, ggplot2::aes(Sample, Abundance, fill = .fill)) +",
      "  ggplot2::geom_col(position = 'fill') +",
      sprintf("  ggplot2::labs(x = 'Muestra', y = 'Abundancia relativa', fill = '%s') +",
              meta$rank),
      "  ggplot2::theme_minimal() +",
      "  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))"
    )
  }
  out
}

# ---- UI helper ----------------------------------------------------------

# Devuelve el HTML para el panel "Reproducir en R" — usable dentro de
# cualquier card_body o tagList. Convención: el módulo registra los outputs
# `r_code` y `download_r` con `r_code_handlers()` (abajo).
r_code_panel <- function(ns, code_id = "r_code", dl_id = "download_r") {
  tags$div(class = "andera-result-section",
    tags$span(class = "andera-eyebrow",
              tagList(bsicons::bs_icon("code-slash"), " Reproducir en R")),
    tags$pre(class = "andera-r-code",
      verbatimTextOutput(ns(code_id))
    ),
    tags$div(class = "andera-actions",
      downloadButton(ns(dl_id), "Descargar .R",
                     class = "btn-outline-secondary")
    )
  )
}

# Conecta el reactive de código con renderText + downloadHandler.
# Llamar dentro de moduleServer pasando el `output` local del módulo.
r_code_handlers <- function(output, code_reactive, prefix,
                             code_id = "r_code", dl_id = "download_r") {
  output[[code_id]] <- shiny::renderText({
    txt <- code_reactive()
    if (is.null(txt) || length(txt) == 0L)
      "# Ejecuta el análisis para ver el código equivalente."
    else paste(txt, collapse = "\n")
  })

  output[[dl_id]] <- shiny::downloadHandler(
    filename = function() paste0("andera-", prefix, "-", Sys.Date(), ".R"),
    content  = function(file) {
      txt <- code_reactive()
      if (is.null(txt) || length(txt) == 0L) {
        writeLines("# Sin análisis ejecutado.", file)
      } else {
        writeLines(txt, file)
      }
    }
  )
}
