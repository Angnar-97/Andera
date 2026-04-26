# Validaciones de entrada para análisis específicos.
#
# Cada validador devuelve una lista con:
#   - $ok      logical, TRUE si el phyloseq es apto para el análisis
#   - $reason  character (presente si !ok), mensaje legible para el usuario
#
# Patrón de uso en módulos:
#   v <- validate_for_grouping(physeq(), input$grouping, min_n_per_group = 3)
#   if (!v$ok) {
#     shinyalert::shinyalert(title = "Configuración inválida",
#                            text = v$reason, type = "warning")
#     return()
#   }

# ---------------------------------------------------------------------------
# Comprobaciones atómicas
# ---------------------------------------------------------------------------

# TRUE si la otu_table contiene únicamente conteos enteros no negativos.
otu_is_integer_counts <- function(ps) {
  otu <- tryCatch(
    as(phyloseq::otu_table(ps), "matrix"),
    error = function(e) NULL
  )
  if (is.null(otu) || !length(otu)) return(FALSE)
  if (any(is.na(otu))) return(FALSE)
  if (any(otu < 0))    return(FALSE)
  # Tolerancia mínima por errores numéricos
  all(abs(otu - round(otu)) < .Machine$double.eps^0.5)
}

# ---------------------------------------------------------------------------
# Validadores por análisis
# ---------------------------------------------------------------------------

# Índices Chao1 / ACE / Observed: requieren conteos enteros (singletons,
# doubletons). Con proporciones devuelven valores sin sentido biológico.
validate_for_count_indices <- function(ps) {
  if (is.null(ps)) {
    return(list(ok = FALSE, reason = "No hay phyloseq cargado."))
  }
  if (!otu_is_integer_counts(ps)) {
    return(list(
      ok = FALSE,
      reason = paste(
        "Los índices Chao1, ACE y Observed requieren conteos enteros",
        "(singletons / doubletons). El phyloseq cargado contiene valores no",
        "enteros — probablemente proporciones o abundancias relativas. Usa",
        "Shannon, Simpson o InvSimpson, que sí funcionan con proporciones."
      )
    ))
  }
  list(ok = TRUE)
}

# Distancias UniFrac (no ponderada y ponderada): requieren árbol filogenético
# en el slot phy_tree. has_tree() ya hace el chequeo en utils_phyloseq.R.
validate_for_unifrac <- function(ps) {
  if (is.null(ps)) {
    return(list(ok = FALSE, reason = "No hay phyloseq cargado."))
  }
  if (!has_tree(ps)) {
    return(list(
      ok = FALSE,
      reason = paste(
        "UniFrac (ponderada o no) requiere un árbol filogenético en el",
        "slot phy_tree. El phyloseq cargado no lo tiene. Usa Bray–Curtis o",
        "Jaccard si no dispones de árbol."
      )
    ))
  }
  list(ok = TRUE)
}

# Variable de agrupación válida para tests multivariantes (PERMANOVA, ANOSIM,
# Mantel). Comprueba existencia, ≥2 niveles, y mínimo de muestras por grupo.
validate_for_grouping <- function(ps, var, min_n_per_group = 3) {
  if (is.null(ps)) {
    return(list(ok = FALSE, reason = "No hay phyloseq cargado."))
  }
  if (is.null(var) || !nzchar(var)) {
    return(list(ok = FALSE, reason = "Falta seleccionar una variable de agrupación."))
  }
  sd <- as(phyloseq::sample_data(ps), "data.frame")
  if (!var %in% colnames(sd)) {
    return(list(
      ok = FALSE,
      reason = sprintf("La variable '%s' no existe en sample_data.", var)
    ))
  }

  groups <- sd[[var]]
  groups <- groups[!is.na(groups)]
  if (!length(groups)) {
    return(list(
      ok = FALSE,
      reason = sprintf("La variable '%s' no tiene valores no-NA.", var)
    ))
  }
  tab <- table(groups)
  if (length(tab) < 2) {
    return(list(
      ok = FALSE,
      reason = sprintf(
        "La variable '%s' tiene un único nivel (%s); se necesitan al menos 2 grupos.",
        var, names(tab)[1]
      )
    ))
  }
  if (min(tab) < min_n_per_group) {
    smallest <- names(tab)[which.min(tab)]
    return(list(
      ok = FALSE,
      reason = sprintf(
        "El grupo '%s' tiene %d muestra(s); se necesita un mínimo de %d por grupo. Usa otra variable o filtra los grupos pequeños.",
        smallest, min(tab), min_n_per_group
      )
    ))
  }
  list(ok = TRUE, n_per_group = as.list(tab))
}

# Combinador: ejecuta varios validadores y devuelve el primero que falle.
# Uso: first_invalid(list(validate_for_grouping(ps, var), validate_for_unifrac(ps)))
first_invalid <- function(validations) {
  for (v in validations) if (!isTRUE(v$ok)) return(v)
  list(ok = TRUE)
}
