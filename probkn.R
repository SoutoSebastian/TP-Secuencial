# ============================================================================
# PROGRAMA: ProbKn
# ----------------------------------------------------------------------------
# OBJETIVO
#   Calcular la probabilidad
#
#                      P(Kn > c*sigma)
#
#   donde Kn es la longitud aleatoria del intervalo de confianza usual
#   para la media mu de una poblacion Normal(mu, sigma^2), con mu y sigma
#   desconocidos, para un nivel de confianza 1-alpha.
#
# FUNDAMENTO TEORICO
#   Sean X1,...,Xn iid N(mu, sigma^2), con mu y sigma desconocidos. El
#   intervalo de confianza usual de nivel 1-alpha para mu es:
#
#       Jn = [ Xbar_n +/- a_(n-1) * Sn / sqrt(n) ]
#
#   donde a_(n-1) es el cuantil superior 50*alpha% de una distribucion t
#   de Student con n-1 grados de libertad; equivalentemente:
#
#       a_(n-1) = qt(1 - alpha/2, df = n-1)
#
#   La longitud del intervalo es:
#
#       Kn = 2 * a_(n-1) * Sn / sqrt(n)
#
#   Como:
#
#       n*(n-1)*Kn^2 / (4*a_(n-1)^2*sigma^2) ~ Chi-cuadrado(n-1),
#
#   se obtiene:
#
#       P(Kn > c*sigma)
#         = P[ Chi-cuadrado(n-1)
#              > n*(n-1)*c^2 / (4*a_(n-1)^2) ].
#
#   VALIDACION CON EL EJEMPLO DEL LIBRO
#   alpha:0.05, sample size n=20
# ============================================================================


# ----------------------------------------------------------------------------
# 1) PANTALLA DE BIENVENIDA DEL PROGRAMA
# ----------------------------------------------------------------------------

#' Imprime la pantalla de bienvenida del programa (imita la del original)
print_banner <- function() {
  banner <- c(
    "***********************************************",
    "* This program computes the probability of    *",
    "* K(n) > c*sigma for given alpha and sample   *",
    "* size.                                       *",
    "***********************************************",
    "",
    "Input following parameters:"
  )
  cat(paste(banner, collapse = "\n"), "\n")
}


# ----------------------------------------------------------------------------
# 2) NUCLEO ESTADISTICO: CALCULO DE P(Kn > c*sigma)
# ----------------------------------------------------------------------------

#' Calcula P(Kn > c*sigma)
#'
#' @param alpha nivel de significacion del intervalo (0 < alpha < 1)
#' @param n     tamaño de muestra (entero, n >= 2)
#' @param c     numero positivo o vector de numeros positivos
#'
#' @return un data.frame con c, el punto critico t, el umbral chi-cuadrado
#'         y la probabilidad p = P(Kn > c*sigma)

#ya se definieron alpha, n, y c
prob_kn <- function(alpha, n, c) {
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha debe ser un unico valor en (0,1)")
  }
  if (length(n) != 1 || !is.finite(n) || n < 2 || n != as.integer(n)) {
    stop("n debe ser un entero mayor o igual que 2")
  }
  if (length(c) < 1 || any(!is.finite(c)) || any(c <= 0)) {
    stop("c debe contener uno o mas valores positivos")
  } #se chequea que los valores dados sean válidos

  n <- as.integer(n)
  
  a_n_minus_1 <- qt(1 - alpha / 2, df = n - 1) #cálculo del cuantil a_n-1
  
  chi_threshold <- n * (n - 1) * c^2 / (4 * a_n_minus_1^2) #definimos el umbral
  
  p <- pchisq(chi_threshold, df = n - 1, lower.tail = FALSE) #calculamos la probabilidad 
                                                             #de la cola derecha 
                                                             #que equivale a P(Kn>c*sigma)

  data.frame(
    c = as.numeric(c),
    a_n_minus_1 = rep(a_n_minus_1, length(c)),
    chi_threshold = as.numeric(chi_threshold),
    p = as.numeric(p)
  )#dataframe final 
}


# ----------------------------------------------------------------------------
# 3) GENERACION DE LA TABLA DE PROBABILIDADES
# ----------------------------------------------------------------------------

#' Genera la tabla de ProbKn para un conjunto de valores de c
#' @param alpha    nivel de significacion
#' @param n        tamaño de muestra
#' @param c_values valores de c a evaluar
#'
#' @return data.frame devuelto por prob_kn()
prob_kn_table <- function(alpha, n, c_values = seq(0.44, 1.48, by = 0.04)) {
  c_values <- round(c_values, 10)
  prob_kn(alpha = alpha, n = n, c = c_values)
}


# ----------------------------------------------------------------------------
# 4) FORMATO DE LOS VALORES Y DE LA TABLA
# ----------------------------------------------------------------------------

format_book_decimal <- function(x, digits) {
  texto <- sprintf(paste0("%.", digits, "f"), x)
  texto <- sub("^0\\.", ".", texto)
  texto <- sub("^-0\\.", "-.", texto)
  texto
}

#' @param tabla data.frame devuelto por prob_kn() o prob_kn_table()
#'
#' @return vector character, una linea por fila de salida
format_probkn_table <- function(tabla) {
  if (!all(c("c", "p") %in% names(tabla))) {
    stop("tabla debe contener las columnas c y p")
  }

  total <- nrow(tabla)
  if (total < 1) return(character(0))

  n_blocks <- 3L
  n_rows <- ceiling(total / n_blocks)
  block_width <- 13L
  lines <- character(n_rows)

  for (row in seq_len(n_rows)) {
    parts <- character(n_blocks)

    for (block in seq_len(n_blocks)) {
      idx <- row + (block - 1L) * n_rows

      if (idx <= total) {
        c_text <- format_book_decimal(tabla$c[idx], 2)
        p_text <- format_book_decimal(tabla$p[idx], 3)
        parts[block] <- sprintf("%5s %5s", c_text, p_text)
      } else {
        parts[block] <- strrep(" ", 11L)
      }
    }

    lines[row] <- paste(parts, collapse = strrep(" ", block_width - 11L))
  }

  lines
}

# ----------------------------------------------------------------------------
# 5) IMPRESION DE LOS RESULTADOS EN CONSOLA Y ARCHIVO
# ----------------------------------------------------------------------------

#' @param tabla data.frame devuelto por prob_kn_table()
#' @param alpha nivel de significacion
#' @param n     tamaño de muestra
print_probkn_screen <- function(tabla, alpha, n) {
  out <- c(
    "Probability p = P(K(n)>c*sigma)",
    sprintf("alpha =%s and Sample Size = %d",
            format_book_decimal(alpha, 3), as.integer(n)),
    "---------------------------------",
    "    c     p      c     p      c     p",
    format_probkn_table(tabla)
  )

  cat(paste(out, collapse = "\n"), "\n")
  invisible(out)
}

#' @param tabla data.frame devuelto por prob_kn_table()
#' @param alpha nivel de significacion
#' @param n     tamaño de muestra
#' @param file  ruta de archivo de salida ("" imprime solo en consola)
print_probkn_output <- function(tabla, alpha, n, file = "") {
  out <- c(
    "OUTPUT FROM THE PROGRAM: ProbKn",
    "===============================",
    "Probability, p=P(K(n)>c*sigma)",
    sprintf("alpha =%s and Sample Size = %d",
            format_book_decimal(alpha, 3), as.integer(n)),
    "---------------------------------",
    "    c     p      c     p      c     p",
    format_probkn_table(tabla)
  )

  cat(paste(out, collapse = "\n"), "\n")

  if (nzchar(file)) {
    cat(paste(out, collapse = "\n"), "\n", file = file, append = TRUE)
  }

  invisible(out)
}

# ----------------------------------------------------------------------------
# 6) INGRESO DE PARAMETROS POR CONSOLA
# ----------------------------------------------------------------------------

#' Pide alpha y n por consola, como el programa original
get_probkn_params <- function() {
  repeat {
    alpha <- as.numeric(readline("alpha = "))
    n <- as.numeric(readline("sample size, n = "))

    ans <- readline("Do you want to change the input parameters (Y or N)? ")
    if (toupper(trimws(ans)) != "Y") break
  }

  list(alpha = alpha, n = n)
}

# ----------------------------------------------------------------------------
# 7) EJECUCION INTERACTIVA COMPLETA DEL PROGRAMA
# ----------------------------------------------------------------------------

#' @param out_file ruta del archivo de salida; use "" para no guardar
#' @param c_values valores de c a evaluar
#'
#' @return invisiblemente, la tabla calculada
run_probkn_interactive <- function(
    out_file = "SEQ02.OUT",
    c_values = seq(0.44, 1.48, by = 0.04)) {

  print_banner()
  params <- get_probkn_params()

  tabla <- prob_kn_table(
    alpha = params$alpha,
    n = params$n,
    c_values = c_values
  )

  cat("\n")
  print_probkn_screen(tabla, params$alpha, params$n)

  if (nzchar(out_file)) {
    cat("\n", file = out_file, append = TRUE)
    out <- c(
      "OUTPUT FROM THE PROGRAM: ProbKn",
      "===============================",
      "Probability, p=P(K(n)>c*sigma)",
      sprintf("alpha =%s and Sample Size = %d",
              format_book_decimal(params$alpha, 3), as.integer(params$n)),
      "---------------------------------",
      "    c     p      c     p      c     p",
      format_probkn_table(tabla)
    )
    cat(paste(out, collapse = "\n"), "\n", file = out_file, append = TRUE)
  }

  invisible(tabla)
}
