# ============================================================================
# PROGRAMA: OC  (subprograma dentro de seq03.exe)
# ----------------------------------------------------------------------------
# OBJETIVO
#   Calcular las aproximaciones de Wald para:
#
#     1) La funcion caracteristica operativa (OC):
#
#            L(theta) = P_theta(aceptar H0)
#
#     2) La funcion de numero muestral promedio (ASN):
#
#            E_theta(N)
#
#   asociadas al SPRT para el contraste de hipotesis simples
#
#            H0: theta = theta0    vs    H1: theta = theta1.
#
#   El programa admite las mismas 9 familias de distribuciones que aparecen
#   en el menu de Seq03.exe y que se implementaron en sprt.R.
#
# FUNDAMENTO TEORICO
#   Sea
#
#       Zi = log{ f(Xi | theta1) / f(Xi | theta0) }
#
#   y sean los limites logaritmicos de Wald:
#
#       a = log((1-beta)/alpha)       limite superior
#       b = log(beta/(1-alpha))       limite inferior.
#
#   Para cada valor verdadero de theta, se busca t0 tal que:
#
#       E_theta[ exp(t0*Z1) ] = 1.
#
#   El programa del libro recorre una grilla de valores de t0 y obtiene el
#   theta correspondiente. Luego calcula:
#
#       L(theta) ~= [exp(t0*a)-1] / [exp(t0*a)-exp(t0*b)],  si t0 != 0
#
#       L(theta) ~= a / (a-b),                              si t0 = 0
#
#   y
#
#       E_theta(N) ~= [b*L(theta)+a*(1-L(theta))] / E_theta(Z1),
#                      si E_theta(Z1) != 0
#
#       E_theta(N) ~= -a*b / E_theta(Z1^2),
#                      si E_theta(Z1) = 0.
#
#   Estas formulas ignoran el exceso del log-cociente de verosimilitud por
#   encima o por debajo de los limites en el momento de parada, tal como lo
#   hace la aproximacion desarrollada en el libro.
# ============================================================================


# ----------------------------------------------------------------------------
# 1) VALIDACIONES Y LIMITES DE WALD
# ----------------------------------------------------------------------------

#' Valida los parametros generales del SPRT utilizado por OC
validate_oc_params <- function(dist_id, theta0, theta1, known, alpha, beta) {
  dist_id <- as.character(dist_id)

  if (!dist_id %in% as.character(1:9)) {
    stop("dist_id debe ser un valor entre '1' y '9'")
  }
  if (length(theta0) != 1 || !is.finite(theta0) ||
      length(theta1) != 1 || !is.finite(theta1) || theta0 == theta1) {
    stop("theta0 y theta1 deben ser valores finitos y diferentes")
  }
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha debe estar en (0,1)")
  }
  if (length(beta) != 1 || !is.finite(beta) || beta <= 0 || beta >= 1) {
    stop("beta debe estar en (0,1)")
  }
  if (alpha + beta >= 1) {
    stop("se requiere alpha + beta < 1")
  }

  # Todas las distribuciones salvo la Normal requieren theta positivo.
  if (dist_id != "1" && (theta0 <= 0 || theta1 <= 0)) {
    stop("theta0 y theta1 deben ser positivos para esta distribucion")
  }

  # Bernoulli y Geometric usan theta como probabilidad de exito.
  if (dist_id %in% c("7", "8") &&
      (theta0 >= 1 || theta1 >= 1)) {
    stop("theta0 y theta1 deben estar en (0,1) para esta distribucion")
  }

  # Parametros conocidos de las distribuciones 1, 2, 3, 4 y 9.
  if (dist_id %in% c("1", "2", "3", "4", "9")) {
    if (length(known) != 1 || !is.finite(known) || known <= 0) {
      stop("el parametro conocido debe ser un valor positivo")
    }
  }

  # En Erlang, k debe ser entero positivo.
  if (dist_id == "4" && known != as.integer(known)) {
    stop("k debe ser un entero positivo para la distribucion Erlang")
  }

  invisible(TRUE)
}


#' Calcula los limites logaritmicos de Wald usando la notacion del libro
#'
#' @return lista con a (limite superior) y b (limite inferior)
oc_wald_boundaries <- function(alpha, beta) {
  list(
    a = log((1 - beta) / alpha),
    b = log(beta / (1 - alpha))
  )
}


# ----------------------------------------------------------------------------
# 2) FUNCION OC Y FUNCION ASN A PARTIR DE t0 Y LOS MOMENTOS DE Z1
# ----------------------------------------------------------------------------

#' Aproximacion de Wald para la funcion OC, L(theta)
#'
#' @param t0    valor o vector de valores que satisfacen E(exp(t0*Z1)) = 1
#' @param alpha probabilidad nominal de error de Tipo I
#' @param beta  probabilidad nominal de error de Tipo II
#'
#' @return vector con L(theta) = P_theta(aceptar H0)

#Definimos OC para cada caso
oc_probability <- function(t0, alpha, beta) {
  lim <- oc_wald_boundaries(alpha, beta)
  a <- lim$a
  b <- lim$b

  t0 <- as.numeric(t0)
  out <- numeric(length(t0))
  #caso t = 0
  zero <- abs(t0) < 1e-10

  out[zero] <- a / (a - b)
  
  #caso t!=0
  if (any(!zero)) {
    t <- t0[!zero]
    out[!zero] <- expm1(t * a) / (exp(t * a) - exp(t * b))
  }

  out
}


#' Aproximacion de Wald para la funcion ASN, E_theta(N)
#'
#' @param L       valores de la funcion OC
#' @param mean_z  E_theta(Z1)
#' @param second_z E_theta(Z1^2)
#' @param alpha,beta errores nominales del SPRT
#'
#' @return vector con la aproximacion del numero muestral promedio

#Definimos ASN para cada caso
asn_from_moments <- function(L, mean_z, second_z, alpha, beta) {
  lim <- oc_wald_boundaries(alpha, beta)
  a <- lim$a
  b <- lim$b

  L <- as.numeric(L)
  mean_z <- as.numeric(mean_z)
  second_z <- as.numeric(second_z)
  
  #validacion de los datos
  if (!(length(L) == length(mean_z) && length(L) == length(second_z))) {
    stop("L, mean_z y second_z deben tener la misma longitud")
  }
  if (any(!is.finite(second_z)) || any(second_z <= 0)) {
    stop("E_theta(Z1^2) debe ser positivo y finito")
  }

  out <- numeric(length(L))
  zero_mean <- abs(mean_z) < 1e-9
  
  #caso E_theta(Z_1)=0
  out[zero_mean] <- -a * b / second_z[zero_mean]
  
  #caso E_theta(Z_1)!=0
  if (any(!zero_mean)) {
    out[!zero_mean] <-
      (b * L[!zero_mean] + a * (1 - L[!zero_mean])) /
      mean_z[!zero_mean]
  }

  out
}


# ----------------------------------------------------------------------------
# 3) RELACION ENTRE t0 Y theta, Y MOMENTOS DE Z1 PARA CADA DISTRIBUCION
# ----------------------------------------------------------------------------

#Definiciones para cada distribución


# 1. Normal(mean = theta, varianza conocida = sigma2)
theta_normal_from_t <- function(t0, theta0, theta1, sigma2) {
  (theta0 + theta1 - t0 * (theta1 - theta0)) / 2
}

moments_normal <- function(theta, theta0, theta1, sigma2) {
  delta <- theta1 - theta0
  midpoint <- (theta0 + theta1) / 2
  mean_z <- delta / sigma2 * (theta - midpoint)
  var_z <- delta^2 / sigma2

  list(
    mean = mean_z,
    variance = rep(var_z, length(theta)),
    second = rep(var_z, length(theta)) + mean_z^2
  )
}


# 2. Gamma(scale = theta, shape conocido)
theta_gamma_from_t <- function(t0, theta0, theta1, shape) {
  log_ratio <- log(theta1 / theta0)
  c_coef <- 1 / theta0 - 1 / theta1
  out <- numeric(length(t0))
  zero <- abs(t0) < 1e-10

  out[zero] <- log_ratio / c_coef

  if (any(!zero)) {
    t <- t0[!zero]
    out[!zero] <- -expm1(-t * log_ratio) / (t * c_coef)
  }

  out
}

moments_gamma <- function(theta, theta0, theta1, shape) {
  log_ratio <- log(theta1 / theta0)
  c_coef <- 1 / theta0 - 1 / theta1
  mean_z <- shape * (c_coef * theta - log_ratio)
  var_z <- shape * c_coef^2 * theta^2

  list(mean = mean_z, variance = var_z, second = var_z + mean_z^2)
}


# 3. Weibull(scale = theta, shape conocido)
theta_weibull_from_t <- function(t0, theta0, theta1, shape) {
  log_ratio_power <- shape * log(theta1 / theta0)
  c_coef <- theta0^(-shape) - theta1^(-shape)
  theta_power <- numeric(length(t0))
  zero <- abs(t0) < 1e-10

  theta_power[zero] <- log_ratio_power / c_coef

  if (any(!zero)) {
    t <- t0[!zero]
    theta_power[!zero] <-
      -expm1(-t * log_ratio_power) / (t * c_coef)
  }

  theta_power^(1 / shape)
}

moments_weibull <- function(theta, theta0, theta1, shape) {
  log_ratio_power <- shape * log(theta1 / theta0)
  c_coef <- theta0^(-shape) - theta1^(-shape)
  theta_power <- theta^shape
  mean_z <- c_coef * theta_power - log_ratio_power
  var_z <- c_coef^2 * theta_power^2

  list(mean = mean_z, variance = var_z, second = var_z + mean_z^2)
}


# 4. Erlang(exp.mean = theta, k conocido, entero)
theta_erlang_from_t <- function(t0, theta0, theta1, k) {
  theta_gamma_from_t(t0, theta0, theta1, shape = k)
}

moments_erlang <- function(theta, theta0, theta1, k) {
  moments_gamma(theta, theta0, theta1, shape = k)
}


# 5. Exponential(mean = theta)
theta_exponential_from_t <- function(t0, theta0, theta1, known = NA) {
  theta_gamma_from_t(t0, theta0, theta1, shape = 1)
}

moments_exponential <- function(theta, theta0, theta1, known = NA) {
  moments_gamma(theta, theta0, theta1, shape = 1)
}


# 6. Poisson(mean = theta)
theta_poisson_from_t <- function(t0, theta0, theta1, known = NA) {
  delta <- theta1 - theta0
  log_ratio <- log(theta1 / theta0)
  out <- numeric(length(t0))
  zero <- abs(t0) < 1e-10

  out[zero] <- delta / log_ratio

  if (any(!zero)) {
    t <- t0[!zero]
    out[!zero] <- t * delta / expm1(t * log_ratio)
  }

  out
}

moments_poisson <- function(theta, theta0, theta1, known = NA) {
  delta <- theta1 - theta0
  log_ratio <- log(theta1 / theta0)
  mean_z <- -delta + theta * log_ratio
  var_z <- theta * log_ratio^2

  list(mean = mean_z, variance = var_z, second = var_z + mean_z^2)
}


# 7. Bernoulli(P(Success) = theta)
theta_bernoulli_from_t <- function(t0, theta0, theta1, known = NA) {
  log_num_base <- log((1 - theta0) / (1 - theta1))
  log_den_base <- log(theta1 * (1 - theta0) /
                        (theta0 * (1 - theta1)))
  out <- numeric(length(t0))
  zero <- abs(t0) < 1e-10

  out[zero] <- log_num_base / log_den_base

  if (any(!zero)) {
    t <- t0[!zero]
    out[!zero] <- expm1(t * log_num_base) /
      expm1(t * log_den_base)
  }

  out
}

moments_bernoulli <- function(theta, theta0, theta1, known = NA) {
  log_odds_ratio <- log(theta1 * (1 - theta0) /
                          (theta0 * (1 - theta1)))
  constant <- log((1 - theta1) / (1 - theta0))
  mean_z <- theta * log_odds_ratio + constant
  var_z <- theta * (1 - theta) * log_odds_ratio^2

  list(mean = mean_z, variance = var_z, second = var_z + mean_z^2)
}


# 8. Geometric(P(Success) = theta)
#    X es el numero de fracasos antes del primer exito.
theta_geometric_from_t <- function(t0, theta0, theta1, known = NA) {
  d_coef <- log(theta1 / theta0)
  c_coef <- log((1 - theta1) / (1 - theta0))
  out <- numeric(length(t0))
  zero <- abs(t0) < 1e-10

  out[zero] <- -c_coef / (d_coef - c_coef)

  if (any(!zero)) {
    t <- t0[!zero]
    out[!zero] <- -expm1(t * c_coef) /
      (exp(t * d_coef) - exp(t * c_coef))
  }

  out
}

moments_geometric <- function(theta, theta0, theta1, known = NA) {
  d_coef <- log(theta1 / theta0)
  c_coef <- log((1 - theta1) / (1 - theta0))
  mean_x <- (1 - theta) / theta
  var_x <- (1 - theta) / theta^2
  mean_z <- d_coef + c_coef * mean_x
  var_z <- c_coef^2 * var_x

  list(mean = mean_z, variance = var_z, second = var_z + mean_z^2)
}


# 9. Negative Binomial(mean = theta, k conocido)
#    p = k/(k+theta), X = numero de fracasos antes del k-esimo exito.
theta_negbinom_from_t <- function(t0, theta0, theta1, k) {
  p0 <- k / (k + theta0)
  p1 <- k / (k + theta1)
  q0 <- 1 - p0
  q1 <- 1 - p1

  d_per_k <- log(p1 / p0)
  c_coef <- log(q1 / q0)
  p <- numeric(length(t0))
  zero <- abs(t0) < 1e-10

  p[zero] <- -c_coef / (d_per_k - c_coef)

  if (any(!zero)) {
    t <- t0[!zero]
    p[!zero] <- -expm1(t * c_coef) /
      (exp(t * d_per_k) - exp(t * c_coef))
  }

  k * (1 - p) / p
}

moments_negbinom <- function(theta, theta0, theta1, k) {
  p0 <- k / (k + theta0)
  p1 <- k / (k + theta1)
  q0 <- 1 - p0
  q1 <- 1 - p1

  d_coef <- k * log(p1 / p0)
  c_coef <- log(q1 / q0)
  mean_z <- d_coef + c_coef * theta
  var_x <- theta + theta^2 / k
  var_z <- c_coef^2 * var_x

  list(mean = mean_z, variance = var_z, second = var_z + mean_z^2)
}


# ----------------------------------------------------------------------------
# 4) TABLA DE DISTRIBUCIONES DEL SUBPROGRAMA OC
# ----------------------------------------------------------------------------

OC_DIST_TABLE <- list(
  "1" = list(
    nombre = "Normal",
    theta_label = "mean",
    param_conocido = "variance",
    theta_from_t = theta_normal_from_t,
    moments = moments_normal
  ),
  "2" = list(
    nombre = "Gamma",
    theta_label = "scale",
    param_conocido = "shape",
    theta_from_t = theta_gamma_from_t,
    moments = moments_gamma
  ),
  "3" = list(
    nombre = "Weibull",
    theta_label = "scale",
    param_conocido = "shape",
    theta_from_t = theta_weibull_from_t,
    moments = moments_weibull
  ),
  "4" = list(
    nombre = "Erlang",
    theta_label = "exp.mean",
    param_conocido = "k (integer)",
    theta_from_t = theta_erlang_from_t,
    moments = moments_erlang
  ),
  "5" = list(
    nombre = "Exponential",
    theta_label = "mean",
    param_conocido = NA,
    theta_from_t = theta_exponential_from_t,
    moments = moments_exponential
  ),
  "6" = list(
    nombre = "Poisson",
    theta_label = "mean",
    param_conocido = NA,
    theta_from_t = theta_poisson_from_t,
    moments = moments_poisson
  ),
  "7" = list(
    nombre = "Bernoulli",
    theta_label = "P(Success)",
    param_conocido = NA,
    theta_from_t = theta_bernoulli_from_t,
    moments = moments_bernoulli
  ),
  "8" = list(
    nombre = "Geometric",
    theta_label = "P(Success)",
    param_conocido = NA,
    theta_from_t = theta_geometric_from_t,
    moments = moments_geometric
  ),
  "9" = list(
    nombre = "Negative Binomial",
    theta_label = "mean",
    param_conocido = "k",
    theta_from_t = theta_negbinom_from_t,
    moments = moments_negbinom
  )
)


# ----------------------------------------------------------------------------
# 5) MOTOR GENERAL DE CALCULO DE OC Y ASN
# ----------------------------------------------------------------------------

#' Calcula las aproximaciones OC y ASN para una de las 9 distribuciones
#'
#' @param dist_id     character "1".."9", segun OC_DIST_TABLE
#' @param theta0      valor de theta bajo H0
#' @param theta1      valor de theta bajo H1
#' @param known       parametro conocido; NA cuando no corresponde
#' @param alpha,beta  errores nominales del SPRT
#' @param t0_values   grilla de valores de t0; por defecto, la del libro
#'
#' @return data.frame con t0, theta, OC, ASN y momentos de Z1
compute_oc_asn <- function(
    dist_id,
    theta0,
    theta1,
    known = NA,
    alpha,
    beta,
    t0_values = seq(-2.2, 2.2, by = 0.1)) {

  dist_id <- as.character(dist_id)
  validate_oc_params(dist_id, theta0, theta1, known, alpha, beta)

  if (length(t0_values) < 1 || any(!is.finite(t0_values))) {
    stop("t0_values debe contener uno o mas valores finitos")
  }

  t0_values <- round(as.numeric(t0_values), 10)
  dist <- OC_DIST_TABLE[[dist_id]]

  theta <- dist$theta_from_t(
    t0 = t0_values,
    theta0 = theta0,
    theta1 = theta1,
    known
  )

  if (any(!is.finite(theta))) {
    stop("no se pudieron obtener valores finitos de theta para la grilla dada")
  }

  moments <- dist$moments(
    theta = theta,
    theta0 = theta0,
    theta1 = theta1,
    known
  )

  L <- oc_probability(t0_values, alpha, beta)
  ASN <- asn_from_moments(
    L = L,
    mean_z = moments$mean,
    second_z = moments$second,
    alpha = alpha,
    beta = beta
  )

  data.frame(
    t0 = t0_values,
    theta = as.numeric(theta),
    OC = as.numeric(L),
    ASN = as.numeric(ASN),
    E_Z1 = as.numeric(moments$mean),
    Var_Z1 = as.numeric(moments$variance),
    E_Z1_sq = as.numeric(moments$second)
  )
}


# ----------------------------------------------------------------------------
# 6) FORMATO DE SALIDA
# ----------------------------------------------------------------------------


format_oc_decimal <- function(x, digits) {
  texto <- sprintf(paste0("%.", digits, "f"), x)
  texto <- sub("^0\\.", ".", texto)
  texto <- sub("^-0\\.", "-.", texto)
  texto
}


#' Convierte una fila de resultados en un bloque t0-theta-OC-ASN
format_oc_row <- function(row) {
  sprintf(
    "%5s %6s %5s %6s",
    format_oc_decimal(row$t0, 1),
    format_oc_decimal(row$theta, 2),
    format_oc_decimal(row$OC, 3),
    format_oc_decimal(row$ASN, 1)
  )
}

#' Para la grilla por defecto, el bloque izquierdo contiene t0=-2.2,...,0
#' y el derecho t0=0.1,...,2.2.
format_oc_asn_table <- function(tabla) {
  required <- c("t0", "theta", "OC", "ASN")
  if (!all(required %in% names(tabla))) {
    stop("tabla debe contener las columnas t0, theta, OC y ASN")
  }
  if (nrow(tabla) < 1) return(character(0))

  tabla <- tabla[order(tabla$t0), , drop = FALSE]

  left <- tabla[tabla$t0 <= 0, , drop = FALSE]
  right <- tabla[tabla$t0 > 0, , drop = FALSE]

  # Si la grilla no cruza cero, se divide aproximadamente por la mitad.
  if (nrow(left) == 0 || nrow(right) == 0) {
    cut <- ceiling(nrow(tabla) / 2)
    left <- tabla[seq_len(cut), , drop = FALSE]
    right <- if (cut < nrow(tabla)) {
      tabla[(cut + 1):nrow(tabla), , drop = FALSE]
    } else {
      tabla[FALSE, , drop = FALSE]
    }
  }

  n_lines <- max(nrow(left), nrow(right))
  lines <- character(n_lines)

  for (i in seq_len(n_lines)) {
    left_text <- if (i <= nrow(left)) {
      format_oc_row(left[i, , drop = FALSE])
    } else {
      strrep(" ", 25)
    }

    right_text <- if (i <= nrow(right)) {
      format_oc_row(right[i, , drop = FALSE])
    } else {
      ""
    }

    lines[i] <- paste0(left_text, "   ", right_text)
  }

  lines
}


#' Descripcion del parametro theta para el encabezado de salida
oc_distribution_description <- function(dist_id, known) {
  dist_id <- as.character(dist_id)

  switch(
    dist_id,
    "1" = sprintf("mean of normal distribution with variance = %.1f", known),
    "2" = sprintf("scale of gamma distribution with shape = %.2f", known),
    "3" = sprintf("scale of Weibull distribution with shape = %.2f", known),
    "4" = sprintf("exp.mean of Erlang distribution with k = %d", as.integer(known)),
    "5" = "mean of exponential distribution",
    "6" = "mean of Poisson distribution",
    "7" = "P(Success) in Bernoulli distribution",
    "8" = "P(Success) in geometric distribution",
    "9" = sprintf("mean of negative binomial distribution with k = %.2f", known)
  )
}


#' Imprime la salida del subprograma OC y opcionalmente la guarda
#'
#' @param tabla data.frame devuelto por compute_oc_asn()
#' @param dist_id,theta0,theta1,known,alpha,beta parametros del test
#' @param file ruta del archivo de salida; "" imprime solo en consola
print_oc_output <- function(
    tabla,
    dist_id,
    theta0,
    theta1,
    known = NA,
    alpha,
    beta,
    file = "",
    print_console = TRUE) {

  description <- oc_distribution_description(dist_id, known)

  out <- c(
    sprintf(
      "SPRT: H0: theta = %s versus H1: theta = %s with",
      format_oc_decimal(theta0, 2),
      format_oc_decimal(theta1, 2)
    ),
    sprintf(
      "alpha = %s and beta = %s where theta is the",
      format_oc_decimal(alpha, 2),
      format_oc_decimal(beta, 2)
    ),
    description,
    "Computed OC and ASN Values",
    "--------------------------",
    "   t0  theta    OC    ASN      t0  theta    OC    ASN",
    format_oc_asn_table(tabla)
  )

  if (isTRUE(print_console)) {
    cat(paste(out, collapse = "\n"), "\n")
  }

  if (nzchar(file)) {
    cat(paste(out, collapse = "\n"), "\n", file = file, append = TRUE)
  }

  invisible(out)
}


# ----------------------------------------------------------------------------
# 7) GRAFICOS DE LAS FUNCIONES OC Y ASN
# ----------------------------------------------------------------------------

#' Grafica las funciones OC y ASN de manera conjunta
#'
#' Reproduce conceptualmente las Figuras 3.4.2 y 3.4.5 del libro.
plot_oc_asn <- function(tabla) {
  required <- c("theta", "OC", "ASN")
  if (!all(required %in% names(tabla))) {
    stop("tabla debe contener las columnas theta, OC y ASN")
  }

  ord <- order(tabla$theta)
  x <- tabla$theta[ord]
  oc <- tabla$OC[ord]
  asn <- tabla$ASN[ord]

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2))

  plot(
    x, oc,
    type = "l",
    xlab = expression(theta),
    ylab = expression(L(theta)),
    main = "OC Function"
  )

  plot(
    x, asn,
    type = "l",
    xlab = expression(theta),
    ylab = expression(E[theta](N)),
    main = "ASN Function"
  )

  invisible(NULL)
}


# ----------------------------------------------------------------------------
# 8) INTERFAZ DE CONSOLA
# ----------------------------------------------------------------------------

print_oc_dist_menu <- function() {
  menu <- c(
    "****************************************************",
    "*  ID#  DISTRIBUTION       Theta      Known         *",
    "*                      (Parameter1)  Parameter2     *",
    "*   1   Normal             mean       variance      *",
    "*   2   Gamma              scale        shape       *",
    "*   3   Weibull            scale        shape       *",
    "*   4   Erlang           exp.mean   k (integer)     *",
    "*   5   Exponential         mean                    *",
    "*   6   Poisson             mean                    *",
    "*   7   Bernoulli       P(Success)                  *",
    "*   8   Geometric       P(Success)                  *",
    "*   9   Negative Binomial   mean         k          *",
    "****************************************************"
  )

  cat("\nOC and ASN functions for Wald's SPRT\n")
  cat("SPRT for H0: theta=theta0 versus H1: theta=theta1\n")
  cat(paste(menu, collapse = "\n"), "\n\n")
}


#' Pide los parametros por consola y permite confirmarlos o corregirlos
get_oc_params <- function(dist_id) {
  dist <- OC_DIST_TABLE[[as.character(dist_id)]]

  repeat {
    theta0 <- as.numeric(readline("theta0 = "))
    theta1 <- as.numeric(readline("theta1 = "))
    alpha <- as.numeric(readline("alpha = "))
    beta <- as.numeric(readline("beta = "))

    known <- NA_real_
    if (!is.na(dist$param_conocido)) {
      known <- as.numeric(readline(sprintf(
        "value of %s (known) = ",
        dist$param_conocido
      )))
    }

    # Valida antes de confirmar para detectar entradas imposibles.
    validate_oc_params(dist_id, theta0, theta1, known, alpha, beta)

    ans <- readline("Do you want to change the input parameters (Y or N)? ")
    if (toupper(trimws(ans)) != "Y") break
  }

  list(
    theta0 = theta0,
    theta1 = theta1,
    alpha = alpha,
    beta = beta,
    known = known
  )
}


#' Version interactiva completa del subprograma OC
#'
#' @param out_file   archivo donde se guarda la salida; use "" para no guardar
#' @param t0_values  grilla de t0; por defecto -2.2,...,2.2, como el libro
#' @param make_plots si TRUE, muestra los graficos OC y ASN
#'
#' @return invisiblemente, la tabla calculada
run_oc_interactive <- function(
    out_file = "SEQ03.OUT",
    t0_values = seq(-2.2, 2.2, by = 0.1),
    make_plots = TRUE) {

  print_oc_dist_menu()

  dist_id <- trimws(readline("Input distribution ID#: "))
  if (!dist_id %in% names(OC_DIST_TABLE)) {
    stop("ID de distribucion invalido")
  }

  params <- get_oc_params(dist_id)

  tabla <- compute_oc_asn(
    dist_id = dist_id,
    theta0 = params$theta0,
    theta1 = params$theta1,
    known = params$known,
    alpha = params$alpha,
    beta = params$beta,
    t0_values = t0_values
  )

  cat("\n")
  print_oc_output(
    tabla = tabla,
    dist_id = dist_id,
    theta0 = params$theta0,
    theta1 = params$theta1,
    known = params$known,
    alpha = params$alpha,
    beta = params$beta,
    file = ""
  )

  if (nzchar(out_file)) {
    cat("\n", file = out_file, append = TRUE)
    print_oc_output(
      tabla = tabla,
      dist_id = dist_id,
      theta0 = params$theta0,
      theta1 = params$theta1,
      known = params$known,
      alpha = params$alpha,
      beta = params$beta,
      file = out_file,
      print_console = FALSE
    )
  }

  if (isTRUE(make_plots)) {
    plot_oc_asn(tabla)
  }

  invisible(tabla)
}


# ----------------------------------------------------------------------------
# 9) EJEMPLOS DE VALIDACION DEL LIBRO
# ----------------------------------------------------------------------------
#
# EJEMPLO 1: NORMAL, Figura 3.4.3
# --------------------------------
# normal_libro <- compute_oc_asn(
#   dist_id = "1",
#   theta0 = 0,
#   theta1 = 1,
#   known = 4,
#   alpha = 0.05,
#   beta = 0.01
# )
#
# print_oc_output(
#   normal_libro,
#  dist_id = "1",
#   theta0 = 0,
#   theta1 = 1,
#   known = 4,
#   alpha = 0.05,
#   beta = 0.01
# )
#
# plot_oc_asn(normal_libro)
#
#
# EJEMPLO 2: BERNOULLI, Figura 3.4.6
# -----------------------------------
# bernoulli_libro <- compute_oc_asn(
#   dist_id = "7",
#   theta0 = 0.5,
#   theta1 = 0.7,
#   alpha = 0.05,
#   beta = 0.05
# )
#
# print_oc_output(
#   bernoulli_libro,
#   dist_id = "7",
#   theta0 = 0.5,
#   theta1 = 0.7,
#   alpha = 0.05,
#   beta = 0.05
# )
#
# plot_oc_asn(bernoulli_libro)
# ----------------------------------------------------------------------------
