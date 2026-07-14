# ============================================================================
# PROGRAMA: SPRT  (subprograma dentro de seq03.exe)
# ----------------------------------------------------------------------------
# Basado en "Sequential Methods and their Applications" de Mukhopadhyay y
# de Silva, Capitulo 3, Seccion 3.2 (SPRT de Wald) y formula (3.2.7)
# (version truncada del SPRT).
#
# OBJETIVO
#   Test secuencial de razon de verosimilitudes (SPRT) para el contraste
#   de hipotesis simples
#         H0: theta = theta0     vs     H1: theta = theta1
#   sobre el parametro theta de 9 familias de distribuciones, con la
#   posibilidad de truncar el muestreo en n = K.
#
# FUNDAMENTO TEORICO
#   Sea Zn = sum_{i=1}^n log( f1(xi) / f0(xi) ) el log-cociente de
#   verosimilitud acumulado. Los limites de Wald son:
#       a = log( beta / (1-alpha) )      (limite inferior)
#       b = log( (1-beta) / alpha )      (limite superior)
#   Regla de decision:
#       Zn <= a  ->  se ACEPTA H0
#       Zn >= b  ->  se RECHAZA H0 (se acepta H1)
#       en otro caso, se sigue observando
#
#   Version TRUNCADA (3.2.7): si se llega a n=K sin cruzar ningun limite,
#   se decide segun el signo de Zk: se acepta H1 si Zk >= 0, se acepta H0
#   si Zk < 0 (regla clasica de Wald para SPRT truncado). 
# ============================================================================


# ----------------------------------------------------------------------------
# 1) MOTOR GENERICO DEL SPRT
# ----------------------------------------------------------------------------

#' Calcula los limites de Wald a partir de alpha y beta
wald_boundaries <- function(alpha, beta) {
  list(a = log(beta / (1 - alpha)),      # limite inferior (acepta H0)
       b = log((1 - beta) / alpha))      # limite superior (acepta H1)
}

#' Corre el SPRT sobre un vector de incrementos log f1/f0 ya calculado
#'
#' @param loglr_inc vector de incrementos log(f1(xi)/f0(xi)), uno por dato
#' @param a         limite inferior (de wald_boundaries)
#' @param b         limite superior (de wald_boundaries)
#' @param K         tamaño de truncamiento (NULL = sin truncar; requiere
#'                  que length(loglr_inc) >= K)
#'
#' @return lista con n (tamaño muestral de parada), decision ("H0"/"H1"),
#'         Z (estadistico final), truncated (logico)
run_sprt <- function(loglr_inc, a, b, K = NULL) {
  n_max <- if (is.null(K)) length(loglr_inc) else min(K, length(loglr_inc))
  Z <- 0
  for (n in seq_len(n_max)) {
    Z <- Z + loglr_inc[n]
    if (Z >= b) return(list(n = n, decision = "H1", Z = Z, truncated = FALSE))
    if (Z <= a) return(list(n = n, decision = "H0", Z = Z, truncated = FALSE))
  }
  # Se llego a K sin cruzar ningun limite -> decision por truncamiento.
  
  
  B <- exp(a)
  A <- exp(b)
  limite <- log((A+B)/2)
  decision <- if (Z > limite) "H1" else "H0"
  list(n = n_max, decision = decision, Z = Z, truncated = TRUE)
}


# ----------------------------------------------------------------------------
# 2) LOG-COCIENTES DE VEROSIMILITUD PARA LAS 9 DISTRIBUCIONES DEL MENU
#    (cada funcion es VECTORIZADA: recibe un vector x y devuelve el vector
#    de incrementos log f1(xi)/f0(xi))
# ----------------------------------------------------------------------------

# 1. Normal(mean = theta, varianza conocida = sigma2)
loglr_normal <- function(x, theta0, theta1, sigma2) {
  (theta1 - theta0) / sigma2 * (x - (theta0 + theta1) / 2)
}

# 2. Gamma(scale = theta, shape conocido = k)
loglr_gamma <- function(x, theta0, theta1, shape) {
  -shape * log(theta1 / theta0) + x * (1 / theta0 - 1 / theta1)
}

# 3. Weibull(scale = theta, shape conocido = k)
loglr_weibull <- function(x, theta0, theta1, shape) {
  -shape * log(theta1 / theta0) + x^shape * (1 / theta0^shape - 1 / theta1^shape)
}

# 4. Erlang(exp.mean = theta, k conocido, entero) -> Gamma con shape = k
loglr_erlang <- function(x, theta0, theta1, k) {
  loglr_gamma(x, theta0, theta1, shape = k)
}

# 5. Exponential(mean = theta)
loglr_exponential <- function(x, theta0, theta1) {
  log(theta0 / theta1) + x * (1 / theta0 - 1 / theta1)
}

# 6. Poisson(mean = theta)
loglr_poisson <- function(x, theta0, theta1) {
  -(theta1 - theta0) + x * log(theta1 / theta0)
}

# 7. Bernoulli(P(exito) = theta)
loglr_bernoulli <- function(x, theta0, theta1) {
  x * log(theta1 / theta0) + (1 - x) * log((1 - theta1) / (1 - theta0))
}

# 8. Geometric(P(exito) = theta); x = num. de fracasos antes del 1er exito
loglr_geometric <- function(x, theta0, theta1) {
  log(theta1 / theta0) + x * log((1 - theta1) / (1 - theta0))
}

# 9. NegativeBinomial(mean = theta, k conocido)
#    p = k / (k + theta);  x = num. de fracasos antes del k-esimo exito
loglr_negbinom <- function(x, theta0, theta1, k) {
  p0 <- k / (k + theta0)
  p1 <- k / (k + theta1)
  k * log(p1 / p0) + x * log((1 - p1) / (1 - p0))
}


# ----------------------------------------------------------------------------
# 3) TABLA DE DISTRIBUCIONES: nombre, funcion de log-verosimilitud,
#    funcion de simulacion (recibe n y el theta VERDADERO + el parametro
#    conocido, cuando corresponda) y nombre del parametro conocido
# ----------------------------------------------------------------------------

DIST_TABLE <- list(
  "1" = list(nombre = "Normal",          param_conocido = "variance",
             loglr = function(x, t0, t1, k) loglr_normal(x, t0, t1, k),
             simular = function(n, theta, k) rnorm(n, mean = theta, sd = sqrt(k))),
  "2" = list(nombre = "Gamma",           param_conocido = "shape",
             loglr = function(x, t0, t1, k) loglr_gamma(x, t0, t1, k),
             simular = function(n, theta, k) rgamma(n, shape = k, scale = theta)),
  "3" = list(nombre = "Weibull",         param_conocido = "shape",
             loglr = function(x, t0, t1, k) loglr_weibull(x, t0, t1, k),
             simular = function(n, theta, k) rweibull(n, shape = k, scale = theta)),
  "4" = list(nombre = "Erlang",          param_conocido = "k (integer)",
             loglr = function(x, t0, t1, k) loglr_erlang(x, t0, t1, k),
             simular = function(n, theta, k) rgamma(n, shape = k, scale = theta)),
  "5" = list(nombre = "Exponential",     param_conocido = NA,
             loglr = function(x, t0, t1, k) loglr_exponential(x, t0, t1),
             simular = function(n, theta, k) rexp(n, rate = 1 / theta)),
  "6" = list(nombre = "Poisson",         param_conocido = NA,
             loglr = function(x, t0, t1, k) loglr_poisson(x, t0, t1),
             simular = function(n, theta, k) rpois(n, lambda = theta)),
  "7" = list(nombre = "Bernoulli",       param_conocido = NA,
             loglr = function(x, t0, t1, k) loglr_bernoulli(x, t0, t1),
             simular = function(n, theta, k) rbinom(n, size = 1, prob = theta)),
  "8" = list(nombre = "Geometric",       param_conocido = NA,
             loglr = function(x, t0, t1, k) loglr_geometric(x, t0, t1),
             simular = function(n, theta, k) rgeom(n, prob = theta)),
  "9" = list(nombre = "NegativeBinomial", param_conocido = "k",
             loglr = function(x, t0, t1, k) loglr_negbinom(x, t0, t1, k),
             simular = function(n, theta, k) rnbinom(n, size = k, prob = k / (k + theta)))
)


# ----------------------------------------------------------------------------
# 4) SPRT SOBRE UN CONJUNTO DE DATOS (reales o ya simulados)
# ----------------------------------------------------------------------------

#' Corre el SPRT sobre un vector de datos (reales o simulados) para una
#' de las 9 distribuciones de DIST_TABLE
#'
#' @param dist_id   character, "1".."9" (ver DIST_TABLE)
#' @param theta0,theta1  valores bajo H0 y H1
#' @param known     valor del parametro conocido (NA si no aplica)
#' @param alpha,beta errores de tipo I y II nominales
#' @param data      vector de datos observados secuencialmente
#' @param K         truncamiento (NULL = sin truncar)
sprt_on_data <- function(dist_id, theta0, theta1, known, alpha, beta, data, K = NULL) {
  dist <- DIST_TABLE[[dist_id]]
  if (is.null(dist)) stop("ID de distribucion invalido")

  inc  <- dist$loglr(data, theta0, theta1, known)
  lim  <- wald_boundaries(alpha, beta)
  run_sprt(inc, a = lim$a, b = lim$b, K = K)
}


# ----------------------------------------------------------------------------
# 5) ESTUDIO DE SIMULACION: repite el SPRT n_reps veces bajo un theta
#    "verdadero" elegido por el usuario, y resume tamaño muestral promedio
#    (ASN) y proporcion de veces que se decidio H1 (potencia empirica /
#    alpha empirico, segun si true_theta = theta0 o theta1)
# ----------------------------------------------------------------------------

#' @param dist_id        "1".."9"
#' @param theta0,theta1  hipotesis
#' @param known          parametro conocido (NA si no aplica)
#' @param alpha,beta     errores nominales
#' @param K              truncamiento deseado por el usuario. Si es NULL,
#'                        no hay truncamiento real: el SPRT corre hasta
#'                        que decide (lo cual ocurre con probabilidad 1
#'                        si theta0 != theta1). Para poder generar los
#'                        datos con las funciones vectorizadas de R, se
#'                        genera igual un numero grande de observaciones
#'                        (SAFE_CAP), pero eso NO se reporta como
#'                        truncamiento salvo que efectivamente se agote
#'                        (lo cual seria anormal y se avisa con un warning)
#' @param true_theta     valor real de theta usado para simular los datos
#' @param n_reps         numero de replicaciones
#' @param seed           semilla para reproducibilidad
#' @param store_file     si no es NULL, se escriben ahi los tamaños
#'                        muestrales de cada replicacion (uno por linea)
SAFE_CAP <- 100000L   # tope de seguridad para simulacion SIN truncamiento

sprt_simulation <- function(dist_id, theta0, theta1, known, alpha, beta, K,
                             true_theta, n_reps, seed = NULL, store_file = NULL) {
  dist <- DIST_TABLE[[dist_id]]
  if (is.null(dist)) stop("ID de distribucion invalido")
  if (!is.null(seed)) set.seed(seed)

  truncar <- !is.null(K)
  n_gen   <- if (truncar) K else SAFE_CAP   # cuantos datos generar por replica

  lim <- wald_boundaries(alpha, beta)
  ns         <- integer(n_reps)
  decisiones <- character(n_reps)
  truncados  <- logical(n_reps)

  for (r in seq_len(n_reps)) {
    datos <- dist$simular(n_gen, true_theta, known)
    inc   <- dist$loglr(datos, theta0, theta1, known)
    # Si el usuario pidio truncar, se lo pasamos a run_sprt tal cual.
    # Si NO pidio truncar, no le pasamos K (NULL) para que run_sprt use
    # todo el vector generado (SAFE_CAP) sin considerarlo "truncamiento".
    res   <- run_sprt(inc, a = lim$a, b = lim$b, K = if (truncar) K else NULL)
    ns[r]         <- res$n
    decisiones[r] <- res$decision
    truncados[r]  <- truncar && res$truncated
    if (!truncar && res$n == SAFE_CAP) {
      warning(sprintf(
        "Replica %d no decidio dentro de SAFE_CAP=%d observaciones; revisar theta0/theta1/alpha/beta",
        r, SAFE_CAP))
    }
  }

  if (!is.null(store_file)) {
    writeLines(as.character(ns), con = store_file)
  }

  summary <- summarize_sprt_runs(ns, decisiones, truncados)
  c(summary,
    list(prop_H1    = mean(decisiones == "H1"),
         prop_H0    = mean(decisiones == "H0"),
         n_vals     = ns,
         decisiones = decisiones))
}


#' Resume un conjunto de corridas del SPRT (sea por simulacion o con datos
#' reales): cuenta truncamientos, decisiones a favor de H0, y estadisticas
#' del tamaño muestral (promedio, desvio, min, max).
#'
#' @param ns          vector de tamaños muestrales de parada
#' @param decisiones  vector character ("H0"/"H1"), una por corrida
#' @param truncados   vector logico, TRUE si esa corrida trunco en K
summarize_sprt_runs <- function(ns, decisiones, truncados) {
  list(n_trunc          = sum(truncados),
       n_H0_after_trunc = sum(truncados & decisiones == "H0"),
       n_total          = length(ns),
       n_H0             = sum(decisiones == "H0"),
       ASN              = mean(ns),
       sd_n             = if (length(ns) > 1) sd(ns) else NA_real_,
       min_n            = min(ns),
       max_n            = max(ns))
}


#' Imprime el resultado de un conjunto de corridas del SPRT (simuladas o
#' con datos reales) con el mismo formato "OUTPUT FROM THE PROGRAM: SPRT"
#' del libro. Es generica: la usan tanto el modo simulacion como el modo
#' de datos reales.
#'
#' @param summary     lista devuelta por summarize_sprt_runs()
#' @param dist        entrada de DIST_TABLE
#' @param theta0,theta1,alpha,beta  parametros del test
#' @param K           truncamiento usado (NULL si no se trunco)
#' @param data_line   linea que describe el origen de los datos (ej.
#'                    "Results for simulated data from ... with mean = X"
#'                    o "Results for user-provided (real) data")
#' @param count_label etiqueta para el conteo total de corridas (ej.
#'                    "Number of simulations" o "Number of datasets")
#' @param store_file  si no es NULL, se agrega la linea de archivo guardado
print_sprt_output <- function(summary, dist, theta0, theta1, alpha, beta, K,
                               data_line, count_label, store_file = NULL) {
  out <- c(
    "OUTPUT FROM THE PROGRAM: SPRT",
    "==============================",
    sprintf("SPRT: H0: theta = %.2f versus H1: theta = %.2f at alpha = %.2f and beta = %.2f",
            theta0, theta1, alpha, beta),
    sprintf("where theta is the mean of %s distribution", tolower(dist$nombre)),
    if (!is.null(store_file))
      sprintf("Simulated sample sizes are given in file: %s", store_file) else NULL,
    data_line,
    "-----------------------------------------",
    if (!is.null(K))
      sprintf("No. of SPRTs truncated (with K = %d) = %d", K, summary$n_trunc)
    else
      "No truncation was used (K not specified)",
    if (!is.null(K))
      sprintf("No. of times H0 is accepted after truncation = %d", summary$n_H0_after_trunc)
    else NULL,
    sprintf("%s = %d", count_label, summary$n_total),
    sprintf("Number of times H0 is accepted = %d", summary$n_H0),
    sprintf("Average sample size (nbar) = %.2f", summary$ASN),
    sprintf("Std. dev. of n(s) = %s",
            if (is.na(summary$sd_n)) "NA (only one run)" else sprintf("%.2f", summary$sd_n)),
    sprintf("Minimum sample size (minn) = %d", summary$min_n),
    sprintf("Maximum sample size (maxn) = %d", summary$max_n)
  )
  cat(paste(out, collapse = "\n"), "\n")
}


#' Wrapper de print_sprt_output() para el caso de datos SIMULADOS
print_sprt_simulation_output <- function(res, dist, theta0, theta1, alpha, beta,
                                          K, true_theta, store_file = NULL) {
  data_line <- sprintf("Results for simulated data from %s distribution with mean = %.2f",
                        tolower(dist$nombre), true_theta)
  print_sprt_output(res, dist, theta0, theta1, alpha, beta, K,
                     data_line = data_line, count_label = "Number of simulations",
                     store_file = store_file)
}


# ----------------------------------------------------------------------------
# 6) INTERFAZ DE CONSOLA (imita la pantalla de seq03.exe -> subprograma SPRT)
# ----------------------------------------------------------------------------

print_dist_menu <- function() {
  menu <- c(
    "****************************************************",
    "*  ID#  DISTRIBUTION       Theta      Known         *",
    "*                      (Parameter1)  Parameter2     *",
    "*   1   Normal             mean       variance      *",
    "*   2   Gamma               scale        shape       *",
    "*   3   Weibull             scale        shape       *",
    "*   4   Erlang            exp.mean   k (integer)     *",
    "*   5   Exponential          mean                    *",
    "*   6   Poisson              mean                    *",
    "*   7   Bernoulli        P(Success)                  *",
    "*   8   Geometric        P(Success)                  *",
    "*   9   NegativeBinomial     mean         k           *",
    "****************************************************"
  )
  cat("\nSPRT for H0: theta=theta0 versus H1: theta=theta1\n")
  cat(paste(menu, collapse = "\n"), "\n\n")
}

#' Pide theta0, theta1, alpha, beta (y el parametro conocido si corresponde)
#' con la misma logica de confirmar/corregir (Y/N) que usamos en Power
get_sprt_params <- function(dist_id) {
  dist <- DIST_TABLE[[dist_id]]
  repeat {
    theta0 <- as.numeric(readline("theta0 = "))
    theta1 <- as.numeric(readline("theta1 = "))
    alpha  <- as.numeric(readline("alpha  = "))
    beta   <- as.numeric(readline("beta   = "))

    known <- NA
    if (!is.na(dist$param_conocido)) {
      known <- as.numeric(readline(sprintf("value of %s (known) = ",
                                            dist$param_conocido)))
    }

    ans <- readline("Do you want to change the input parameters (Y or N)? ")
    if (toupper(trimws(ans)) != "Y") break
  }
  list(theta0 = theta0, theta1 = theta1, alpha = alpha, beta = beta, known = known)
}

#' Corre el SPRT sobre UN conjunto de datos reales (sin imprimir nada,
#' solo calcula). El detalle paso a paso (Zn) no se muestra a proposito:
#' la salida final se arma en run_sprt_real_data() con el mismo formato
#' que el modo simulado.
#'
#' @param dist    entrada de DIST_TABLE
#' @param theta0,theta1,known,alpha,beta  parametros del test
#' @param datos   vector de datos observados secuencialmente
#' @param K       truncamiento (NULL = sin truncar)
run_sprt_on_dataset <- function(dist, theta0, theta1, known, alpha, beta, datos, K = NULL) {
  lim <- wald_boundaries(alpha, beta)
  inc <- dist$loglr(datos, theta0, theta1, known)
  run_sprt(inc, a = lim$a, b = lim$b, K = K)
}


#' Pide uno o varios conjuntos de datos reales, corre el SPRT sobre cada
#' uno, y muestra un unico resultado final con el mismo formato que el
#' modo de simulacion (incluyendo si hubo truncamiento).
run_sprt_real_data <- function(dist, p, K = NULL) {
  n_sets <- as.integer(readline("How many datasets do you want to enter? : "))
  if (is.na(n_sets) || n_sets < 1) n_sets <- 1

  resultados <- vector("list", n_sets)
  for (i in seq_len(n_sets)) {
    cat(sprintf("\n--- Dataset %d of %d ---\n", i, n_sets))
    texto <- readline("Enter the data values separated by spaces or commas: ")
    datos <- parse_numeric_list(texto)
    resultados[[i]] <- run_sprt_on_dataset(dist, p$theta0, p$theta1, p$known,
                                            p$alpha, p$beta, datos, K)
  }

  ns         <- sapply(resultados, function(r) r$n)
  decisiones <- sapply(resultados, function(r) r$decision)
  truncados  <- sapply(resultados, function(r) isTRUE(r$truncated))

  summary <- summarize_sprt_runs(ns, decisiones, truncados)
  print_sprt_output(summary, dist, p$theta0, p$theta1, p$alpha, p$beta, K,
                     data_line = "Results for user-provided (real) data",
                     count_label = "Number of datasets")

  invisible(resultados)
}


#' Version interactiva completa del subprograma SPRT
run_sprt_interactive <- function() {
  print_dist_menu()
  dist_id <- trimws(readline("Input distribution ID# : "))
  if (!dist_id %in% names(DIST_TABLE)) stop("ID de distribucion invalida")
  dist <- DIST_TABLE[[dist_id]]

  p <- get_sprt_params(dist_id)

  trunc_ans <- readline("Would you like to truncate the SPRT when sample size reaches K (Y or N)? ")
  K <- NULL
  if (toupper(trimws(trunc_ans)) == "Y") {
    K <- as.integer(readline("Enter a value for K: "))
  }

  modo <- toupper(trimws(readline(
    "Do you want to simulate data or input real data? Type S (Simulation) or D (Real Data): ")))

  if (modo == "D") {
    return(invisible(run_sprt_real_data(dist, p, K)))
  }

  # modo == "S": simulacion
  store_ans <- toupper(trimws(readline("Do you want to store simulated sample sizes (Y or N)? ")))
  store_file <- NULL
  if (store_ans == "Y") {
    store_file <- readline("Input file name to store sample sizes (no more than 9 characters): ")
  }

  n_reps <- as.integer(readline("Number of replications for a simulation? "))
  seed   <- as.integer(readline("Enter a positive integer (<10) to initialize the random number sequence: "))

  true_theta <- as.numeric(readline(sprintf(
    "For the data simulation, input %s of the %s distribution: ",
    if (is.na(dist$param_conocido)) "the true value of theta" else "theta",
    dist$nombre)))

  res <- sprt_simulation(dist_id, p$theta0, p$theta1, p$known, p$alpha, p$beta,
                          K, true_theta, n_reps, seed, store_file)

  cat("\n")
  print_sprt_simulation_output(res, dist, p$theta0, p$theta1, p$alpha, p$beta,
                                K, true_theta, store_file)

  invisible(res)
}


# ----------------------------------------------------------------------------
# 7) UTILIDAD (reutilizada de power.R)
# ----------------------------------------------------------------------------
parse_numeric_list <- function(texto) {
  partes <- strsplit(texto, "[, ]+")[[1]]
  partes <- partes[nzchar(partes)]
  as.numeric(partes)
}



