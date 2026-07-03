# ============================================================================
# PROGRAMA: Power
# ----------------------------------------------------------------------------
# Basado en el programa "Power" incluido en seq02.exe, del libro
#"Sequential Methods and their Applications" de Nitis Mukhopadhyay
# y Basil M. de Silva (Capitulo 2, Seccion 2.5, Figura 2.5.2).
#
# OBJETIVO
#   Calcular la probabilidad de error de Tipo II (beta) dada un tamaño de muestra n 
#   o el tamaño de muestra dado un beta del test MP (mas potente / Neyman-Pearson), 
#   para el contraste de hipotesis simples
#
#         H0: mu = mu0     vs     H1: mu = mu1
#
#   sobre la media de una poblacion Normal(mu, sigma^2), con sigma
#   CONOCIDO, a un nivel de significacion alpha.
#
# FUNDAMENTO TEORICO
#   Sean X1,...,Xn iid N(mu, sigma^2), sigma conocido. El test MP
#   (Neyman-Pearson) para H0: mu=mu0 vs H1: mu=mu1 rechaza H0 cuando:
#
#       si mu1 > mu0:  Xbar > mu0 + z_alpha * sigma/sqrt(n)
#       si mu1 < mu0:  Xbar < mu0 - z_alpha * sigma/sqrt(n)
#
#   donde z_alpha = Phi^{-1}(1-alpha) (cuantil superior alpha de la Normal
#   estandar). La potencia de este test evaluada en mu1 es:
#
#       Power(n) = 1 - Phi( z_alpha - |mu1-mu0| * sqrt(n) / sigma )
#       beta(n)  =     Phi( z_alpha - |mu1-mu0| * sqrt(n) / sigma )
#   
#   y en caso de querer calcular n, se despeja de esa fórmula   
#
#   Estos son exactamente los cálculos que realiza el libro.
#
#
# VALIDACION CON EL EJEMPLO DEL LIBRO (Figura 2.5.2)
#   H0: mu=0 vs H1: mu=1, sigma=2, alpha=0.05, n=16  ->  beta = 0.361
# ============================================================================


#' Imprime la pantalla de bienvenida del programa (imita la del original)
print_banner <- function() {
  banner <- c(
    "*********************************************************",
    "*  This program considers the testing of hypotheses:    *",
    "*                                                        *",
    "*     H0: mu = mu0  against  H1: mu = mu1  at level      *",
    "*     alpha                                              *",
    "*                                                        *",
    "*  where mu is the mean of a normal population and       *",
    "*  sigma is the known standard deviation.                *",
    "*                                                        *",
    "*  The program computes:                                 *",
    "*                                                        *",
    "*  (a) the probability of type II error for a given      *",
    "*      sample size, or                                   *",
    "*  (b) the sample size for a given probability of        *",
    "*      type II error.                                    *",
    "*                                                        *",
    "*********************************************************",
    "",
    "Now, input the parameter values for the above hypotheses."
  )
  cat(paste(banner, collapse = "\n"), "\n\n")
}


#' Calcula beta y la potencia del test MP de tamaño fijo n
#'
#' @param mu0    valor de la media bajo H0
#' @param mu1    valor de la media bajo H1 (mu1 != mu0)
#' @param sigma  desvio estandar poblacional, conocido (sigma > 0)
#' @param alpha  nivel de significacion (0 < alpha < 1)
#' @param n      tamaño de muestra (n > 0, entero)
#'
#' @return una lista con mu0, mu1, sigma, alpha, n, z_alpha, beta, power
power_MP <- function(mu0, mu1, sigma, alpha, n) {
  if (sigma <= 0)              stop("sigma debe ser positivo")
  if (alpha <= 0 || alpha >= 1) stop("alpha debe estar en (0,1)")
  if (n <= 0)                  stop("n debe ser positivo")
  if (mu0 == mu1)               stop("mu0 y mu1 no pueden ser iguales")

  z_alpha  <- qnorm(1 - alpha)
  delta    <- abs(mu1 - mu0)
  beta     <- pnorm(z_alpha - delta * sqrt(n) / sigma)
  potencia <- 1 - beta

  list(mu0 = mu0, mu1 = mu1, sigma = sigma, alpha = alpha, n = n,
       z_alpha = z_alpha, beta = beta, power = potencia)
}


#' Calcula el tamaño de muestra n necesario para lograr un beta dado
#'
#' @param mu0    valor de la media bajo H0
#' @param mu1    valor de la media bajo H1 (mu1 != mu0)
#' @param sigma  desvio estandar poblacional, conocido (sigma > 0)
#' @param alpha  nivel de significacion (0 < alpha < 1)
#' @param beta   probabilidad de error de Tipo II objetivo (0 < beta < 1)
#'
#' @return una lista con mu0, mu1, sigma, alpha, beta, z_alpha, z_beta,
#'         n_exact (valor real, no entero) y n (redondeado hacia arriba)
sample_size <- function(mu0, mu1, sigma, alpha, beta) {
  if (sigma <= 0)               stop("sigma debe ser positivo")
  if (alpha <= 0 || alpha >= 1) stop("alpha debe estar en (0,1)")
  if (beta <= 0 || beta >= 1)   stop("beta debe estar en (0,1)")
  if (mu0 == mu1)               stop("mu0 y mu1 no pueden ser iguales")

  z_alpha <- qnorm(1 - alpha)
  z_beta  <- qnorm(beta)          # cuantil de cola inferior (negativo si beta<0.5)
  delta   <- abs(mu1 - mu0)

  n_exact <- ((sigma * (z_alpha - z_beta)) / delta)^2
  n       <- ceiling(n_exact)     # el tamaño de muestra debe ser entero

  list(mu0 = mu0, mu1 = mu1, sigma = sigma, alpha = alpha, beta = beta,
       z_alpha = z_alpha, z_beta = z_beta, n_exact = n_exact, n = n)
}


#' Imprime el resultado con un formato similar al de SEQ02.OUT
#'
#' @param res   lista devuelta por power_MP()
#' @param file  ruta de archivo de salida ("" imprime solo en consola)
print_power_output <- function(res, file = "") {
  out <- c(
    "OUTPUT FROM THE PROGRAM: Power",
    "==============================",
    sprintf("For the Test H0: mu = %s against H1: mu = %s at alpha = %.3f",
            format(res$mu0), format(res$mu1), res$alpha),
    sprintf("where mu is the mean of normal with sigma = %.1f", res$sigma),
    sprintf("If n = %d then beta = %.3f  (power = %.3f)",
            res$n, res$beta, res$power)
  )
  cat(paste(out, collapse = "\n"), "\n")
  if (nzchar(file)) {
    cat(paste(out, collapse = "\n"), "\n", file = file, append = TRUE)
  }
}


#' Imprime el resultado de sample_size() con formato similar al de SEQ02.OUT
#'
#' @param res   lista devuelta por sample_size()
#' @param file  ruta de archivo de salida ("" imprime solo en consola)
print_sample_size_output <- function(res, file = "") {
  out <- c(
    "OUTPUT FROM THE PROGRAM: Power",
    "==============================",
    sprintf("For the Test H0: mu = %s against H1: mu = %s at alpha = %.3f",
            format(res$mu0), format(res$mu1), res$alpha),
    sprintf("where mu is the mean of normal with sigma = %.1f", res$sigma),
    sprintf("If beta = %.3f then n = %d  (n exacto = %.3f)",
            res$beta, res$n, res$n_exact)
  )
  cat(paste(out, collapse = "\n"), "\n")
  if (nzchar(file)) {
    cat(paste(out, collapse = "\n"), "\n", file = file, append = TRUE)
  }
}


#' Pide mu0, mu1, sigma, alpha y permite confirmar/corregir (Y/N),
#' tal como lo hace el programa original.
get_hypothesis_params <- function() {
  repeat {
    mu0   <- as.numeric(readline("mu0 = "))
    mu1   <- as.numeric(readline("mu1 = "))
    sigma <- as.numeric(readline("sigma = "))
    alpha <- as.numeric(readline("alpha = "))

    ans <- readline("Do you want to change the input parameters (Y or N)? ")
    if (toupper(trimws(ans)) != "Y") break
  }
  list(mu0 = mu0, mu1 = mu1, sigma = sigma, alpha = alpha)
}


#' Muestra el segundo menu del programa (opciones 1, 2 o 3)
print_options_menu <- function() {
  cat("\nTo compute probability of the type II error (beta), input the sample size (n)\n")
  cat("OR\n")
  cat("To compute the sample size (n), input the type II error probability (beta).\n\n")

  menu <- c(
    "****************************************************",
    "*  Choose one of the following options:             *",
    "*                                                    *",
    "*  1 : Only one computation of beta or n             *",
    "*  2 : Compute beta values for a set of n values     *",
    "*  3 : Compute n values for a set of beta values     *",
    "*                                                    *",
    "****************************************************"
  )
  cat(paste(menu, collapse = "\n"), "\n\n")
}


#' Interpreta un valor ingresado por el usuario como beta (si esta en (0,1))
#' o como n (si es >= 1), y hace el calculo correspondiente. Imprime el
#' resultado con el mismo formato que el programa original.
compute_beta_or_n <- function(params, value, file = "") {
  if (value > 0 && value < 1) {
    res <- sample_size(params$mu0, params$mu1, params$sigma, params$alpha, value)
    linea <- sprintf("When beta = %.3f sample size, n=%d", value, res$n)
  } else {
    n   <- as.integer(round(value))
    res <- power_MP(params$mu0, params$mu1, params$sigma, params$alpha, n)
    linea <- sprintf("When n = %d type II error probability, beta=%.3f", n, res$beta)
  }
  cat(linea, "\n")
  if (nzchar(file)) cat(linea, "\n", file = file, append = TRUE)
  invisible(res)
}


#' Parsea una lista de numeros separados por espacios y/o comas
parse_numeric_list <- function(texto) {
  partes <- strsplit(texto, "[, ]+")[[1]]
  partes <- partes[nzchar(partes)]
  as.numeric(partes)
}


#' Version interactiva completa (imita el uso por consola de seq02.exe)
run_power_interactive <- function(out_file = "SEQ02.OUT") {
  print_banner()

  params <- get_hypothesis_params()

  print_options_menu()
  opcion <- as.integer(readline("Enter your option : "))

  if (opcion == 1) {
    value <- as.numeric(readline("Enter a value for beta or n : "))
    compute_beta_or_n(params, value, file = out_file)

  } else if (opcion == 2) {
    texto <- readline("Enter the set of n values (separados por espacio o coma): ")
    n_vals <- parse_numeric_list(texto)
    for (n in n_vals) compute_beta_or_n(params, n, file = out_file)

  } else if (opcion == 3) {
    texto <- readline("Enter the set of beta values (separados por espacio o coma): ")
    beta_vals <- parse_numeric_list(texto)
    for (b in beta_vals) compute_beta_or_n(params, b, file = out_file)

  } else {
    stop("Opcion invalida: debe ser 1, 2 o 3")
  }

  invisible(NULL)
}


