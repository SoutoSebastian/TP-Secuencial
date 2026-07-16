stein_two_stage <- function(alpha, d, m, mu, sigma2, repeticiones) {
  
  # Cuantil t de Student usando la muestra piloto
  a_m <- qt(1 - alpha / 2, df = m - 1)
  
  # Cuantil normal para el caso ideal de varianza conocida
  a <- qnorm(1 - alpha / 2)
  
  # Tamaño muestral óptimo si sigma^2 fuera conocida
  C <- (a^2 * sigma2) / d^2
  
  # Vectores para guardar resultados de las simulaciones
  tamaños <- numeric(repeticiones)
  coberturas <- logical(repeticiones)
  
  for (i in 1:repeticiones) {
    
    # Primera etapa: muestra piloto
    muestra_piloto <- rnorm(
      n = m,
      mean = mu,
      sd = sqrt(sigma2)
    )
    
    # Varianza de la muestra piloto
    S2_m <- var(muestra_piloto)
    
    # Tamaño muestral total de Stein
    N <- max(
      m,
      floor((a_m^2 * S2_m) / d^2) + 1
    )
    
    # Segunda etapa
    if (N > m) {
      
      muestra_segunda <- rnorm(
        n = N - m,
        mean = mu,
        sd = sqrt(sigma2)
      )
      
      muestra_final <- c(
        muestra_piloto,
        muestra_segunda
      )
      
    } else {
      
      muestra_final <- muestra_piloto
    }
    
    # Intervalo final de longitud fija 2d
    media_final <- mean(muestra_final)
    
    limite_inferior <- media_final - d
    limite_superior <- media_final + d
    
    # Guardamos el tamaño muestral final
    tamaños[i] <- N
    
    # Verificamos si el intervalo contiene a mu
    coberturas[i] <- (
      limite_inferior <= mu &&
        mu <= limite_superior
    )
  }
  
  # Resúmenes
  cobertura <- mean(coberturas)
  
  error_estandar_cobertura <- sqrt(
    cobertura * (1 - cobertura) / repeticiones
  )
  
  resultados <- list(
    C = C,
    promedio_N = mean(tamaños),
    sd_N = sd(tamaños),
    minimo_N = min(tamaños),
    maximo_N = max(tamaños),
    cobertura = cobertura,
    se_cobertura = error_estandar_cobertura
  )
  
  return(resultados)
}


# ========================
# INPUTS INTERACTIVOS
# ========================

cat("STEIN TWO-STAGE PROCEDURE\n")
cat("=========================\n\n")

alpha <- as.numeric(
  readline("Ingrese alpha: ")
)

d <- as.numeric(
  readline("Ingrese d (semiancho del intervalo): ")
)

m <- as.integer(
  readline("Ingrese el tamaño de la muestra piloto m: ")
)

mu <- as.numeric(
  readline("Ingrese la media mu de la distribución normal: ")
)

sigma2 <- as.numeric(
  readline("Ingrese la varianza sigma^2: ")
)

repeticiones <- as.integer(
  readline("Ingrese el número de repeticiones: ")
)


# ========================
# EJECUCIÓN
# ========================

resultado <- stein_two_stage(
  alpha = alpha,
  d = d,
  m = m,
  mu = mu,
  sigma2 = sigma2,
  repeticiones = repeticiones
)


# ========================
# OUTPUT ESTILO LIBRO
# ========================

confianza <- 100 * (1 - alpha)
ancho <- 2 * d

cat("\n")
cat("      Fixed-Width Confidence Interval\n")
cat("          Using Two-Stage Procedure\n")
cat("          =========================\n\n")

cat(
  sprintf(
    "Simulation study for %.0f%% confidence interval for\n",
    confianza
  )
)

cat("the normal mean. Data were generated from normal\n")

cat(
  sprintf(
    "distribution with mean = %.2f and variance = %.2f\n\n",
    mu,
    sigma2
  )
)

cat(
  sprintf(
    "Number of simulation replications      = %d\n",
    repeticiones
  )
)

cat(
  sprintf(
    "Width of the confidence interval (2d)  = %.2f\n\n",
    ancho
  )
)

cat(
  sprintf(
    "Initial sample size                    = %d\n\n",
    m
  )
)

cat(
  sprintf(
    "     Optimal sample size (C)            = %7.2f\n\n",
    resultado$C
  )
)

cat(
  sprintf(
    "     Average sample size (n_bar)        = %7.2f\n",
    resultado$promedio_N
  )
)

cat(
  sprintf(
    "     Std. dev. of n (s)                 = %7.2f\n",
    resultado$sd_N
  )
)

cat(
  sprintf(
    "     Minimum sample size (min n)        = %7d\n",
    resultado$minimo_N
  )
)

cat(
  sprintf(
    "     Maximum sample size (max n)        = %7d\n\n",
    resultado$maximo_N
  )
)

cat(
  sprintf(
    "     Coverage probability (p_bar)       = %7.3f\n",
    resultado$cobertura
  )
)

cat(
  sprintf(
    "     Standard error of p_bar            = %7.3f\n",
    resultado$se_cobertura
  )
)

