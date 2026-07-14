stein_two_stage <- function(alpha, d, m, mu, sigma2, repeticiones) {
  
  #Cuantil t de Student usando la muestra piloto
  a_m <- qt(1 - alpha / 2, df = m - 1)
  
  #Vectores para guardar los resultados
  tamaños <- numeric(repeticiones)
  coberturas <- logical(repeticiones)
  
  for (i in 1:repeticiones) {
    
    #Primera etapa: muetra piloto
    
    muestra_piloto <- rnorm(
      n = m,
      mean = mu,
      sd = sqrt(sigma2)
    )
    
    #Varianza de la muestra piloto
    S2_m <- var(muestra_piloto)
    
    #Tamaño muestral total de Stein
    N <- max(
      m,
      floor((a_m^2 * S2_m) / d^2) + 1
    )
    
    #Segunda etapa
    
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
    
    #Intervalo final
    
    media_final <- mean(muestra_final)
    
    limite_inferior <- media_final - d
    limite_superior <- media_final + d
    
    #Guardamos N
    tamanios[i] <- N
    
    #El intervalo contiene a mu?
    coberturas[i] <- (
      limite_inferior <= mu &&
        mu <= limite_superior
    )
  }
  
  #RES
  
  resultados <- list(
    promedio_N = mean(tamaños),
    sd_N = sd(tamaños),
    minimo_N = min(tamaños),
    maximo_N = max(tamaños),
    cobertura = mean(coberturas)
  )
  
  return(resultados)
}