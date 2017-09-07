################### INIT #####################
Ns <- 24
A <- 27 #cm^2
lam <- .0127 #um to CM
ilimitden <- .86 # mA/cm^2 # to amps/cm^2
RHa <- 1
RHc <- 1
Tt <- 353.15 # K
Pa <- 3 # barm
Pc <- 5 # bar

################ FUNCTIONS ###################

randomChrom <- function(){
  chrom <- c(
    runif(1, min = -1.19969, max = -0.8532), # X1   1 
    runif(1, min = .001, max = .005), #X2    2
    runif(1, min = 3.6e-05, max = 9.8e-05), #X3     3
    runif(1, min = -2.60e-04, max = -9.54e-05), #X4    4
    runif(1, min = 10, max = 24), # h    5
    runif(1, min = .0001, max = .0008), # Rc    6
    runif(1, min = .0136, max = .5 # Bv    7
    ))
  return(chrom)
}
Fitness <- function(actual, model){
  
  SSE <- 0
  for (i in 2:23){
    SSE <- SSE + (Vstack(actual, i) - Vstack(model, i))^2
  }
  return (SSE)
}
Vstack <- function(chrom, i){
  Phsat <- (2.95e-02 * (Tt - 273.15) - 9.19e-05 * ((Tt - 273.15)^2) + 1.44e-07 * ((Tt - 273.15)^3) - 2.18)
  Phsat <- 10^Phsat
  ## Get PH2
  Ph2a <- (RHa * Phsat) / Pa
  Ph2b <- (1.635 * (i/A))/(Tt^1.334)  
  Ph2 <- (.5 * RHa) * Phsat * ((1 / (Ph2a * exp(Ph2b))) -1)
  ## Get PO2
  PO2a <- (RHc * Phsat) / Pc
  PO2b <- (4.192 * (i/A))/(Tt^1.334)  
  PO2 <-  RHc * Phsat * ((1 / (PO2a * exp(PO2b))) -1)
  ## Use PH2 and PO2 to get Enernst
  Enernst <- 1.229 - 0.85 * .001 *( Tt - 298.15) + 4.3085e-05 * Tt * log(Ph2 * sqrt(PO2))
  ## Get Pm
  Pmnum <- 181.6* (1 + (.03 * (i/A)) + .062*(Tt/303)*(i/A)^2.5)
  Pmdenom <- ((chrom[5] - .634 - 3*(i/A)) * exp(4.18 * ((Tt-303)/Tt))) 
  Pm <- Pmnum / Pmdenom
  ## Get Nohm
  Rm <- (Pm * lam) / A
  Rc <- chrom[6]
  Nohm <- i * (Rm + Rc)
  ## Get Nact
  Co <- PO2 / (5.08e06 * exp(-498/Tt))
  Nact <- -(chrom[1] + (chrom[2]* Tt) + (chrom[3] * Tt * log(Co)) + (chrom[4] * Tt * log(i)))
  # Get Nconc
  iden <- (i/A)
  Nconc <- -chrom[7]*log(1 - (iden/ilimitden))
  return (Ns * (Enernst - Nact - Nohm - Nconc))
}
Vstack_actual <- c(-0.944957,0.00301801,7.401e-05,-1.88e-04, 23, .0001, 0.02914489)
Crossover <- function(inDF){
  ## SELECT THE PARENTS
  pop <- df
  selectors <- sample(1:35, 16, replace=F)
  
  if (pop[selectors[1],8] < pop[selectors[2],8]){
    a <- pop[selectors[1],]
  } else { a <- pop[selectors[2],] }
  
  if (pop[selectors[3],8] < pop[selectors[4],8]){
    b <- pop[selectors[3],]
  } else { b <- pop[selectors[4],] }
  
  if (pop[selectors[5],8] < pop[selectors[6],8]){
    c <- pop[selectors[5],]
  } else { c <- pop[selectors[6],] }
  
  if (pop[selectors[7],8] < pop[selectors[8],8]){
    d <- pop[selectors[7],]
  } else { d <- pop[selectors[8],] }
  
  if (pop[selectors[9],8] < pop[selectors[10],8]){
    ee <- pop[selectors[9],]
  } else { ee <- pop[selectors[10],] }
  
  if (pop[selectors[11],8] < pop[selectors[12],8]){
    ff <- pop[selectors[11],]
  } else { ff <- pop[selectors[12],] }
  
  if (pop[selectors[13],8] < pop[selectors[14],8]){
    gg <- pop[selectors[13],]
  } else { gg <- pop[selectors[14],] }
  
  if (pop[selectors[15],8] < pop[selectors[16],8]){
    hh <- pop[selectors[15],]
  } else { hh <- pop[selectors[16],] }
  
  if (a[,8] < b[,8]){
    semi1 <- a
  } else { semi1 <- b }
  
  if (c[,8] < d[,8]){
    semi2 <- c
  } else { semi2 <- d }
  
  if (ee[,8] < ff[,8]){
    semi3 <- ee
  } else { semi3 <- ff }
  
  if (gg[,8] < hh[,8]){
    semi4 <- gg
  } else { semi4 <- hh }
  
  if (semi1[,8] < semi3[,8]){
    a <- semi1
  } else { a <- semi2 }
  
  if (semi2[,8] < semi4[,8]){
    b <- semi2
  } else { b <- semi4 }
  
  ## WE have our parents A and B
  ## randomly select up to 5 crossover points
  selectors <- sample(1:7, runif(1, min=1, max=7), replace=F)
  ## make child
  child <- a
  ## cross
  for(i in 1:length(selectors)){
    child[selectors[i]] <- b[selectors[i]]
  }

  ## Chance for selected gene mutation 11%
  mut <- runif(1, 0, 1)
  if (mut < .11){
    chrom <- randomChrom()
    for(i in 1:length(selectors)){
      child[selectors[i]] <- chrom[selectors[i]]
    }
  }
  return (child)
}
getRandPop <- function(){
  population <- list()
  for (m in 1:50){
    population <- c(population, list(randomChrom()))  
  }
  return (population)
}

################## MAIN BODY ################

## Create initial pop
population <- getRandPop()

for (index in 1:100){
  
  ## Evaluate the fitness of each individual
  df <- data.frame(matrix(unlist(population), nrow=50, byrow=T))
  df$fit <- 0
  colnames(df) <- c('x1', 'x2', 'x3', 'x4', 'h', 'Rc', 'Bv', 'SSE')
  for(j in 1:50){
    df[j,8] <- Fitness(Vstack_actual, unlist(population[j]))  
  }
  df <- df[order(df$SSE),]
  elite <- df[1,]
  if (elite[8] < 2 ){
    print('Vector Found')
    print(elite[1:8])
    break
  }
  ## Create a new population
  population <- NULL
  population <- list()
  for (m in 1:30){
    population <- c(population, list(Crossover(df)))  
  }
  ## Add the elite to new generation
  population <- c(population, as.list(as.data.frame(t(df[1:19,]))))
  population <- c(population, as.list(as.data.frame(t(df[40,]))))
}
  
##########  END MAIN BODY ###########



