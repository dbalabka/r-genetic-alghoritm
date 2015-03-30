#importFrom(graphics, "plot") 
setClassUnion("numericOrNA", members = c("numeric", "logical"))

my_ga <- setClass(
  Class="my_ga",
  slots=c(
    fitness="function",
    fitnessValues="numericOrNA",
    fitnessValue="numeric",
    pcrossover="numeric", 
    pmutation="numeric",
    popSize="numeric",
    step="numeric",
    min="numeric",
    max="numeric",
    population="matrix",
    maxiter="numeric",
    iter="numeric",
    monitor="function",
    maxfitness="numeric",
    eps="numeric",
    keepBest="logical",
    solution="matrix",
    run="numeric",
    maxRun="numeric",
    elitism="numeric",
    summary="matrix",
    seed="numeric",
    type="character"
  ),
  prototype=list(
    pcrossover = 0.8, 
    pmutation = 0.1,
    popSize = 50,
    step = 0,
    maxiter = 100,
    iter = 0,
    maxfitness = Inf,
    eps=sqrt(.Machine$double.eps),
    keepBest=TRUE,
    run=0,
    seed=NULL,
    type="real-valued"
  )
);

setMethod (
  "initialize",
  "my_ga",
  function (.Object, min, max, fitness, popSize, monitor, maxfitness, maxiter, maxRun ,elitism, seed, ...) {
    if (missing(min) | missing(max)) {
      stop("Minimum and maximum range of real number must be set");
    }
    .Object@min <- min;
    .Object@max <- max;
    
    if (missing(fitness)) {
      stop("Fitness function must provided");
    }
    .Object@fitness <- fitness;
    
    if (!missing(popSize)) {
      if (!is.numeric(popSize) | popSize <= 0) {
        stop("Population size must be number bigger than 0");
      }
      .Object@popSize <- popSize;
    }
    .Object@fitnessValues <- rep(NA, .Object@popSize)
    
    if (!missing(monitor)) {
      if (!is.function(monitor)) {
        stop("monitor argument must be function")
      }
      .Object@monitor <- monitor
    }
    
    if (!missing(maxfitness)) {
      if (!is.numeric(maxfitness)) {
        stop("maxfitness argument must be numeric")
      }
      .Object@maxfitness <- maxfitness
    }
    
    if (!missing(maxiter)) {
      if (!is.numeric(maxiter)) {
        stop("maxiter argument must be numeric")
      }
      .Object@maxiter <- maxiter
    }
    
    if (!missing(maxRun)) {
      if (!is.numeric(maxRun)) {
        stop("max argument must be numeric")
      }
      .Object@maxRun <- maxRun
    } else {
      .Object@maxRun <- .Object@maxiter
    }
    
    if (!missing(elitism)) {
      if (!is.numeric(elitism)) {
        stop("elitism argument must be numeric")
      }
      .Object@elitism <- elitism
    } else {
      .Object@elitism <- base::max(1, round(.Object@popSize*0.05));
    }
    
    if (!missing(seed)) {
      .Object@seed <- seed
    }
    
    .Object@summary <- matrix(as.double(NA), nrow = .Object@maxiter, ncol = 6)
    colnames(.Object@summary) <- names(my_ga.summary(.Object, x=rnorm(10)))
    
    return(my_ga.run(.Object));
  }
)

setGeneric(
  "my_ga.run", 
  function(.Object, ...) standardGeneric("my_ga.run")
)
setMethod(
  "my_ga.run",
  signature = "my_ga",
  definition = function(.Object, ...){
    
    if (!is.null(.Object@seed)) {
      set.seed(.Object@seed);
    }
    
    # initial population
    .Object@population <- my_ga.generatePopulation(.Object)
    
    # start iterations
    for (iter in seq_len(.Object@maxiter)) {
      # evalute fitness function (when needed) 
      for (i in seq_len(.Object@popSize)) {
        if (is.na(.Object@fitnessValues[i])) { 
          .Object@fitnessValues[i] <- .Object@fitness(.Object@population[i,], ...) 
        }
      }
      
      .Object@summary[iter,] <- my_ga.summary(.Object);
      
      .Object@iter = iter;
      
#       if(.Object@keepBest) {
#         .Object@bestSolution[[iter]] <- unique(.Object@population[Fitness == fitnessSummary[iter,1],,drop=FALSE])
#       }
      
      if (is.function(.Object@monitor)) {
        .Object@monitor(.Object)
      }
      
      # check stopping criteria
      if(iter > 1) { 
        if (.Object@summary[iter,1] > .Object@summary[iter-1,1] + .Object@eps) {
          .Object@run <- 1 
        } else {
          .Object@run <- .Object@run + 1 
        }
      }
      if (.Object@run >= .Object@maxRun) {
        break  
      }
      if (max(.Object@fitnessValues, na.rm = TRUE) >= .Object@maxfitness) {
        break
      }
      if (.Object@iter == .Object@maxiter) { 
        break 
      }

      # make popultaion copy for Еlitist selection
      ord <- order(.Object@fitnessValues, decreasing = TRUE)
      PopSorted <- .Object@population[ord,,drop=FALSE]
      FitnessSorted <- .Object@fitnessValues[ord]
      
      # selection
      sel <- my_ga.selection(.Object)
      .Object@population <- sel$population
      .Object@fitnessValues <- sel$fitnessValues
      
      # crossover
      if (.Object@pcrossover > 0) { 
        nmating <- floor(.Object@popSize/2)
        mating <- matrix(sample(1:(2*nmating), size = (2*nmating)), ncol = 2)
        for (i in seq_len(nmating)) {
          if (.Object@pcrossover > runif(1)) {
            parents <- mating[i,]
            Crossover <- my_ga.crossover(.Object, parents)
            .Object@population[parents,] <- Crossover$children
            .Object@fitnessValues[parents] <- Crossover$fitness
          }
        }
      }
      
      # mutation
      if (.Object@pmutation > 0) { 
        for (i in seq_len(.Object@popSize)) { 
          if (.Object@pmutation > runif(1)) { 
            Mutation <- my_ga.mutation(.Object, i)
            .Object@population[i,] <- Mutation
            .Object@fitnessValues[i] <- NA
          }
        }
      }

      # Еlitist selection is strategy that guarantees that the solution quality 
      # obtained by the GA will not decrease from one generation to the next.
      if(.Object@elitism > 0) { 
        ord <- order(.Object@fitnessValues, na.last = TRUE)
        u <- which(!duplicated(PopSorted, margin = 1))
        .Object@population[ord[1:.Object@elitism],] <- PopSorted[u[1:.Object@elitism],]
        .Object@fitnessValues[ord[1:.Object@elitism]] <- FitnessSorted[u[1:.Object@elitism]]
      }
    }

    # in case of premature convergence remove NA from summary fitness evalutations
#     .Object@summary <- na.exclude(.Object@summary)
#     attr(.Object@summary, "na.action") <- NULL
    
    # get solution(s)
    .Object@fitnessValue <- max(.Object@fitnessValues, na.rm = TRUE)
    valueAt <- which(.Object@fitnessValues == .Object@fitnessValue)
    solution <- .Object@population[valueAt,,drop=FALSE]
    if (nrow(solution) > 1) { 
      # find unique solutions to precision given by default tolerance
      eps <- .Object@eps
      solution <- unique(round(solution/eps)*eps, margin = 1)
    }
    #colnames(solution) <- parNames(object)
    .Object@solution <- solution

#     if(.Object@keepBest)
#       .Object@bestSolution <- .Object@bestSolution[!sapply(.Object@bestSolution, is.null)]  
    
    return(.Object)
  }
);

setGeneric(
  "my_ga.generatePopulation", 
  function(.Object, ...) standardGeneric("my_ga.generatePopulation")
)
setMethod(
  "my_ga.generatePopulation",
  signature = "my_ga",
  definition = function(.Object, ...){
    # Generate a random population of size popSize in the range [min, max]  
    min <- .Object@min
    max <- .Object@max
    nvars <- length(min)
    population <- matrix(as.double(NA), nrow = .Object@popSize, ncol = nvars)
    for (j in 1:nvars) {
      population[,j] <- runif(.Object@popSize, min[j], max[j]) 
    }
    return(population)
  }
);

# setGeneric(
#   "my_ga.fitnessEvaluation", 
#   function(.Object) standardGeneric("my_ga.fitnessEvaluation")
# )
# setMethod(
#   "my_ga.fitnessEvaluation",
#   signature = "my_ga",
#   definition = function(.Object){
#   }
# );
# 
# setGeneric(
#   "my_ga.convergenceCheck", 
#   function(.Object) standardGeneric("my_ga.convergenceCheck")
# )
# setMethod(
#   "my_ga.convergenceCheck",
#   signature = "my_ga",
#   definition = function(.Object){
#   }
# );

setGeneric(
  "my_ga.selection", 
  function(.Object, ...) standardGeneric("my_ga.selection")
)
setMethod(
  "my_ga.selection",
  signature = "my_ga",
  definition = function(.Object, ...){
    # Fitness proportional selection with fitness linear scaling  
    f <- .Object@fitnessValues
    fmin <- min(f, na.rm = TRUE)
    if (fmin < 0) { # will fix all negative numbers
      f <- f - fmin
      fmin <- min(f, na.rm = TRUE) 
    }
    fave <- mean(f, na.rm = TRUE)
    fmax <- max(f, na.rm = TRUE)
    sfactor <- 2 # scaling factor
    # transform f -> f' = a*f + b such that
    if (fmin > (sfactor * fave - fmax) / (sfactor - 1)) {
      # ave(f) = ave(f')
      # 2*ave(f') = max(f')
      delta <- fmax - fave
      a <- (sfactor - 1.0)*fave/delta
      b <- fave * (fmax - sfactor*fave)/delta 
    } else { 
      # ave(f) = ave(f')
      # min(f') = 0
      delta <- fave - fmin
      a <- fave/delta
      b <- -1*fmin*fave/delta 
    }
    fscaled <- a*f + b
    prob <- abs(fscaled)/sum(abs(fscaled))
    sel <- sample(1:.Object@popSize, size = .Object@popSize, 
                  prob = pmin(pmax(0, prob), 1, na.rm = TRUE), 
                  replace = TRUE)
    out <- list(population = .Object@population[sel,,drop=FALSE],
                fitnessValues = .Object@fitnessValues[sel])
    return(out)
  }
);

setGeneric(
  "my_ga.crossover", 
  function(.Object, parents, ...) standardGeneric("my_ga.crossover")
)
setMethod(
  "my_ga.crossover",
  signature = "my_ga",
  definition = function (.Object, parents, ...) {
    # Local arithmetic crossover
    # This crossover alghoritm choose number that are somewhere between two parent 
    # according to choosen random weigth
    parents <- .Object@population[parents,,drop = FALSE]
    n <- ncol(parents)
    children <- matrix(as.double(NA), nrow = 2, ncol = n)
    a <- runif(n) # weight factor
    children[1,] <- a*parents[1,] + (1-a)*parents[2,]
    children[2,] <- a*parents[2,] + (1-a)*parents[1,]
    out <- list(children = children, fitness = rep(NA,2))
    return(out)
  }
);


setGeneric(
  "my_ga.mutation", 
  function(.Object, parent, ...) standardGeneric("my_ga.mutation")
)
setMethod(
  "my_ga.mutation",
  signature = "my_ga",
  definition = function (.Object, parent, ...) {
    # Uniform random mutation. 
    # Replace specified population value with random number
    mutate <- parent <- as.vector(.Object@population[parent,])
    n <- length(parent)
    j <- sample(1:n, size = 1)
    mutate[j] <- runif(1, .Object@min[j], .Object@max[j])
    return(mutate)
  }
);

setGeneric(
  "my_ga.summary", 
  function(.Object, x = .Object@fitnessValues, ...) standardGeneric("my_ga.summary")
)
setMethod(
  "my_ga.summary",
  signature = "my_ga",
  definition = function (.Object, x = .Object@fitnessValues, ...) {
    # compute summary for each step
    x <- na.exclude(as.vector(x))
    q <- fivenum(x)
    c(max = q[5], mean = mean(x), q3 = q[4], median = q[3], q1 = q[2], min = q[1])
  }
);

setMethod("plot", signature="my_ga", plot.ga)

summary.my_ga = function (object, digits = getOption("digits"), ...) {    
  dotargs <- list(...)
  if(is.null(dotargs$head)) dotargs$head <- 10
  if(is.null(dotargs$tail)) dotargs$tail <- 1
  if(is.null(dotargs$chead)) dotargs$chead <- 20
  if(is.null(dotargs$ctail)) dotargs$ctail <- 1
  
  cat("+-----------------------------------+\n")
  cat("|         Genetic Algorithm         |\n")
  cat("+-----------------------------------+\n\n")
  cat("GA settings: \n")
  cat(paste("Type                  = ", object@type, "\n"))
  cat(paste("Population size       = ", object@popSize, "\n"))
  cat(paste("Number of generations = ", object@maxiter, "\n"))
  cat(paste("Elitism               = ", object@elitism, "\n"))
  cat(paste("Crossover probability = ", format(object@pcrossover, digits = digits), "\n"))
  cat(paste("Mutation probability  = ", format(object@pmutation, digits = digits), "\n"))
  
  if (object@type == "real-valued") { 
    cat(paste("Search domain \n"))
    domain <- rbind(object@min, object@max)
    rownames(domain) <- c("Min", "Max")
    if(ncol(domain) == length(object@min)) {
      nvars <- ncol(object@population)
      varnames <- paste("x", 1:nvars, sep = "")
      colnames(domain) <- varnames
    }
    print(domain, digits = digits)
  }
  
  cat("\nGA results: \n")
  cat(paste("Iterations             =", format(object@iter, digits = digits), "\n"))
  cat(paste("Fitness function value =", format(object@fitnessValue, digits = digits), "\n"))
  if (nrow(object@solution) > 1) { 
    cat(paste("Solutions              = \n")) 
  } else { 
    cat(paste("Solution               = \n")) 
  }
  print(object@solution)
  
  invisible()
}
setMethod(
  "summary", 
  signature="my_ga", 
  summary.my_ga
)