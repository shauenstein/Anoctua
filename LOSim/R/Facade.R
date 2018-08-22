#' @useDynLib LOSim
#' @importFrom Rcpp evalCpp
#' @export
runSimulation <- function(environment, iterations, parameters, randomSeed = NULL, 
                          output = "raw", #  "raw", "observed", "rawplus")
                          TimeAtObservation = NULL, maxDuration = NULL,
                          upperleft = NULL, environmentResolution = 1,
                          startTime, startDayOfYear, 
                          solarLatitude, solarLongitude){
  
  # start timing 
  ptm <- proc.time()
  
  # translate output to integers
  if(output == "raw"){
    output <- 1
  }else{
    if(output == "observed"){
      output <- 2
    }else{
      if(output == "rawplus"){
        stop("Both simulated and observed output is currently not implemented.")
        output <- 3
      }else{
        stop("output must be defined as one of the following: raw, observed, rawplus")
      }
    }
  }
  
  # check settings and input
  if(output != 1 & is.null(TimeAtObservation)){
    stop("To obtain simulated observations please provide observation times")
  }
  # TimeAtObservation cannot be NULL
  if(is.null(TimeAtObservation)) TimeAtObservation <- c(-999)
  
  # upperleft cannot be NULL
  if(is.null(upperleft)) stop("Please provide the dimensions of the environmental input.")
  
  # check if parameters is a matrix with 10 columns (11 if ABCcalibration)
  if(output != 1 & NCOL(parameters) != 12){
    stop("make sure to provide all 12 parameters if you want to simulate observed dispersal.")
  }
  if(output == 1 & NCOL(parameters) != 10){
    stop("make sure to provide all 10 parameters if you simply simulate dispersal.")
  }
  
  # translate spatial coordinates
  parameters[,1] <- (parameters[,1] - upperleft[1]) / environmentResolution
  parameters[,2] <- -(parameters[,2] - upperleft[2]) / environmentResolution
  # translate observation error
  if(output != 1) parameters[,11] <- parameters[,11] / environmentResolution
  if(output != 1) parameters[,12] <- parameters[,12] / environmentResolution 
  # translate maximum effort
  parameters[,9] <- parameters[9] / environmentResolution
  # translate maxDuration
  if(is.null(maxDuration) | output %in% c(2,3)){
    maxDuration <- -999
  }else{
    maxDuration <- maxDuration * 24 * 60
  }
  
  # environment values, read from file or matrix
  if(class(environment) == "character"){
    filename <- environment
    environment <- matrix(0, 0, 0) # 0x0 matrix
    if(!file.exists(filename)){
      stop("Provide correct path to environment input.")
    }
  }else{
    filename <- "empty"
    if(class(environment) != "matrix"){
      stop("Provide environment as character (read from file) or as numeric matrix.")
    } 
  }
  
  
  # set random seed
  if(is.null(randomSeed)) randomSeed <- sample(10000, size = 1)
  
  # RCPP function call
  out <- callModel(environment, iterations, maxDuration, parameters, 
                   randomSeed, output, TimeAtObservation, 
                   upperleft, environmentResolution, filename,
                   startTime, startDayOfYear, 
                   solarLatitude, solarLongitude)
  
  # stop timing
  print(proc.time() - ptm)
  
  return(out)
  
}