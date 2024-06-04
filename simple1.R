### Generating simple example data #############################################

### Packages, functions, settings ##############################################

  # Seed
  set.seed(8219)

  # Packages
  library(markovchain)
  library(tidyverse)
  library(dtms)
  
  # Functions
  source("Functions/functions_simulation.R")


### Model setup ################################################################

  model <- list(time_steps = 0:19,
                transient  = c("A","B"),
                absorbing  = "X",
                probs      = list(A=c(A=0.1,B=0.89,X=0.01),
                                  B=c(A=0.495,B=0.495,X=0.01)),
                gen_duration = T,
                which_duration = "A",
                diff_duration = list(A=c(A=0.79,B=-0.79,X=0.1)),
                interpolation_duration = list(A="linear"),
                gen_age = T,
                which_age = c("A","B"),
                diff_age = list(A=c(A=0,B=0,X=0.1),
                B=c(A=0,B=0,X=0.2)),
                interpolation_age = list(A="sigmoid",B="sigmoid"),
                sample_size=1000,
                replications=1,
                initial_distr=c(0.5,0.5))
  
### Generate ###################################################################
  
  # Generate transition probabilities
  simmodel <- generate_sim(model)
  
  # Transient states
  transient_states <- levels(interaction(model$which_duration,1:(max(model$time_steps)+1),sep=""))
  transient_states <- c(transient_states,model$transient[!model$transient%in%model$which_duration])
  
  # Model
  simdtms <- dtms( transient = transient_states,
                   absorbing = model$absorbing,
                   timescale = model$time_steps)
  
  # Settings
  sample_size <- model$sample_size
  nlength <- length(model$time_steps)
  
  # Transition matrix
  Tm <- dtms_matrix(dtms=simdtms,
    probs=simmodel
  )
  
  class(Tm) <- "matrix"
  
  # Setup markovchain package
  mcsim <- new("markovchain", 
               states = rownames(Tm),
               transitionMatrix = as.matrix(Tm),
               name = "mcsim")
  simdata <- data.frame(matrix(nrow=0,ncol=nlength))
  
  # Get names of starting states right i
  starting_states <- paste0(model$which_duration,"1_0")
  other_states <- model$transient[!model$transient%in%model$which_duration]
  starting_states <- c(starting_states,paste0(other_states,"_0"))

  # Generate individual trajectories
  for(i in 1:sample_size) {
    
    # Get initial state
    initial_state <- sample(starting_states,1,prob=model$initial_distr)
    simseq <- rmarkovchain(n = nlength-1, object = mcsim, t0 = initial_state,include.t0=T)  
    names(simseq) <- names(simdata)
    simdata[i,] <- simseq
    
  }  
  
  # Simplify
  simdata <- simplifydata(simdata)

  # Rename variables
  namesdata <- paste("T",0:(dim(simdata)[2]-1),sep="_")
  names(simdata) <- namesdata
  
  # Censor
  simdata <- censoring(data=simdata,probleft=0.25,probright=0.25)
  
  # Add a few gaps
  simdata <- gaps(data=simdata,prob=0.025)
  
  # Drop dead obs
  simdata <- drop_dead(data=simdata,absorbing="X")
  
  # Add IDs
  simdata$ID <- 1:sample_size
  
  # Reshape
  simdata <- simdata %>% pivot_longer(cols=starts_with("T_"),
                                      names_prefix="T_",
                                      names_to="time",
                                      values_to="from")
  
  # Change time variable to numeric
  class(simdata$time) <- "numeric"
  
  # Drop missing values
  simdata <- na.omit(simdata)

  # Rename
  names(simdata) <- c("id","time","state")

  # Save (1)
  write.csv(simdata,
            file="Data/simple_example.csv",row.names=F)
  
  # Save (2)
  save(simdata,file="Data/simple.rda")
  