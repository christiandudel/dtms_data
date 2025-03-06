### Generate work-trajectory example ###########################################

### Packages, functions, settings ##############################################

  # Seed
  set.seed(8219)
  
  # Packages
  library(markovchain)
  library(tidyverse)
  
  # Functions
  source("Functions/functions_simulation.R")


### Load transition matrices ###################################################

  # Women
  women <- read.csv(file="Input/Pmat_f_2009.csv",row.names=1)
  women <- as.matrix(women)
  states <- rownames(women)
  colnames(women) <- states
  
  # Men
  men <- read.csv(file="Input/Pmat_m_2009.csv",row.names=1)
  men <- as.matrix(men)
  colnames(men) <- states
  
  # Row-orientation, not column orientation
  men <- t(men)
  women <- t(women)
  
  
### Simulate data ##############################################################
  
  # Setup markovchain package
  mensim <- new("markovchain", 
               states = states,
               transitionMatrix = men,
               name = "mensim")
  
  womensim <- new("markovchain", 
                states = states,
                transitionMatrix = women,
                name = "womensim")
  
  # Further sum setup
  sample_size <- 5000
  nlength <- length(50:99)
  
  # Empty data frame
  simdata <- data.frame(matrix(nrow=0,ncol=nlength-1))
  gendervar <- numeric(0)
  
  # Get names of starting states right i
  starting_states <- c("50::Employed","50::Inactive","50::Retired")

  # Generate individual trajectories
  for(i in 1:sample_size) {
    
    # Gender
    gender <- sample(0:1,1)
    
    # Get initial state
    if(gender==0) initial_state <- sample(starting_states,1,prob=c(0.8,0.15,0.05))
    if(gender==1) initial_state <- sample(starting_states,1,prob=c(0.65,0.30,0.05))
    
    # Simulate trajectory
    if(gender==0) simseq <- rmarkovchain(n = nlength-1, object = mensim, t0 = initial_state,include.t0=T)  
    if(gender==1) simseq <- rmarkovchain(n = nlength-1, object = womensim, t0 = initial_state,include.t0=T)
    
    # Put gender and trajectory in object
    gendervar <- c(gendervar,gender)
    simdata <- rbind(simdata,simseq)
    
  }  
  
### Edit simulated data ########################################################
  
  # Variable names
  names(simdata) <- paste("T",50:99,sep="_")

  # Function to change state names
  simplify_other <- function(x,sep="_") {
    res <- lapply(x,function(y) str_remove_all(string=paste(y,collapse=""),pattern="[:digit:][:digit:]::"))  
    res <- unlist(res)
    return(res)
  }
  
  # Change state names
  simdata <- apply(simdata,2,function(x) simplify_other(x,sep="::"))
  
  # Back to data frame
  simdata <- as.data.frame(simdata)
  
  # Censor
  simdata <- censoring(data=simdata)
  
  # Add a few gaps
  simdata <- gaps(data=simdata,prob=0.025)
  
  # Drop dead obs
  # simdata <- drop_dead(data=simdata,absorbing="X")
  
  # Add IDs
  simdata$ID <- 1:sample_size  
  
  # Add Gender
  simdata$gender <- gendervar

  # Reshape
  simdata <- simdata %>% pivot_longer(cols=starts_with("T_"),
                                      names_prefix="T_",
                                      names_to="age",
                                      values_to="state")
  
  # Change time variable to numeric
  class(simdata$age) <- "numeric"

  # Rename
  names(simdata) <- c("ID","Gender","Age","State")
  
  # Change state names
  simdata$State[simdata$State%in%"Employed"] <- "Working"
  simdata$State[simdata$State%in%"Inactive"] <- "Non-working"
  
### Save #######################################################################
  
  # Save
  write.csv(simdata,
            file="Output/work.csv",row.names=F)  
  
  # Save (2)
  save(simdata,file="Output/work.rda")
  