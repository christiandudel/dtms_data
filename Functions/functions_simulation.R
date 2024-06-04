### For simulations ############################################################
  
  # Helper functions to simplify data with duration dependence
  # General idea is: STATENAME%D[sep]%T
  # Where STATENAME is the name of the state, e.g., 'healthy'
  # %D is the duration, so %D=1,2,...
  # %T is the time, so %T=0,1,2,...
  # [sep] is the separator between duration %D and time %T
  
  # Lowest level: 
  # Removes numbers from string and thus turns
  # STATENAME%D into STATENAME
  paster <- function(x) {
    res <- str_remove_all(x[1],"[:digit:]")
    return(res)
  }
  
  # Splits [sep]%T from STATENAME%D[sep]%T, 
  # which gives STATENAME%D and passes it to the lowest level 
  # function 'paster'
  simplify <- function(data,sep="_") {
    res <- lapply(str_split(data,pattern=sep),paster)
    res <- unlist(res)
    return(res)
  }
  
  # Applies above functions to all variables in a data.frame
  # This turns a duration-dependent state space STATENAME%D[sep]%T
  # into a duration-independent state space STATENAME[sep]%T
  simplifydata <- function(data,sep="_") {
    res <- apply(data,2,simplify,sep=sep)
    res <- as.data.frame(res)
    return(res)
  }
  
  
### Generate trans probs in data frame, basic ##################################
  
  basic_sim <- function(transient,  # Character vector of transient states
                        absorbing,  # Character vector of absorbing states
                        time_steps, # Time steps of the process
                        probs){     # Transition probabilities as named list
    
    # Packages
    require(tidyverse)
    
    # Get full state space
    state_space <- c(transient,absorbing)
    
    # Number of states
    n_transient <- length(transient)
    
    # Create data frame
    model_frame <- expand.grid(from=state_space,
                               to=state_space,
                               time=time_steps,
                               stringsAsFactors=F)
    
    # Drop absorbing states as starting states
    model_frame <- model_frame %>% filter(!from%in%absorbing)
    
    # Transition probabilities empty
    model_frame$P <- NA
    
    # Place probabilities (no generated duration dependence)
    for(i in transient) {
      for(j in state_space) {
        model_frame <- model_frame %>% mutate(P=ifelse(from==i & to==j,probs[[i]][j],P))
      }
    }
 
    # Everybody dies
    model_frame <- model_frame %>% mutate(
      P=ifelse(time==max(time_steps) & to==absorbing[1],1,P),
      P=ifelse(time==max(time_steps) & to!=absorbing[1],0,P) 
    )
    
    # Rename values
    model_frame <- model_frame %>% mutate(
      from=ifelse(from%in%transient,paste(from,time,sep="_"),from),
      to  =ifelse(to%in%transient,paste(to,time+1,sep="_"),to)
    )
    
    # Return
    return(model_frame)
    
  }
  
  # # Test:
  # model_test <- basic_sim(time_steps = 0:10,
  #                            transient  = c("A","B"),
  #                            absorbing  = "X",
  #                            probs      = list(A=c(A=0.495,B=0.495,X=0.01),
  #                                              B=c(A=0.495,B=0.495,X=0.01)))
  
  
### Generate transition probabilities in data frame ############################

  # Things to be aware of:
  # - Changes applied by age- or duration-dependence will not necessarily be applied 1-to-1, as probabilities have to sum to 1
  # - After maximum duration is reached, there is a forced transition to either another state or death, with rescaling such that probabilities sum to 1

  # Function
  generate_sim <- function(sim){ 
                           
    # Packages
    require(tidyverse)
    
    # Assign
    transient <- sim$transient
    absorbing <- sim$absorbing
    time_steps <- sim$time_steps
    probs <- sim$probs
    gen_duration <- sim$gen_duration
    gen_age <- sim$gen_age
    which_duration <- sim$which_duration
    diff_duration <- sim$diff_duration
    interpolation_duration <- sim$interpolation_duration
    which_age <- sim$which_age
    diff_age <- sim$diff_age
    interpolation_age <- sim$interpolation_age
    
    # Get full state space
    if(gen_duration) {
      dependent_states <- levels(interaction(which_duration,1:(max(time_steps)+1),sep=""))
      independent_states <- transient[!(transient%in%which_duration)]
      state_space_full <- c(dependent_states,independent_states,absorbing)
      state_space_red <- c(transient,absorbing)
    } else state_space_full <- state_space_red <- c(transient,absorbing)
    
    # Number of states
    n_transient <- length(transient)
    
    # Create data frame
    model_frame <- expand.grid(from=state_space_full,
                               to=state_space_full,
                               time=time_steps,
                               stringsAsFactors=F)
    
    # Drop absorbing states as starting states
    model_frame <- model_frame %>% filter(!from%in%absorbing)
    
    # Transition probabilities empty
    model_frame$P <- 0
    
    # Duration dependence: select states which are not duration dependent
    take <- rep(T,n_transient)
    if(gen_duration) {
      donttake <- which(transient%in%which_duration)
      take[donttake] <- F
    } 
    
    # Place probabilities (no generated duration dependence)
    for(i in transient[take]) {
      for(j in state_space_red) {
        p <- probs[[i]][j]
        if(j %in% which_duration) {
          j_to <- paste0(j,"1")
        } else j_to <- j
        model_frame <- model_frame %>% mutate(P=ifelse(from==i & to==j_to,p,P))
      }
    }
    
    # Generated duration dependence
    for(i in which_duration) {
      for(j in state_space_red) {
        
        # Get duration states
        i_states <- paste0(i,1:(max(time_steps)))
        if(i==j) it_states <- paste0(i,2:(max(time_steps)+1))
        
        # Get starting probability
        first <- probs[[i]][j]
        
        # Is j also duration dependent? Not used if i==j
        if(j %in% which_duration & j!=i) {
          j_to <- paste0(j,"1")
        } else j_to <- j
        
        # Check if any dependence for transition to state j, if no dependence move on to next j
        difference <- diff_duration[[i]][j]
        if(difference==0) {
          model_frame <- model_frame %>% mutate(P=ifelse(from%in%i_states & to==j_to,first,P))
          next
        } 

        # Interpolate
        nsteps <- length(i_states)-1
        
        # Linear interpolation
        if(str_detect(interpolation_duration[[i]],"linear")) {
          values <- seq(first,first+difference,length.out=nsteps+1)
        }
        
        # Sigmoid interpolation
        if(str_detect(interpolation_duration[[i]],"sigmoid")) {
          values <- seq(-15,15,
                        length.out=nsteps)
          values <- plogis(values)*difference
          values <- values+first
          values <- c(first,values) 
        }
        
        # Random values
        if(str_detect(interpolation_duration[[i]],"random")) {
          low <- min(first,first+difference)
          high <- max(first,first+difference)
          values <- c(first,
                      runif(min=low,
                            max=high,
                            n=nsteps))
        }
        
        # Switching regime
        if(str_detect(interpolation_duration[[i]],"switch")) {
          where <- str_sub(interpolation_duration[[i]],start=7)
          where <- as.numeric(where)
          values <- c(rep(first,where),
                      rep(first+difference,
                          (nsteps+1)-where))
        }

        # Place values
        for(time in 1:length(values)) {
          i_from <- i_states[time]
          if(j==i) {
            j_to <- it_states[time]} else if(j %in% which_duration & j!=i) {
            j_to <- paste0(j,"1")
          } else j_to <- j
          p <- values[time]
          model_frame <- model_frame %>% mutate(P=ifelse(from==i_from & to==j_to,p,P))
        }
        
        # Last state only dies
        last_state <- last(it_states)
        model_frame <- model_frame %>% mutate(P=ifelse(from==last_state & to=="X",1,P),
                                              P=ifelse(from==last_state & to!="X",0,P))
        
        
      }
    }
    
    # Everybody dies
    model_frame <- model_frame %>% mutate(
      P=ifelse(time==max(time_steps) & to==absorbing[1],1,P),
      P=ifelse(time==max(time_steps) & to!=absorbing[1],0,P) 
    )
  
    # # Age effects
    if(gen_age) {
      
      # Check all age-dependent states
      for(i in which_age) {
        
        # All receiving states
        for(j in state_space_red) {
          
          # Get starting probability
          first <- probs[[i]][j]
          
          # Check if any dependence for transition to state j, if no age-dependence move on to next j
          difference <- diff_age[[i]][j]
          if(difference==0) {
            next
          }
          
          # Interpolate between all values
          nsteps <- length(time_steps)-1
          
          # Linear interpolation
          if(str_detect(interpolation_age[[i]],"linear")) {
            values <- seq(0,difference,length.out=nsteps)
          }
          
          # Sigmoid interpolation
          if(str_detect(interpolation_age[[i]],"sigmoid")) {
            values <- seq(-15,15,
                          length.out=nsteps-1)
            values <- plogis(values)*difference
            values <- values
            values <- c(0,values) 
          }
          
          # Simple states
          model_frame$simple_from <- str_match(model_frame$from,i)
          model_frame$simple_to <- str_match(model_frame$to,j)
          
          # Loop over time
          lasttime <- max(time_steps)-1
          for(timeloop in 0:lasttime) {
            
            if(!j%in% which_duration) {
              i_from <- i
              j_to <- j
              model_frame <- model_frame %>% mutate(P=ifelse(simple_from%in%i_from & 
                                                               simple_to%in%j_to &
                                                               time==timeloop,P+values[timeloop+1],P))
            }
            
            if(j %in% which_duration & i!=j) {
              j_to <- paste0(j,1) 
              model_frame <- model_frame %>% mutate(P=ifelse(simple_from%in%i_from & 
                                                           to%in%j_to &
                                                           time==timeloop,P+values[timeloop+1],P))
            }
            
            if(j %in% which_duration & i==j) {
              
              for(duration in 1:(max(time_steps)-1)) {
              i_from <- paste0(i,duration)
              j_to <- paste0(j,duration+1)
              model_frame <- model_frame %>% mutate(P=ifelse(from%in%i_from & 
                                                               to%in%j_to &
                                                               time==timeloop,P+values[timeloop+1],P))
              }
            }
          }
          
          # Remove simple states
          model_frame <- model_frame %>% select(!c("simple_from","simple_to"))
        }
      }
      
    }
    
    # No values below 0 (can happen of you go wild with values)
    if(any(model_frame$P<0)) {
      warning("Created transition probabilities below 0, readjusting")
      model_frame$P[model_frame$P<0] <- 0
    }
    
    # Rescaling 
    model_frame <- model_frame %>% group_by(from,time) %>% mutate(P=P/sum(P)) %>% ungroup
    class(model_frame) <- "data.frame"

    # Rename values
    if(gen_duration) changestates <- c(dependent_states,independent_states) else changestates <- transient
    model_frame <- model_frame %>% mutate(
      from=ifelse(from%in%changestates,paste(from,time,sep="_"),from),
      to  =ifelse(to%in%changestates,paste(to,time+1,sep="_"),to)
    )
    
    # Return
    return(model_frame)
    
  }
  
  # Old test cases
  # # Test case 1
  # model_test <- generate_sim(time_steps = 0:10,
  #                            transient  = c("A","B"),
  #                            absorbing  = "X",
  #                            probs      = list(A=c(A=0.495,B=0.495,X=0.01),
  #                                              B=c(A=0.495,B=0.495,X=0.01)))
  # 
  # # Test case 2
  # model_test <- generate_sim(time_steps = 0:10,
  #               transient  = c("A","B"),
  #               absorbing  = "X",
  #               probs      = list(A=c(A=0.495,B=0.495,X=0.01),
  #                                B=c(A=0.495,B=0.495,X=0.01)),
  #               gen_duration = T,
  #               diff_duration = list(A=c(A=-0.1,B=0,X=0.1)),
  #               which_duration = "A",
  #               interpolation_duration = list(A="random"))
  # test <- paste(paste0("A",0:10),"0",sep="_")
  # test2 <- paste(paste0("A",0:10),"1",sep="_")
  # model_test %>% filter(from%in%test & to%in%test2)
  # model_test %>% filter(from=="A10_10")   
  # 
  # # Test case 3
  # model_test <- generate_sim(time_steps = 0:100,
  #                            transient  = c("A","B"),
  #                            absorbing  = "X",
  #                            probs      = list(A=c(A=0.1,B=0.89,X=0.01),
  #                                              B=c(A=0.495,B=0.495,X=0.01)),
  #                            gen_duration = T,
  #                            which_duration = "A",
  #                            diff_duration = list(A=c(A=0.79,B=-0.79,X=0.1)),
  #                            interpolation_duration = list(A="linear"))
  # model_test %>% filter(from=="A100_0")   
  # 
  # # Test case 3
  # model_test <- generate_sim(time_steps = 0:100,
  #                            transient  = c("A","B"),
  #                            absorbing  = "X",
  #                            probs      = list(A=c(A=0.1,B=0.89,X=0.01),
  #                                              B=c(A=0.495,B=0.495,X=0.01)),
  #                            gen_duration = T,
  #                            which_duration = "A",
  #                            diff_duration = list(A=c(A=0.79,B=-0.79,X=0.1)),
  #                            interpolation_duration = list(A="linear"),
  #                            gen_age = T,
  #                            which_age = c("A","B"),
  #                            diff_age = list(A=c(A=0,B=0,X=0.1),
  #                                            B=c(A=0,B=0,X=0.2)),
  #                            interpolation_age = list(A="linear",B="linear"))
  # 
  # model_test %>% filter(from=="A1_0" & to=="B_1")   
  # model_test %>% filter(from=="A100_0" & to=="B_1")   
  # 

  
### Complete sim objects if not all entries ####################################
  
  complete_sim <- function(sim) {
    
    simnames <- names(sim)
    allnames <- c("transient","absorbing","time_steps","probs","gen_duration",
                  "which_duration","diff_duration","interpolation_duration",
                  "gen_age","which_age","diff_age","interpolation_age")
    which_missing <- allnames[!allnames%in%simnames]
    n_missing <- length(which_missing)
    add <- vector("list",n_missing)
    names(add) <- which_missing
    sim <- c(sim,add)
    return(sim)
    
  }
  
### Functions to be applied on simulated data ##################################
  
  # Absorbing states only once
  drop_dead <- function(data,absorbing){
    result <- apply(data,1,
                    FUN=na_after,
                    absorbing=absorbing,
                    simplify=T)  
    result <- as.data.frame(t(result))
    return(result)
  }

  # Helper for previous  
  na_after <- function(x,absorbing) {
    check <- which(x==absorbing)
    if(length(check)>1) x[check[-1]] <- NA
    return(x)
  }
  
  # Missings in between/gaps
  gaps <- function(data,prob=0.005) {
    result <- apply(data,1,
                    FUN=na_gap,
                    prob=prob,
                    simplify=T)  
    result <- as.data.frame(t(result))
    return(result)
  }
  
  # Helper for previous
  na_gap <- function(x,prob) {
    nlength <- length(x)
    drop <- runif(n=nlength)
    drop <- which(drop<prob)
    x[drop] <- NA
    return(x)
  }
  
  # Left/right censoring (uninformative)
  censoring <- function(data,probleft=0.5,probright=0.5) {
    result <- apply(data,1,
                    FUN=na_cens,
                    probleft=probleft,
                    probright=probright,
                    simplify=T)  
    result <- as.data.frame(t(result))
    return(result)
  }
  
  # Function for above
  na_cens <- function(x,probleft,probright) {
    nlength <- length(x)
    leftdrop <- runif(1)<probleft
    rightdrop <- runif(1)<probright
    places <- sort(sample(1:nlength,2))
    
    if(leftdrop&!rightdrop) drop <- 1:places[1]
    if(!leftdrop&rightdrop) drop <- places[2]:nlength
    if(leftdrop&rightdrop) drop <- c(1:places[1],places[2]:nlength)
    if(leftdrop|rightdrop) x[drop] <- NA
    return(x)
  }
  

### Function for variance ######################################################
  
  sim_var <- function(matrix) {
    
    # Get values
    values <- as.numeric(colnames(matrix))
    
    # Drop NAs (e.g., TOTAL column)
    matrix <- matrix[,!is.na(values)]
    values <- na.omit(values)
    
    # Calculate
    result <- apply(matrix,1,
                    function(x) sum( (x*(values-weighted.mean(x=values,w=x))^2)/sum(x) ))
    
    # Return
    return(result)
    
  }