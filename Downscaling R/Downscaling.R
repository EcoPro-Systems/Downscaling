Downscaling=function(path,freq,BGL){
  # Check if the package is installed; if not, install it
  # Define the package name you want to load
  package_name <- "data.table"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "akima"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "parallel"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "doParallel"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "foreach"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "pryr"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  ##################################################################################################################################################
  # Load data from CSV files
  Model=as.matrix(data.table::fread(paste0(path,"/Model_data.csv")))
  Obs=as.matrix(data.table::fread(paste0(path,"/Obs_data.csv")))
  
  N=nrow(Obs)
  T_o=ncol(Obs[,-c(1:2)])
  T_all=ncol(Model[,-c(1:2)])
  
  # Source a custom function 'Stand_downscaling' from 'Stand_downscaling.R'
  source(paste0(path,"/Stand_downscaling.R"))
  Stand_SD=Stand_downscaling(Obs,Model,freq)
  Stand_SD=cbind(Obs[,1:2],Stand_SD)
  
  # Write the results to CSV files
  fwrite(Stand_SD[,-(c(1:T_o)+2)],paste0(path,"Stand_Downscaled.csv"))
  fwrite(Stand_SD,paste0(path,"Stand_Downscaled_training.csv"))
  
  
  if(BGL=="on"){
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    # Initialize matrices to store results
    Mu1=matrix(NA, nrow = N, ncol = T_all)
    Mu2=matrix(NA, nrow = N, ncol = T_all)
    
    # Create time vectors 't' and 't_new'
    t=rep(1:freq, times = T_o/freq)
    t_new=rep(1:freq, times = T_all/freq)
    
    # Calculate means 'X' and 'X_new'
    X=colMeans(Model[, c(1:T_o)+2])
    X_new=colMeans(Model[,-c(1:2)])
    
    # Create data frames
    data=data.frame(cbind(X,t))
    names(data)=c("X","t")
    
    data_new=data.frame(cbind(X_new,t_new))
    names(data_new)=c("X","t")
    
    cat("Mean estimation starts!\n")
    
    # Loop through locations (N)
    for (n in 1:N) {
      Y1=Stand_SD[n, c(1:T_o)+2]
      Y2=Obs[n, -c(1:2)]
      
      # Fit linear regression models
      reg_1=lm(Y1 ~ X+t+X*t,data=data)
      reg_2=lm(Y2 ~ X+t+X*t,data=data)
      
      Mu1[n,]=predict(reg_1, data_new)
      Mu2[n,]=predict(reg_2, data_new)
      
      if (n %% 10000 == 1) {
        cat(paste("Location", n, "\n"))
      }
    }
    
    cat("Mean estimation ends!\n")
    
    # Write the results to CSV files
    fwrite(Mu1,paste0(path,"Mu1.csv"))
    fwrite(Mu2,paste0(path,"Mu2.csv"))
    
    # Create a command to run an external script
    #items=c(paste0(path,"run_BGL.sh"),path_matlabRuntime,path)
    #items = source("BGL")
    # Wrap each item with double quotes and collapse them into a single string
    #command=paste0("\"", paste(items, collapse = "\" \""), "\"")
    source("/Users/lihancheng/Desktop/Downscaling Projec/Matlab2R/BGL.R")
    # Execute the external script using the 'system' function
    #system(print(command))
  }
  
  
}


