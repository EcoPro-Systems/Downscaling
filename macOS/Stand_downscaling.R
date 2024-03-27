#This is a function to perform Standard downscaling
#The function takes below inputs and returns the estimated trend as a matrix. Each row is corresponds to the location and each column corresponds to the respective time point.
#Inputs::::
#Obs:This is the matrix of observational data
#Model:This is the matrix of model data
#Fre:Frequency of your data

Stand_downscaling=function(Obs,Model,Freq){
  
  T_o=ncol(Obs)-2 #Observational time points
  T_all=ncol(Model)-2 #Total time points
  ############################################ Calculate climatology at fine pixels
  #Take the average 
  T_o_sequence=rep(1:Freq,times=T_o/Freq)
  T_all_sequence=rep(1:Freq,times=T_all/Freq)
  
  Model_current=Model[,c(1:T_o)+2] #Model projections for the current period
  
  ############################################## Correct for model bias (Subtract model climatology from model data
  
  #Subtract model mean
  Detrend.Model=Model[,-c(1:2)]
  Obs_average=matrix(NA,nrow(Obs),Freq)
  for(i in 1:Freq){
    Detrend.Model[,which(T_all_sequence==i)]=Detrend.Model[,which(T_all_sequence==i)]-rowMeans(Model_current[,which(T_o_sequence==i)])
    Obs_average[,i]=rowMeans(Obs[,which(T_o_sequence==i)+2])
    }
  
  x0=Obs[,1]#list of fine pixel longitudes
  y0=Obs[,2]#list of fine pixel latitudes
  
  #interpolation
  x=as.vector(Model[,1]) # list of coarse grid longitudes (including coarse coral grids and nearest neighbors)
  y=as.vector(Model[,2]) # list of coarse grid latitudes 
  
  W=as.data.frame(Detrend.Model) # De-trended values
  
  library(foreach)
  library(doParallel)
  
  cores=detectCores()-2
  parallelCluster = makePSOCKcluster(cores, master=nsl(Sys.info()['nodename']),outfile='')
  
  setDefaultCluster(parallelCluster)
  doParallel::registerDoParallel(parallelCluster)  
  
  Stand_SD=foreach::foreach(i=1:ncol(W), .combine = 'cbind', .inorder = T) %dopar% {
    suppressWarnings(akima::interpp(x=x,y=y,z=as.vector(W[,i]),xo=x0,yo=y0,linear = FALSE, extrap = TRUE)$z)
  }
  
  parallel::stopCluster(parallelCluster)
  
  Stand_SD=Stand_SD+Obs_average[,T_all_sequence]
  return(Stand_SD)
}