#!/usr/bin/Rscript

# Extract command line arguments
args=commandArgs(trailingOnly = TRUE)

path=args[1]#Path to your directory
file=args[2]#getwd()#file name
Freq=eval(parse(text = args[3])) #Frequency of your data. If it is monthly frequency is 12. If it is daily, it is 365.
obs_year=eval(parse(text = args[4])) #Stat year of the observational period
downscale_year=eval(parse(text = args[5])) #Till the your you plan to perform downscaling

# Check if the package is installed; if not, install it
# Define the package name you want to load
package_name <- "ncdf4"

# Check if the package is installed; if not, install it
if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
  install.packages(package_name, repos = "http://cran.us.r-project.org")
  library(package_name, character.only = TRUE)
}

package_name <- "data.table"

# Check if the package is installed; if not, install it
if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
  install.packages(package_name, repos = "http://cran.us.r-project.org")
  library(package_name, character.only = TRUE)
}


print("Start processing netcdf file!")
#Read netcdf file
nc_file=ncdf4::nc_open(paste0(path,file))

lon=ncvar_get(nc_file,'lon')-360
lat=ncvar_get(nc_file,'lat')
time=ncvar_get(nc_file,'time')
years=rep(1970:2099,each=Freq)
tos=ncvar_get(nc_file,'sst')

print("Details of netcdf file:")
print("Dimenstions:")
print(dim(tos))
print("Longitudes range:")
range(lon)
print("Latitudes range:")
range(lat)
print("Time range:")
print(range(years))

tos=tos[,,which(years>=obs_year & downscale_year>=years)]

cat("\n")

Model_data=cbind(expand.grid(lon,lat))

for(i in 1:dim(tos)[3]){
  Model_data=cbind(Model_data,as.vector(tos[,,i]))
}

Model_data=na.omit(Model_data)
colnames(Model_data)=c("lon","lat",years[which(years>=obs_year & downscale_year>=years)])

print("Dimenstions of the Model_data.csv file:")
dim(Model_data)
fwrite(Model_data,paste0(path,"Model_data.csv"))

print("Model_data.csv file created!")
