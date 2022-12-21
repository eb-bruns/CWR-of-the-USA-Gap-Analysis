###
# Emily Beckman Bruns testing Dan Carver's modeling method
# 20 December 2022
# R version 4.2.2
###

### 
# Primary script for running CWRNA. This should be the only location where users 
# have to edit information.
# 20200414
# dan.carver@carverd.com
### 

pacman::p_load(tidyverse,sp,raster,rgdal,tmap,devtools,
               randomForest,rgeos,VSURF,modelr,maxnet,
               pROC,dismo,redlistr,fasterize,devtools,DT,
               terra,beepr,sf)
# EBB: install earlier version of spatialEco -- kernal density doesn't work
#      sometimes when using new version (2.0)
#devtools::install_version("spatialEco", version = "1.3.1", repos = "http://cran.us.r-project.org")
library(spatialEco)

#install_github("ccsosa/GapAnalysis")
#devtools::install_github("DFJL/SamplingUtil")
#devtools::install_github("hunzikp/velox")
#devtools::install_github("valentinitnelav/geobuffer")
library(geobuffer)
library(velox)
library(SamplingUtil)
tmap::tmap_mode("view")


# set all standard directories
base_dir <<- "/Users/emily/Desktop/*work/NorthAm-CWR"
gap_dir <<- paste0(base_dir , "/gap_analysis")
par_dir <<- paste0(base_dir , "/parameters")
occ_dir <<- paste0(par_dir, "/occurrenceData")
#temp_dir <<- paste0(base_dir , "/TEMP")
repo_dir <<- "/Users/emily/Documents/GitHub/CWR-of-the-USA-Gap-Analysis"


# set name of the run version 
run_version <<- paste0("temp-",Sys.Date())

# set adjustable parameters 
numPoints <<- 2000 # maximum number of points used in model (subSampleCountry.R)
bufferDist <<- 50000 # used to define buffer distance in gBuffer.r
set.seed(1234)

# load the sources scripts
# get list of scripts in repo
source.files = list.files(repo_dir, ".[rR]$", full.names = TRUE, recursive = T)
source.files = source.files[ !grepl("dataBaseTransform", source.files) ]
source.files = source.files[ !grepl("test", source.files) ]
source.files = source.files[ !grepl("lineal", source.files) ]
source.files = source.files[ !grepl("summaryMarkdown", source.files) ]
# EBB: removing more here that I don't want to run right now
source.files = source.files[ !grepl("gatherAllCode", source.files) ]
# load everything
lapply(source.files, source)
#for(i in 1:length(source.files)){
#  cat(i,"\n")
#  source(source.files[i])
#}


## set all primary file sources...

# worldclim variables
bioVars <<- read_rds(paste0(par_dir,"/bioLayer_2.5/parameters.RDS"))

# country boundaries; downloaded from:
#   https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/
countrySHP <<- rgdal::readOGR(paste0(par_dir,"/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp"),verbose = FALSE)

# country shapefile clipped to target countries
naSHP <<- countrySHP[(countrySHP@data$ADMIN == "United States of America" |
                        countrySHP@data$ADMIN == "Canada" |
                        countrySHP@data$ADMIN == "Mexico"), ]
#naSHP <<- rgdal::readOGR(paste0(par_dir,"/allUSAArea/NorthAmerica_AllUSA.shp"), verbose = FALSE)
# excluding pacific territories- runs near all species faster 
#naSHP <<- readOGR(paste0(par_dir,"/northAmericaArea/northAmericaArea.shp"),verbose = FALSE)

# global ecoregions & protected areas; downloaded from:
#   https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/WTLNRG
ecoReg <<- rgdal::readOGR(paste0(par_dir,"/ecoregions/tnc_terr_ecoregions.shp"),verbose = FALSE)
proArea <<- raster::raster(paste0(par_dir,"/protectedAreas/wdpa_reclass.tif"))
# running with the PAUD dataset 
#proArea <<- raster::raster(paste0(par_dir, "/protectedAreas/PAUDrasters/allAreas.tif"))

# target taxa occurrence point data
# EBB: compiled using this repo: https://github.com/eb-bruns/SDBG_CWR-trees-gap-analysis
#      can read in as one file like this...
#occData <- data.table::fread(paste0(base_dir,"/occurrence_data2019_05_29/combinedOccurance2020-07-21a.csv"),header = TRUE)
#occData <<- occData[,2:ncol(occData)]
# ...or as taxon-level files and combine:
taxonFiles <- list.files(paste0(occ_dir,"/taxon_edited_points"),pattern = ".csv",full.names = TRUE)
taxonDfs <- lapply(taxonFiles,read.csv,header = TRUE,na.strings = c("","NA"),colClasses = "character")
occData <<- Reduce(bind_rows, taxonDfs)
# as needed: filter based on flagging columns from 5-refine_occurrence_data.R
occData <- occData %>%
  mutate(
    .cen = as.logical(.cen),
    #.urb = as.logical(.urb),
    .inst = as.logical(.inst),
    #.con = as.logical(.con),
    .outl = as.logical(.outl),
    #.gtsnative = as.logical(.gtsnative),
    #.rlnative = as.logical(.rlnative),
    .yr1950 = as.logical(.yr1950),
    .yr1980 = as.logical(.yr1980),
    .yrna = as.logical(.yrna)
  ) %>%
  # select or deselect these filters as desired:
  filter(
    #database == "Ex_situ" |
    (.cen & .inst & .outl &
       #.con &
       #.urb & .yr1950 & .yr1980 & .yrna &
       #(.gtsnative | is.na(.gtsnative)) &
       #(.rlnative  | is.na(.rlnative)) &
       #(.rlintroduced | is.na(.rlintroduced)) &
       basisOfRecord != "FOSSIL_SPECIMEN" & 
       basisOfRecord != "LIVING_SPECIMEN" &
       establishmentMeans != "INTRODUCED" & 
       establishmentMeans != "MANAGED" &
       establishmentMeans != "CULTIVATED" #&
       #latlong_countryCode %in% c("US","CA","MX")
     )) %>%
  ## remove specific ids
  filter(
    !(id == "id00048470")
  ) 
  # edit data to match format needed for GapAnalysis
  mutate(type = recode(database,
                       "Ex_situ" = "G",
                       .default = "H")) %>%
  rename(taxon = taxon_name_accepted,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         databaseSource = database) %>%
  dplyr::select(-genus) %>%
  separate(col="taxon",into="genus",sep=" ",extra="drop",remove=F) #%>%
#dplyr::select(genus,taxon,latitude,longitude,type)
#occData$species <- gsub(" ","_",occData$species)

# list of target taxa
spList <- unique(occData$taxon)

# reference names for the raster data. RDS store files as layer.1, layer.2 etc;
# name them at later point in the workflow as the full names have a lot of 
# spaces/special characters that will make indexing annoying.
#layerDescription <<- read.csv(paste0(par_dir,"/layerDesrciptions.csv"))

# reference file for what states specific species were native too; used as a
# data quality check, basically only accepting location data from states that
# were reference to contain that species in GRIN. If you have more confidence 
# in the quality of your input dataset you can skip this step.
#statesData <<- read.csv(paste0(par_dir,"/statePerTaxon/CWRofUSA_nativeareas_2020_1_30.csv"))

# states shapefile, used for filtering occurrences by state;not needed 
# if you don't do the previous step. 
#statesSpObject <<- readRDS(paste0(par_dir,"/statesByCountry/gadmCanUsaMex_sp.rds"))


## final setup...

# set loop at genus level
genera <- unique(occData$genus)
#testGen <- genera[c(32, 25, 47 )]

# select all species at the genus level and apply master script to run process
#beepr::beep_on_error(
#  for(i in genera){
i <- genera[5]

    t2a <- Sys.time()
    genus <<- i
    if (!file.exists(paste0(gap_dir,"/summaryDocs"))) 
      {dir.create(paste0(gap_dir,"/summaryDocs"),recursive=T)}
    allSpec <- occData %>% dplyr::filter(genus == i)
    # generate a folder within the gap analysis
    folder <- paste0(occ_dir, "/",i)
    if (!file.exists(folder)) {dir.create(paste0(folder),recursive=T)}
    # test for genus level folder
    genFolder <- paste0(gap_dir, "/", i)
    if (!file.exists(genFolder)) {dir.create(paste0(genFolder),recursive=T)}
    write.csv(allSpec, paste0(folder, "/", "raw",i,".csv"), row.names = FALSE)
    genusOcc <<- read.csv(paste0(folder, "/", "raw",i,".csv"))
    speciesList <<- sort(unique(allSpec$taxon))
    write.csv(x = speciesList, file = paste0(gap_dir,'/', genus, "/", 'speciesList.csv'))

    speciesList <- speciesList[speciesList %in% spList]
  if(!is.na(speciesList[1])){
    #calls the main function 
    result_main = lapply(speciesList[1: length(speciesList)], main_run)
  }
#}  
#)

    
    
# pull all maxent evaluation metrics into one file...

# compile all taxon-level csv results files, then save as one output
compile_results <- function(file_name){
  files <- list.files(path = base_dir, pattern = file_name, 
                      full.names = TRUE,  recursive = TRUE)
  dfs <- lapply(files, read.csv, header=TRUE)
  results <- Reduce(bind_rows, dfs)
  print(nrow(results))
  write.csv(results, file.path(base_dir,"gap_analysis","summaryDocs",file_name))
}

compile_results("eval_metrics.csv")















# function for filtering list based on character values
include <- function (theList, toMatch){
  matches <- unique (grep(paste(toMatch,collapse="|"),
                          theList, value=TRUE))
  return(matches)
}
# 
eval <- include(theList = eval, toMatch = spList)
eval <-fcs[grepl(pattern = "insitu",x = eval)]
eval <- fcs[grepl(pattern = "test20200203", x = eval)]
fcs
    for(i in 1:length(fcs)){
      if(i == 1){
        all <- read.csv(fcs[i])
      }else{
        all <- rbind(all, read.csv(fcs[i]))
      }
      
  
eval <- list.files(path = base_dir,pattern = "_Run20200203_2020-08-18.html", full.names = TRUE,  recursive = TRUE)
folder <- "F:/nrelD/cwrNA/runSummaries/speciesLevelHTML"
fcs <- include(theList = fcs, toMatch = spList)
for(i in fcs){
  ## currently set to only include the trouble shooting I've been working on today
    file.copy(i, folder)
}

### test to see what species are included
ot <- sort(unique(occData$taxon))
ot <- data.frame(taxon = sort(unique(occData$taxon)), summaryDoc = NA)
n = 1 
for(i in ot$taxon){
   a <- include(fcs, i)
   if(length(a) > 0){
     ot$summaryDoc[n] <- a
   }
   n = n+1
}

reRun <- rmSpec[!rmSpec %in% ot$taxon]

# Ficus 

# interestingly no real issue with Elymus lanceolatus subsp. lanceolatus, Helianthus praecox subsp. praecox, 
# , Ipomoea violacea

# compile a list of all fcs data for troublesome species 
fcs <- list.files(path = base_dir,pattern = "summary.csv", full.names = TRUE,  recursive = TRUE)
# function for flitering list based on character values
include <- function (theList, toMatch){
  matches <- unique (grep(paste(toMatch,collapse="|"),
                          theList, value=TRUE))
  return(matches)
}

fcs <- include(theList = fcs, toMatch = spList)
fcs <-fcs[grepl(pattern = "insitu",x = fcs)]
fcs <- fcs[grepl(pattern = "test20200203", x = fcs)]
fcs
for(i in 1:length(fcs)){
  if(i == 1){
    all <- read.csv(fcs[i])
  }else{
    all <- rbind(all, read.csv(fcs[i]))
  }
}
write.csv(x = all, file = "F:/nrelD/cwrNA/troubleshooting/srsInsituWDPAreruns20200813.csv")
