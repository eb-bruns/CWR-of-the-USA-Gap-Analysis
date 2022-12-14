###
# Emily Beckman Bruns testing Dan Carver's modeling method
# 13 December 2022
# R version 4.2.2
# Pulled over variable definitions from run_lineal.R so everything is in
#  one place for testing.
###

###
# Primary function to call all functions that are a part of the modeling method. 
# dan.carver@carver.com 
# 20200414
### 


## pulled from run_lineal.R...


pacman::p_load(tidyverse,sp,raster,rgdal,tmap,devtools,
               randomForest,rgeos,VSURF,modelr,maxnet,
               pROC,dismo,redlistr,fasterize,devtools,DT,
               terra,beepr)
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
run_version <<- "temp20221213"

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
  # EBB: removing more here
source.files = source.files[ !grepl("master", source.files) ]
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
  # can read in as one file like this...
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
    .urb = as.logical(.urb),
    .inst = as.logical(.inst),
    .con = as.logical(.con),
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
    (.cen & .inst & .con & .outl &
       #.urb & .yr1950 & .yr1980 & .yrna &
       #(.gtsnative | is.na(.gtsnative)) &
       #(.rlnative  | is.na(.rlnative)) &
       #(.rlintroduced | is.na(.rlintroduced)) &
       basisOfRecord != "FOSSIL_SPECIMEN" & 
       basisOfRecord != "LIVING_SPECIMEN" &
       establishmentMeans != "INTRODUCED" & 
       establishmentMeans != "MANAGED" &
       establishmentMeans != "CULTIVATED" &
       latlong_countryCode %in% c("US","CA","MX"))) %>%
# edit data to match format needed for GapAnalysis
  mutate(type = recode(database,
                       "Ex_situ" = "G",
                       .default = "H")) %>%
  rename(taxon = taxon_name_accepted,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         databaseSource = database) #%>%
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
genera <- sort(unique(occData$genus))
#testGen <- genera[c(32, 25, 47 )]

# select all species at the genus level and apply master script to run process
#beepr::beep_on_error(
#  for(i in testGen){
i <- genera[1]
    t2a <- Sys.time()
    genus <<- i 
    if (!file.exists(paste0(gap_dir,"/summaryDocs"))) 
      {dir.create(paste0(gap_dir,"/summaryDocs"),recursive=T)}
    allSpec <- occData %>%
      dplyr::filter(genus == i)
    # generate a folder within the gap analysis
    folder <- paste0(occ_dir, "/",i)
    if (!file.exists(folder)) {dir.create(paste0(folder),recursive=T)}
    # test for genus level folder.
    genFolder <- paste0(gap_dir, "/", i)
    if (!file.exists(genFolder)) {dir.create(paste0(genFolder),recursive=T)}
    write.csv(allSpec, paste0(folder, "/", "raw",i,".csv"), row.names = FALSE)
    genusOcc <<- read.csv(paste0(folder, "/", "raw",i,".csv"))
    speciesList <<- sort(unique(allSpec$taxon))
    write.csv(x = speciesList, file = paste0(gap_dir,'/', genus, "/", 'speciesList.csv'))
    
    speciesList <- speciesList[speciesList %in% spList]
    #if(!is.na(speciesList[1])){
      #calls the master function 
      #result_master = lapply(speciesList[1: length(speciesList)], master_run)
    #}
    
    

### main function for modeling...
    
#master_run <- function(species){
  #species <<- species
  species <<- spList[9]
  print(paste0("the process for ", species, " has begun."))
  
  # build a datframe that captures the total run time for a process.
  time_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(time_df) <- c("functionUsed", "runTime")
  
  # start time for totaling run time at the end
  startTime <- Sys.time()
  
  t1a <- Sys.time()
  cat("...creating directories\n")
  create_sp_dirs(species)
  time_df <- rbind(time_df, data.frame(functionUsed="create_sp_dirs", 
                                       runTime=difftime(Sys.time(), t1a, 
                                                        units='secs')))
  # sp_dir <<- primary directory for this run
  
  # if(!file.exists(paste0(sp_dir, '/sdm.rds'))){
  #   print('been modeled, moving on')
    # run <- "Run20200203"
    # try(
    #   rmarkdown::render(paste0(repo_dir, "/summaryMarkdown/singleSpeciesSummary.rmd"),
    #                     output_file =  paste(species, '_', run,'_' , Sys.Date(), ".html", sep=''),
    #                     output_dir = paste0(gap_dir,"/", genus,"/summaryDocs"))
    # )
  #
  # }else{
  t1a <- Sys.time()
  cat("...developing raw data\n")
  developRaw(species)
  time_df <- rbind(time_df, data.frame(functionUsed="developing raw data",  runTime=difftime(Sys.time(), t1a, units='secs')))
  # rawData <<- listing of all data for the species
  #
    t1a <- Sys.time()
    cat("...generating counts csv\n")
    developCounts(species)
    time_df <- rbind(time_df, data.frame(functionUsed="testLatLong",  runTime=difftime(Sys.time(), t1a, units='secs')))
    ## dataThin <<- raw data with complete lat long and only necessary

    t1a <- Sys.time()
    cat("...conducting SRSex assessment\n")
    srs_exsitu(species)
    time_df <- rbind(time_df, data.frame(functionUsed="srs_exsitu",  runTime=difftime(Sys.time(), t1a, units='secs')))

    t1a <- Sys.time()
    cat("...creating spatial points data frame\n")
    spPoints(species)
    time_df <- rbind(time_df, data.frame(functionUsed="spPoints",  runTime=difftime(Sys.time(), t1a, units='secs')))
    # # spPoint <<- all usable occurence data for NA
    
    #if(class(spPoint) == "character"){
    #  run <- "Run20200203"
    #  rmarkdown::render(paste0(repo_dir, "/summaryMarkdown/singleSpeciesSummary.rmd"),
    #                    output_file =  paste(species, '_', run,'_' , Sys.Date(), ".html", sep=''),
    #                    output_dir = paste0(gap_dir,"/", genus,"/summaryDocs"))
    #}else{
      t1a <- Sys.time()
      cat("...adding North American Points to Counts.csv\n")
      try(addNorthAmericanCounts(species))
      time_df <- rbind(time_df, data.frame(functionUsed="addNorthAmericanCounts",  runTime=difftime(Sys.time(), t1a, units='secs')))
      #
      #
      t1a <- Sys.time()
      cat("...srsIn before the modeling process\n")
      srs_insitu_preModel(species)
      time_df <- rbind(time_df, data.frame(functionUsed="srs_insitu_preModel",  runTime=difftime(Sys.time(), t1a, units='secs')))
      
      #if(class(spPoint) == 'character' ){
      #  print("there is not locational data for this species this is the end of the modeling process")
      #}else{
        cat("...extracting country values to points and remove duplicate lat long\n")
        countryCheck(species)
        time_df <- rbind(time_df, data.frame(functionUsed="countryCheck",  runTime=difftime(Sys.time(), t1a, units='secs')))
        # xyData <<- coords for elements in spPoint object
        # cleanPoints <<- spatial points layer were all points are on land.
           #
        #if(class(cleanPoints) == "character"){
        #  print('no usable records for the modeling method')
        #}else{
          t1a <- Sys.time()
          cat("...spatial sample to under 2000 tota\n")
          subSampleByCountry(species)
          time_df <- rbind(time_df, data.frame(functionUsed="sampling",  runTime=difftime(Sys.time(), t1a, units='secs')))
          # cleanPoints <<- if spatail thining is applied the cleanPoint var is redeclared.
          # code needs to be ran to this location inorder to produce htmls
          
            t1a <- Sys.time()
            cat("...generate native area shp\n")
            nat_area_shp(species)
            time_df <- rbind(time_df, data.frame(functionUsed="nat_area_shp",  runTime=difftime(Sys.time(), t1a, units='secs')))
            ##nativeArea <<- ecoregions clipped to countries with species present

                    # t1a <- Sys.time()
                    # cat("...conducting EOO and AOO assessment\n")
                    # eooAoo(species)
                    # time_df <- rbind(time_df, data.frame(functionUsed="eooAoo",  runTime=difftime(Sys.time(), t1a, units='secs')))

                    #if(class(nativeArea)== "character"){
                    #print("not enough points for modelings")
                    #  }else{
                      if(nrow(cleanPoints)<= 3){
                      print("not enough points for modelings")
                        }else{

                        t1a <- Sys.time()
                        cat("...generate ga50Raster\n")
                       create_buffers(species)
                        time_df <- rbind(time_df, data.frame(functionUsed="create_buffers",  runTime=difftime(Sys.time(), t1a, units='secs')))

                        t1a <- Sys.time()
                        cat("...generate background and extract raster data to background and presence data\n")
                        generateModelingData(species)
                        time_df <- rbind(time_df, data.frame(functionUsed="generateModelingData",  runTime=difftime(Sys.time(), t1a, units='secs')))
                        # bioValues <<- Presence and background points with predictor data attached

                        t1a <- Sys.time()
                        cat("...perform variable selection and correlation\n")
                        varaibleSelection(species)
                        time_df <- rbind(time_df, data.frame(functionUsed="varaibleSelection",  runTime=difftime(Sys.time(), t1a, units='secs')))
                        # bioValues <<- redefined if any NA are present in predictor datasets.
                        # variblesToModel <<- a listing of variable names used to modeling

                        t1a <- Sys.time()
                        cat("...perform maxent model\n")
                        runMaxnet(species)
                        time_df <- rbind(time_df, data.frame(functionUsed="runMaxnet",  runTime=difftime(Sys.time(), t1a, units='secs')))
                        # sdm_results <<- output of the Maxnet modeling process

          if(!file.exists(paste0(sp_dir, "/sdm.rds"))){
          print("the model did not successfully run")
            }else{

                          t1a <- Sys.time()
                          evaluate_sdm_function(species)
                          time_df <- rbind(time_df, data.frame(functionUsed="evaluate_sdm_function",  runTime=difftime(Sys.time(), t1a, units='secs')))
                          ## thrshold <<- threshold raster

                          #if(length(unique(values(thrshold)))==2){
                          #print("the model did not successfully run")
                          #  }else{

                            t1a <- Sys.time()
                            cat("...create a mess map based on top predictor \n")
                            messMap(species)
                            time_df <- rbind(time_df, data.frame(functionUsed="messMap",  runTime=difftime(Sys.time(), t1a, units='secs')))

                            t1a <- Sys.time()
                            cat("...create a kernal density map of the sample points  \n")
                            kernalDensity(species)
                            time_df <- rbind(time_df, data.frame(functionUsed="kernalDensity",  runTime=difftime(Sys.time(), t1a, units='secs')))
          #
#                            ### start of the gap analysis metrics
#                            t1a <- Sys.time()
#                            cat("...conducting GRSex assessment\n")
#                            grs_exsitu(species)
#                            time_df <- rbind(time_df, data.frame(functionUsed="grs_exsitu",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#              t1a <- Sys.time()
#              cat("...conducting ERSex assessment\n")
#              try(ers_exsitu(species))
#              time_df <- rbind(time_df, data.frame(functionUsed="ers_exsitu",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#              t1a <- Sys.time()
#              cat("...conducting  fcsex assessment\n")
#              fcs_exsitu(species)
#              time_df <- rbind(time_df, data.frame(functionUsed="fcs_exsitu",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#                            cat("...Ex situ assessment is complete\n")
#
#          t1a <- Sys.time()
#          cat("...conducting srs In assessment\n")
#          try(srs_insitu(species))
#          time_df <- rbind(time_df, data.frame(functionUsed="srs_insitu",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#          t1a <- Sys.time()
#          cat("...conducting grs In assessment\n")
#          try(insitu_grs(species))
#          time_df <- rbind(time_df, data.frame(functionUsed="insitu_grs",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#          t1a <- Sys.time()
#          cat("...conducting ers In assessment\n")
#          try(ers_insitu(species))
#          time_df <- rbind(time_df, data.frame(functionUsed="ers_insitu",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#          t1a <- Sys.time()
#          cat("...conducting fcs In assessment\n")
#          try(fcs_insitu(species))
#          time_df <- rbind(time_df, data.frame(functionUsed="fcs_insitu",  runTime=difftime(Sys.time(), t1a, units='secs')))
#                  #
#                                cat("...in situ assessment is complete\n")
#
#          t1a <- Sys.time()
#          cat("...conducting fcs combined assessment\n")
#          fcs_combine(species)
#          time_df <- rbind(time_df, data.frame(functionUsed="fcs_combine",  runTime=difftime(Sys.time(), t1a, units='secs')))
#
#          t1a <- Sys.time()
#          run <- "Run20200915"
#          rmarkdown::render(paste0(repo_dir, "/summaryMarkdown/singleSpeciesSummary.rmd"),
#                            output_file =  paste(species, '_', run,'_' , Sys.Date(), ".html", sep=''),
#                            output_dir = paste0(gap_dir,"/", genus,"/summaryDocs"))
#          
#          
          time_df <- rbind(time_df, data.frame(functionUsed="sinlge species Summary",  runTime=difftime(Sys.time(), t1a, units='secs')))
          time_df$Minutes <- time_df$runTime/60
          write.csv(x = time_df, file = paste0(sp_dir, "/time_df.csv"))
          print(paste0("the process for ", species, " has ended."))

          
          # remove all global variables that are specific to the species
          #try(globalVars <- c(xyData, cleanPoints, sp_dir,rawData,dataThin, nativeArea, spPoint,pa_spp,
          #                    ecoValsAllPointsLen,pa_spp_area, gBufferRas1, gBufferRas_area,
          #                    sp_counts,bioValues,sdm,evaluate_table,sdm_results, thrshold,
          #                    evaluate_table_f, crossValDir,rastersToModel, data_train, bioValues,
          #                    variblesToModel,species))
          #
          #try(rm(globalVars))
          beepr::beep(1)
          
    }
    
    
 
      }
     }
    }
  }
 }
}
   }
# #}




