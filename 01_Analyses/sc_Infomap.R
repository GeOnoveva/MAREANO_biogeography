# README
# Script to run the analyses and visualize the results
# ___________________________________________________________
          
          library(tidyverse)
          library(sp)
          library(qdapTools)
          library(openxlsx)
          
# 
# Functions                                       ----
# ____________________________________________________
          
          source('01_Analyses/sc_func_Infomap.R')
          
# 
# Load
# ____________________________________________________________________
          
          # person  <- 'Genoveva'
          person  <- 'Ruben'  
          
          if (person == 'Ruben') {comp_path <- "00_Data/"}
          if (person == 'Genoveva') {comp_path <- "U:\\Mareano\\VIDEOLAB\\VIDEO DATA\\200m_scale_species_by_sample\\Data_Delivery_2024\\"}
          
          if (person == 'Genoveva') {source('sc_func_Infomap.R')}
          if (person == 'Ruben') {source('01_Analyses/sc_func_Infomap.R')}
          
          spp_dens <- read.csv(paste0(comp_path,"species_densities.csv")) %>% as.data.frame
          samp_info_cs <- read.csv(paste0(comp_path,"sample_info.csv")) %>% as.data.frame
          
# 
# Station filter
# ____________________________________________________________________
          spp_dens <- spp_dens %>% 
            mutate(VL = SampID %l% data.frame(samp_info_cs$SampID2, samp_info_cs$VL))
          
          otu_orig_cs <- spp_dens %>% filter(!is.na(SampID))
          
# 
# Data cleaning
# ___________________________________________________________________
          
          if (person == 'Ruben') {PathTaxonary <- "00_Data/Taxonary.xlsx"}
          # if (person == 'Genoveva') {PathTaxonary <- "data"}
          Taxonary<-read.xlsx(PathTaxonary , sheet = 1)
          
          # filter small
          otu_new_cs <- otu_orig_cs %>% filter(clean_taxonomy %l% data.frame(Taxonary$Reference_List, Taxonary$basal_area_cm2)>4|
                                                 is.na(clean_taxonomy %l% data.frame(Taxonary$Reference_List, Taxonary$basal_area_cm2)))
          
          # filter non benthic
          otu_new_cs <- otu_new_cs %>% filter(clean_taxonomy %l% data.frame(Taxonary$Reference_List, Taxonary$Ecosystem_section=="Benthic"))
          
          # filter non organism
          otu_new_cs <- otu_new_cs %>% filter(clean_taxonomy %l% data.frame(Taxonary$Reference_List, Taxonary$Object_type=="Organism"))
          
          
          
          
# 
# Load data                                              ----
# ___________________________________________________________

          # La funcion actual requiere el siguiente formato
          # Columna 1: unidades geograficas
          # Columna 2: unidades taxonomicas operacionales
          # Columna 3: weight of the occurrences (abundance, density...)
          
          # Read the file
          da <- otu_new_cs
          names(da)<- c("SampID", "clean_taxonomy",  "count","density_n100m2", "VL")
          head(da)
          da <- data.frame(da) %>% 
            select(SampID, clean_taxonomy, density_n100m2)
          
          
# 
# Write input infomap                      ----
# _____________________________________________

          if(person=='Genoveva'){
            path_write_Input <-  "C:/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/"
            # Directorio del ejecutable de Infomap
            path.info <- "/mnt/c/R_projects/Hortal_collaboration/MAREANO_biogeography/infomap-ubuntu/"
            # Directorio donde esta el input
            path.in <- "/mnt/c/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/"
            # Directorio donde quieres guardar el output
            path.out <- "/mnt/c/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/"
            # Run the function to get the 
            dir.info.out <- 'C:/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/Ip_Benthos.tree'
            dir.edges <- 'C:/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/edge_Benthos'
          }
          
          
          if(person=='Ruben'){
            # Input
            path_write_Input <-  "C:/Users/Ruben/temporal/Benthos_Norway/00_Data/"
            # Directorio del ejecutable de Infomap
            path.info <- "/mnt/c/Users/Ruben/temporal/Benthos_Norway/infomap-ubuntu/"
            # Directorio donde esta el input
            path.in <- "/mnt/c/Users/Ruben/temporal/Benthos_Norway/00_Data/"
            # Directorio donde quieres guardar el output
            path.out <- "/mnt/c/Users/Ruben/temporal/Benthos_Norway/01_Analyses/"  
            # Run the function to get the 
            dir.info.out <- 'C:/Users/Ruben/temporal/Benthos_Norway/01_Analyses/Ip_Benthos.tree'
            dir.edges <- 'C:/Users/Ruben/temporal/Benthos_Norway/00_Data/edge_Benthos'
          }
          
          
          Input_name <- 'Ip_Benthos.net'
          edges_file <- 'edge_Benthos'
          # Create input in the same folder
          write.input.info(da, path_write_Input, Input_name, edges_file, weighted = T)
          # En la funcion he puesto que te imprima las partes en las que consta el input para que te hagas una idea

# 
# Run Infomap                              ----
# _____________________________________________
  
          # -N X = numero de analisis
          # -s = seed
          # empty = hierarchical // -2 = two-levels
          N_trials <- 100 # Yo suelo poner 100
          seed <- 1
          # Si quieres hacer una biorregionalizacion a un solo nivel o explorar la existencia de multiples niveles jerarquicos
          Hierarchical <- T
          if (Hierarchical == T){
            Type <- ''
          } 
          if (Hierarchical == F){
            Type <- '-2'
          } 
          
          parameters <- paste("-N", N_trials, "-s", seed, Type)
          
          # Run 
          par.n <- run.infomap (path.info, path.in, Input_name, path.out, parameters)
         
          # 
          # First exploration 
          # :::::::::::::::::::::::::::::::::::
          
          # partition
          head(par.n [[1]])
          # Quality of the partition
          par.n [[2]]
          # Summary of the results
          par.n [[3]]
  
          
# 
# Grid
# _________________________________________________
          
          library(sf)
          gr <- samp_info_cs %>% 
            distinct(SampID, x_coord, y_coord) %>% 
            st_as_sf(coords = c("x_coord","y_coord"))
          
          
# 
# Explore saved output                          ----
# __________________________________________________
          
          # :::::::::::::: 
          # Only the table
          # :::::::::::::: 
          db_output <- info.out(dir.info.out)
          db_output
          
          
  # 
  # Tamaño regiones
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

          db_sit <- db_output %>% filter(Name %in% da$SampID)
          db_spp <- db_output %>% filter(Name %in% da$clean_taxonomy)
          
          # Que especies o sitios estan solas
          db_output %>% filter(Module_3 %in% unique(db_spp$Module_3)[!unique(db_spp$Module_3) %in% db_sit$Module_3])
          db_output %>% filter(Module_3 %in% unique(db_sit$Module_3)[!unique(db_sit$Module_3) %in% db_spp$Module_3])
          
          da %>% filter(clean_taxonomy == 'Poraniomorpha sp.') %>% tally()
          da %>% filter(clean_taxonomy %in% 'Colossendeis sp.') %>% tally()
          
          # Frecuencia de los sitios
          db_sit %>% count(Module_3) %>% arrange(n) %>% filter(n<=10)
          db_sit %>% count(Module_3) %>% arrange(n) %>% filter(n>10)
          
          db_sit %>% count(Module_2) %>% arrange(n) %>% filter(n<=10)
          db_sit %>% count(Module_2) %>% arrange(n) %>% filter(n>10)
          
          db_sit %>% count(Module_1) %>% arrange(n) %>% filter(n<=10)
          db_sit %>% count(Module_1) %>% arrange(n) %>% filter(n>10)
          
          gr_3 <- gr %>% left_join(db_sit, by=c('SampID'='Name')) %>% filter(!is.na(Module_3)) # %>%  group_by(Module_3) %>% summarise()
          gr_2 <- gr %>% left_join(db_sit, by=c('SampID'='Name')) %>% filter(!is.na(Module_2)) # %>%  group_by(Module_2) %>% summarise()
          gr_1 <- gr %>% left_join(db_sit, by=c('SampID'='Name')) %>% filter(!is.na(Module_1)) # %>%  group_by(Module_1) %>% summarise()
          
          gr_1 <- gr %>% left_join(db_sit, by=c('SampID'='Name')) %>% filter(!is.na(Module_1)) 
          
          plot(gr_1[8])
          
          
          
          
          
# Muchas biorregiones tienen solo un sitio
          # 
          # 
          
          
          # ::::::::::::::::::::::::::::::::
          # Table plus metrics for the roles
          # ::::::::::::::::::::::::::::::::
          Output <- func_Link_Info(dir.info.out, dir.edges, Hierarchical_level = 'last_mod', show.add.metrics = T, ret='all', links_as_cell_spp = T)
          # Esto genera un output de 4 archivos
          names(Output)
          # Output of Infomap
          info <- Output[[which(names(Output)=='Infomap')]]
          # Links and associated modules
          db_Links <- Output[[which(names(Output)=='links')]]
          # Indicative values of the species of each Module
          db_IndVal <- Output[[which(names(Output)=='db_IndVal')]]
          # Metrics to estimate the biogeographical roles
          db_Metrics <- Output[[which(names(Output)=='db_Metrics')]]
          db_Metrics <- db_Metrics %>%  
            dplyr::select(cell, Module, Rel_Species_Richness, Biota_Overlap, Endemicity, Rel_Occupancy) 
          
          
# plot in GIS
          library(sf)

          db_Spatial <- db_Metrics %>% mutate(x_coord = cell %l% data.frame(samp_info_cs$SampID2, samp_info_cs$x_coord)) %>%
            mutate(y_coord = cell %l% data.frame(samp_info_cs$SampID2, samp_info_cs$y_coord)) 

          db_Spatial <- db_Spatial %>% st_as_sf(coords = c("x_coord","y_coord"))
          st_write(db_Spatial, ("infomap_result.shp"))
          
          plot(db_Spatial)
        
          