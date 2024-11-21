# README
# Script to run the analyses and visualize the results
# ___________________________________________________________
          
          library(tidyverse)
          
# 
# Functions                                       ----
# ____________________________________________________
          
          source('sc_func_Infomap.R')
          
# 
# Load data                                              ----
# ___________________________________________________________

          # La funcion actual requiere el siguiente formato
          # Columna 1: unidades geograficas
          # Columna 2: unidades taxonomicas operacionales
          # Columna 3: weight of the occurrences (abundance, density...)
          
          # Read the file
          da <- otu_orig_cs
          names(da)<- c("SampID", "clean_taxonomy",  "count","density_n100m2", "VL")
          head(da)
          da <- data.frame(da)
# 
# Write input infomap                      ----
# _____________________________________________

          path_write_Input <-  "C:/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/"
          Input_name <- 'Ip_Benthos.net'
          edges_file <- 'edge_Benthos'
          # Create input in the same folder
          write.input.info(da, path_write_Input, Input_name, edges_file, weighted = T)
          # En la funcion he puesto que te imprima las partes en las que consta el input para que te hagas una idea

# 
# Run Infomap                              ----
# _____________________________________________
  
          # Directorio del ejecutable de Infomap
          path.info <- "/mnt/c/R_projects/Hortal_collaboration/MAREANO_biogeography/infomap-ubuntu/"
          # Directorio donde esta el input
          path.in <- "/mnt/c/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/"
          # Directorio donde quieres guardar el output
          path.out <- "/mnt/c/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/"
          
          # -N X = numero de analisis
          # -s = seed
          # empty = hierarchical // -2 = two-levels
          N_trials <- 10 # Yo suelo poner 100
          seed <- 1
          # Si quieres hacer una biorregionalizacion a un solo nivel o explorar la existencia de multiples niveles jerarquicos
          Hierarchical <- T
          if (Hierarchical == T){
            Type <- ''
          } 
          if (Hierarchical == F){
            Type <- '-2'
          } 
          
          parameters <- paste("-N", N_trials, "-s", 1, Type)
          
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
# Explore saved output                          ----
# __________________________________________________
          
          # Run the function to get the 
          dir.info.out <- 'C:/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/Ip_Benthos.tree'
          dir.edges <- 'C:/R_projects/Hortal_collaboration/MAREANO_biogeography/regionalization/edge_Benthos'
            
          # :::::::::::::: 
          # Only the table
          # :::::::::::::: 
          db_output <- info.out(dir.info.out)
          db_output
          
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
        
          