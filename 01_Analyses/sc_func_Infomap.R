# README                                    ----
# scripts with all functions
# ______________________________________________

# 
# Execute Infomap                          ----
# ______________________________________________

          run.infomap <- function(path.info, path.in, name, path.out, parameters){
            # Run Infomap
            system(paste("bash -c", shQuote(paste("cd ", path.info ," && ./Infomap ", path.in, name , " ", path.out, " ", parameters, sep=""),type="cmd")),wait=T)
            # Read partition
            out.win<-gsub("/mnt/c","C:",path.out)
            par<-read.table(paste(out.win, gsub('.net', '.tree', name), sep=""))
            # Compression codelength
            load_quality <- readLines(paste(out.win, gsub('.net', '.tree', name), sep=""))
            quality <- load_quality[grep('# relative codelength', load_quality)]
            quality <- gsub('# relative codelength savings ', '', quality)
            quality <- as.numeric(gsub('%', '', quality))
            # Store results
            res <- list(par, # partition
                        quality, # Quality of the partition. Range [0-1]
                        load_quality[1:8] # Summary of the results
            )
            return(res)
          }

# 
# Write Infomap Input     ----
# ________________________________________________
# This function transform a long         


          write.input.info<-function(mat, path, file, edges_file, weighted=F){
            # Identify unique gegraphical units
            id.cell <- unique(mat[,1]) %>% pull
            # Identify unique OTUs
            id.sps <- unique(mat[,2])	%>% pull
            # Occurrences
            if(weighted==F){
              link.node<-mat[,1:2]  
            }
            if(weighted==T){
              link.node<-mat[,1:3]  
            }
            # Create a database to relate the real names with the names in Infomap (a vector from 1 to the number of nodes)
            vert<-cbind(c(1:(length(id.cell)+length(id.sps))),c(id.cell,id.sps))
            # Substitute the real names by Infomap ones in the links
            link.node[,1] <- as.numeric(vert[match(link.node[,1], vert[,2]),1])
            link.node[,2] <- as.numeric(vert[match(link.node[,2], vert[,2]),1])
            # Store the four elements of the input
            input<-list(cbind("*Vertices",nrow(vert)),vert,cbind("*Edges",nrow(link.node)),link.node)
            # Write the input
            write.table(input[[1]], paste(path,file,sep=""),row.names=F,col.names=F,quote=F )
            write.table(input[[2]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=c(2))
            write.table(input[[3]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=F )
            write.table(input[[4]], paste(path,file,sep=""), append= T,row.names=F,col.names=F,quote=F)
            
            # Write edges for Roles
            write.table(input[[4]], paste(path,edges_file,sep=""), row.names=F,col.names=F,quote=F )
            
            # Print for visualizing
            print(head(input[[1]]))
            print(head(input[[2]]))
            print(tail(input[[2]]))
            print(head(input[[3]]))
            print(head(input[[4]]))
          }


# 
# info.output                                          ----
# _________________________________________________________
# Function to transform Infomap output (.tree) into a table
          
          # Example
          # x <- 'C:/Users/Ruben/temporal/Benthos_Norway/01_Analyses/Ip_Benthos.tree'
          
          info.out <- function(x){
            
            library(tidyverse)
            
            ## 
            ## Create a processed data.frame of the Infomap partition     ####
            ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
            # read.table skips the first lines with comments from Infomap output, providing a table with four variables
            par <- read.table(x) %>% 
              # Set the names
              # Classif: The classification each node belongs to. Originally called "path" in Infomap (https://www.mapequation.org/infomap/)
              # Flow:
              # Name: Original name of the node
              # Network_Name: Numerical name of the node assigned when creating Infomap Input.
              setNames(c("Classif", "Flow", "Name", "Network_Name")) %>% 
              
              ### 
              ### Extract hierarchical modules in distinct variables (non-used for estimating the four biodiversity-metrics; see next section) ####
            ### 
            # Create a reverse copy to delete the first element (ID od the node in the cluster) and generate the hierarchical levels.
            mutate (Level_mirror = stringi::stri_reverse(Classif)) %>% 
              # Remove the the ID of each node in each cluster
              mutate (Level_mirror_ss = substr(x = Level_mirror,
                                               start = str_locate(Level_mirror, pattern = ":")[,1]+1,
                                               stop = nchar(Level_mirror))) %>%
              # Revert to original order
              mutate (Level = stringi::stri_reverse(Level_mirror_ss)) %>%
              # Remove temporal and unnecesary variables
              dplyr::select(-Level_mirror, -Level_mirror_ss) %>%
              # Put the hierarchical levels in independent columns
              # The code will create as many columns as the maximum number of hierarchical levels of the Infomap partition.
              splitstackshape::cSplit(., c("Level"), sep = ":", type.convert = F) %>% 
              tibble() # Convert to tibble to rename variables in the future (see below db_Module)
            
            
            ###  
            ### Extract hierarchical modules in a format to study biodiversity metrics and plots   ####
            ### 
            
            # Indicate maximum hierarchical levels observed
            N_max_levels = max(stringr::str_count(par$Classif, ":")) 
            # Empty list to include the hierarchical levels of each node
            sto_modules <- vector(mode = "list", length = nrow(par))
            
            for (r in 1:nrow(par)){
              # Position of the character ":" indicating the nestedness in the Infomap output
              pos <- as.vector(gregexpr("\\:", par$Classif[r])[[1]])
              
              # Get the hierarchical levels for that node
              sto_row <- c()
              for (p in 1:length(pos)){
                sto_row[length(sto_row)+1] <- substr(par$Classif[r], start=1, stop=pos[p]-1)
              }
              
              # For nodes with less hierarchical levels than the maximum, set to 0 until the hierarchical levels are complete. 
              while(length(sto_row)<N_max_levels){
                sto_row[length(sto_row)+1] <- paste0(sto_row[length(sto_row)], ":0")
              }
              
              # Store the hierarchical levels of the node in the list
              sto_modules[[r]] <- sto_row
            }
            
            ##  
            ## Merge Infomap output with processed hierarchical levels             ####
            ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
            # Create a table with the hierarchical modules of each node
            Modules_assign <- do.call(rbind, sto_modules) %>% 
              as.data.frame() %>% # To use setNames in the following command line
              setNames(paste0("Module_", 1:length(sto_modules[[1]]))) %>% 
              tibble()
            
            # Merge the modules and Infomap output
            db_Module <- par %>% 
              bind_cols(Modules_assign) 
            
            # Rename the variables modules
            n_pos <- grep("Module_", colnames(db_Module))
            if (length(n_pos)>1){
              db_Module[,n_pos] <- apply(db_Module[,n_pos],2,function(x)gsub(":","_",x))  
            }
            
            return(db_Module)
          }
          
# 
# Roles                                                ----
# _________________________________________________________
          
          # # # Files for example
          # dir.info.out <- 'C:/Users/Ruben/temporal/Benthos_Norway/01_Analyses/Ip_Benthos.tree'
          # dir.edges <- 'C:/Users/Ruben/temporal/Benthos_Norway/00_Data/edge_Benthos'
          
          func_Link_Info <- function(dir.info.out, dir.edges, Hierarchical_level, show.add.metrics, ret, links_as_cell_spp){
            
            library(tidyverse)
            
            # Load files
            info <- info.out (dir.info.out)  # output infomap
            links <- read.table(dir.edges, header=F) %>% tibble() # Links with numeric names used in the input of Infomap (called Edges in mapequation; https://www.mapequation.org/infomap/)
            links <- links[,1:2]
            
            ## 
            ## Select the Name and the desired Module              ####
            ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            
            if (Hierarchical_level %in% c('last_mod', colnames(info)[grep('Module_', colnames(info))])) {
              
              # Focus on the last hierarchical level (default)
              if (Hierarchical_level %in% c('last_mod')){
                # Get the names of columns that match the pattern "Module_"
                module_cols <- colnames(info) %>%  grep("^Module_", ., value = TRUE)
                # Select column A and the last column that matches the pattern "Module_"
                info <- info %>% 
                  dplyr::select(Name, Network_Name, all_of(tail(module_cols, 1))) %>% 
                  setNames(c('Name', 'Network_Name', 'Module'))
              } 
              
              # Focus on a different hierarchical level than the last one.
              if (Hierarchical_level %in% colnames(info)[grep('Module_', colnames(info))]){
                info <- info[,c('Name', Hierarchical_level)]
                info <- info %>%  setNames(c('Name', 'Network_Name', 'Module'))
              } 
              
            } else {
              print('Indicate which hierarchical level:')
              print(colnames(info)[grep('Module_', colnames(info))])
              break()
            } 
            
            
            ## 
            ## Table with links and modules assigned to each node                 ####
            ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            if (links_as_cell_spp==F){
              links <- cbind(links[,2], links[,1])
            } 
            
            db_links <- links %>% 
              tibble() %>% 
              setNames(c('cell_Network', 'spp_Network')) %>% 
              left_join(info, by=c('cell_Network' = 'Network_Name')) %>% 
              rename(Module_cell = Module,
                     cell = Name) %>% 
              left_join(info, by=c('spp_Network' = 'Network_Name')) %>% 
              rename(Module_spp = Module,
                     spp = Name) %>% 
              dplyr::select(cell, spp, Module_cell, Module_spp)
            
            
            # Check if any column has NA
            if (any(is.na(db_links))==T){
              print('There can be no unassigned nodes')
              break()
            }
            
            ## 
            ## Calculate Metrics                     ####
            ## ::::::::::::::::::::::::::::::::::::::::::
            
            if (ret %in% c('db_Metrics', 'all')) {
              
              ### 
              ### Regional Species Pool (RSP)        ####
              ### :::::::::::::::::::::::::::::::::::::::
              
              db_RSP <- db_links %>% 
                # To store only RSP of modules with grid cells
                filter(Module_cell == Module_spp) %>% 
                # To calculate RSP
                distinct(spp, Module_spp, Module_cell) %>% 
                group_by(Module_cell) %>% 
                summarise(RSP = n())
              
              ### 
              ### Species Richness in grid cell        ####
              ### :::::::::::::::::::::::::::::::::::::::
              
              # Total species richness
              db_S_cell <- db_links %>% 
                group_by(cell) %>% 
                summarise(S_Total = n()) 
              
              # Species richness of regional species
              db_S_reg <- db_links %>% 
                filter(Module_cell==Module_spp) %>% 
                group_by(cell) %>% 
                summarise(S_reg = n()) 
              
              # Species richness of non-regional species
              db_S_non_reg <- db_links %>% 
                filter(Module_cell!=Module_spp) %>% 
                group_by(cell) %>% 
                summarise(S_non_reg = n()) 
              
              # Data.frame with the three metrics of species richness (total, regional and non-regional)
              db_S_all <- db_S_cell %>% 
                left_join(db_S_reg, by='cell') %>% 
                mutate(S_reg = if_else(is.na(S_reg), 0, S_reg)) %>% 
                left_join(db_S_non_reg, by='cell') %>% 
                mutate(S_non_reg = if_else(is.na(S_non_reg), 0, S_non_reg)) 
              
              # Check if the calculation is right
              if( sum(db_S_all$S_Total != db_S_all$S_reg+db_S_all$S_non_reg)!=0){print('Error when calculating richness'); break()}
              
              
              ### 
              ### Occupancy in Bioregion             ####
              ### :::::::::::::::::::::::::::::::::::::::
              
              db_Occ <- db_links %>% 
                # To store only RSP of modules with grid cells
                filter(Module_cell == Module_spp) %>% 
                # To calculate occ in their module
                group_by(spp) %>% 
                summarise(Occ_Bior = n())
              
              ### 
              ### Total distribution range             ####
              ### :::::::::::::::::::::::::::::::::::::::::
              
              db_total_range <- db_links %>% 
                group_by(spp) %>% 
                summarise(Total_Range = n())
              
              ### 
              ### db to calculate metrics                       ####
              ### ::::::::::::::::::::::::::::::::::::::::::::::::::
              db_all <- db_links %>% 
                # To check values
                unite (ID, 'cell', 'spp', sep = '_', remove = F) %>% 
                # Merge
                left_join(db_RSP, by=c('Module_cell')) %>% 
                left_join(db_S_all, by=c('cell')) %>% 
                left_join(db_Occ, by=c('spp')) %>% 
                left_join(db_total_range, by=c('spp'))  
              
              #### 
              #### (1) Overlap of biota                     ####
              #### :::::::::::::::::::::::::::::::::::::::::
              
              db_Overlap <- db_all %>% 
                mutate(Overlap = S_non_reg / S_Total) %>% 
                distinct(cell, Overlap)
              
              #### 
              #### (2) Endemicity                           ####
              #### :::::::::::::::::::::::::::::::::::::::::
              
              db_Endemicity <- db_all %>% 
                # Calculate Endemicity per species as the number of occurrences in their bioregion / total occurrences
                mutate(Endemicity_spp = Occ_Bior/Total_Range) %>% 
                # Calculate mean value per grid cell (only for regional spp)
                filter(Module_spp == Module_cell) %>% 
                distinct(cell, spp, .keep_all = T) %>% 
                group_by(cell) %>% 
                mutate(Endemicity_cell = median(Endemicity_spp)) %>% 
                ungroup() %>% 
                distinct(cell, Endemicity_cell)
              
              #### 
              #### (3) Relative richness (z-score)              ####
              #### :::::::::::::::::::::::::::::::::::::::::::::
              
              db_Relative_S <- db_all %>% 
                distinct(cell, Module_cell, S_reg, S_non_reg) %>% 
                group_by(Module_cell) %>% 
                # Mean reg and non-reg
                mutate(S_reg_Mean     = mean(S_reg),
                       S_non_reg_Mean = mean(S_non_reg)) %>% 
                # SD reg and non-reg
                mutate(S_reg_SD     = sd(S_reg),
                       S_non_reg_SD = sd(S_non_reg)) %>% 
                # z reg and non-reg
                mutate(S_reg_z     = (S_reg     - S_reg_Mean)     / S_reg_SD,
                       S_non_reg_z = (S_non_reg - S_non_reg_Mean) / S_non_reg_SD) %>% 
                # centered reg and non-reg
                mutate(S_reg_centred     = (S_reg-S_reg_Mean),
                       S_non_reg_centred = (S_non_reg-S_non_reg_Mean)) %>% 
                ungroup() %>% 
                dplyr::select(cell, S_reg_z, S_non_reg_z, S_reg_centred, S_non_reg_centred)
              
              
              #### 
              #### (4) Relative Occupancy (z-score)              ####
              #### ::::::::::::::::::::::::::::::::::::::::::::::::::
              
              db_Relative_Occ <- db_all %>% 
                # To avoid modules with only species
                filter(Module_cell==Module_spp) %>% 
                # One value per species
                distinct(spp, Module_spp, Occ_Bior) %>% 
                # Relative Richness of Regional Species per grid cell
                group_by(Module_spp) %>% 
                mutate(Occ_Bior_Mean_spp = mean(Occ_Bior)) %>% 
                mutate(Occ_Bior_SD_spp = sd(Occ_Bior)) %>% 
                mutate(Occ_Bior_z_spp = (Occ_Bior-Occ_Bior_Mean_spp)/Occ_Bior_SD_spp) %>% 
                mutate(Occ_Bior_centred = (Occ_Bior-Occ_Bior_Mean_spp)) %>% 
                ungroup() %>% 
                # Ensure one value per species
                distinct(spp, Module_spp, Occ_Bior, Occ_Bior_Mean_spp, Occ_Bior_SD_spp, Occ_Bior_z_spp, Occ_Bior_centred) %>% 
                # Till here the code only shows data for species, so merge values to db_all 
                # and calculate the median of the occurrence values per grid cells
                full_join(db_all, by=c('spp', 'Occ_Bior', 'Module_spp')) %>% 
                # Median of regional species per grid cell
                filter(Module_cell==Module_spp) %>% 
                group_by(cell) %>% 
                summarise(Occ_Bior_z_cell = median(Occ_Bior_z_spp))
              
              
              ## 
              ## Merge all metrics                        ####
              ## :::::::::::::::::::::::::::::::::::::::::::::
              
              db_Metrics <- db_Overlap %>% 
                full_join(db_Endemicity, by='cell') %>% 
                full_join(db_Relative_S, by='cell') %>% 
                full_join(db_Relative_Occ, by='cell') %>% 
                dplyr::select(cell, Overlap, Endemicity_cell, Occ_Bior_z_cell, S_reg_z, S_non_reg_z, S_reg_centred, S_non_reg_centred) %>% 
                left_join(db_all[,c('cell', 'Module_cell', 'RSP')], by='cell') %>% 
                rename(Module = Module_cell) %>% 
                rename(Rel_Species_Richness = S_reg_z, 
                       Biota_Overlap = Overlap, 
                       Endemicity = Endemicity_cell, 
                       Rel_Occupancy = Occ_Bior_z_cell) %>% 
                distinct()
              db_Metrics 
              
              
              if(show.add.metrics==F){
                db_Metrics <- db_Metrics %>%  
                  dplyr::select(cell, Module, S_reg_z, Overlap, Endemicity_cell, Occ_Bior_z_cell) %>% 
                  rename(Rel_Species_Richness = S_reg_z, 
                         Biota_Overlap = Overlap, 
                         Endemicity = Endemicity_cell, 
                         Rel_Occupancy = Occ_Bior_z_cell)
              }
            }
            
            
            
  # 
  # Indicative species
  # ___________________________________________
            
            # Extent bioregion
            db_extent_bioregion <- db_all %>% 
              distinct(cell, Module_cell) %>% 
              count(Module_cell)
            
            # Calculate proportion spp in Bior 
            db_IndVal <- db_all %>% 
              # pre-select desired variables
              distinct(spp, Module_spp, Occ_Bior, Total_Range) %>% 
              # Merge extent
              left_join(db_extent_bioregion, by=c('Module_spp'='Module_cell')) %>% 
              rename(Extent_Bior = n) %>% 
              # Calculate Relative Occ
              mutate(Rel_Occ = Occ_Bior/Extent_Bior) %>% 
              # Calculate Relative Endemicity
              mutate(Endem = Occ_Bior/Total_Range) %>% 
              # Calculate IndVal
              mutate(IndVal = Rel_Occ * Endem) %>% 
              # select desired variables
              dplyr::select(spp, Module_spp, IndVal)
            
            
            
            if (ret %in% c('links')) { return(db_links) } 
            if (ret == 'db_Metrics') { return(db_Metrics) } 
            if (ret == 'info') { return(info) } 
            if (ret == 'all') { 
              list_return <- list(db_links, info, db_Metrics, db_IndVal)
              names(list_return) <- c('links', 'Infomap', 'db_Metrics', 'db_IndVal')
              return(list_return)
            }
          }
          