marviddata <- "U:\\Mareano\\VIDEOLAB\\VIDEO DATA\\200m_scale_species_by_sample\\marvid_data"

gisdata <- "U:\\Mareano\\VIDEOLAB\\VIDEO DATA\\200m_scale_species_by_sample\\comparison_study"

ftp_url <- "ftp://ftp2.ngu.no/toGeno/EnvLayersNiN2022/all"

ftp_filepath <- paste0(strsplit(ftp_url, "\\//")[[1]][1], "//havforsk:Havf8776@", 
                       strsplit(ftp_url, "\\//")[[1]][2])#insert the credentials into the url

PathTaxonary <- "U:\\Mareano\\VIDEOLAB\\VIDEO DATA\\200m_scale_species_by_sample\\Data_Delivery_2024"