# BSH_colours.R
## colours & codes for BSHs

BSH_name <- c("Intertidal coarse sediment","Intertidal sand and muddy sand","Intertidal mud","Intertidal mixed sediments","Intertidal seagrass beds","Subtidal coarse sediment","Subtidal sand","Subtidal mud","Subtidal mixed sediment","Maerl beds","Subtidal seagrass beds","Intertidal rock","Intertidal biogenic reefs: Sabellaria alveolata","Intertidal biogenic reefs: mussel beds","Infralittoral rock","Circalittoral rock","Subtidal biogenic reefs: Sabellaria spp.","Subtidal biogenic reefs: mussel beds","Sea caves","Saltmarsh")
BSH_code <- c("A2.1","A2.2","A2.3","A2.4","A2.61","A5.1","A5.2","A5.3","A5.4","A5.51","A5.53","A1","A2.71","SF_SH_5","A3","A4","A5.61","SF_SH_6","H8330","A2.5")

R <- c(255,255,255,195,76,255,255,229,221,230,153,255,255,255,170,0,255,255,115,56)
G <- c(128,255,189,255,255,187,255,197,255,152,214,85,127,170,255,255,127,211,115,168)
B <- c(0,76,64,76,195,153,128,115,153,0,221,0,0,0,0,255,0,127,0,0)

BSH_codes <- as_tibble(data.frame(BSH_name,BSH_code,R,G,B))
rm(BSH_name,BSH_code,R,G,B)
