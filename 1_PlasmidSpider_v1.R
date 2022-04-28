# This scripts makes the table of all 9,172 plasmid-associated genes across
# 4,445 bacterial samples. It takes forever to run, FYI.
# 
# You need "stringInput.csv" in your local folder. Get it from Figshare doi:
# https://doi.org/10.6084/m9.figshare.19674027.v1

# install.packages("genbankr")  # update.packages(ask=F)
library(genbankr)
# install.packages("igraph") # v1.0.0
library(igraph)
# BiocManager::install(c("STRINGdb"), ask=F) # , version="3.8") #  
# see https://bioconductor.org/packages/devel/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf
#STRINGdb$methods()              # To list all the methods available.
#STRINGdb$help("get_graph")      # To visualize their documentation.
library(STRINGdb)     # activate the STRINGdb library # eg species_id=9606 is Homo sapiens 
#install.packages("VennDiagram")
library(VennDiagram)
#install.packages("rentrez")
library(rentrez)
#install.packages("tidyverse")
library(tidyverse)
# install.packages("dplyr")
library(dplyr)
# install.packages("stringr")
library(stringr)
#BiocManager::install("GenomicRanges")
library(GenomicRanges)

input <- data.frame(read.csv("stringInput.csv")) # species, proteins, String_ID
input$Species <- gsub("/",  replacement="-", input$Species)
input$interactions <- rep(0, length(input$Species))

#get plasmid accession IDs from GenBank file
#   STRING_id   compact_name
#      511145   Escherichia coli K12 MG1655 # done
#      316385   Escherichia coli K12 DH10B # K12 done already # ignore
#      316407   Escherichia coli K12 W3110 # K12 done already # ignore
#      155864   Escherichia coli O157H7
#      199310   Escherichia coli CFT073
#      362663   Escherichia coli 536
#      469008   Escherichia coli BL21DE3
#      481805   Escherichia coli ATCC8739
#         573   Klebsiella pneumoniae ?? 
#         571   Klebsiella oxytoca ?? 
#     1028307   Enterobacter aerogenes KCTC2190
speciesID = c(511145, 155864, 199310, 362663, 469008, 481805, 573, 571, 1028307)
id_name = c("K12", "O157H7", "CFT073", "536", "BL21DE3", "ATCC8739",
            "Klebsiella_p", "Klebsiella_o", "Enterobacter_a")
table_plasmid <- data.frame(Plasmid=character(), Chrom=character(), Genes=integer(), 
                            Interacting_genes=integer(), Interactions=integer(),
                            Total_Genes=integer(), Total_Interactions=integer())

#--- function1 --- function for plasmids 
getplasmidgenes <- function(GBAccession_ID){ #  function plotstuffs
    plasmid <- readGenBank(GBAccession(GBAccession_ID), partial=TRUE)
    # cds(plasmid)$product # gene long names
    # cds(plasmid)$gene    # gene short names
    plasmid_genes <- unique(sort(tolower(as.vector(na.omit( cds(plasmid)$gene))))) 
    plasmid_genes <- gsub("[()]",  replacement="", plasmid_genes) # remove ( + )
    # unique protein-coding genes all uppercase sorted #  length(plasmid_genes)
    return(tolower(plasmid_genes))
} # end function 

accessID = c("EU935739", #  pEK499
             "EU935740", #  pEK204
             "EU935738", #  pEK516
             "NZ_HG941719.1", # pEC958 2
            "CP009231.1", # pCA14 # 5
             "NC_013655.1", # pSE15
             "NC_020271.1",  # pJE186 
             "NZ_CP072324.1",  # pOXA-48 EC-JS426 plasmid pOXA-48_EC-JS426
             "JN626286" # pOXA-48  Kleb pneumo
)  # end vector # check in https://www.ncbi.nlm.nih.gov/nuccore

#accessID = c("EU935740")
#input id for species ("K12", "O157H7", "CFT073"...) (from input.csv)
o <- which(input$ID %in% speciesID)

for (jj in 1:length(accessID)){
  access_ID = accessID[jj] # select plasmid 
  
########### get chromosomal genes  # 
 for (k in 1:4445){           # change k  1:4445 -> 4423 valid chromosome

     if((k != 22)&&(k != 25)&&(k != 51)&&(k != 99)&&(k != 235)&&(k != 246)
        &&(k !=266)&&(k !=345)&&(k !=349)&&(k !=376) # 10 total 
     &&(k !=403)&&(k !=407)&&(k !=416)&&(k !=419)&&(k !=485)&&(k !=529)&&(k !=546) # 7
     &&(k !=559)&&(k !=622)&&(k !=691)&&(k != 1502)&&(k != 1873)&&(k != 1876)
     &&(k != 1916)&&(k != 2275)&&(k !=4377)&&(k !=4378)# 5 => 22 not work
     ){  #  if
    
  string_db <- STRINGdb$new(version="11", species=input$ID[k], # new STRINGdb object
                            score_threshold=400, input_directory="")
  # string_db$proteins$preferred_name contains the list of genes, eg "DR97_1"

  mapped   <- string_db$map(data.frame(gene_name = string_db$proteins$preferred_name),
                           'gene_name', removeUnmappedRows=T) 
  # mapped is a table with gene_name and STRING_id
  
  links_all1 = string_db$get_interactions(string_db$mp(string_db$proteins$preferred_name)) 
  links_all <- links_all1[!duplicated(links_all1[,c('from','to')]),]
  # links_all is a table with the STRING_ids of the 
  #                                  connected proteins (cols 1 & 2) & their score (col 3)
  # links_all does not contain duplicates (unlike links_all1)
  
  # Now create a table for all genes with STRING_ids and combined score
  from_vector <- c()			# empty vector for gene names rather than STRING_ids
  to_vector <- c()				# empty vector for gene names rather than STRING_ids  
  for (i in 1:length(links_all$from)){		# for each STRING_id, extract the gene name
    from_vector[i] <- subset(mapped, STRING_id == links_all$from[i])$gene_name
    to_vector[i] <- subset(mapped, STRING_id == links_all$to[i])$gene_name }
  links_all$from_ID <- tolower(from_vector)
  links_all$to_ID <- tolower(to_vector) 
  # str(links_all) # now: from, to, combined_score, from_ID, to_ID
  #write.csv(links_all, paste("OUTPUT/", input$Species[k], "_chrom_data.csv", sep=""))
  input$interactions[k] <- length(links_all$combined_score)
  
  ######## get plasmid's gene interactions # select plasmid-associated interactions
  p1 <- getplasmidgenes(access_ID)
  str(p1) # 84 genes, single vector of names
  mapped_p <- data.frame() # mapped_pEK499$STRING_id # has from/to data
  for (m in 1:length(p1)){	# get plasmid genes in K12
  mapped_p <- rbind(mapped_p, data.frame(subset(links_all, (to_ID==p1[m]) | (from_ID==p1[m])))) }
  str(mapped_p) # 684 for K12 and pEK499
  #write.csv(mapped_p, paste("OUTPUT/", input$Species[k], "_", access_ID, ".csv", sep=""))

  # add row to plasmid table
  table_plasmid[length(table_plasmid$Plasmid)+1,] <- c(access_ID, input$Species[k], 
                                                       length(p1), # number of plasmid genes
                length(unique(c(mapped_p$from_ID, mapped_p$to_ID))), length(mapped_p$from),
                length(mapped[,1]),input$interactions[k])
  str(table_plasmid)

     } # end if
 } # end for
} # end for plasmids
# str(table_plasmid) # check

write.csv( table_plasmid, "plasmids_species.csv")

#--- function2 --- get plasmid information by species name as input #modified
getchrom.plasmid <- function(chrom_name){ 
  dt <- subset(table_plasmid,Chrom == chrom_name)
  names(dt)[names(dt) == "Plasmid"] <- chrom_name
  dt <- subset(dt, select = -c(Chrom))
  dt } 

for (j in 1:length(id_name)){  #e.g. Five plasmids information for six species
   id_name[j] = unique(table_plasmid$Chrom)[j]
   print(getchrom.plasmid(id_name[j]))  } # end loop

#--------------------------------------------

#--- function3 --- get intersection of a chromosome across all plasmids 
jaccard <- function(a) {  # a is dataframe (e.g.only data of E.coli K12)
 #one <- subset(table_plasmid,Chrom == "Escherichia_coli_str._K-12_substr._MG1655_")
  intersection <- as.integer(a$Interacting_genes)
  chrom.total.genes <- as.integer(a$Total_Genes)
  plasmid.genes <- as.integer(a$Genes)
  union = chrom.total.genes + plasmid.genes - intersection
  return (intersection/union) #value of intersection between chromosome and plasmids 
}                             # (e.g.K12 vs pEK499, K12 vs pEK204, K12 vs pEK516)

chrom <- unique(table_plasmid$Chrom)
#for (i in 1:length(chrom)) {
#  print(round(jaccard(subset(table_plasmid,Chrom == chrom[i])),3)) }# end for

#-----add jaccard index of each chrom and plasmid into dataframe
table_jaccard <- data.frame(Chrom=character()) #create a dataframe with an empty column
    #create dataframe with column number of length(plasmid),
    #change column name to accessID
table_jaccard <- cbind(table_jaccard, 
        setNames(data.frame(matrix(ncol = length(accessID), nrow = 0)),accessID))
#bind 1st column with other columns with Plasmid ID
for (i in 1:length(chrom)) { #put in jaccard index
  table_jaccard[length(table_jaccard$Chrom)+1,] <- c(chrom[i],
    format(round(jaccard(subset(table_plasmid,Chrom == chrom[i])),3), nsmall=3)) }
write.csv( table_jaccard, "plasmids_species_jac.csv")

#-----table of interaction#
table_interaction <- c()
table_interaction <- data.frame(Chrom=character()) # dataframe with  empty column
m <- setNames(data.frame(matrix(ncol = length(accessID), nrow = 0)),accessID)
          #create dataframe with column number of length(plasmid),
          #change column name to accessID
table_interaction <- cbind(table_interaction,m)
#bind the first column with other columns with Plasmid ID
for (i in 1:length(id_name)) { #add in interactions
  id_name[i] = unique(table_plasmid$Chrom)[i]
  table_interaction[length(table_interaction$Chrom)+1,] <- c(id_name[i],
                                            getchrom.plasmid(id_name[i])[,4]) }

for (i in 1:length(id_name)) { #add in interactions
  id_name[i] = unique(table_plasmid$Chrom)[i]
  table_interaction[length(table_interaction$Chrom)+1,] <- c(id_name[i],
                                  getchrom.plasmid(id_name[i])$Interactions) }
write.csv( table_interaction, "plasmids_species_interactions.csv")

#-----table of interacting genes
table_interacting_genes <- data.frame(Chrom=character()) # dataframe with an empty column
m <- setNames(data.frame(matrix(ncol = length(accessID), nrow = 0)),accessID)
       #create dataframe with column number of length(plasmid),
       #change column name to accessID
table_interacting_genes <- cbind(table_interacting_genes,m) 
          # bind the first column with other columns with Plasmid ID
  
for (i in 1:length(id_name)) { #add in interacting genes
  id_name[i] = unique(table_plasmid$Chrom)[i]
  table_interacting_genes[length(table_interacting_genes$Chrom)+1,] <- c(id_name[i],
                                                getchrom.plasmid(id_name[i])[,3])
}

for (i in 1:length(id_name)) { #add in interacting genes
  id_name[i] = unique(table_plasmid$Chrom)[i]
  table_interacting_genes[length(table_interacting_genes$Chrom)+1,] <- c(id_name[i],
                                    getchrom.plasmid(id_name[i])$Interacting_genes)
}
write.csv( table_interacting_genes, "plasmids_species_genes.csv")


#--- function4 --- get plasmid metadata: (1) plasmid nickname, (2) species, (3) strain, (4) length, (5) number of genes
#table <- setDT(table, keep.rownames = TRUE)[] #convert rowname to the first column

getdata <- function(GBAccession_ID){
    plasmid <- readGenBank(GBAccession(GBAccession_ID), partial=TRUE) #read in information from plasmid accession number
    w <- elementMetadata(sources(plasmid)) #get metadata(type,organism,mol_type,strain,db_xref,plasmid,loctype,temp_grouping_id)
    name <- w$plasmid #plasmid name
    species <- w$organism #species
    species <- gsub(" ",  replacement="_", species)
    strain <- w$strain #strain
    genes <- unique(sort(tolower(as.vector(na.omit(cds(plasmid)$gene))))) #remove redundant genes
    genes <- gsub("[()]",  replacement="", genes) 
    genes <- length(genes)    #number of unique plasmid genes
    info <- seqinfo(plasmid)  #get plasmid information
    df_plasmid <- as.data.frame(info)
    length <- df_plasmid$seqlengths #plasmid length
    metadata <- c(name,species,strain,genes,length) #combine all data
}

for (i in 1:length(accessID)) { #show metadata of plasmids
  print(getdata(accessID[i])) }

#---- table shows with plasmid metadata 
mt <- c("name","species","strain","genes","length")
table_metadata <- data.frame(Accession=character()) #create a dataframe with an empty column
m <- setNames(data.frame(matrix(ncol = length(mt), nrow = 0)),mt)
      #create dataframe with column number of length(mt),
      #change colnames with variables in mt
table_metadata <- cbind(table_metadata,m) #combine the first column with other variables
  
for (i in 1:length(accessID)) { #get plasmids metadata
  table_metadata[length(table_metadata$Accession)+1,] <- c(accessID[i], getdata(accessID[i]))
} #  end for

write.csv( t(table_metadata), "plasmids_species_metadata.csv")

quit()
y
