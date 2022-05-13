plasmids = read.csv("Plasmid_PPI_data_PMNLE_corrected.csv")
FP = read.csv("approximatePMNLE_fullPPIN.csv")
PC = read.csv("approximatePMNLE_purely_chromosomal.csv")


ccs_FP = NULL
for(j in 1:length(FP$cConnectedComponents)){
  numberCCs <- eval(parse(text=as.character(FP$cConnectedComponents[[j]])))
  ccs_FP <- c(ccs_FP, numberCCs[26])
}
min(ccs_FP[is.na(ccs_FP)==FALSE])


ccs_PC = NULL
for(j in 1:length(PC$cConnectedComponents)){
  numberCCs <- eval(parse(text=as.character(PC$cConnectedComponents[[j]])))
  ccs_PC <- c(ccs_PC, numberCCs[26])
}
min(ccs_PC[is.na(ccs_PC)==FALSE])


FP$ccs <- ccs_FP
FP$STRING.bacteria.name[is.na(FP$ccs)==FALSE&FP$ccs == 0]


PC$ccs <- ccs_PC

Edges_FP = NULL
for(j in 1:length(FP$cEdges)){
  numberEdges <- eval(parse(text=as.character(FP$cEdges[[j]])))
  Edges_FP <- c(Edges_FP, numberEdges[26])
}
min(Edges_FP[is.na(Edges_FP)==FALSE])


Edges_PC = NULL
for(j in 1:length(PC$cEdges)){
  numberEdges <- eval(parse(text=as.character(PC$cEdges[[j]])))
  Edges_PC <- c(Edges_PC, numberEdges[26])
}
min(Edges_PC[is.na(Edges_PC)==FALSE])


FP$Edges <- Edges_FP

PC$Edges <- Edges_PC



Vertices_FP = NULL
for(j in 1:length(FP$cVertices)){
  numberVertices <- eval(parse(text=as.character(FP$cVertices[[j]])))
  Vertices_FP <- c(Vertices_FP, numberVertices[26])
}
min(Vertices_FP[is.na(Vertices_FP)==FALSE])


Vertices_PC = NULL
for(j in 1:length(PC$cVertices)){
  numberVertices <- eval(parse(text=as.character(PC$cVertices[[j]])))
  Vertices_PC <- c(Vertices_PC, numberVertices[26])
}
min(Vertices_PC[is.na(Vertices_PC)==FALSE])


FP$Vertices <- Vertices_FP

PC$Vertices <- Vertices_PC



Triangles_FP = NULL
for(j in 1:length(FP$cTriangles)){
  numberTriangles <- eval(parse(text=as.character(FP$cTriangles[[j]])))
  Triangles_FP <- c(Triangles_FP, numberTriangles[26])
}
min(Triangles_FP[is.na(Triangles_FP)==FALSE])


Triangles_PC = NULL
for(j in 1:length(PC$cTriangles)){
  numberTriangles <- eval(parse(text=as.character(PC$cTriangles[[j]])))
  Triangles_PC <- c(Triangles_PC, numberTriangles[26])
}
min(Triangles_PC[is.na(Triangles_PC)==FALSE])


FP$Triangles <- Triangles_FP

PC$Triangles <- Triangles_PC


NonTrivialLoops_FP = NULL
for(j in 1:length(FP$cNonTrivialLoops)){
  numberNonTrivialLoops <- eval(parse(text=as.character(FP$cNonTrivialLoops[[j]])))
  NonTrivialLoops_FP <- c(NonTrivialLoops_FP, numberNonTrivialLoops[26])
}
min(NonTrivialLoops_FP[is.na(NonTrivialLoops_FP)==FALSE])


NonTrivialLoops_PC = NULL
for(j in 1:length(PC$cNonTrivialLoops)){
  numberNonTrivialLoops <- eval(parse(text=as.character(PC$cNonTrivialLoops[[j]])))
  NonTrivialLoops_PC <- c(NonTrivialLoops_PC, numberNonTrivialLoops[26])
}
min(NonTrivialLoops_PC[is.na(NonTrivialLoops_PC)==FALSE])


FP$NonTrivialLoops <- NonTrivialLoops_FP

PC$NonTrivialLoops <- NonTrivialLoops_PC

chromosomal_PMNLE <- NULL
chromosomal_triangles <- NULL
chromosomal_loops <- NULL
chromosomal_connected_components <- NULL
for (i in 1:nrow(plasmids)) {
  j = grep( plasmids$Sample[i], PC$STRING.bacteria.name)
  if(length(j) == 0){
        correct_PMNLE_value = NA;
        correct_number_of_triangles = NA;
        correct_number_of_loops = NA;        
        correct_connected_components = NA;
  } else {
     if(length(j) == 1){
           correct_PMNLE_value = PC$aPMNLE[j];
           correct_number_of_triangles = PC$Triangles[j];
           correct_number_of_loops = PC$NonTrivialLoops[j];
           correct_connected_components = PC$ccs[j];
     } else {
         for (k in j){
           if ( plasmids$Sample[i] == PC$STRING.bacteria.name[k]){
            correct_PMNLE_value = PC$aPMNLE[k];
            correct_number_of_triangles = PC$Triangles[k];
            correct_number_of_loops = PC$NonTrivialLoops[k];
            correct_connected_components = PC$ccs[k];
          }
         }
     }
  }
  chromosomal_PMNLE  <- c(chromosomal_PMNLE, correct_PMNLE_value)
  chromosomal_triangles <- c(chromosomal_triangles, correct_number_of_triangles)
  chromosomal_loops <- c(chromosomal_loops, correct_number_of_loops)
  chromosomal_connected_components <- c(chromosomal_connected_components, correct_connected_components)
}

plasmids$aPMNLE = chromosomal_PMNLE
plasmids$cTri = chromosomal_triangles
plasmids$loops = chromosomal_loops
plasmids$ccs = chromosomal_connected_components

PMNLE_on_full_PPIN <- NULL
triangles_on_full_PPIN <- NULL
loops_on_full_PPIN <- NULL
connected_components_on_full_PPIN <- NULL
for (i in 1:nrow(plasmids)) {
  j = grep( plasmids$Sample[i], FP$STRING.bacteria.name)
  if(length(j) == 0){
        correct_PMNLE_value = NA;
        correct_number_of_triangles = NA;
        correct_number_of_loops = NA;        
        correct_connected_components = NA;
  } else {
     if(length(j) == 1){
           correct_PMNLE_value = FP$aPMNLE[j];
           correct_number_of_triangles = FP$Triangles[j];
           correct_number_of_loops = FP$NonTrivialLoops[j];
           correct_connected_components = FP$ccs[j];
     } else {
         for (k in j){
           if ( plasmids$Sample[i] == FP$STRING.bacteria.name[k]){
            correct_PMNLE_value = FP$aPMNLE[k];
            correct_number_of_triangles = FP$Triangles[k];
            correct_number_of_loops = FP$NonTrivialLoops[k];
            correct_connected_components = FP$ccs[k];
          }
         }
     }
  }
  PMNLE_on_full_PPIN  <- c(PMNLE_on_full_PPIN, correct_PMNLE_value)
  triangles_on_full_PPIN <- c(triangles_on_full_PPIN, correct_number_of_triangles)
  loops_on_full_PPIN <- c(loops_on_full_PPIN, correct_number_of_loops)
  connected_components_on_full_PPIN <- c(connected_components_on_full_PPIN, correct_connected_components)
}

plasmids$aPMNLE_all = PMNLE_on_full_PPIN 
plasmids$cTri_all = triangles_on_full_PPIN
plasmids$loops_all = loops_on_full_PPIN
plasmids$ccs_all = connected_components_on_full_PPIN

write.csv(plasmids,"Plasmid_PPI_data_PMNLE_corrected_twice.csv")
