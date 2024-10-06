encoded_data_for_predictor <- function(df){
  
  proteome_encoded<-df
  
  ### simplify the predictor variables. 
  ## simple struc --> 6 levels in total: loop, ntermini, ctermini, alphahelix, strand, 310helix
  ## structure before + after = 6 levels each (start/end, ntermini/ctermini, loop, alphahelix, strand, 310helix)
  ## align_to_centre (position) should be 0, 1, 2, 3 and 4+ for the rest
  ## aling_to_centre for termini should be a sepereate column 
  ## add length but simplified: short, medium, long. 
  
  ## recode align_to_centre for termini
  proteome_encoded$align_to_center_termini<-""
  
  for (i in 1:nrow(proteome_encoded)){
    if (proteome_encoded[i,]$simple_struc == "ntermini"){
      if (proteome_encoded[i,]$align_to_center %in% c("-1","-2","-3")){
        proteome_encoded[i,]$align_to_center_termini<-proteome_encoded[i,]$align_to_center
      }else{
        proteome_encoded[i,]$align_to_center_termini<-"-4+"
      }
    }else if (proteome_encoded[i,]$simple_struc == "ctermini"){
      if (proteome_encoded[i,]$align_to_center %in% c("1","2","3")){
        proteome_encoded[i,]$align_to_center_termini<-proteome_encoded[i,]$align_to_center
      }else{
        proteome_encoded[i,]$align_to_center_termini<-"4+"
      }
    }
  }
  
  ## recode align_to_centre for structure/loop
  for (i in 1:nrow(proteome_encoded)){
    if (proteome_encoded[i,]$align_to_center %in% c(-1:-3)){
      proteome_encoded[i,]$align_to_center<-proteome_encoded[i,]$align_to_center
    }else if(proteome_encoded[i,]$align_to_center %in% c(-4:-400)){
      proteome_encoded[i,]$align_to_center<-"-4+"
    }else if (proteome_encoded[i,]$align_to_center %in% c(1:3)){
      proteome_encoded[i,]$align_to_center<-proteome_encoded[i,]$align_to_center
    }else if(proteome_encoded[i,]$align_to_center %in% c(4:400)){
      proteome_encoded[i,]$align_to_center<-"4+"
    }
  }
  
  ## now remove the align_to_centre for n and c termini....
  for (i in 1:nrow(proteome_encoded)){
    if (proteome_encoded[i,]$simple_struc %in% c("ntermini", "ctermini")){
      proteome_encoded[i,]$align_to_center<-""
    }
  }
  
  ### simplify the lenght
  ## helix: short: 4-8, medium: 9-15, long: 16<
  ## strand short: 2-4, medium: 5-6, long: 7<
  ## 310helix: medium: all lenghts
  ## loop: short: 1-3, medium: 4-5, long: 6< 
  ## ntermini: short: 1-3, long: 4+
  ## ctermini: short: 1-3, long: 4+
  
  proteome_encoded$length<-""
  
  for (i in 1:nrow(proteome_encoded)){
    if (proteome_encoded[i,]$simple_struc == "Strand"){
      if (proteome_encoded[i,]$element_lenght %in% c(2:4)){
        proteome_encoded[i,]$length<-"short"
      }else if (proteome_encoded[i,]$element_lenght %in% c(5:6)){
        proteome_encoded[i,]$length<-"medium"
      }else{
        proteome_encoded[i,]$length<-"long"
      }
    }else if (proteome_encoded[i,]$simple_struc == "AlphaHelix"){
      if (proteome_encoded[i,]$element_lenght %in% c(4:8)){
        proteome_encoded[i,]$length<-"short"
      }else if (proteome_encoded[i,]$element_lenght %in% c(9:15)){
        proteome_encoded[i,]$length<-"medium"
      }else{
        proteome_encoded[i,]$length<-"long"
      }
    }else if (proteome_encoded[i,]$simple_struc == "loop"){
      if (proteome_encoded[i,]$element_lenght %in% c(1:3)){
        proteome_encoded[i,]$length<-"short"
      }else if (proteome_encoded[i,]$element_lenght %in% c(4:5)){
        proteome_encoded[i,]$length<-"medium"
      }else{
        proteome_encoded[i,]$length<-"long"
      }
    }else if (proteome_encoded[i,]$simple_struc == "310Helix"){
      proteome_encoded[i,]$length<-"medium"
    }else if (proteome_encoded[i,]$simple_struc %in% c("ntermini", "ctermini")){
      if (proteome_encoded[i,]$element_lenght %in% c(1:3)){
        proteome_encoded[i,]$length<-"short"
      }else{
        proteome_encoded[i,]$length<-"long"
      }
    }
  }
  
  ## change the terminiso both ntermini/ctermini are called termini
  # can be distinguished by the align_to_centre column. 
  
  for (i in 1:nrow(proteome_encoded)){
    if (proteome_encoded[i,]$simple_struc %in% c("ntermini", "ctermini")){
      proteome_encoded[i,]$simple_struc <- "termini"
    }
  }
  
  proteome_encoded<-proteome_encoded[,c("pdb_name", "Pos",
                                    "resid", "simple_struc","structure_before", "structure_after",
                                    "align_to_center", "neighbour_n", "neighbour_c", 
                                    "align_to_center_termini", "length")]
  return(list(proteome_encoded = proteome_encoded))
  
}
