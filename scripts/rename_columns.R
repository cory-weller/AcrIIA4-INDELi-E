rename_columns <- function(df_encoded, df){
  
  one_hot_encoded_df<-as.data.frame(df_encoded)
  one_hot_encoded_df <- cbind(df$pdb_name, df$Pos, one_hot_encoded_df)
  
  colnames(one_hot_encoded_df)[c(1:2)] <- c("domain", "Pos")
  
  ## this might differ between proteins, hence we need a "if" statement
  
  # Create a lookup table 
  column_mapping <- c("align_to_center-1" = "align_to_center_min_1",
                      "align_to_center-2" = "align_to_center_min_2",
                      "align_to_center-3" = "align_to_center_min_3",
                      "align_to_center-4+" = "align_to_center_min_4plus",
                      "align_to_center_termini-1" = "align_to_center_termini_min_1",
                      "align_to_center_termini-2" = "align_to_center_termini_min_2",
                      "align_to_center_termini-3" = "align_to_center_termini_min_3",
                      "align_to_center_termini-4+" = "align_to_center_termini_min_4plus",
                      "align_to_center4+" = "align_to_center4plus",
                      "align_to_center_termini4+" = "align_to_center_termini4plus")
  
  for (i in seq_along(colnames(one_hot_encoded_df))) {
    old_col_name <- colnames(one_hot_encoded_df)[i]
    if (old_col_name %in% names(column_mapping)) {
      new_col_name <- column_mapping[old_col_name]
      colnames(one_hot_encoded_df)[i] <- new_col_name
    }
  }
  
  
  one_hot_encoded_df[,-c(1)]<-lapply(one_hot_encoded_df[,-c(1)], as.numeric)
  return(list(one_hot_encoded_df = one_hot_encoded_df))
  
  
}
  
