
# --------------------------------------------------------
# --- Get correction vector from single-cell reference ---
# --------------------------------------------------------

get_rnaCorrection_vector = function(sc.reference, cellType.var){
  
  # Get total number of counts per cell:
  total_genes_cell = colSums(sc.reference$data)
  
  # Get cell-types:
  unique_cts = unique(sc.reference$metadata[, cellType.var])
  
  # Get correction vector:
  correction_vec = rep(0, length(unique_cts))
  names(correction_vec) = unique_cts
  for(ct in unique_cts){
    ct_cells = colnames(sc.reference$metadata)[sc.reference$metadata[,cellType.var]==ct]
    correction_vec[ct] = mean(total_genes_cell[ct_cells])
  }
  return(correction_vec)
}




# -------------------------------------
# --- Correct estimated proportions ---
# -------------------------------------

correct_fractions_mRNABias = function(proportions_matrix, correction_vector){
  c_matrix = diag(correction_vector[rownames(proportions_matrix)])
  c_matrix_inv = matlib::Inverse(c_matrix)
  
  num = c_matrix_inv %*% as.matrix(proportions_matrix)
  
  corrected_fractions = round(num / colSums(num), 3)
  rownames(corrected_fractions) = rownames(proportions_matrix)
  colnames(corrected_fractions) = colnames(proportions_matrix)
  
  return(corrected_fractions)
}





# ----------------------------------------------------------
# --- Correct reference, prior to estimating proportions ---
# ----------------------------------------------------------

correct_reference_mRNABias = function(sig_matrix, correction_vector){
  sig_matrix = t(t(sig_matrix)*correction_vector[colnames(sig_matrix)])
  return(sig_matrix)
}
