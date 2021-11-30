
# ------------------------------------
# --- single-cell RNAseq Reference ---
# ------------------------------------

# Read scRNAseq reference data and metadata files.
# scRNAseq_ref_path  --> string with path to the CSV file containing the single-cell data counts. The first column in the file must be
#                        the names of the genes, while the first row must be the names of the cells.
# scRNAseq_meta_path --> string with path to the CSV file containing the single-cell metadata. The first column in the file must be
#                        the names of the cells, while the first row must be the names of the metadata variables.
# Returns a list with two data.frames: counts data ($data) and metadata ($metadata)
read_scRNAseq_reference = function(scRNAseq_ref_path, scRNAseq_meta_path){
  scRef = list()
  scRef$data = read.csv(scRNAseq_ref_path, row.names=1)
  scRef$metadata = read.csv(scRNAseq_meta_path, row.names=1)
  
  return(scRef)
}





# -----------------
# --- Bulk Data ---
# -----------------

# Read Bulk data and, optionally, the cell-type proportions' ground truth.
# bulk_data_file    --> string with path to the CSV file containing the bulk data counts. The first column in the file must be
#                       the names of the genes, while the first row must be the names of the cells.
# ground_truth_file --> OPTIONAL, defaults to NULL. string with path to the .csv file with the ground_truth proportions.
read_bulk_data = function(bulk_data_file, ground_truth_file=NULL){
  bulkData = list()
  bulkData$data = read.csv(bulk_data_file, row.names=1)
  if(!is.null(ground_truth_file)) bulkData$ground_truth = read.csv(ground_truth_file, row.names=1)
  
  return(bulkData)
}


