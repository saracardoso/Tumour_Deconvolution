
# ---
# - CRC - NICs
# ---

# 0. Store data.csv also in data.txt
data = read.csv('./Data/bulk/CRC/NICs/data/data.csv', row.names=1)
write.table(data, './Data/bulk/CRC/NICs/data/data.txt', col.names=NA, row.names=TRUE,
            quote=FALSE, sep='\t')


# 1. Proportions of bulk data:
# 1.1. Read original file
bulk_cells_counts_hyperion = read.table('./Data/bulk/CRC/NICs/original_files/cell_count_perPhenotype_perSample.txt', header=T)

# 1.2. Get phenotypes per cell-type:
all_phenotypes = unique(bulk_cells_counts_hyperion$phenotype)
cancer_cells = c(grep('_tum', all_phenotypes, value=TRUE), 'Tum')
stromal_cells = c(grep('essels', all_phenotypes, value=TRUE), grep('ibroblasts', all_phenotypes, value=TRUE))
Bcells = 'Bcells'
cd4_tcells = c("CD39+_CD8-_Tcells", "CD4+_Tcells", "CD57+_CD8-_Tcells", "CD8-_Tcells", "Intra_CD8-_Tcells")
regulatory_tcells = c("Regulatory_Tcells", "ICOS+_Regulatory_Tcells")
cd8_tcells = c("CD8+_Tcell", "CD57+_CD8+_Tcells", "Intra_CD39+_CD8+_Tcells", "Intra_CD8+_Tcells", "Intra_GZMB+_CD8+_Tcell")
proliferative_tcells = 'Prol_Tcells'
NK_cells = c('NK_cells', "CD7+_CD3-_cells", "CD56+_D2-40+_cells")
monomacro = c("TGFb+_monocytes", "HLA-DR+_CD163+_macrophages", "HLA-DR+_macrophages", "HLA-DR+_monocytes",
              "VISTA+_monocytes", "Macrophages_undefined", "Monocytes", "CD11c+_macrophages",
              "CD45RO_undefined")
other_cells = c("Apoptotic_cells", "VISTA+_CD31+_CD38+_cells", "TGFb+_CD31+_CD38+_cells", "CD31+_CD38+_cells",
                "Granulocytes", "TGFb+_granulocytes")

# 1.3. Get cell-type counts per sample:
n_cell_types = 10
n_samples = length(unique(bulk_cells_counts_hyperion$sample))
cell_counts = as.data.frame(matrix(rep(0, n_cell_types*n_samples), nrow=n_samples))
rownames(cell_counts) = unique(bulk_cells_counts_hyperion$sample)
colnames(cell_counts) = c('Cancer cells', 'Stromal cells', 'Bcells', 'CD4 Tcells', 'Regulatory CD4 Tcells', 'CD8 Tcells',
                          'Proliferative Tcells', 'NK cells', 'Macro/mono Lineage', 'Other cells')
first_counts = bulk_cells_counts_hyperion[,c('sample', 'phenotype', 'cells.count.mean')]
for(samp in rownames(cell_counts)){
  cell_counts[samp, 'Cancer cells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%cancer_cells, 'cells.count.mean'])
  cell_counts[samp, 'Stromal cells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%stromal_cells, 'cells.count.mean'])
  cell_counts[samp, 'Bcells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%Bcells, 'cells.count.mean'])
  cell_counts[samp, 'CD4 Tcells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%cd4_tcells, 'cells.count.mean'])
  cell_counts[samp, 'Regulatory CD4 Tcells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%regulatory_tcells, 'cells.count.mean'])
  cell_counts[samp, 'CD8 Tcells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%cd8_tcells, 'cells.count.mean'])
  cell_counts[samp, 'Proliferative Tcells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%proliferative_tcells, 'cells.count.mean'])
  cell_counts[samp, 'NK cells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%NK_cells, 'cells.count.mean'])
  cell_counts[samp, 'Macro/mono Lineage'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%monomacro, 'cells.count.mean'])
  cell_counts[samp, 'Other cells'] = sum(first_counts[first_counts$sample==samp & first_counts$phenotype%in%other_cells, 'cells.count.mean'])
}
cell_counts = rbind(cell_counts, cell_counts['NIC12',])
cell_counts = rbind(cell_counts, cell_counts['NIC22',])
rownames(cell_counts)[c(22,23)] = c('NIC12nova', 'NIC22nova')
rownames(cell_counts) = gsub('C0', 'C', rownames(cell_counts))

# 1.4. Get cell-type proportions per sample:
# 1.4.1. Including other cells in the calculation:
cell_proportions = cell_counts
total_sample_counts = rowSums(cell_counts)
for(samp in names(total_sample_counts)) cell_proportions[samp, ] = cell_proportions[samp, ] / total_sample_counts[samp]
# 1.4.2. Not including other cells in the calculation:
cell_proportions_woOthers = cell_counts[,colnames(cell_counts)!='Other cells']
total_sample_counts = rowSums(cell_counts[,colnames(cell_counts)!='Other cells'])
for(samp in names(total_sample_counts)) cell_proportions_woOthers[samp, ] = cell_proportions_woOthers[samp, ] / total_sample_counts[samp]


# 1.5. Save counts and proportions in csv format:
write.csv(t(cell_counts), './Data/bulk/CRC/NICs/data/new_metadata/cell_counts.csv')
write.csv(t(cell_proportions), './Data/bulk/CRC/NICs/data/new_metadata/cell_proportions.csv')
write.csv(t(cell_proportions_woOthers), './Data/bulk/CRC/NICs/data/new_metadata/cell_proportions_woOthers.csv')



