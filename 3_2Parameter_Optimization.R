code_dir = './code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)



# ---
# - Set the Train and Test Data
# ---

# 1. Read Bulk Data:
CRC_bulk = read_bulk_data('./Data/bulk/CRC/NICs/data/data.csv', ground_truth_file='./Data/bulk/CRC/NICs/data/cell_proportions.csv')

# 2. Randomly choose train (~70% of samples = 16 samples) and test samples (~30% of samples = 7 samples):
train_samples = sample(colnames(CRC_bulk$ground_truth), 16)
test_samples = colnames(CRC_bulk$ground_truth)[is.na(match(colnames(CRC_bulk$ground_truth), train_samples))]

# 2. Create train data:
train = list()
train$data = CRC_bulk$data[, train_samples]
train$metadata = CRC_bulk$ground_truth[, train_samples]

# 3. Create test data:
test = list()
test$data = CRC_bulk$data[, test_samples]
test$metadata = CRC_bulk$ground_truth[, test_samples]

# 4. Save the data:
saveRDS(train, './Data/bulk/CRC/NICs/train_test_sets/train.Rdata')
saveRDS(test, './Data/bulk/CRC/NICs/train_test_sets/test.Rdata')

write.table(train$data, './Data/bulk/CRC/NICs/train_test_sets/train.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
write.table(test$data, './Data/bulk/CRC/NICs/train_test_sets/test.txt', col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
