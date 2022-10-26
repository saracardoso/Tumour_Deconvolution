# Deconvolution of bulk RNAseq colorectal cancer samples

This project contains the code developed to benchmark tumour deconvolution methods that use scRNAseq data as reference, for colorectal cancer (CRC).

To deconvolute the CRC samples, we used the scRNAseq tumour samples from this [colorectal cancer (CRC) atlas](https://github.com/saracardoso/CRC_ATLAS).


## Methods Tested

| Method          | Language | 		                                            |
|-----------------|----------|------------------------------------------------------|
| AutoGeneS       | python   | [Link](https://github.com/theislab/AutoGeneS)        |
| BisqueRNA       | R        | [Link](https://github.com/cran/BisqueRNA)            |
| BSeqSC          | R        | [Link](https://github.com/shenorrLab/bseqsc)         |
| CIBERSORTx      | web      | [Link](https://cibersortx.stanford.edu/)             |
| DigitalDLSorter | R        | [Link](https://github.com/diegommcc/digitalDLSorteR) |
| DWLS            | R        | [Link](https://github.com/dtsoucas/DWLS)             |
| MOMF            | R        | [Link](https://github.com/sqsun/MOMF)                |
| MuSiC           | R        | [Link](https://github.com/xuranw/MuSiC)              |
| Scaden          | python   | [Link](https://github.com/KevinMenden/scaden)        |
| SCDC            | R        | [Link](https://github.com/meichendong/SCDC)          |


## Main Results

+ Overall, *CIBERSORTx*, *DigitalDLSorter* and *Scaden* are the best methods;

+ The best methods to predict each cell-type are:

	+ Cancer cells: *Scaden*
	
	+ Stromal cells: *BisqueRNA*

	+ Macro/mono lineage cells: *BSeqSC*

	+ B-cells: *Scaden*

	+ CD4 T-cells: *CIBERSORTx*

	+ Regulatory CD4 T-cells: *Scaden*

	+ CD8 T-cells: *DigitalDLSorter*

	+ Proliferative T-cells: *MuSiC_wGrouping*

	+ NK cells: *Scaden*

+ RNA content bias correction did not improve predictions


<!-- ## 4. How to reference this work -->
