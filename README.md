# VCFArray

`VCFArray` is a Bioconductor package that represents VCF files as objects
 derived from the `DelayedArray` package and `DelayedArray` class. 
It converts data entries from VCF file into a DelayedArray-derived 
data structure. The backend VCF file could either be saved on-disk 
locally or remote as online resources. Data entries that could be 
extracted include the fixed data fields (REF, ALT, QUAL, FILTER), 
information field (e.g., AA, AF...), and the individual format field
 (e.g., GT, DP...). The array data generated from fixed/information 
fields are one-dimensional `VCFArray`, with the dimension being the 
length of the variants. The array data generated from individual 
"FORMAT" field are always returned with the first dimension being 
"variants" and the second dimension being "samples". This feature 
is consistent with the assay data saved in `SummarizedExperiment`, 
and makes the `VCFArray` package interoperable with other established 
_Bioconductor_ data infrastructure.

