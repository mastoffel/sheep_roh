library(readr)
library(data.table)

correct_vcf <- function(vcf_file) {
      
meta <- readr::read_delim(vcf_file, delim = "\t", n_max = 5, 
                          col_names = FALSE)[[1]]
# read genotypes, data.table does this correctly
# i.e. jumping over the first 5 lines
gt <- data.table::fread(vcf_file, skip = 5)
# correct ref and alt columns
gt$REF <- 0
gt$ALT <- 1
# write out meta 
readr::write_lines(meta, path = vcf_file)
# append genotypes to make file complete
data.table::fwrite(gt, file = vcf_file, append = TRUE, 
       col.names = TRUE, sep = "\t")

}
