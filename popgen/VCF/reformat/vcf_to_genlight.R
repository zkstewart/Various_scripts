# vcf_to_genlight.R

library(vcfR)
library(argparser)

# Parse arguments
p <- arg_parser("Convert a VCF into genlight format for use by other R popgen software")

p <- add_argument(p, "v", help="VCF file input", type = "character")
p <- add_argument(p, "o", help="Output prefix for genlight file; will append '.gl.rds' as the file suffix automatically", type = "character")

args <- parse_args(p)

# Correct the output file name if needed
if (! endsWith(args$v, ".gl.rds"))
{
  args$o = paste0(args$o, ".gl.rds")
}

# Validate file paths
if (! file.exists(args$v))
{
  stop(paste0(
    "Input VCF file '", args$v, "' does not exist or is not a file!"
  ))
}
if (file.exists(args$o))
{
  stop(paste0(
    "Output VCF file '", args$o, "' already exists and will not be overwritten!"
  ))
}

# Load VCF as a vcfR object and convert to genlight
vcf <- read.vcfR(args$v, verbose = FALSE)
vcfgl <- vcfR2genlight(vcf)

# Save output object
saveRDS(vcfgl, file = args$o)
