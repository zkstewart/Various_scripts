library(argparser)

# Establish parser
p <- arg_parser("Plot PLINK2 missingness data")

# Add arguments
p <- add_argument(p, "s", help="smiss file", type = "character")
p <- add_argument(p, "v", help="vmiss file", type = "character")
p <- add_argument(p, "o", help="output prefix", type = "character")

args <- parse_args(p)

# Load data
samplemiss <- read.table(file=args$s, header=TRUE)
variantmiss <- read.table(file=args$v, header=TRUE)

# Output sample missingness plot
pdf(paste0(args$o, "_sample_miss.pdf")) #indicates pdf format and gives title to file
hist(samplemiss[,6],main="Histogram sample missingness") #selects column 6, names header of file
dev.off()

# Output variant missingness plot
pdf(paste0(args$o, "_variant_miss.pdf"))
hist(variantmiss[,5],main="Histogram variant missingness")
dev.off()
