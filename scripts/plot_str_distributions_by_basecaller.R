require(ggplot2)
require(tidyr)
require(argparse)

load <- function(filename) {
    basecaller <- strsplit(filename, ".", fixed=T)[[1]][2]
    d <- read.table(filename, header=T)
    d$basecaller <- basecaller
    return(d)
}

parser <- ArgumentParser()
parser$add_argument("-o", "--output", required=TRUE)
parser$add_argument("-m", "--maximum-length", type="integer", required=TRUE)
parser$add_argument("input", nargs='+')
args <- parser$parse_args()

# load data
data <- do.call(rbind, lapply(args$input, load))

reshaped <- data %>% tidyr::gather(method, count, graphaligner, simplecount, tg, strique)
bw = args$maximum_length / 50
p <- ggplot(reshaped, aes(count, fill=strand)) + 
  geom_histogram(binwidth=bw, position="identity", alpha=0.5) + 
  xlim(0, args$maximum_length) +
  theme_bw(base_size=18) + 
  facet_grid(basecaller ~ method, scales="free")
  
ggsave(args$output, p, height=20, width=20)
