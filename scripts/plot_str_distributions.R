require(ggplot2)
require(tidyr)
require(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", required=TRUE)
parser$add_argument("-o", "--output", required=TRUE)
args <- parser$parse_args()

data <- read.table(args$input, header=T)
reshaped <- data %>% tidyr::gather(method, count, strique, graphaligner, strscore)
p <- ggplot(reshaped, aes(count, fill=strand)) + 
  geom_histogram(binwidth=1, position="identity", alpha=0.5) + 
  theme_bw() + 
  facet_grid(method ~ ., scales="free")
ggsave(args$output, p, height=20, width=15)