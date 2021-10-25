require(ggplot2)
require(tidyr)
require(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", required=TRUE)
parser$add_argument("-o", "--output", required=TRUE)
parser$add_argument("-m", "--maximum-length", type="integer", required=TRUE)
args <- parser$parse_args()

data <- read.table(args$input, header=T)
reshaped <- data %>% tidyr::gather(method, count, graphaligner, simplecount)
p <- ggplot(reshaped, aes(count, fill=strand)) + 
  geom_histogram(binwidth=1, position="identity", alpha=0.5) + 
  xlim(0, args$maximum_length) +
  theme_bw(base_size=18) + 
  facet_grid(. ~ method, scales="free")
  
ggsave(args$output, p, height=10, width=20)
