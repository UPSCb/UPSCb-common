args <- commandArgs(trailingOnly=T)

if (length(args) < 3) {
    cat('Usage: Rscript plotVCFQuals.R <input.vcf> <path_to_folder> <title>\n')
    stop('Too few arguments.')
}

vcf_file = args[1]
out_path = args[2]
plot_title = args[3]

library(stringr)
library(ggplot2)
library(LSD)
library(VariantAnnotation)

# Only load quality and depth
vcfparam <- ScanVcfParam(fixed="QUAL", info="DP", geno=NA)
vcf <- readVcf(vcf_file, 'Potri', vcfparam)

dfr <- data.frame(Quality=fixed(vcf)$QUAL, Depth=info(vcf)$DP)

summary(dfr$Quality)
summary(dfr$Depth)

# Plot a histogram of the quality scores
theme_set(theme_bw(base_size=12))
cat('Plotting quality histogram\n')
pdf(file.path(out_path, paste(plot_title, '.raw_snp_quals.pdf', sep='')),
    height=4, width=6)
ggplot(dfr, aes(x = Quality)) + geom_histogram(binwidth=0.01) +
    scale_x_log10() + ggtitle(paste(plot_title, ': SNP quality score distribution',
        sep=''))
dev.off()

# Add pseudocounts
dfr$Depth = dfr$Depth + 1
dfr$Quality = dfr$Quality + 1

# Plot a comparison of quality scores vs sequencing depth
cat('Plotting quality vs depth\n')
pdf(file.path(out_path, paste(plot_title, '.raw_qual_vs_depth.pdf', sep='')),
    height=4, width=6)
ggplot(dfr, aes(x = Depth, y = Quality)) + stat_binhex() +
    scale_x_log10() + scale_y_log10() +
    xlab('Depth + 1') + ylab('Quality + 1')
dev.off()
