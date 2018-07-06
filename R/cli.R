#' @export
align10x <- function() {
    option_list = list(
        optparse::make_option(c("-s", "--sample"), type="character",
                              default=NA_character_,
                              help="sample"),
        optparse::make_option(c("-o", "--output.dir"), type="character",
                              default=NA_character_,
                              help="sample"),
        optparse::make_option(c("-r", "--reference"), type="character",
                              default=NA_character_,
                              help="Reference genome FASTA file."),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=detectCores(),
                              help="set the number of cores")
    )
    parser = optparse::OptionParser(
                           "Rscript -e 'library(methods);gscars::align10x()' [options] [input_dirs..]",
                           description=c("Align 10X FASTQ preserving barcode information.\n"),
                           epilogue=c(
                               "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
                               "Michigan Center for Translational Pathology (c) 2018\n"),
                           option_list=option_list
                       )
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    .align10x(opt$options$sample, opt$args, opt$options$output.dir, opt$options$reference, opt$options$cores)
}

#' @export
bcl10x <- function() {
    option_list = list(
        optparse::make_option(c("-r", "--run.dir"), type="character",
                              default=NA_character_,
                              help="Illumina run directory."),
        optparse::make_option(c("-b", "--bcl.sheet"), type="character",
                              default=NA_character_,
                              help="Simple BCL Sheet"),
        optparse::make_option(c("-o", "--output.dir"), type="character",
                              default=NA_character_,
                              help="Outpute directory"),
        optparse::make_option(c("-j", "--read.1.len"), type="integer",
                              default=150,
                              help="Read 1 length"),
        optparse::make_option(c("-k", "--read.2.len"), type="integer",
                              default=150,
                              help="Read 2 length"),
        optparse::make_option(c("-i", "--read.i.len"), type="integer",
                              default=8,
                              help="Index Read length"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=detectCores(),
                              help="set the number of cores")
    )
    parser = optparse::OptionParser(
                           "Rscript -e 'library(methods);gscars::align10x()' [options] [input_dirs..]",
                           description=c("Align 10X FASTQ preserving barcode information.\n"),
                           epilogue=c(
                               "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
                               "Michigan Center for Translational Pathology (c) 2018\n"),
                           option_list=option_list
                       )
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    .bcl10x(opt$options$run.dir, opt$options$bcl.sheet, opt$options$output.dir,
            opt$options$read.1.len, opt$options$read.2.len, opt$options$read.i.len, opt$options$cores)
}
