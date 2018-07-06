.bcl10x <- function(run.dir, bcl.sheet, output.dir, read.1.len, read.2.len, read.i.len, cores) {
    ## input directory
    wd <- getwd()
    dir.create(output.dir, recursive=TRUE)
    output.dir <- normalizePath(output.dir)
    sample.sheet.file <- file.path(output.dir, "SampleSheet.csv")
    setwd(output.dir)
    ## create SampleSheet
    tmp <- fread(bcl.sheet)
    chromium.index <- fread(system.file("extdata/chromium-dna-sample-indexes-plate.csv", package="gscars"), header=FALSE)
    chromium.index <- melt(chromium.index, id.vars="V1")[,.(Index=V1,Barcode=value)]
    setkey(tmp, Index)
    setkey(chromium.index, Index)
    tmp <- chromium.index[tmp,.(Lane,Sample_ID=Index,Sample_Name=Library,Index=str_sub(Barcode, 1, read.i.len),
                                Sample_Project=Flowcell)]
    sample.sheet <- c("[Header]","EMFileVersion,4","[Reads]", read.1.len, read.2.len, "[Data]", capture.output(fwrite(tmp)))
    writeLines(sample.sheet, sample.sheet.file)
    ## BCL command script
    base.mask <- sprintf("Y%d,I%d,Y%d", read.1.len, read.i.len, read.2.len)
    bcl.cmd <- c(
        sprintf("--use-bases-mask=%s -p%d -r%d -w%d -R %s --output-dir=%s --sample-sheet=%s",
                base.mask, cores, 8, 8, run.dir, output.dir, sample.sheet.file),
        sprintf("--minimum-trimmed-read-length=%d", read.i.len),
        sprintf("--mask-short-adapter-reads=%d", read.i.len),
        "--create-fastq-for-index-reads",
        "--ignore-missing-positions",
        "--ignore-missing-controls",
        "--ignore-missing-filter",
        "--ignore-missing-bcls"
    )
    ##
    system2("bcl2fastq", bcl.cmd)
    setwd(wd)
}

.align10x <- function(sample, input.dirs, output.dir, reference, cores, clean=TRUE) {
    wd <- getwd()
    
    input.dirs <- normalizePath(input.dirs)
    fq.unpairs <- list(
        "1"=list.files(input.dirs, sprintf("^%s_.*R1", sample), full.names=TRUE),
        "2"=list.files(input.dirs, sprintf("^%s_.*R2", sample), full.names=TRUE)
    )
    

    run.dir <- file.path(output.dir, sample)
    dir.create(run.dir, recursive=TRUE)
    run.dir <- normalizePath(run.dir)
    setwd(run.dir)
    dir.create("tmp")
    
    { ## CPU-bound
        mclapply(names(fq.unpairs), function(read) {
            fqs <- fq.unpairs[[read]]
            system2("bash", args=c("-c", shQuote(sprintf("zcat %s | tee >(wc -l > read_%s.counts) > merge_%s.fq ",
                                                         paste(fqs, collapse=" "), read, read))))
        }, mc.cores=2)
        stopifnot(readLines("read_1.counts") == readLines("read_2.counts"))
    }

    { ## IO-bound no need to parallelize
        dir.create("chunks")
        n.lines <- as.numeric(readLines("read_1.counts"))
        n.reads <- n.lines / 4
        n.reads.per.chunk <- ceiling(n.reads/cores)
        n.lines.per.chunk <- n.reads.per.chunk * 4L
        system2("split", sprintf("-a 3 -d -l %d merge_1.fq chunks/merge_1_", n.lines.per.chunk))
        system2("split", sprintf("-a 3 -d -l %d merge_2.fq chunks/merge_2_", n.lines.per.chunk))
    }

    { ## count barcodes
        wl_fn_gz <- system.file("extdata/barcode-whitelist-4mwafeb2016-10x.txt.gz", package="gscars")
        system2("zcat", sprintf("%s > barcode-whitelist-4mwafeb2016-10x.txt", wl_fn_gz))
        countBarcodes("merge_1.fq", "barcode_counts.bin", wl_fn="barcode-whitelist-4mwafeb2016-10x.txt")
    }
    
    { ## CPU-bound
        chunk.pairs <- split(list.files("chunks", "merge", full.names = TRUE),
                             str_match(list.files("chunks", "merge"), "_([^_]*)$")[,2])
        mclapply(chunk.pairs, function(chunk) {
            fq1 <- chunk[1]
            fq2 <- chunk[2]
            ofq1 <- str_replace(fq1, "chunks/merge_", "chunks/process_")
            ofq2 <- str_replace(fq2, "chunks/merge_", "chunks/process_")
            preprocessFastq("barcode_counts.bin", fq1, fq2, ofq1, ofq2)
        }, mc.cores=cores)
    }

    { ## CPU-bound alignment
        pchunk.pairs <- split(list.files("chunks", "process", full.names = TRUE),
                              str_match(list.files("chunks", "process"), "_([^_]*)$")[,2])
        mclapply(names(pchunk.pairs), function(pchunk.name) {
            pchunk <- pchunk.pairs[[pchunk.name]]
            pfq1 <- pchunk[1]
            pfq2 <- pchunk[2]
            ofn <- file.path(dirname(pfq1), paste0(pchunk.name, ".bam"))
            system2("bwa", sprintf("mem -C -t %s %s %s %s | samtools view -bS - > %s",
                                   8, reference,
                                   pfq1, pfq2, ofn))
        }, mc.cores=ceiling(cores / 8))
    }

    { ## CPU/IO
        bams = list.files("chunks", "*.bam", full.names = TRUE)
        mclapply(bams, function(bam) {
            obam <- str_replace(bam, ".bam", "_sort.bam")
            system2("novosort", sprintf("--markduplicates --keeptags -c%s --tmpdir tmp -o %s %s", 8, obam, bam))
        }, mc.cores=ceiling(cores / 8))
    }
    system2("novosort", sprintf("-i --markduplicates --tmpdir tmp -c%d -r '@RG\tID:%s\tSM:%s' -o final.bam chunks/*sort.bam",
                                cores, sample, sample))

    ##
    system2("extract-sv-reads", sprintf("--threads %d final.bam -s splitters.bam -d discordants.bam", cores))
    system2("sambamba", sprintf("index -t%d splitters.bam", cores))
    system2("sambamba", sprintf("index -t%d discordants.bam", cores))

    if (clean) {
        system2("rm", "-rf tmp chunks merge_1.fq merge_2.fq read_1.counts read_2.counts")
        system2("rm", "barcode_counts.bin barcode-whitelist-4mwafeb2016-10x.txt")
    }
    
    setwd(wd)
}
