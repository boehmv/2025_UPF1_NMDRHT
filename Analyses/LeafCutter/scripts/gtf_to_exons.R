#!/usr/bin/env Rscript

require(dplyr)
require(data.table)

args=commandArgs(trailingOnly = T)

if (length(args)<2) {
  stop("Usage is: ./gtf_to_exons.R input.gtf.gz output.txt.gz")
}

cat("Reading in ",args[1],"\n")
gtf=fread(cmd=paste("zcat <", args[1]), data.table = F, col.names=c("chr","source","feature","start","end","a","strand","b","dat"))

cat("Processing...\n")
gtf = gtf %>% filter( feature=="exon" )

gn_where=regexpr("gene_name \"[^ ]+\"" , gtf$dat) # find gene_names in dat
gn_where=gn_where + 11 # ignore "gene_name" label
attr(gn_where,"match.length")=attr(gn_where,"match.length") - 11- 1 # cutoff trailing quote mark

gtf$gene_name=regmatches(gtf$dat, gn_where )

if( any( gtf$gene_name== "" ) ){
  cat("Warning: there are empty 'gene_name' attributes, using 'gene_id' for them\n")
  gi_where=regexpr("gene_id \"[^ ]+\"" , gtf$dat) # find gene_ids in dat
  gi_where=gi_where + 9 # ignore "gene_id" label
  attr(gi_where,"match.length")=attr(gi_where,"match.length") - 9- 1 # cutoff trailing quote mark
  gtf$gene_id=regmatches(gtf$dat, gi_where )
  gtf$gene_name[ gtf$gene_name == "" ] <- gtf$gene_id[ gtf$gene_name == "" ]
  gtf=select (gtf,-gene_id )
}

# Get Ensembl_ID
gid_where=regexpr("gene_id \"[^ ]+\"" , gtf$dat) # find gene_ids in dat
gid_where=gid_where + 9 # ignore "gene_name" label
attr(gid_where,"match.length")=attr(gid_where,"match.length") - 9- 1 # cutoff trailing quote mark
gtf$gene_id=regmatches(gtf$dat, gid_where ) 

#gtf$gene=foreach(s=strsplit(gtf$dat," "), .combine=c) %dopar% { s[which(s=="gene_name")+1] }
#gtf$gene=substr(gtf$gene, 1, nchar(gtf$gene)-1)

gtf = gtf %>% select( chr, start, end, strand, gene_name, gene_id ) %>% distinct()

cat("Saving exons to ",args[2],"\n")
gz=gzfile(args[2],"w")
write.table(gtf, gz, row.names = F, quote=F, sep="\t")
close(gz)
