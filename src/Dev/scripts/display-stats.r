#!/usr/bin/Rscript
cmds = commandArgs(trailingOnly=TRUE)

if (length(cmds) < 2)
{
  print('must be given 2 stats files')
  quit(status=1)
}

read_stats_file = cmds[1]
contig_stats_file = cmds[2]

# read stats, don't convert to factors
r = read.table(read_stats_file, stringsAsFactors=FALSE, header=TRUE, row.names=1)
c = read.table(contig_stats_file, stringsAsFactors=FALSE, header=TRUE)

pdf('output.pdf')

# histogram of contig lengths
hist(c$'len', breaks=10000, xlim=c(1,500), ylim=c(0,10000), main='', xlab='contig length')

# histogram of number of chunks
hist(r$'num.chunks', breaks=1000, xlim=c(0,500), main='', xlab='number of chunks per read')

# chunks per read vs read length
plot(r$'len', r$'num.chunks', xlim=c(1000,10000), ylim=c(0,500), pch='.', cex=1, xlab='read length', ylab='number of chunks')

# histogram of contig coverage
c$'covg' = c$'bp.chunks' / c$'len'
hist(c$'covg', breaks=100000, xlim=c(0,100), main='', xlab='contig coverage')

# contig coverage vs length
plot(c$'len', c$'covg', pch='.', cex=1, xlim=c(0,500), ylim=c(0,500), xlab='contig length', ylab='contig coverage')

#c$'mut.dens' = c$'num.mut' / c$'len'
# number of mutations vs contig length
plot(c$'len', c$'num.mut', pch='.', cex=1, xlab='contig length', ylab='number of mutations', xlim=c(0,500), ylim=c(0,100))

# bp in mutations vs contig length
plot(c$'len', c$'bp.mut', pch='.', cex=1, xlab='contig length', ylab='bp in mutations', xlim=c(0,500), ylim=c(0,500))

# histogram of contig end coverage
end_covg = c(c$'covg.left', c$'covg.right')
hist(end_covg, breaks=0:max(end_covg), main='', xlab='end coverage', xlim=c(0,100), ylim=c(0,100000))

dev.off()
