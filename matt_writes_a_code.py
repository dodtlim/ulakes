import numpy
import vcf
import distance
import csv

block_size = 10000
size_cut   = 2e7

vcf_reader = vcf.Reader(open('SNPs_only.kivu.vcf.gz', 'rb'))

contig_lengths = {}
for contig in vcf_reader.contigs :
    length = vcf_reader.contigs[contig].length
    contig_lengths[contig] = length




#plot(sorted(contig_lengths.values()))
#axhline(size_cut,color='black')
#semilogy()

bigcontigs = []
mysamples_df = {}
mysamples_dxy = {}
for name in vcf_reader.samples :
    mycontigs = {}
    for contig in vcf_reader.contigs :
        if contig_lengths[contig] > size_cut :
            if not bigcontigs.__contains__(contig) :
                bigcontigs.append(contig)
            mycontigs[contig] = { 'length' : length, 'snps' : numpy.zeros(contig_lengths[contig]/block_size+1) }
    mysamples_df[name] = mycontigs
    mysamples_dxy[name] = mycontigs

for n,obs in enumerate(vcf_reader) :
    if n%100000 == 0 : print n
    if n == 100 : break
    if not bigcontigs.__contains__(obs.CHROM) : continue
    if obs.is_indel : continue
    samp1 = obs.samples[0].gt_bases
    samp2 = obs.samples[1].gt_bases
    if samp1 is None or samp2 is None : continue
    samp1 = ''.join(sorted(samp1))
    samp2 = ''.join(sorted(samp2))
    sampdist = distance.hamming(samp1,samp2)
    mysamples_df[s.sample][obs.CHROM]['snps'][s.site.POS/block_size] \
        = mysamples_df[s.sample][obs.CHROM]['snps'][s.site.POS/block_size] + sampdist

