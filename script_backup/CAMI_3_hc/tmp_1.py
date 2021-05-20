import os

bin_id_file = '/Users/songweizhi/Desktop/refined_MAGs_prefix.txt'

for each in open(bin_id_file):
    bin_id = each.strip()
    add_prefix_cmd = 'blastn -query /srv/scratch/z5039045/MarkerMAG_wd/CAMI_high/CAMI_high_refined_bins_renamed/%s.fasta -db combined_refs.fna -out %s_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12' % (bin_id, bin_id)
    print(add_prefix_cmd)


# /srv/scratch/z5039045/MarkerMAG_wd/CAMI_high/blast_between_bin_and_ref/combined_refs.fna