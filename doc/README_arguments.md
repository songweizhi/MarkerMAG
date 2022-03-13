
### Arguments


+ `skip_16s_qc`     
  Skip checking if there are non-16S sequences at the end of provided 16S rRNA gene sequences. Specify only if you are confident with the quality of your 16S rRNA gene sequences.


+ `min_16s_len`    
  16S rRNA gene sequences shorter than specified length (default: 1200) will be ignored.
  

+ `max_16s_div`    
   Maximum genetic divergence (%) of 16S rRNA genes that allow to be linked to the same MAG (default : 1).


+ `keep_ctg_end_16s`      
   16S rRNA gene (fragment) sequences located at MAG contig ends will be removed by default, 
   as false positives might be introduced by aligning of reads to the conserved regions of these sequences.


+ `mismatch`     
   Maximum mismatch percentage for reads alignment, (default: 2)


+ `aln_len`   
   Minimum alignment length for reads alignment, (default: 45)


+ `aln_pct`   
   Minimum read alignment percentage, (default: 35)


+ `min_link`  
   Minimum number of linkages to report, (default: 9)
