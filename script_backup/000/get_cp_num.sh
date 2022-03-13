
###################################################### MBARC26 #####################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
module load hmmer/3.3
module load blast+/2.11.0
module load seqtk/20190219
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/get_cp_num
read_r1='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_R1.fasta'
read_r2='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_R2.fasta'
read_16s='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_16S_reads.fasta'
cp_mags='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0830_iden99_MarkerMAG_wd/MBARC26_0830_iden99_rd1_wd/input_MAGs/Refined_refined_bins_renamed_combined.fa'
cp_mag_gff='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0830_iden99_MarkerMAG_wd/MBARC26_0830_iden99_rd1_wd/input_MAGs/combined_barrnap.gff'
identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0830_iden99_MarkerMAG_wd/MBARC26_0830_iden99_linkages_by_genome.txt'
provided_cp_num='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/get_cp_num/MBARC26_ref_16S_cp_num.txt'
mag_gc_bias_folder='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/get_cp_num/MAGs_2_50_GC_bias'
mag_cov_gc_txt='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/get_cp_num/MAGs_2_50_MAG_depth_GC_content.txt'
seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0830_iden99_MarkerMAG_wd/MBARC26_SILVA138_polished_polished_min1200bp_c99.0.fasta'
#seq_16s_c999='/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_polished.fasta'
vxtractor_pl='/srv/scratch/z5039045/Softwares/vxtractor/vxtractor.pl'
silver_order_seq='/srv/scratch/z5039045/MarkerMAG_wd/SILVA_138.1_one_seq_per_order.fasta'
vxtractor_hmm_bac='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/bacteria'
vxtractor_hmm_arc='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/archaea'

#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p MBARC26_0906 -both_pair_mapped -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p MBARC26_0906_AllReads -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p MBARC26_0906_AllReads_35_15 -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p MBARC26_0906_AllReads_45_05_2 -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p MBARC26_0906_AllReads_49_1 -ignore_lowest_pct 49 -ignore_highest_pct 1 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc

python3 get_cp_num.py -p MBARC26_0911 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
python3 get_cp_num.py -p MBARC26_0911_2 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc

python3 get_cp_num.py -p MBARC26_0911_3 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref_16s_cp_num $provided_cp_num -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc


######################################################## GI iden99 ########################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
module load hmmer/3.3
module load blast+/2.11.0
module load seqtk/20190219
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num
read_r1='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_R1.fasta'
read_r2='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_R2.fasta'
read_16s='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/3_GI_Matam16S_wd/3_GI_16S_reads.fasta'
cp_mags='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99_MarkerMAG_wd/GI_0830_iden99_rd1_wd/input_MAGs/GI_refined_bins_complete50_contain5_combined.fa'
cp_mag_gff='//srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99_MarkerMAG_wd/GI_0830_iden99_rd1_wd/input_MAGs/combined_barrnap.gff'
identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99_MarkerMAG_wd/GI_0830_iden99_linkages_by_genome.txt'
provided_cp_num='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num/provided_mag_16s_cp_num.txt'
mag_gc_bias_folder='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num/MAGs_2_50_GC_bias'
mag_cov_gc_txt='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num/MAGs_2_50_MAG_depth_GC_content.txt'
seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99_MarkerMAG_wd/GI_128_16S_0.999_polished_min1200bp_c99.0.fasta'
#seq_16s_c999='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_128_16S_0.999.fasta'
vxtractor_pl='/srv/scratch/z5039045/Softwares/vxtractor/vxtractor.pl'
silver_order_seq='/srv/scratch/z5039045/MarkerMAG_wd/SILVA_138.1_one_seq_per_order.fasta'
vxtractor_hmm_bac='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/bacteria'
vxtractor_hmm_arc='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/archaea'

#python3 get_cp_num.py -p GI_iden99_45_5_allow_singleton -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_iden99_45_5_allow_singleton_all_16S -map_to_all_16s -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_iden99_45_5_only_paired -both_pair_mapped -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_iden99_45_5_only_paired_uniq -by_uniq_aln -both_pair_mapped -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_iden99_45_5_only_paired_overall_16S -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_iden99_45_5_only_paired_overall_16S_uniq -map_to_all_16s -by_uniq_aln -both_pair_mapped -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_16S_cp_num_iden99_59_1 -both_pair_mapped -map_to_all_16s -ignore_lowest_pct 59 -ignore_highest_pct 1 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_16S_cp_num_iden99_59_1_uniq -by_uniq_aln -both_pair_mapped -map_to_all_16s -ignore_lowest_pct 59 -ignore_highest_pct 1 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_16S_cp_num_iden99_59_1_linked -both_pair_mapped -ignore_lowest_pct 59 -ignore_highest_pct 1 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p GI_0831_49_1_map_to_all_16s -map_to_all_16s -rm_mis_alignments -both_pair_mapped -ignore_lowest_pct 49 -ignore_highest_pct 1 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15 -both_pair_mapped -map_to_all_16s -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15 -both_pair_mapped -map_to_all_16s -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_linked -both_pair_mapped -map_to_all_16s -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_20 -both_pair_mapped -map_to_all_16s -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_20_linked -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_linked -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_linked_mapc10 -mapc 10 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_linked_mapc20 -mapc 20 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_linked_mapc30 -mapc 30 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_linked_mapc50 -mapc 50 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc0 -mapc 0 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc10 -mapc 10 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc20 -mapc 20 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc30 -mapc 30 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc40 -mapc 40 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc50 -mapc 50 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99_MarkerMAG_wd/GI_128_16S_0.999_polished_min1200bp_c99.0_no_3361.fasta'
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0901_35_15_15_mapc10_no_3361 -mapc 10 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0902_35_15_mapc10_linked -mapc 10 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0902_35_15_mapc10_all -mapc 10 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0906_AllReads -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0906_AllReads_35_15 -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0906_AllReads_45_05 -ignore_lowest_pct 45 -ignore_highest_pct 05 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0906_AllReads_49_01 -ignore_lowest_pct 49 -ignore_highest_pct 1 -force -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc


######################################################## GI iden99.5 ########################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
module load hmmer/3.3
module load blast+/2.11.0
module load seqtk/20190219
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num
read_r1='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_R1.fasta'
read_r2='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_R2.fasta'
read_16s='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/3_GI_Matam16S_wd/3_GI_16S_reads.fasta'
cp_mags='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99.5_MarkerMAG_wd/GI_0830_iden99.5_rd1_wd/input_MAGs/GI_refined_bins_complete50_contain5_combined.fa'
cp_mag_gff='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99.5_MarkerMAG_wd/GI_0830_iden99.5_rd1_wd/input_MAGs/combined_barrnap.gff'
identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99.5_MarkerMAG_wd/GI_0830_iden99.5_linkages_by_genome.txt'
provided_cp_num='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num/provided_mag_16s_cp_num.txt'
mag_gc_bias_folder='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num/MAGs_2_50_GC_bias'
mag_cov_gc_txt='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/get_cp_num/MAGs_2_50_MAG_depth_GC_content.txt'
seq_16s_c995='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0830_iden99.5_MarkerMAG_wd/GI_128_16S_0.999_polished_min1200bp_c99.5.fasta'
seq_16s_c999='/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_128_16S_0.999.fasta'
vxtractor_pl='/srv/scratch/z5039045/Softwares/vxtractor/vxtractor.pl'
silver_order_seq='/srv/scratch/z5039045/MarkerMAG_wd/SILVA_138.1_one_seq_per_order.fasta'
vxtractor_hmm_bac='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/bacteria'
vxtractor_hmm_arc='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/archaea'

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0902_iden995_35_15_mapc10 -mapc 10 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c995 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p GI_0902_iden995_25_25_mapc10 -mapc 10 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 25 -ignore_highest_pct 25 -r1 $read_r1 -r2 $read_r2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c995 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder


####################################################### Oral 99% #######################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
module load hmmer/3.3
module load blast+/2.11.0
module load seqtk/20190219
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num
cp_mags='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0828_iden99_MarkerMAG_wd/Oral_0828_iden99_rd1_wd/input_MAGs/3_Oral_refined_MAGs_complete50_contain5_combined.fa'
cp_mag_gff='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0828_iden99_MarkerMAG_wd/Oral_0828_iden99_rd1_wd/input_MAGs/combined_barrnap.gff'
identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0828_iden99_MarkerMAG_wd/Oral_0828_iden99_linkages_by_genome.txt'
read_16s='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_16S_reads.fasta'
provided_cp_num='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/provided_mag_16s_cp_num.txt'
mag_gc_bias_folder='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/MAGs_2_50_GC_bias'
mag_cov_gc_txt='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/MAGs_2_50_MAG_depth_GC_content.txt'
seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0828_iden99_MarkerMAG_wd/CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.0.fa'
#seq_16s_c999='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_138_16S_0.999.polished_min1200.fa'
vxtractor_pl='/srv/scratch/z5039045/Softwares/vxtractor/vxtractor.pl'
silver_order_seq='/srv/scratch/z5039045/MarkerMAG_wd/SILVA_138.1_one_seq_per_order.fasta'
vxtractor_hmm_bac='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/bacteria'
vxtractor_hmm_arc='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/archaea'

#python3 get_cp_num.py -p Oral_0831_49_1 -map_to_all_16s -ignore_lowest_pct 49 -ignore_highest_pct 1 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p Oral_0831_map_to_all_45_05 -map_to_all_16s -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p Oral_0831_map_to_all_40_10 -map_to_all_16s -ignore_lowest_pct 40 -ignore_highest_pct 10 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 get_cp_num.py -p Oral_0831_map_to_all_35_15 -map_to_all_16s -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder

#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_45 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_45_linked -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_linked_mapc0 -mapc 0 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_linked_mapc20 -mapc 20 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_linked_mapc30 -mapc 30 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_linked_mapc50 -mapc 50 -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_mapc0 -mapc 0 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder

#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_mapc10 -mapc 10 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_mapc20 -mapc 20 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_mapc30 -mapc 30 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_mapc40 -mapc 40 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder
#python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0901_35_15_mapc50 -mapc 50 -map_to_all_16s -both_pair_mapped -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -unclustered_marker $seq_16s_c999 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0906_AllReads -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0906_AllReads_35_15 -ignore_lowest_pct 35 -ignore_highest_pct 15 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0906_AllReads_45_05 -ignore_lowest_pct 45 -ignore_highest_pct 05 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Oral_0906_AllReads_49_01 -ignore_lowest_pct 49 -ignore_highest_pct 1 -force -r1 ../Oral_5_25_R1.fasta -r2 ../Oral_5_25_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc


####################################################### Oral 99% (why that bad) #######################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
module load hmmer/3.3
module load blast+/2.11.0
module load seqtk/20190219
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd
cp_mags='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/input_MAGs/3_Oral_refined_MAGs_complete50_contain5_combined.fa'
cp_mag_gff='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/input_MAGs/combined_barrnap.gff'
identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_linkages_by_genome.txt'
read_1='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_5_25_R1.fasta'
read_2='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_5_25_R2.fasta'
#read_16s='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_16S_reads.fasta'
provided_cp_num='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/provided_mag_16s_cp_num.txt'
mag_gc_bias_folder='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_get_16S_cp_num_wd/Oral_0909_iden990_GC_bias'
mag_cov_gc_txt='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_get_16S_cp_num_wd/Oral_0909_iden990_MAG_depth_GC_content.txt'
seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.0.fa'
vxtractor_pl='/srv/scratch/z5039045/Softwares/vxtractor/vxtractor.pl'
silver_order_seq='/srv/scratch/z5039045/MarkerMAG_wd/SILVA_138.1_one_seq_per_order.fasta'
vxtractor_hmm_bac='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/bacteria'
vxtractor_hmm_arc='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/archaea'

python3 get_cp_num.py -p test_1_25_25 -force -r1 $read_1 -r2 $read_2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc -sam_16s_sorted_by_read /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/Oral_0909_iden990_input_reads_to_16S_sorted.sam
python3 get_cp_num.py -p test_2_45_05 -ignore_lowest_pct 45 -ignore_highest_pct 5 -force -r1 $read_1 -r2 $read_2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc
python3 get_cp_num.py -p test_2_25_25 -force -r1 $read_1 -r2 $read_2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc -sam_16s_sorted_by_read /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/Oral_0909_iden990_input_reads_to_16S_sorted.sam
python3 get_cp_num2.py -p test_3_25_25 -force -r1 $read_1 -r2 $read_2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc -sam_16s_sorted_by_read /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/Oral_0909_iden990_input_reads_to_16S_sorted.sam
python3 get_cp_num.py -p test_22_25_25 -subsample_pct 100 -force -r1 $read_1 -r2 $read_2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref $provided_cp_num -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc -sam_16s_sorted_by_read /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/Oral_0909_iden990_input_reads_to_16S_sorted.sam

python3 get_cp_num.py -p test_sub20 -subsample_pct 20 -force -r1 $read_1 -r2 $read_2 -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -ref_16s_cp_num $provided_cp_num -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc -sam_16s_sorted_by_read /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/Oral_0909_iden990_MarkerMAG_wd/Oral_0909_iden990_rd1_wd/Oral_0909_iden990_input_reads_to_16S_sorted.sam


####################################################### Kelp 99% #######################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
module load hmmer/3.3
module load blast+/2.11.0
module load seqtk/20190219
cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp/get_cp_num
cp_mags='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_0907_MarkerMAG_wd/Kelp_0907_rd1_wd/input_MAGs/BH_ER_050417_refined_bins_complete50_contain5_combined.fa'
cp_mag_gff='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_0907_MarkerMAG_wd/Kelp_0907_rd1_wd/input_MAGs/combined_barrnap.gff'
identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_0907_MarkerMAG_wd/Kelp_0907_linkages_by_genome.txt'
read_16s='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_SILVA138_id99_Matam16S_wd/Kelp_SILVA138_id99_16S_reads.fasta'
seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_0907_MarkerMAG_wd/Kelp_SILVA138_id99_assembled_16S_uclust_0.999_polished_min1200bp_c99.0.fasta'
vxtractor_pl='/srv/scratch/z5039045/Softwares/vxtractor/vxtractor.pl'
silver_order_seq='/srv/scratch/z5039045/MarkerMAG_wd/SILVA_138.1_one_seq_per_order.fasta'
vxtractor_hmm_bac='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/bacteria'
vxtractor_hmm_arc='/srv/scratch/z5039045/Softwares/vxtractor/HMMs/SSU/archaea'
#provided_cp_num='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/provided_mag_16s_cp_num.txt'
#mag_gc_bias_folder='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/MAGs_2_50_GC_bias'
#mag_cov_gc_txt='/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/MAGs_2_50_MAG_depth_GC_content.txt'

python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Kelp_0907 -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 ../Kelp_R1.fasta -r2 ../Kelp_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc






identified_linkages='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/get_cp_num/fake_linkages.txt'
seq_16s_c990='/srv/scratch/z5039045/MarkerMAG_wd/Kelp/get_cp_num/fake_Matam_16S.fasta'
python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Kelp_0907_fake -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 ../Kelp_R1.fasta -r2 ../Kelp_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc


python3 /srv/scratch/z5039045/MarkerMAG_wd/get_cp_num.py -p Kelp_0907_fake -ignore_lowest_pct 25 -ignore_highest_pct 25 -force -r1 ../Kelp_R1.fasta -r2 ../Kelp_R2.fasta -cp_mags $cp_mags -cp_mag_gff $cp_mag_gff -linkages $identified_linkages -marker $seq_16s_c990 -r16s $read_16s -t 16 -vxtractor $vxtractor_pl -silva_order_refs $silver_order_seq -hmm_bac $vxtractor_hmm_bac -hmm_arc $vxtractor_hmm_arc -mag_cov_gc $mag_cov_gc_txt -mag_gc_bias $mag_gc_bias_folder















########################################################################################################################

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10

cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/gc_bias
bowtie2-build --quiet --threads 16 -f combined_prefixed_MAGs.fa combined_prefixed_MAGs
bowtie2 -x combined_prefixed_MAGs -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S combined_prefixed_MAGs_global_report_one.sam -p 16 -f --xeq --no-unal -N 1 -L 30

python3 filter_sam.py -in combined_prefixed_MAGs_global_report_one.sam -out combined_prefixed_MAGs_global_report_one_mis0_50bp.sam -mm 0 -aln 50
samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis0_50bp.sam
samtools depth -a combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted_depth.txt

python3 filter_sam.py -in combined_prefixed_MAGs_global_report_one.sam -out combined_prefixed_MAGs_global_report_one_mis0_100bp.sam -mm 0 -aln 100
samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis0_100bp.sam
samtools depth -a combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted_depth.txt

python3 filter_sam.py -in combined_prefixed_MAGs_global_report_one.sam -out combined_prefixed_MAGs_global_report_one_mis2_50bp.sam -mm 2 -aln 50
samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis2_50bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis2_50bp.sam
samtools depth -a combined_prefixed_MAGs_global_report_one_mis2_50bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis2_50bp_sorted_depth.txt

python3 filter_sam.py -in combined_prefixed_MAGs_global_report_one.sam -out combined_prefixed_MAGs_global_report_one_mis2_75bp.sam -mm 2 -aln 75
samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis2_75bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis2_75bp.sam
samtools depth -a combined_prefixed_MAGs_global_report_one_mis2_75bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis2_75bp_sorted_depth.txt

python3 filter_sam.py -in combined_prefixed_MAGs_global_report_one.sam -out combined_prefixed_MAGs_global_report_one_mis2_100bp.sam -mm 2 -aln 100
samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis2_100bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis2_100bp.sam
samtools depth -a combined_prefixed_MAGs_global_report_one_mis2_100bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis2_100bp_sorted_depth.txt

# subset depth file
head -100000 /Users/songweizhi/Desktop/combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted_depth.txt > /Users/songweizhi/Desktop/combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted_depth_subset.txt
head -5000000 /Users/songweizhi/Desktop/combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted_depth.txt > /Users/songweizhi/Desktop/combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted_depth_subset.txt


module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/gc_bias
bowtie2-build --quiet --threads 16 -f linked_16s.fa linked_16s
bowtie2 -x linked_16s -U MBARC26_SILVA138_id99_16S_reads.fasta -S linked_16S_local_report_all.sam -p 16 -f --xeq --local --all --no-unal -N 1 -L 30
samtools sort -n -O sam --threads 12 -o linked_16S_local_report_all_sorted_by_read.sam linked_16S_local_report_all.sam


python3 filter_sam.py -in Map_to_linked_16S_linked_16s.sam -out Map_to_linked_16S_linked_16s_mis2_50bp.sam -mm 2 -aln 50

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num_3/Map_to_linked_16S_get_16S_cp_num_wd
bowtie2 -x Map_to_linked_16S_linked_16s -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S linked_16S_only_paired.sam -p 16 -f --local --xeq --no-unal -N 1 -L 30

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num_3
bowtie2 -x CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5 -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired.sam -p 16 -f --local --xeq --no-unal -N 1 -L 30

python3 filter_sam.py -in linked_16S_only_paired.sam -out linked_16S_only_paired_mis2_50bp.sam -mm 2 -aln 50
python3 filter_sam.py -in CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired.sam -out CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired_mis2_50bp.sam -mm 2 -aln 50


cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num_3/mapping
bowtie2-build --quiet --threads 16 -f all_16s.fa all_16s
bowtie2-build --quiet --threads 16 -f all_16s_iden99.fa all_16s_iden99
bowtie2-build --quiet --threads 16 -f linked_16s.fa linked_16s
bowtie2 -x all_16s -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S all_16s_only_paired.sam -p 16 -f --local --xeq --no-unal -N 1 -L 30
bowtie2 -x all_16s_iden99 -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S all_16s_iden99_only_paired.sam -p 16 -f --local --xeq --no-unal -N 1 -L 30
bowtie2 -x linked_16s -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S linked_16S_only_paired.sam -p 16 -f --local --xeq --no-unal -N 1 -L 30

python3 filter_sam.py -in linked_16S_only_paired.sam -out linked_16S_only_paired_mis2_50bp.sam -mm 2 -aln 50
python3 filter_sam.py -in all_16s_only_paired.sam -out all_16s_only_paired_mis2_50bp.sam -mm 2 -aln 50
python3 filter_sam.py -in all_16s_iden99_only_paired.sam -out all_16s_iden99_only_paired_mis2_50bp.sam -mm 2 -aln 50

samtools sort -O sam --threads 4 -o linked_16s_only_paired_mis2_50bp_only_paired_sorted.sam linked_16s_only_paired_mis2_50bp_only_paired.sam
samtools sort -O sam --threads 4 -o all_16s_only_paired_mis2_50bp_only_paired_sorted.sam all_16s_only_paired_mis2_50bp_only_paired.sam
samtools sort -O sam --threads 4 -o all_16s_iden99_only_paired_mis2_50bp_only_paired_sorted.sam all_16s_iden99_only_paired_mis2_50bp_only_paired.sam


cd /Users/songweizhi/Desktop/111
samtools sort -O sam --threads 4 -o CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired_mis2_50bp_sorted.sam CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired_mis2_50bp.sam
samtools sort -O sam --threads 4 -o CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired_mis2_50bp_only_paired_sorted.sam CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.5_only_paired_mis2_50bp_only_paired.sam
samtools sort -O sam --threads 4 -o combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis0_50bp.sam
samtools sort -O sam --threads 4 -o combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis0_50bp.sam

module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
module load samtools/1.10

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/mapping
bowtie2-build --quiet --threads 16 -f 16s_iden99.fa 16s_iden990
bowtie2 -x 16s_iden990 -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S 16s_iden990_only_paired.sam -p 16 -f --local --xeq --no-unal
python3 filter_sam.py -in 16s_iden990_only_paired.sam -out 16s_iden990_only_paired_mis2_50bp.sam -mm 2 -aln 50
python3 keep_both_mapped.py -in 16s_iden990_only_paired_mis2_50bp.sam -out 16s_iden990_only_paired_mis2_50bp_both_mapped.sam
samtools sort -O sam --threads 16 -o 16s_iden990_only_paired_mis2_50bp_both_mapped_sorted.sam 16s_iden990_only_paired_mis2_50bp_both_mapped.sam

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num/mapping
bowtie2-build --quiet --threads 16 -f 16s_iden99.9.fa 16s_iden999
bowtie2 -x 16s_iden999 -1 CAMI_Oral_16S_reads_R1.fasta -2 CAMI_Oral_16S_reads_R2.fasta -S 16s_iden999_only_paired.sam -p 16 -f --local --xeq --no-unal
python3 filter_sam.py -in 16s_iden999_only_paired.sam -out 16s_iden999_only_paired_mis2_50bp.sam -mm 2 -aln 50
python3 keep_both_mapped.py -in 16s_iden999_only_paired_mis2_50bp.sam -out 16s_iden999_only_paired_mis2_50bp_both_mapped.sam
samtools sort -O sam --threads 16 -o 16s_iden999_only_paired_mis2_50bp_both_mapped_sorted.sam 16s_iden999_only_paired_mis2_50bp_both_mapped.sam


BioSAK top_16S_hits -p 16S_vs_GTDB -q 16S_sequences.fa -r GTDB_ssu_all_r95.fna -t 6

BioSAK top_16S_hits -p 16S_vs_SILVA -q MBARC26_SILVA138_polished_polished_min1200bp_c99.0.fasta -r /srv/scratch/z5039045/MarkerMAG_wd/SILVA_ref_16S_order/SILVA_138.1_one_seq_per_order.fasta -t 6
