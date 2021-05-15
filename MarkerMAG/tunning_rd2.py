

class MappingRecord:

    def __init__(self):

        #################### overall ####################

        self.qualified_reads = False

        #################### round 1 16s ####################

        self.r1_16s_ref_dict = dict()
        self.r2_16s_ref_dict = dict()

        self.r1_16s_refs_lowest_mismatch = None
        self.r2_16s_refs_lowest_mismatch = None

        self.r1_16s_refs_no_ignored = dict()
        self.r2_16s_refs_no_ignored = dict()
        self.shared_16s_refs_no_ignored = dict()

        self.both_mapped_to_16s = False

        #################### round 1 ctg ####################

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()

        self.r1_ctg_refs_lowest_mismatch = None
        self.r2_ctg_refs_lowest_mismatch = None

        self.r1_ctg_refs_no_ignored = dict()
        self.r2_ctg_refs_no_ignored = dict()
        self.shared_ctg_refs_no_ignored = dict()

        self.matched_to_ctg = False

        #################### round 2 ####################

        self.r1_ctg_ref_dict_rd2 = dict()
        self.r2_ctg_ref_dict_rd2 = dict()

        self.r1_ctg_refs_lowest_mismatch_rd2 = None
        self.r2_ctg_refs_lowest_mismatch_rd2 = None

        #################################################


        self.r1_longest_clp_cigar = ''
        self.r1_longest_clp_falg = ''
        self.r2_longest_clp_cigar = ''
        self.r2_longest_clp_falg = ''

        self.consider_r1_unmapped_mate = False
        self.consider_r1_clipping_part = False
        self.consider_r2_unmapped_mate = False
        self.consider_r2_clipping_part = False

        self.consider_round_2 = False

        self.r1_filtered_refs = set()
        self.r2_filtered_refs = set()
        self.unmapped_r1_refs = set()
        self.unmapped_r2_refs = set()
        self.clipping_r1_refs = set()
        self.clipping_r2_refs = set()
        self.unmapped_r1_refs_with_pos = set()
        self.unmapped_r2_refs_with_pos = set()
        self.clipping_r1_refs_with_pos = set()
        self.clipping_r2_refs_with_pos = set()


rd1_unlinked_mags_sam_bowtie_reformat = '/Users/songweizhi/Desktop/tunning_rd2/round_1_unlinked_gnm_bowtie.sam'

round_2_ctg_end_seq_len_dict = {}
round_2_MappingRecord_dict = {}
for each_line in open(rd1_unlinked_mags_sam_bowtie_reformat):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('@'):
        rd2_ref_id = ''
        rd2_ref_len = 0
        for each_element in each_line_split:
            if each_element.startswith('SN:'):
                rd2_ref_id = each_element[3:]
            if each_element.startswith('LN:'):
                rd2_ref_len = int(each_element[3:])
        round_2_ctg_end_seq_len_dict[rd2_ref_id] = rd2_ref_len
    else:
        cigar = each_line_split[5]
        if cigar != '*':
            read_id = each_line_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_line_split[2]
            ref_pos = int(each_line_split[3])

            if read_id_base not in round_2_MappingRecord_dict:
                round_2_MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if ref_id not in round_2_MappingRecord_dict[read_id_base].r1_ctg_ref_dict_rd2:
                    round_2_MappingRecord_dict[read_id_base].r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                else:
                    round_2_MappingRecord_dict[read_id_base].r1_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar

            if read_strand == '2':
                if ref_id not in round_2_MappingRecord_dict[read_id_base].r2_ctg_ref_dict_rd2:
                    round_2_MappingRecord_dict[read_id_base].r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                else:
                    round_2_MappingRecord_dict[read_id_base].r2_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar

