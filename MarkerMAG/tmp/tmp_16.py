
def r12_16s_ref_dict_to_str(r2_16s_ref_dict):

    top_dict_list = []
    for sub_dict_key in r2_16s_ref_dict:
        sub_dict_value = r2_16s_ref_dict[sub_dict_key]
        bottom_dict_list = []
        for bottom_dict_key in sub_dict_value:
            bottom_dict_value = sub_dict_value[bottom_dict_key]
            bottom_dict_list.append('%s:%s' % (bottom_dict_key, bottom_dict_value))
        bottom_dict_str = '%s' % (';'.join(bottom_dict_list))
        top_dict_list.append('%s:::%s' % (sub_dict_key, bottom_dict_str))

    top_dict_str = ';;;'.join(top_dict_list)

    return top_dict_str


def get_r12_16s_ref_dict_from_str(r2_16s_ref_dict_str):

    r2_16s_ref_dict_str_split = r2_16s_ref_dict_str.split(';;;')

    ref_dict_from_str = {}
    for each_sub_dict_str in r2_16s_ref_dict_str_split:
        each_sub_dict_str_split = each_sub_dict_str.split(':::')
        sub_dict_key = each_sub_dict_str_split[0]
        sub_dict_str = each_sub_dict_str_split[1]
        sub_dict_str_split = sub_dict_str.split(';')
        bottom_dict = dict()
        for each_bottom_dict_str in sub_dict_str_split:
            bottom_dict_str_split = each_bottom_dict_str.split(':')
            bottom_dict_key = int(bottom_dict_str_split[0])
            bottom_dict_value = bottom_dict_str_split[1]
            bottom_dict[bottom_dict_key] = bottom_dict_value
        ref_dict_from_str[sub_dict_key] = bottom_dict

    return ref_dict_from_str


def no_ignored_dict_to_str(r2_16s_refs_no_ignored_dict):
    top_dict_list = []
    for each_ref in r2_16s_refs_no_ignored_dict:
        each_ref_cigar_list = r2_16s_refs_no_ignored_dict[each_ref]
        each_ref_cigar_str = ','.join(each_ref_cigar_list)
        top_dict_list.append('%s:%s' % (each_ref, each_ref_cigar_str))
    top_dict_str = ';'.join(top_dict_list)

    return top_dict_str


def get_no_ignored_dict_from_str(r2_16s_refs_no_ignored_dict_str):
    r2_16s_refs_no_ignored_dict_str_split = r2_16s_refs_no_ignored_dict_str.split(';')
    no_ignored_dict_from_str = {}
    for each_sub_dict_str in r2_16s_refs_no_ignored_dict_str_split:
        each_sub_dict_str_split = each_sub_dict_str.split(':')
        ref_id = each_sub_dict_str_split[0]
        cigar_list_str = each_sub_dict_str_split[1]
        cigar_list = cigar_list_str.split(',')
        no_ignored_dict_from_str[ref_id] = cigar_list

    return no_ignored_dict_from_str


'''
MBARC26_891084	r2_16s_ref_dict	{'FA_m1': {426: '123='}, 'PS_m1': {406: '61S62='}, 'EV_m1': {396: '61S19=1X24=1X14=3S'}, 'EC_m3': {495: '65S15=1X24=1X17='}, 'EC_m6': {357: '65S15=1X24=1X17='}, 'AS_m1': {409: '65S15=1X24=1X17='}, 'AS_m2': {495: '65S15=1X24=1X17='}, 'SB_m1': {227: '65S15=1X24=1X17='}, 'AS_m3': {493: '65S15=1X24=1X17='}, 'SB_m2': {409: '65S15=1X24=1X17='}, 'EC_m2': {257: '65S15=1X24=1X17='}, 'CG_m3': {367: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'CG_m5': {366: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'CG_m2': {304: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'CG_m1': {456: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'OU_m1': {452: '47S12=2X3=1X4=1X10=1X24=1X4=1X10=2S'}, 'TR_m1': {419: '47S11=3X3=1X4=1X36=1X3=1X9=3S'}, 'SR_m1': {444: '47S12=3X2=1X15=1X24=1X7=1X6=3S'}, 'HB_m1': {416: '47S10=1X1=2X3=1X4=1X40=1X2=1X6=3S'}, 'FP_m2': {346: '45S5=1X6=1X2=1X8=1X10=1X24=1X4=1X2=1X7=2S'}, 'FP_m1': {434: '46S4=1X6=1X2=1X8=1X10=1X24=1X4=1X2=1X7=2S'}, 'MS_m1': {452: '66S14=1X24=1X7=1X7=2S'}, 'DA_m5': {573: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DA_m1': {269: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DM_m7': {466: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DM_m6': {464: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DM_m4': {411: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'TC_m2': {394: '47S14=1X1=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DM_m2': {464: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DM_m5': {464: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DA_m3': {333: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DA_m2': {412: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DM_m1': {46: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DA_m4': {401: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'CA_m1': {509: '72S8=1X24=1X2=2X11=2S'}, 'HT_m1': {453: '42S4=1X10=2X1=1X3=1X4=2X9=1X24=1X3=14S'}, 'HT_m2': {453: '42S4=1X10=2X1=1X3=1X4=2X9=1X24=1X3=14S'}, 'HT_m3': {453: '42S4=1X10=2X1=1X3=1X4=2X9=1X24=1X3=14S'}, 'CP_m4': {458: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m6': {265: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m3': {475: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m5': {458: '71S9=1X24=1X3=2X2=1X7=2S'}, 'SP_m1': {506: '71S9=1X24=1X3=14S'}, 'CP_m1': {456: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m2': {475: '71S9=1X24=1X3=2X2=1X7=2S'}}

MBARC26_891084	r1_16s_refs_no_ignored	{}
MBARC26_891084	r2_16s_refs_no_ignored	{}
MBARC26_891084	shared_16s_refs_no_ignored	{'FA_m1': ['125=', '123=']}

'''

# r2_16s_ref_dict = {'FA_m1': {426: '123=', 416: '61S62='}, 'PS_m1': {406: '61S62='}, 'EV_m1': {406: '61S62=', 396: '61S19=1X24=1X14=3S'}, 'EC_m3': {495: '65S15=1X24=1X17='}, 'EC_m6': {357: '65S15=1X24=1X17='}, 'AS_m1': {409: '65S15=1X24=1X17='}, 'AS_m2': {495: '65S15=1X24=1X17='}, 'SB_m1': {227: '65S15=1X24=1X17='}, 'AS_m3': {493: '65S15=1X24=1X17='}, 'SB_m2': {409: '65S15=1X24=1X17='}, 'EC_m2': {257: '65S15=1X24=1X17='}, 'CG_m3': {367: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'CG_m5': {366: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'CG_m2': {304: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'CG_m1': {456: '47S12=2X3=1X15=1X24=1X7=1X6=3S'}, 'OU_m1': {452: '47S12=2X3=1X4=1X10=1X24=1X4=1X10=2S'}, 'TR_m1': {419: '47S11=3X3=1X4=1X36=1X3=1X9=3S'}, 'SR_m1': {444: '47S12=3X2=1X15=1X24=1X7=1X6=3S'}, 'HB_m1': {416: '47S10=1X1=2X3=1X4=1X40=1X2=1X6=3S'}, 'FP_m2': {346: '45S5=1X6=1X2=1X8=1X10=1X24=1X4=1X2=1X7=2S'}, 'FP_m1': {434: '46S4=1X6=1X2=1X8=1X10=1X24=1X4=1X2=1X7=2S'}, 'MS_m1': {452: '66S14=1X24=1X7=1X7=2S'}, 'DA_m5': {573: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DA_m1': {269: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DM_m7': {466: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DM_m6': {464: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DM_m4': {411: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'TC_m2': {394: '47S14=1X1=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DM_m2': {464: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DM_m5': {464: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DA_m3': {333: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DA_m2': {412: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'DM_m1': {46: '47S10=2X4=2X4=1X10=1X24=1X4=1X9=3S'}, 'DA_m4': {401: '47S11=1X4=2X4=1X10=1X24=1X4=1X2=1X6=3S'}, 'CA_m1': {509: '72S8=1X24=1X2=2X11=2S'}, 'HT_m1': {453: '42S4=1X10=2X1=1X3=1X4=2X9=1X24=1X3=14S'}, 'HT_m2': {453: '42S4=1X10=2X1=1X3=1X4=2X9=1X24=1X3=14S'}, 'HT_m3': {453: '42S4=1X10=2X1=1X3=1X4=2X9=1X24=1X3=14S'}, 'CP_m4': {458: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m6': {265: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m3': {475: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m5': {458: '71S9=1X24=1X3=2X2=1X7=2S'}, 'SP_m1': {458: '71S9=1X24=1X3=2X2=1X7=2S', 506: '71S9=1X24=1X3=14S'}, 'CP_m1': {456: '71S9=1X24=1X3=2X2=1X7=2S'}, 'CP_m2': {506: '71S9=1X24=1X3=14S', 475: '71S9=1X24=1X3=2X2=1X7=2S'}}
# r2_16s_ref_dict_str = r12_16s_ref_dict_to_str(r2_16s_ref_dict)
# ref_dict_from_str = r12_16s_ref_str_to_dict(r2_16s_ref_dict_str)
# print(r2_16s_ref_dict_str)
# print(ref_dict_from_str)
# print(ref_dict_from_str == r2_16s_ref_dict)

r2_16s_refs_no_ignored_dict = {'DG_m4': ['69=34S', '69=34S', '69=34S'], 'DG_m7': ['68=35S'], 'DG_m6': ['68=35S'], 'DG_m5': ['54=49S']}



r2_16s_refs_no_ignored_dict_str = no_ignored_dict_to_str(r2_16s_refs_no_ignored_dict)
print(r2_16s_refs_no_ignored_dict_str)
r2_16s_refs_no_ignored_dict_str = 'DG_m4:69=34S,69=34S,69=34S;DG_m7:68=35S;DG_m6:68=35S;DG_m5:54=49S'
r2_16s_refs_no_ignored_dict_str = ''

r2_16s_refs_no_ignored_dict_from_str = get_no_ignored_dict_from_str(r2_16s_refs_no_ignored_dict_str)
print(r2_16s_refs_no_ignored_dict_from_str)
print(r2_16s_refs_no_ignored_dict_from_str == r2_16s_refs_no_ignored_dict)