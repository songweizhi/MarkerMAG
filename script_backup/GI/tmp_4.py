import os
from Bio import SeqIO


###################################################### file in/out #####################################################

wd = '/Users/songweizhi/Desktop'

# file in
drep_cdb_file           = '%s/Cdb.csv'                         % wd
ref_to_strain_file      = '%s/ref_to_strain.txt'               % wd


########################################################################################################################

# get ref_to_strain_dict
ref_to_strain_dict = {}
for ref in open(ref_to_strain_file):
    ref_split = ref.strip().split('\t')
    ref_to_strain_dict[ref_split[0]] = ref_split[1]


cluster_to_ref_dict = {}
ref_to_cluster_dict = {}
for each_ref in open(drep_cdb_file):
    if not each_ref.startswith('genome,secondary_cluster'):
        each_ref_split = each_ref.strip().split(',')
        ref_file_name = each_ref_split[0]
        ref_file_name_no_ext = '.'.join(ref_file_name.split('.')[:-1])
        ref_cluster = 'C' + each_ref_split[1]
        ref_to_cluster_dict[ref_file_name_no_ext] = ref_cluster
        if ref_cluster not in cluster_to_ref_dict:
            cluster_to_ref_dict[ref_cluster] = [ref_file_name_no_ext]
        else:
            cluster_to_ref_dict[ref_cluster].append(ref_file_name_no_ext)

print(cluster_to_ref_dict)
print(ref_to_cluster_dict)

for each_cluster in cluster_to_ref_dict:
    for each_gnm in cluster_to_ref_dict[each_cluster]:
        #print('%s\t%s\t%s' % (each_cluster, each_gnm, ref_to_strain_dict[each_gnm]))
        print('%s\t%s' % (each_cluster, ref_to_strain_dict[each_gnm]))
    print()


'''
C29_1	CP008713.1 Candidatus Arthromitus sp. SFB-mouse-NL
C29_1	AP012209.1 Candidatus Arthromitus sp. SFB-mouse-Yit DNA
C29_1	AP012202.1 Candidatus Arthromitus sp. SFB-mouse-Japan DNA
C30_0	AP012210.1 Candidatus Arthromitus sp. SFB-rat-Yit DNA

C34_1	CP000673.1 Clostridium kluyveri DSM 555
C34_1	AP009049.1 Clostridium kluyveri NBRC 12016 DNA
C34_2	CP018335.1 Clostridium kluyveri strain JZZ

C36_1	CP001666.1 Clostridium ljungdahlii DSM 13528
C36_1	CP006763.1 Clostridium autoethanogenum DSM 10061
C36_1	CP012395.1 Clostridium autoethanogenum DSM 10061

C40_1	AP014696.1 Clostridium botulinum DNA, strain: 111
C40_1	CP000726.1 Clostridium botulinum A str. ATCC 19397
C40_1	AM412317.1 Clostridium botulinum A str. ATCC 3502
C40_1	CP001083.1 Clostridium botulinum Ba4 str. 657
C40_1	CP014148.1 Clostridium botulinum strain CDC_67190
C40_1	CP013247.1 Clostridium botulinum strain CDC_53174
C40_1	CP006907.1 Clostridium botulinum CDC_297
C40_1	CP000939.1 Clostridium botulinum B1 str. Okra
C40_1	CP013246.1 Clostridium botulinum strain CDC_69094
C40_1	FR773526.1 Clostridium botulinum H04402 065 sequence
C40_1	CP001581.1 Clostridium botulinum A2 str. Kyoto
C40_1	CP000728.1 Clostridium botulinum F str. Langeland
C40_1	CP006908.1 Clostridium botulinum CDC_1436
C40_1	CP000727.1 Clostridium botulinum A str. Hall
C40_1	CP002011.1 Clostridium botulinum F str. 230613
C40_2	CP000962.1 Clostridium botulinum A3 str. Loch Maree

C44_1	CP000721.1 Clostridium beijerinckii NCIMB 8052
C44_1	CP011966.1 Clostridium pasteurianum NRRL B-598
C44_1	CP006777.1 Clostridium beijerinckii ATCC 35702
C44_1	CP010086.2 Clostridium beijerinckii strain NCIMB 14988

C47_1	CP010520.1 Clostridium botulinum strain NCTC 8266
C47_1	CP010521.1 Clostridium botulinum strain NCTC 8550
C47_1	CP001078.1 Clostridium botulinum E3 str. Alaska E43

C47_2	CP001056.1 Clostridium botulinum B str. Eklund 17B
C47_2	FR745875.1 Clostridium botulinum B str. Eklund 17B(NRP)
C47_2	CP006903.1 Clostridium botulinum 202F

C52_1	CP009267.1 Clostridium pasteurianum DSM 525 = ATCC 6013
C52_1	CP013019.1 Clostridium pasteurianum strain M150B
C52_1	CP009268.1 Clostridium pasteurianum DSM 525 = ATCC 6013
C52_1	CP013018.1 Clostridium pasteurianum DSM 525 = ATCC 6013
C53_0	CP003261.1 Clostridium pasteurianum BC1


C75_1	CP014065.1 Achromobacter xylosoxidans strain FDAARGOS_162
C75_1	CP012046.1 Achromobacter xylosoxidans strain MN001
C75_2	CP014028.1 Achromobacter xylosoxidans strain FDAARGOS_150





'''