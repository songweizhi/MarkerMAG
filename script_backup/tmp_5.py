import os

full_name_file = '/Users/songweizhi/Desktop/full_name_sorted.txt'
Kelp_sample_id = '/Users/songweizhi/Desktop/Kelp_sample_id.txt'
full_name_file = 'full_name_sorted.txt'
Kelp_sample_id = 'Kelp_sample_id.txt'


sample_file_dict = {}
for file in open(full_name_file):
    file = file.strip()
    sample_id = file.split('_')[0]
    if sample_id not in sample_file_dict:
        if 'pairedForward' in file:
            sample_file_dict[sample_id] = [[file],[]]
        if 'pairedReverse' in file:
            sample_file_dict[sample_id] = [[],[file]]
    else:
        if 'pairedForward' in file:
            sample_file_dict[sample_id][0].append(file)
        if 'pairedReverse' in file:
            sample_file_dict[sample_id][1].append(file)

for sample in sample_file_dict:
    current_sample_r1 = sorted(sample_file_dict[sample][0])
    current_sample_r2 = sorted(sample_file_dict[sample][1])
    current_sample_r1_no_gz = [i[:-3] for i in current_sample_r1]
    current_sample_r2_no_gz = [i[:-3] for i in current_sample_r2]

    for r1_gz in current_sample_r1:
        os.system('gunzip %s' % r1_gz)
    for r2_gz in current_sample_r2:
        os.system('gunzip %s' % r2_gz)

    os.system('cat %s > %s_R1.fastq' % (' '.join(current_sample_r1_no_gz), sample))
    os.system('cat %s > %s_R2.fastq' % (' '.join(current_sample_r2_no_gz), sample))

    for r1_gz in current_sample_r1_no_gz:
        os.system('rm %s' % r1_gz)
    for r2_gz in current_sample_r2_no_gz:
        os.system('rm %s' % r2_gz)



co_assembly_dict = {}
for co_assembly in open(Kelp_sample_id):
    co_assembly_split = co_assembly.strip().split('\t')
    if co_assembly_split[0] not in co_assembly_dict:
        co_assembly_dict[co_assembly_split[0]] = [co_assembly_split[1]]
    else:
        co_assembly_dict[co_assembly_split[0]].append(co_assembly_split[1])

for each_co_assembly in co_assembly_dict:
    if len(co_assembly_dict[each_co_assembly]) == 1:
        os.system('mv %s_R1.fastq %s_R1.fastq' % (co_assembly_dict[each_co_assembly][0], each_co_assembly))
        os.system('mv %s_R2.fastq %s_R2.fastq ' % (co_assembly_dict[each_co_assembly][0], each_co_assembly))
    else:
        os.system('cat %s > %s_R1.fastq' % (' '.join([(i + '_R1.fastq') for i in co_assembly_dict[each_co_assembly]]), each_co_assembly))
        os.system('cat %s > %s_R2.fastq' % (' '.join([(i + '_R2.fastq') for i in co_assembly_dict[each_co_assembly]]), each_co_assembly))
    for each_id in co_assembly_dict[each_co_assembly]:
        os.system('rm %s_R1.fastq' % each_id)
        os.system('rm %s_R2.fastq' % each_id)

    os.system('fq2fa %s_R1.fastq %s_R1.fasta' % (each_co_assembly, each_co_assembly))
    os.system('fq2fa %s_R2.fastq %s_R2.fasta' % (each_co_assembly, each_co_assembly))

    os.system('rm %s_R1.fastq' % each_co_assembly)
    os.system('rm %s_R2.fastq' % each_co_assembly)

