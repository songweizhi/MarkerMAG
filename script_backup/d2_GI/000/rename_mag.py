import os
import glob
import shutil
import argparse


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def rename_mag(args):

    output_prefix   = args['p']
    mag_folder      = args['m']
    mag_ext         = args['mx']
    ref_folder      = args['r']
    ref_ext         = args['rx']
    ani_cutoff      = args['c']
    num_thread      = args['t']


    ################################################### define file name ###################################################

    link_wd = '%s_link_mag_ref_wd' % output_prefix

    mag_folder_no_path = mag_folder
    if '/' in mag_folder:
        mag_folder_no_path = mag_folder.split('/')[-1]

    ref_folder_no_path = ref_folder
    if '/' in ref_folder:
        ref_folder_no_path = ref_folder.split('/')[-1]

    mag_list_file       = '%s/%s_%s_file_list.txt'      % (link_wd, output_prefix, mag_folder_no_path)
    ref_list_file       = '%s/%s_%s_file_list.txt'      % (link_wd, output_prefix, ref_folder_no_path)
    fastANI_output      = '%s/%s_fastANI.txt'           % (link_wd, output_prefix)
    ref_to_mag_table    = '%s/%s_ref_to_mag_table.txt'  % (link_wd, output_prefix)
    renamed_to_file     = '%s/%s_renamed_to.txt'        % (link_wd, output_prefix)
    renamed_mag_folder  = '%s/%s_renamed'               % (link_wd, mag_folder_no_path )

    force_create_folder(link_wd)


    ############################################### prepare files for fastANI ##############################################

    # get file list
    mag_file_re = '%s/*.%s' % (mag_folder, mag_ext)
    ref_file_re = '%s/*.%s' % (ref_folder, ref_ext)
    mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
    ref_file_list = [os.path.basename(file_name) for file_name in glob.glob(ref_file_re)]

    mag_list_file_handle = open(mag_list_file, 'w')
    for mag_file in mag_file_list:
        mag_list_file_handle.write('%s/%s\n' % (mag_folder, mag_file))
    mag_list_file_handle.close()

    ref_list_file_handle = open(ref_list_file, 'w')
    for ref_file in ref_file_list:
        ref_list_file_handle.write('%s/%s\n' % (ref_folder, ref_file))
    ref_list_file_handle.close()

    fastANI_cmd = 'fastANI --refList %s --ql %s -o %s --minFrag 1 -t %s' % (ref_list_file, mag_list_file, fastANI_output, num_thread)
    os.system(fastANI_cmd)


    ################################################# parse fastANI output #################################################

    os.mkdir(renamed_mag_folder)

    ref_to_query_dict = {}
    query_to_ref_dict = {}
    for line in open(fastANI_output):

        line_split = line.strip().split('\t')
        query_name = line_split[0]
        if '/' in query_name:
            query_name = query_name.split('/')[-1]

        ref_name   = line_split[1]
        if '/' in ref_name:
            ref_name = ref_name.split('/')[-1]

        ani_value  = float(line_split[2])

        if ani_value >= ani_cutoff:

            if ref_name not in ref_to_query_dict:
                ref_to_query_dict[ref_name] = [query_name]
            else:
                ref_to_query_dict[ref_name].append(query_name)

            if query_name not in query_to_ref_dict:
                query_to_ref_dict[query_name] = [ref_name]
            else:
                query_to_ref_dict[query_name].append(ref_name)


    ref_to_mag_table_handle = open(ref_to_mag_table, 'w')
    renamed_to_file_handle = open(renamed_to_file, 'w')
    ref_list_match_to_multiple_mag = []
    for ref in ref_to_query_dict:

        ref_no_ext                  = '.'.join(ref.split('.')[:-1])
        current_query_list          = sorted(ref_to_query_dict[ref])
        current_query_list_no_ext   = ['.'.join(i.split('.')[:-1]) for i in current_query_list]

        ref_to_mag_table_handle.write('%s\t%s\n' % (ref_no_ext, ','.join(current_query_list_no_ext)))

        # rename mags
        if len(current_query_list) == 1:
            rename_cmds = 'cp %s/%s %s/%s' % (mag_folder, current_query_list[0], renamed_mag_folder, ref)
            renamed_to_file_handle.write('%s\t%s\n' % ('.'.join(current_query_list[0].split('.')[:-1]), ref_no_ext))
            os.system(rename_cmds)
        elif len(current_query_list) > 1:
            rename_index = 1
            for each_query in current_query_list:
                rename_cmds = 'cp %s/%s %s/%s_%s.%s' % (mag_folder,  each_query, renamed_mag_folder,  ref_no_ext, rename_index, ref_ext)
                renamed_to_file_handle.write('%s\t%s_%s\n' % ('.'.join(each_query.split('.')[:-1]), ref_no_ext, rename_index))
                os.system(rename_cmds)
                rename_index += 1
            ref_list_match_to_multiple_mag.append(ref_no_ext)

    ref_to_mag_table_handle.close()
    renamed_to_file_handle.close()


    ######################################################## report ########################################################

    mags_matched_to_multi_refs_list = []
    for mag in query_to_ref_dict:
        if len(query_to_ref_dict[mag]) >= 2:
            mags_matched_to_multi_refs_list.append(mag)
    if len(mags_matched_to_multi_refs_list) == 0:
        print('Good! No MAG matched to multiple reference genomes.')
    else:
        print('Caution! Found MAG(s) matched to multiple reference genomes.')


    ref_list_not_match_to_mag = []
    for each_ref in ref_file_list:
        if each_ref not in ref_to_query_dict:
            ref_list_not_match_to_mag.append(each_ref)
    if len(ref_list_not_match_to_mag) == 0:
        print('Good! All reference(s) have matched MAGs')
    else:
        print('Caution! Found reference(s) do not have matched MAG(%s): %s' % (len(ref_list_not_match_to_mag), ','.join(ref_list_not_match_to_mag)))


    if len(ref_list_match_to_multiple_mag) == 0:
        print('Good! No reference matched to multiple MAGs')
    else:
        print('Caution! Found reference(s) matched to multiple MAG(%s): %s' % (len(ref_list_match_to_multiple_mag), ','.join(sorted(ref_list_match_to_multiple_mag))))


if __name__ == '__main__':

    rename_mag_parser = argparse.ArgumentParser()

    rename_mag_parser.add_argument('-p',             required=True,                              help='output prefix')
    rename_mag_parser.add_argument('-m',             required=True,  default=None,               help='MAG folder')
    rename_mag_parser.add_argument('-mx',            required=False, default='fasta',            help='MAG file extension, default: fasta')
    rename_mag_parser.add_argument('-r',             required=True,  default=None,               help='REF folder')
    rename_mag_parser.add_argument('-rx',            required=False, default='fasta',            help='REF file extension, default: fasta')
    rename_mag_parser.add_argument('-c',             required=False, type=float, default=99,     help='ANI cutoff, default: 99')
    rename_mag_parser.add_argument('-t',             required=False, type=int, default=1,        help='number of threads, default: 1')

    args = vars(rename_mag_parser.parse_args())

    rename_mag(args)


'''

module load python/3.7.3
source ~/mypython3env/bin/activate
module load fastani/1.1
cd /srv/scratch/z5039045/test
python3 rename_mag.py -p Test -m Refined_refined_bins -mx fasta -r selected_genomes_renamed_no_plasmid -rx fna -c 99 -t 12

'''
