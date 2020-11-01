

subsample_pcts = '1.000,5,10,25,50.5,75'


def str_to_um_list(nums_str):

    subsample_pct_list = []
    for pct_value in [str(float(i)) for i in subsample_pcts.split(',')]:
        if pct_value[-2:] == '.0':
            subsample_pct_list.append(int(float(pct_value)))
        else:
            subsample_pct_list.append(float(pct_value))

    return subsample_pct_list






print(str_to_um_list(subsample_pcts))
# subsample_pct_list = []
# for each_pct in subsample_pcts.split(','):
#     if isinstance()

