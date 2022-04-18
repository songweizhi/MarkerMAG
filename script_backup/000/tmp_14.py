import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_scatter_plot(num_list_1, num_list_2, png_file):

    if (len(num_list_1) > 0) and (len(num_list_2) > 0):

        num_arrary_1 = np.array(num_list_1)
        num_arrary_2 = np.array(num_list_2)

        fig = plt.figure(figsize=(6, 6))
        plt.margins(0)

        plt.scatter(num_arrary_1, num_arrary_2)
        plt.xlabel("Estimated copy number", fontsize=15)
        plt.ylabel("User provided copy number", fontsize=15)

        # set axis range
        max_value = max([max(num_list_1), max(num_list_2)])
        plt.xlim(0, round(max_value + 1))
        plt.ylim(0, round(max_value + 1))

        # add fit line
        coeffs = np.polyfit(num_arrary_1, num_arrary_2, 1)
        slope = coeffs[0]
        intercept = coeffs[1]
        plt.plot(num_arrary_1, slope * num_arrary_1 + intercept)

        # get R-squared value
        p = np.poly1d(coeffs)
        yhat = p(num_arrary_1)
        ybar = np.sum(num_arrary_2)/len(num_arrary_2)
        ssreg = np.sum((yhat-ybar)**2)
        sstot = np.sum((num_arrary_2 - ybar)**2)
        r_squared = ssreg/sstot

        if intercept >= 0:
            title_str = 'y = %sx + %s, R2 = %s' % (float("{0:.4f}".format(slope)), float("{0:.4f}".format(intercept)), float("{0:.4f}".format(r_squared)))
        else:
            title_str = 'y = %sx - %s, R2 = %s' % (float("{0:.4f}".format(slope)), abs(float("{0:.4f}".format(intercept))), float("{0:.4f}".format(r_squared)))

        plt.title(title_str, fontsize=16)

        # save plot
        plt.tight_layout()
        plt.savefig(png_file, format='svg')
        plt.close()


num_list_1 = []
num_list_2 = []
for each in open('/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/Figures/Figure 2 copy_number/GI_0912_copy_num_by_MAG.txt'):
    if not each.startswith('MAG'):
        each_split = each.strip().split('\t')
        num_list_1.append(float(each_split[3]))
        num_list_2.append(float(each_split[1]))
get_scatter_plot(num_list_1, num_list_2, '/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/Figures/Figure 2 copy_number/GI_0912_copy_num_by_MAG.svg')


num_list_1 = []
num_list_2 = []
for each in open('/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/Figures/Figure 2 copy_number/Oral_0913_copy_num_by_MAG.txt'):
    if not each.startswith('MAG'):
        each_split = each.strip().split('\t')
        num_list_1.append(float(each_split[3]))
        num_list_2.append(float(each_split[1]))
get_scatter_plot(num_list_1, num_list_2, '/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/Figures/Figure 2 copy_number/Oral_0913_copy_num_by_MAG_2.svg')


num_list_1 = []
num_list_2 = []
for each in open('/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/Figures/Figure 2 copy_number/MBARC26_0912_copy_num_by_MAG.txt'):
    if not each.startswith('MAG'):
        each_split = each.strip().split('\t')
        num_list_1.append(float(each_split[3]))
        num_list_2.append(float(each_split[2]))
get_scatter_plot(num_list_1, num_list_2, '/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/Figures/Figure 2 copy_number/MBARC26_0912_copy_num_by_MAG_norm_2.svg')
