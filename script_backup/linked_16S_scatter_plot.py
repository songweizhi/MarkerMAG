import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


value_list_x    = [1,2,3,4,5,3,2,5]
value_list_y    = [1,2,3,4,5,3,2,5]
pwd_figure      = '/Users/songweizhi/Desktop/multi_subplots.png'
axis_label_x    = 'Number of 16S rRNA gene in genome'
axis_label_y    = 'Linked 16S rRNA gene'


plt.scatter(value_list_x, value_list_y, c="b")
plt.xlabel(axis_label_x)
plt.ylabel(axis_label_y)
plt.tight_layout()
plt.savefig(pwd_figure)
plt.close()
