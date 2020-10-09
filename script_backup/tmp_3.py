import os


cmd_sankey_intersection = 'Rscript /Users/songweizhi/PycharmProjects/MarkerMAG/MarkerMAG/get_sankey_plot.R -f /Users/songweizhi/Desktop/EvenDepth3x_stats_intersection.txt -x 300 -y 600'


os.system(cmd_sankey_intersection)

count = len(open('/Users/songweizhi/Desktop/EvenDepth3x_stats_intersection.txt').readlines())
print(count)