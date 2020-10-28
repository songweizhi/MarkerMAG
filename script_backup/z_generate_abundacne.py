import random


genome_id_list = ['s1', 's2', 's3', 's4', 's5', 'g1', 'g2', 'g3', 'g4', 'g5', 'f1', 'f2', 'f3', 'f4', 'f5', 'o1', 'o2', 'o3', 'o4', 'o5', 'c1', 'c2', 'c3', 'p1', 'p2', 'p3', 'p4', 'p5']

# 25x       47.58
# 1x        1.903
# 8.33x     15.86 x 28 = 444


random_lol = []
for i in range(0,3):
    randomlist = random.sample(range(1, 100), 28)
    randomlist_sum = sum(randomlist)
    divid_factor = randomlist_sum/444
    randomlist_norm = [round(i/divid_factor) for i in randomlist]

    depth_file = '/Users/songweizhi/Desktop/depth_rep%s.txt' % (i + 1)
    depth_file_handle = open(depth_file, 'w')
    for (genome, depth) in zip(genome_id_list, randomlist_norm):
        depth_file_handle.write('%s\t%s\n' % (genome, depth))
    depth_file_handle.close()

    print(randomlist_norm)
