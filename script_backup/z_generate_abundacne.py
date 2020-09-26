import random


genome_id_list = ['s1', 's2', 's3', 's4', 's5', 'g1', 'g2', 'g3', 'g4', 'g5', 'f1', 'f2', 'f3', 'f4', 'f5', 'o1', 'o2', 'o3', 'o4', 'o5', 'c1', 'c2', 'c3', 'p1', 'p2', 'p3', 'p4', 'p5']
randomlist = random.sample(range(10, 30), 5)


random_lol = []
for i in range(0,3):
    randomlist = random.sample(range(1, 100), 28)
    randomlist_pct = [float("{0:.4f}".format(i/sum(randomlist))) for i in randomlist]
    random_lol.append(randomlist_pct)


print('id\trep_1\trep_2\trep_3')
for (genome, x, y, z) in zip(genome_id_list, random_lol[0], random_lol[1], random_lol[2]):
    print('%s\t%s\t%s\t%s' % (genome, x, y, z))



