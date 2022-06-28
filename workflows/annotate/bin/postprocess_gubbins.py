

cnt = 0
with open('gubbins.recombination_predictions.gff', 'r') as file:
    for line in file:
        if not line[0] == '#':
            # 'taxa="   e j  d  b  c a  k  Reference f"'
            l = line.strip().split('\t')
            taxa = l[-1].split(';')[2].split(' ')
            n_taxa = len([i for i in taxa if i][1:])
            if n_taxa > 2:
               cnt += 1
               print(line)