from pathlib import Path
import screed


# 
with screed.open('coding.faa') as file:
    cnt, folder = 0, 0

    for line in file:
        name = line.name.split(' # ')[0]  # eg 1_3 ie contig 1 ORF 3 etc
        seq = line.sequence
        # AlphaFold2 does not like "*"
        if seq[-1] == '*':
            seq = seq[:-1]

        if cnt % 10000 == 0:
            folder += 1

            p = Path(f'{folder}')
            _ = p.mkdir(exist_ok=False)

        
        with open(p / f'{name}.fasta', 'w+') as out:
            out.write(f'>{name}\n{seq}\n')

        cnt += 1
