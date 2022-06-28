import screed


fp = 'sequences.fna'

cnt = 0
with screed.open(fp) as file:
    for line in file:
        cnt += 1
        with open(f'primate_{cnt}.fna', 'w+') as out:
            out.write(f'>{line.name}\n{line.sequence}\n')


# Manually restructure


