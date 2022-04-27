nextflow.enable.dsl = 2


/*
To prepare an input file:


from pathlib import Path
from uuid import uuid4


files = Path('genomes').glob('*.fa')
d = {}
for i in files:
    name = uuid4().__str__().split('-')[0]
    if name in d:
        raise ValueError('Dublicate key generated')
    d[name] = i.resolve()

with open('input.csv', 'w+') as out:
    for k, v in d.items():
        out.write(f'{k},{v}\n')
*/


process coding {
    /*
    MAGs with -p single and -gc 11:

    > By default, Prodigal tries genetic code 11. -- https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#draft-genomes

    https://github.com/hyattpd/Prodigal/issues/57
    */
    container 'nanozoo/prodigal:2.6.3--2769024'
    // Some genomes might have been corrupted during download
    errorStrategy 'ignore'
    publishDir "${params.results}/coding", mode: 'copy', overwrite: true
    // Fix cpus to force 1:1 parallelization
    cpus 1

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path("${name}.faa"))

    """
    prodigal -p single -gc 11 -q -i ${genome} -a ${name}.faa > /dev/null
    """
}


process cluster {
    /*
    On cluster modes:

    https://github.com/soedinglab/MMseqs2/wiki#how-to-set-the-right-alignment-coverage-to-cluster
    */
    container 'nanozoo/prodigal:13.45111--810c0f2'
    publishDir "${params.results}", mode: 'copy', overwrite: true

    input:
        path(frames)

    output:
        path('cluster_rep_seq.fasta')

    """
    cat ${frames} > frames.faa
    mmseqs easy-cluster frames.faa cluster tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0
    """
}


workflow {
    genomes = channel.fromPath(params.genomes)
                     .splitCsv(header: false)
                     .map{ row -> [row[0], row[1]] }
                     // .view()
    coding(genomes)
    frames = coding.out.map{ it -> it[1] }.collect()
    cluster(frames)
}