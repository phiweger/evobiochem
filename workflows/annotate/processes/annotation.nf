process coding {
    /*
    MAGs with -p single and -gc 11:

    > By default, Prodigal tries genetic code 11. -- https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#draft-genomes

    https://github.com/hyattpd/Prodigal/issues/57
    */
    // container 'nanozoo/prodigal:2.6.3--2769024'
    // Some genomes might have been corrupted during download
    // errorStrategy 'ignore'
    publishDir "${params.results}/coding", mode: 'copy', overwrite: true
    // Fix cpus to force 1:1 parallelization
    cpus 1

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path("${name}.coding.fna"), path("${name}.faa"))

    script:
    
    if( params.test )
        """
        prodigal -p meta -g 11 -q -i ${genome} -a ${name}.faa -d ${name}.coding.fna > /dev/null
        """
    else
        """
        prodigal -p single -g 11 -q -i ${genome} -a ${name}.faa -d ${name}.coding.fna > /dev/null
        """
}

