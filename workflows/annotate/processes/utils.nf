process rename {
    cache 'lenient'  // TODO: not sure why this is needed, file sys changes?

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), path("${name}.fna"), path("${name}.json"))

    script:
    """
    rename.py --name ${name} --genome ${genome}
    """

}


process checksum {
    cache 'lenient'  // path and size-based
    // cache 'deep'  // content-based but slow

    input:
        tuple(val(name), path(genome))

    output:
        tuple(val(name), env(checksum), path(genome))

    shell:
    '''
    checksum=$(md5sum !{genome} | awk '{ printf $1 }')
    '''
    // awk print will insert newline, printf does not
    // stackoverflow.com/questions/2021982
    // alternative: ... | tr -d '\n'
    // stackoverflow.com/questions/3134791
}
