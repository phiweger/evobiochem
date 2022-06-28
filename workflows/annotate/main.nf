nextflow.enable.dsl = 2

include { coding; coding as coding_ref } from './processes/annotation.nf'
include { rename; rename as rename_ref } from './processes/annotation.nf'

/*
Run:

from pathlib import Path

p = Path('more_genomes')

cnt = 0
with open('more_input.csv', 'w+') as out:
    for i in p.glob('*.fna'):
        cnt += 1
        out.write(f's{cnt}\t{i.resolve()}\n')


nextflow run ~/Dropbox/repos/evobiochem/workflows/annotate/main.nf \
    --genomes input.csv \
    --ref GCF_000008685.2_ASM868v2_genomic.fna \
    --ptable recombination/PVT.3SEQ.2017.700 \
    --dereplicate \
    -resume

Test w/ primate data:

TODO: automate ptable creation or download from osf
TODO: Add neisseria seqs

nextflow run ~/Dropbox/repos/evobiochem/workflows/annotate/main.nf \
    --ptable recombination/PVT.3SEQ.2017.700 \
    --test \
    -resume
*/


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



/*
// Specify number of CPUs, memory
// https://bioinformatics.stackexchange.com/questions/8712/autodetect-max-number-of-cores-and-pass-as-an-argument-in-nextflow
params.requestedCPU = 80
params.requestedMem = 50

maxcpus = Runtime.runtime.availableProcessors()
requestedCpus = Math.round((params.requestedCPU * maxcpus) / 100)

maxmem = Runtime.runtime.totalMemory() / 20482048
requestedMem = (Math.round((params.requestedMem * maxmem) / 100))


process resources {
    echo true
    cpus = requestedCpus
    memory = requestedMem.GB

    script:
        """
        echo $requestedCpus
        echo $requestedMem
        """
}
*/


process cluster {
    /*
    On cluster modes:

    https://github.com/soedinglab/MMseqs2/wiki#how-to-set-the-right-alignment-coverage-to-cluster
    */
    // debug true
    container 'nanozoo/mmseqs2:13.45111--810c0f2'
    publishDir "${params.results}", mode: 'copy', overwrite: true
    // Give all resources to this process
    // cpus = requestedCpus
    // memory = requestedMem.GB
    memory = params.maxmem
    cpus = params.maxcpu

    input:
        path(frames)

    output:
        path('cluster_rep_seq.fasta')

    """
    # echo ${task.cpus}
    # echo ${task.memory}
    cat ${frames} > frames.faa
    mmseqs easy-cluster frames.faa cluster tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0
    """
}


// process prokka {
//     /*
//     Only annotation
    
//     sed '/^##FASTA/Q' prokka.gff > nosequence.gff
//     */
//     publishDir "${params.results}/annotation", mode: 'copy', overwrite: true
//     container 'nanozoo/prokka:1.14.6--c99ff65'
//     // container 'staphb/prokka:latest'
//     cpus 1

//     input:
//         tuple(val(name), path(genome))

//     output:
//         tuple(val(name), path("annotation/${name}.gff"))

//     """
//     gunzip -c ${genome} > tmp
//     prokka --mincontiglen 2000 --addgenes --locustag ${name} --cpus ${task.cpus} --prefix ${name} --outdir annotation tmp
//     rm tmp
//     """
// }



/*
mmseqs easy-search -h

 # Search multiple FASTA against FASTA (like BLASTP, TBLASTN, BLASTX, BLASTN --search-type 3, TBLASTX --search-type 2)
 mmseqs easy-search examples/QUERY.fasta examples/QUERY.fasta examples/DB.fasta result.m8 tmp
*/
process search {
    /*
    --search-type INT                Search type 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment [0]

    > A translated nucleotide/nucleotide (TBLASTX) search can be trigged using the flag --search-type 2

    -- https://github.com/soedinglab/mmseqs2/wiki

    // tblastx does 6-frame translation, which is not what we want; we want the nucleotide seq translated, and these proteins aligned

    */
    input:
        path(reference)
        path(proteins)

    output:
        path('aln.m8')

    """
    mmseqs easy-search ${proteins} ${reference} aln.m8 tmp --min-seq-id 0.8 -c 0.95 --max-accept 1 --threads 8 --format-output "query,target,qseq,tseq" --search-type 3
    """
    // --translation-table 11 --search-type 1
    // --search-type 3
    // --cov-mode 4 reduced the number of seqs in the test alignment (primates)
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


process prepare_msa {
    input:
        path(aln)

    output:
        path('4msa/*.fna'), emit: prep_na
        path('4msa/*.faa'), emit: prep_aa

    script:
    """
    4msa.py --aln ${aln} --outdir 4msa
    """
}


process msa {
    input:
        tuple(val(name), path(aa))

    output:
        tuple(val(name), path("${name}.msa.faa"))

    script:
    """
    linsi ${aa} > ${name}.msa.faa
    """
}


/*
TODO: nextalign as alternative
*/
process codon_aln {
    input:
        tuple(val(name), path(msa), path(coding))

    output:
        tuple(val(name), path("${name}.codon.fna"))


    """
    pal2nal.pl ${msa} ${coding} -output fasta > ${name}.codon.fna
    """
}


/*
pal2nal.pl msa.faa test.fna -output fasta > codon.fna
FastTree codon.fna > tree.nwk
hyphy slac --alignment codon.fna --tree tree.nwk --output hyphy.json
*/
process tree {
    input:
        tuple(val(name), path(codon_aln))

    output:
        tuple(val(name), path(codon_aln), path("${name}.tree.nwk"))

    """
    FastTree -nt ${codon_aln} > ${name}.tree.nwk
    """
}


/*
Contig names will lead to errors for the .gbk file, but not other formats. This
is bc/ Genbank enforces < 37 characters in a name and prokka calls tbl2asn which
enforces this.

https://github.com/tseemann/prokka/issues/249
*/
process annotate {
    publishDir "${params.results}/hyphy/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/prokka:1.14.6--c99ff65'
    cpus 8

    input:
        tuple(val(name), path(genome))
    
    output:
        path("${name}.nosequence.gff")
    
    """
    if [[ ${genome} =~ \\.gz\$ ]]
        then
            gunzip -c ${genome} > tmp
        else
            cat ${genome} > tmp
    fi

    prokka --cpus ${task.cpus} --mincontiglen 2000 --outdir annotation --prefix ${name} tmp
    sed '/^##FASTA/Q' annotation/${name}.gff > ${name}.nosequence.gff
    """
}


/*
snippy --outdir 'a' --ctgs '/Users/phi/tmp/evobiochem/workflow/genomes/GCA_000171755.2_ASM17175v2_genomic.fna' --ref GCF_000008685.2_ASM868v2_genomic.fna --cpus 8
snippy --outdir 'b' --ctgs '/Users/phi/tmp/evobiochem/workflow/genomes/GCA_000172315.2_ASM17231v2_genomic.fna' --ref GCF_000008685.2_ASM868v2_genomic.fna --cpus 8
snippy-core --ref 'a/ref.fa' a b c d e f g h i j k

snippy-multi ${spec} --ref ${ref} --cpus 8 > runme.sh
sh ./runme.sh
snippy-clean_full_aln core.full.aln > clean.full.aln
*/
process snippy {
    container 'staphb/snippy:4.6.0-SC2'
    
    input:
        tuple(path(ref), val(name), path(genome))

    output:
        path("${name}")

    """
    snippy --outdir '${name}' --ctgs ${genome} --ref ${ref} --cpus 8
    """
}


process core {
    container 'staphb/snippy:4.6.0-SC2'

    input:
        tuple(val(name), path(ref))
        path(variants)

    output:
        path('clean.full.aln')

    """
    snippy-core --ref ${ref} ${variants}
    snippy-clean_full_aln core.full.aln > clean.full.aln
    """
}


/*
Lots of recombination in borrelia

- https://pubmed.ncbi.nlm.nih.gov/21890743/
- VlsE https://www.cell.com/fulltext/S0092-8674(00)80206-8

process recombination {
    container 'nanozoo/gubbins:3.2.1--d9bcf34'

    input:
        path(core_aln)

    output:
        path('gubbins.filtered_polymorphic_sites.fasta')

    """
    run_gubbins.py --filter-percentage 50 --threads 8 -p gubbins clean.full.aln
    """
}
*/


/*
Can fail bc/:

- < 3 sequences in alignment either bc/ accessory gene or after deduplication
- "ERROR: The alignment does not contain any polymorphic columns; analysis halted."

The output is optional; most likely nothing is found:

- no breakpoints found (only header is returned in x.3s.rec.csv)

Lots of recombination (50) [100%] 50 of 50, cached: 10, failed: 18

We could do this recursively but maybe it is too much trouble for 2 iterations:

https://github.com/nextflow-io/nextflow/discussions/2521
*/
process recombination {
    container 'nanozoo/3seq:1.9.0--dbd8fce'
    memory '1 GB'
    // echo true
    errorStrategy 'ignore'

    input:
        tuple(val(name), path(aln), path(ptable))

    output:
        tuple(val(name), path(aln), env(REC), path("${name}.3s.rec.csv"))

    """
    echo 'y' | 3seq -f ${aln} -d -id ${name} -L100 -ptable ${ptable}

    if [ "\$(wc -l < ${name}.3s.rec.csv)" != 1 ]
    then
        REC=1
    else
        REC=0
    fi
    """
}


// process recombination {
//     container 'nanozoo/3seq:1.9.0--dbd8fce'
//     memory '1 GB'
//     // echo true
//     errorStrategy 'ignore'

//     input:
//         tuple(val(name), path(aln), path(ptable))

//     output:
//         tuple(val(name), path(aln), path("${name}.3s.rec.csv"), optional: true)

//     """
//     echo 'y' | 3seq -f ${aln} -d -id ${name} -L100 -ptable ${ptable}

//     if [ "\$(wc -l < ${name}.3s.rec.csv)" != 1 ]; then rm ${name}.3s.rec.csv; fi
//     """
// }


// TODO: disentange, replace "tree" process
// process tree2 {
//     input:
//         path(no_recomb)

//     output:
//         path('clean.core.tree')

//     """
//     snp-sites -c ${no_recomb} > clean.core.aln
//     FastTree -gtr -nt clean.core.aln > clean.core.tree
//     """
// }


process bar {
    echo true

    input:
        tuple(val(name), path(aln), val(rec), path(recombination))

    script:
    """
    if (( ${rec} == 1 ))
    then
        cat ${recombination} | head -n2
    fi
    """
}


process hyphy_slac {
    publishDir "${params.results}/hyphy/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/hyphy:2.5.39--92c6c14'
    errorStrategy 'ignore'
    cpus 8

    input:
        tuple(val(name), path(codon_aln), path(tree))

    output:
        tuple(val(name), path("${name}__slac"))        

    """
    hyphy CPU=${task.cpus} slac --alignment ${codon_aln} --tree ${tree} --output ${name}.hyphy.slac.json
    cp ${name}.hyphy.slac.json ${name}__slac
    """
}


/*
CPU=8 from https://github.com/veg/hyphy/issues/955
*/
process hyphy_meme {
    publishDir "${params.results}/hyphy/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/hyphy:2.5.39--92c6c14'
    errorStrategy 'ignore'
    cpus 8

    input:
        tuple(val(name), path(codon_aln), path(tree))

    output:
        tuple(val(name), path("${name}__meme"))        

    """
    hyphy CPU=${task.cpus} meme --alignment ${codon_aln} --tree ${tree} --output ${name}.hyphy.meme.json
    cp ${name}.hyphy.meme.json ${name}__meme
    """
}


/*
FUBAR is described in Murrell et al. which is intended to supersede (owing to its remarkable speed and statistical performance), previous REL, SLAC, and and FEL methods (although note SLAC and FEL may still be used). -- https://stevenweaver.github.io/hyphy-site/tutorials/current-release-tutorial/
*/
process hyphy_fubar {
    publishDir "${params.results}/hyphy/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/hyphy:2.5.39--92c6c14'
    errorStrategy 'ignore'
    cpus 8

    input:
        tuple(val(name), path(codon_aln), path(tree))

    output:
        tuple(val(name), path("${name}__fubar"))        

    """
    hyphy CPU=${task.cpus} fubar --alignment ${codon_aln} --tree ${tree} --output ${name}.hyphy.fubar.json
    cp ${name}.hyphy.fubar.json ${name}__fubar
    """
}


/*
> MEME is described in Murrell et al. and is our default recommendation for finding individual sites under selection. It is MUCH slower than FUBAR, however, so there's room for both. -- https://stevenweaver.github.io/hyphy-site/tutorials/current-release-tutorial/
*/
process hyphy_busted {
    publishDir "${params.results}/hyphy/${name}", mode: 'copy', overwrite: true
    container 'nanozoo/hyphy:2.5.39--92c6c14'
    errorStrategy 'ignore'
    cpus 8

    input:
        tuple(val(name), path(codon_aln), path(tree))

    output:
        tuple(val(name), path("${name}__busted"))        

    """
    hyphy CPU=${task.cpus} busted --alignment ${codon_aln} --tree ${tree} --output ${name}.hyphy.busted.json
    cp ${name}.hyphy.busted.json ${name}__busted
    """
}


// process foo {
//     echo true
    
//     input:
//         val(x)
    
//     when:
//         x == 'EMPTY'

//     script:
//     '''
//     echo hello
//     '''
// }


process parse_test {
    publishDir "${params.results}/hyphy/${name}", mode: 'copy', overwrite: true
    // echo true

    input:
        tuple(val(name), path(fubar), path(meme), path(busted), path(msa))

    output:
        tuple(val(name), path("${name}.hyphy.parsed.json"))

    script:
    """
    parse_hyphy.py --selection ${fubar} ${meme} ${busted} --msa ${msa} --out ${name}.hyphy.parsed.json
    """
}


// TODO: Does not work, use python API
// process create_db {
//     output:
//         path('results.db')

//     """
//     sqlite3 results.db
//     """
// }


/*
> [K]eeping identical sequences in the reference alignment can lead to incorrect rearrangement reconstruction down the road. -- http://hyphy.org/w/index.php/IgSCUEAL#Remove_duplicates

More processes downstream will fail (recombination, ...) bc/ we have fewer
MSAs (less duplicates and thus less sequence clusters and more orphan seqs)
*/
process dereplicate {
    cpus 8

    input:
        path(genomes)

    output:
        path('outdir/*')

    """
    # script from Docker container, not .../bin/
    dereplicator.py --threads ${task.cpus} --threshold 0.005 . outdir
    """
}


process annotate_structure {
    container 'nanozoo/foldvis:d91a017--eb4c1ab'
    input:
        tuple(path(fold), path(pfam))

//   Pfam-A.hmm    aln.faa            result.txt
// HmmPy.py                  Pfam-A.hmm.dat    aln.stk            test.tsv
// InteracDome_v0.3-confident.tsv        Pfam-A.hmm.h3f    aln2.stk           test2.tsv
// InteracDome_v0.3-representable.tsv    Pfam-A.hmm.h3i    anvio_parser.py        ynjC
// InteracDome_v0.3-representableNR.tsv  Pfam-A.hmm.h3m    rcsb_pdb_1HSO.dom.tsv
// NDM-1                     Pfam-A.hmm.h3p


    output:
        path('domains.csv')

    script:
    """
    annotate_structure.py --model ${fold} --pfam . --out domains.csv
    """
}


workflow {
    // TODO: probably easiest wf is fold af and __all__ dmasif then filter later/ dynamically in the wf

    // resources()

    // TODO: create results.db, pass fp to subsequent processes, use name
    // db = create_db()
    // as unique ID -- workflow.runName for example
    // https://www.nextflow.io/docs/latest/metadata.html
    // https://docs.python.org/3/library/sqlite3.html
    // onComplete turn db into json or whatever
    // no nosql alternative
    // https://news.ycombinator.com/item?id=27490361
    // sqlite concurrent write
    // https://stackoverflow.com/questions/10325683/can-i-read-and-write-to-a-sqlite-database-concurrently-from-multiple-connections
    // https://www.sqlite.org/quickstart.html
    // TODO: or just dump stuff in the results folder for now?

    // TODO; params file in json etc.
    // -params-file run_42.yaml 

    if ( params.test ) {

        ref = channel.fromPath("${workflow.projectDir}/data/human.fna")
                     .map{ x -> ['ref', x] }

        fold = channel.fromPath("${workflow.projectDir}/data/alphafold2/test_08df6_unrelaxed_rank_1_model_3.pdb")
        
        surface = channel.fromPath("${workflow.projectDir}/data/dmasif/transferrin_binding.pdb")

    } else {

        ref = channel.fromPath(params.ref)
                     .map{ x -> ['ref', x] }

    }

    ptable = channel.fromPath(params.ptable)
    pfam   = channel.fromPath("${params.pfam}/Pfam-A.hmm*").collect().toList()
    /*
    The Pfam database index has several files as components (.h3f, .h3i, ...);
    to link them into the working directory for this process, we need to
    explicitely state that we want them all. We could create one channel for
    each file and then combine them, but this is cumbersome. So we just pass
    the entire folder with a wildcard as a list, and then the process can just
    assume that Pfam.hmm* is present. Of course, this means we cannot rename
    the Pfam database files, but why should we.
    */
    
    if ( params.test ) {
        annotate_structure(fold.combine(pfam))
    }


    // TODO: Add path to alphafold; we can then have a script do
    // all filtering/ window averaging/ ... across all generated results
    // (selection, conservation, ...); really, DMASIF should be there, too;
    // Counterargument: Leave wf as is, publish, then sell 3D stuff separate
    
    if ( params.test ) {

        genomes = channel.fromPath("${workflow.projectDir}/data/primates/*.fna")
                         .map { x -> [x.baseName, x] }

    } else {

        genomes = channel.fromPath(params.genomes)
                         .splitCsv(header: false, sep:'\t')
                         .map{ row -> [row[0], row[1]] }
                         // .view()
    }
    // spec = channel.fromPath(params.genomes)

    // rename(checksum(genomes))

    rename(genomes)
    
    // Dereplicate genomes?
    if ( params.dereplicate ) {
        
        original = rename.out.map { name, fp, log -> [fp] }.collect()
        dereplicate(original)
        rn = dereplicate.out.flatten().map { fp -> [fp.baseName, fp] }
    
    } else {
        
        rn = rename.out.map { name, fp, log -> [name, fp] }
    
    }
    
    coding(rn)
    

    // frames = coding.out.map{ it -> it[1] }.collect()
    // cluster(frames)

    // prokka(genomes)

    rn_ref = rename_ref(ref).map { name, fp, log -> [name, fp] }
    annotate(rn_ref)

    coding_ref(rn_ref)

    search(
        coding_ref.out.map { name, fna, faa -> fna }, 
        coding.out.map { name, fna, faa -> fna }.collect()
    )
    // create_db(ref)

    prepare_msa(search.out)
    
    // TODO: Remove take 3
    x = prepare_msa.out.prep_aa.flatten()
                               .map { it -> [it.baseName, it] }
                               .take(50)
    
    // TODO: inject test MSA (transferrin, ...) here 
    msa(x)
    // TODO: conservation script

    // na .. nucleic acid, aa .. amino acid
    ch_na = prepare_msa.out.prep_na.flatten().map { it -> [it.baseName, it] }
    codon_aln(msa.out.join(ch_na))
    // TODO: clean (remove indels relative to reference)
    
    recombination(codon_aln.out.combine(ptable))
    no_recomb = recombination.out.map { name, aln, flag, recomb -> [name, aln] }

    // bar(recombination.out)
    // recombination.out | ifEmpty { 'EMPTY' } | foo

    tree(no_recomb)
    

    // TODO: remove duplicate seqs or does hyphy offer to do this?
        
    // pervasive diversification
    h1 = hyphy_fubar(tree.out)
    // episodic diversification
    h2 = hyphy_meme(tree.out)
    // gene level
    h3 = hyphy_busted(tree.out)
    parse_test(h1.join(h2).join(h3).join(msa.out))
    



    // TODO: Optional; if there is a channel w/ proteins and binding provided,
    // then integrate this; that way the wf can actually do both, without
    // forcing the user to.



    // println workflow.scriptFile
    // println workflow.projectDir
    

    // TODO: merge json results using jq
    // https://stackoverflow.com/questions/29636331/merging-json-files-using-a-bash-script
    // jq -s 'reduce .[] as $item ({}; . * $item)' json_files/*
    // TODO: Or, just have the parse script do this; that way we can collect
    // and not have to invoke the parse script x times; if no result, put
    // None in there

    // codon_aln.out.view()
    // codon_aln.out.join(hyphy_slac.out).view()

    // TODO: parse hyphy result into annotation track
    // note: we need the codon alignment to remove positions with gaps in the
    // reference

    // TODO: execute when empty https://nextflow-io.github.io/patterns/index.html#_problem_23
    // optional input https://github.com/nextflow-io/patterns/blob/master/optional-input.nf
    // codon_aln.out.view()

    // selection(tree(codon_aln.out))

    // annotate(rn_ref)

    // ch_snippy = annotate.out.combine(genomes)
    // ch_snippy = rn_ref.map { name, fp -> fp }.combine(genomes)
    // ch_snippy.view()
    // snippy(ch_snippy)

    // core(rn_ref, snippy.out.collect())
    // tree2(recombination(core.out))


// https://stackoverflow.com/questions/68317153/nextflow-how-do-you-pass-an-output-multiple-files-from-the-publishdir-to-the
// https://www.nextflow.io/docs/edge/dsl2.html#process-named-output


}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
}

// TODO: search for restricted HMM set through proteins

