# gene_loader.py
import bisect
import HTSeq
from Gene_Site_Iter_Graph_v1_dev import Gene

def createGenes(annotation, aType, qGene):
    """
    Parse an annotation file and construct an array of Gene objects,
    organized by chromosome, for later splice-site analysis.

    Expects the annotation file to contain gene entries with at least the following:
    - Chromosome identifier
    - Strand
    - Transcriptional start and end positions
    - Gene identifier

    For each line in the annotation describing a gene, a Gene object is created
    and placed in a chromosome-indexed list. If a specific gene is queried (via qGene),
    that Gene object is also returned as QUERY_gene.

    Assumes the gene list is sorted by ascending transcriptional start site.

    Parameters
    ----------
    annotation : str
        Path to the GFF3/GTF annotation file.
    aType : str
        The annotation feature type to filter for (e.g., 'gene').
        Currently unused, retained for compatibility.
    qGene : str
        If not 'All', limits parsing to a specific gene name.

    Returns
    -------
    chrom_index : list of str
        List of chromosome names seen in the annotation.
    gene2D_array : list of list of Gene
        A 2D list where each sublist contains Gene objects for a given chromosome.
    QUERY_gene : Gene or None
        The Gene object matching qGene, if found; otherwise None.
    """
    chrom_index = []
    gene2D_array = []
    QUERY_gene = None
    geneCounter = 0

    for line in HTSeq.GFF_Reader(annotation):
        if line.type == 'gene':
            GeneName = line.name
            chrom = line.iv.chrom
            GeneLeft = line.iv.start
            GeneRight = line.iv.end
            strand = line.iv.strand

            if chrom not in chrom_index:
                chrom_index.append(chrom)
                gene2D_array.append([])

            if qGene == 'All':
                geneCounter += 1
                bisect.insort(
                    gene2D_array[chrom_index.index(chrom)],
                    Gene(
                        chromosome=chrom,
                        name=GeneName,
                        leftPos=GeneLeft,
                        rightPos=GeneRight,
                        readNums=None,
                        samples=1,
                        strand=strand,
                        source=''
                    )
                )
            elif GeneName == qGene:
                print('Query Gene found')
                geneCounter += 1
                QUERY_gene = Gene(
                    chromosome=chrom,
                    name=GeneName,
                    leftPos=GeneLeft,
                    rightPos=GeneRight,
                    readNums=None,
                    samples=1,
                    strand=strand,
                    source=''
                )
                gene2D_array[chrom_index.index(chrom)].append(QUERY_gene)

    print(f"{geneCounter} Genes created in {len(chrom_index)} bins")

    NA_gene = Gene(chromosome = None,
                    name = 'NA',
                    leftPos = -1,
                    rightPos = -1,
                    readNums = None,
                    samples = 1,
                    strand = None,
                    source = '')

    return chrom_index, gene2D_array, QUERY_gene, NA_gene