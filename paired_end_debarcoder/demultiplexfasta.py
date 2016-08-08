# Code for Demultiplexing paired-end Fasta files.
#
#
#
#

import click

from Bio import SeqIO
from collections import defaultdict
from cytoolz.dicttoolz import assoc
from cytoolz.functoolz import thread_first
from .barcodes import fastadataSchema, check_barcode, process_barcodefile


###################################################################
## Helper Functions

def reversecomplement(s):
    "reverse complement a DNA strand"
    compdict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    srev = s.upper()[::-1] #
    return "".join([compdict[char] for char in srev])

def fasta_to_dict(fasta):
    "Provide Sequence data as a dictionary"
    fastaf, fastar =  fasta
    return {"forward_id":       fastaf.id,
            "forward_desc":     fastaf.description,
            "forward_sequence": str(fastaf.seq),
            "reverse_id":       fastar.id,
            "reverse_desc":     fastar.description,
            "reverse_sequence": str(fastar.seq)}

def truncate_by_size(fastadict, trimsize_forward, trimsize_reverse):
    "subset sequence and indicate if short"
    fseq = fastadict['forward_sequence']
    rseq = fastadict['reverse_sequence']
    tooshort = False
    if len(fseq) < trimsize_forward:
        tooshort= True
    if len(rseq) < trimsize_reverse:
        tooshort= True

    return thread_first(fastadict,
                        (assoc, "tooshort", tooshort),
                        (assoc, "forward_sequence", fseq[:trimsize_forward]),
                        (assoc, "reverse_sequence", rseq[:trimsize_reverse]))


###################################################################
## Public Functions

@click.command()
@click.option('--forward_fasta', type=click.File('r'), prompt=True,help="name of the fasta forward file")
@click.option('--reverse_fasta', type=click.File('r'), prompt=True, help="name of the fasta reverse file")
@click.option('--barcodefile', type=click.Path(exists=True), prompt=True,help="name of the barcode file")
@click.option('--barcodelength', type=click.INT, prompt=True, help="how long is the barcode")
@click.option('--outfile', type=click.File('w'), prompt=True, help="output fasta file")
@click.option('--logfile', type=click.File('w'), prompt=True, help="output log file")
@click.option('--max_mismatches', type=click.INT, default=1, help="maximum difference between sequence and barcode")
@click.option('--trimsize_forward',type=click.INT, default=1000)
@click.option('--trimsize_reverse',type=click.INT, default=1000)
@click.option('--includeshort/--no-includeshort', default=False)
@click.option('--spacersequence', default="N")
@click.option('--sampleindex', type=click.INT, default=1)
@click.option('--keepunassigned/--no-keepunassigned', default=True)
def demultiplex(forward_fasta, reverse_fasta, barcodefile, barcodelength, outfile,logfile, max_mismatches,
                trimsize_forward, trimsize_reverse, includeshort, spacersequence, sampleindex, keepunassigned):
    """
    Demultiplexing paired Fasta files with a barcode file.

    This is intended for use with MiSeq reads where the paired ends have
    barcodes inside of the Illlumina sequencing primers. The full barcode is
    a concatenation of barcodes on  the forward and reverse read. A barcode file
    must be provided that gives the barcode-to-sample relationship.

    Fasta files are expected to have been pretrimmed for quality using, for example,
    seqtk.

    This script will:

    1. check the barcode allowing for a minimal distance to the supplied barcode file.
    2. truncate forward and reverse to specified values
    3. concatenate the forward and reverse together
    4. trim by size
    5. output fasta to specified file

    """

    # get the barcode and fasta data
    barcodes   = process_barcodefile(barcodefile, barcodelength)
    fastas     = zip(SeqIO.parse(forward_fasta, 'fasta'), SeqIO.parse(reverse_fasta,'fasta'))
    fastadicts = (fasta_to_dict(fasta) for fasta in fastas)

    # get barcode information
    fastabarcodes = (check_barcode(fastadict,
                                    barcodedict=barcodes,
                                    barcodelength=barcodelength,
                                    maxdistance=max_mismatches)
                     for fastadict in fastadicts)

    #filter sizes, and reverse complement
    fastasizetruncated = (truncate_by_size(fastadict,
                                           trimsize_forward=trimsize_forward,
                                           trimsize_reverse=trimsize_reverse)
                          for fastadict in fastabarcodes)

    # validate data before progressing
    fastadata = (fastadataSchema.validate(d) for d in fastasizetruncated)

    #iterate through and keep relevant data
    tooshortcount = 0
    badbarcodecount = 0
    errorcount = 0
    count = 0
    samplecounts = defaultdict(int)

    for result in fastadata:
        #sampledata
        forward_id     = result['forward_id']
        forward_desc   = result["forward_desc"]
        forward_seq    = result["forward_sequence"]
        reverse_id     = result["reverse_id"]
        reverse_desc   = result["reverse_desc"]
        reverse_seq    = result["reverse_sequence"]
        sample         = result["sample"]
        barcode        = result["barcode"]
        brcd_dist      = result["barcode_distance"]
        tooshort       = result["tooshort"]
        spacermismatch = result['spacermismatch']

        #accounting
        count += 1
        samplecounts[sample] += 1
        if not sample: badbarcodecount += 1
        if tooshort: tooshortcount += 1

        #write sample
        def writesample(forward_seq=forward_seq,
                        reverse_seq=reverse_seq,
                        sample=sample,forward_id=forward_id, count=count, barcode=barcode, brcd_dist=brcd_dist):


            #combine the forward and reverse sequence
            allseq = forward_seq +  spacersequence + reversecomplement(reverse_seq)

            # write out sequences
            if sample is None:
                sample = "Unassigned"

            fastaheader = "{}.{}.{:06d} barcode:{} barcodemismatches:{} spacermismatch: {}".format(
                sample, forward_id, count, barcode, brcd_dist, str(spacermismatch))

            outfile.write(">{}\n{}\n".format(fastaheader,allseq))

        def shouldwritesample(sample=sample,includeshort=includeshort,tooshort=tooshort,
                              brcd_dist=brcd_dist,max_mismatches=max_mismatches):
            "encapsulate sequence-writing logic in a function"

            # Only use sequences samples that have a sample
            if not sample:
                if keepunassigned:
                    return True
                else:
                    return False

            # Ignore short sequences if the flag is false
            if includeshort is False and tooshort is True:
                return False

            # Ignore sequences with barcode mismatches above the threshold
            if brcd_dist > max_mismatches:
                return False

            return True

        shouldwrite = shouldwritesample()

        if shouldwrite == True:
            writesample()

    # write out log information
    logfile.write("""
       Barcode File: {}
       Sequenced Processed: {}
       Samples Below the Length Cutoff: {}
       Samples Unassigned due to Barcodes: {}

       """.format(barcodefile, count, tooshortcount, badbarcodecount))

    for sam, cnt in samplecounts.items():
        logfile.write("Observed Counts for Sample {}: {}\n".format(sam,cnt))

    print("Finished Demultiplexing")