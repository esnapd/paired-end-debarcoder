# Code for Demultiplexing paired-end FastQ files.
#
#
#
#

import click
import os

from Bio import SeqIO
from collections import defaultdict
from itertools import zip_longest
from .barcodes import process_barcodefile, hamdist


###################################################################
## Public Functions

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


@click.command()
@click.option('--forward_fastq', type=click.Path(exists=True), prompt=True,help="name of the fastq forward file")
@click.option('--reverse_fastq', type=click.Path(exists=True), prompt=True, help="name of the fastq reverse file")
@click.option('--barcodefile', type=click.Path(exists=True), prompt=True,help="name of the barcode file")
@click.option('--barcodelength', type=click.INT, prompt=True, help="how long is the barcode")
@click.option('--max_mismatches', type=click.INT, default=1, help="maximum difference between sequence and barcode")
@click.option('--outdirectory', type=click.Path(exists=False, dir_okay=True), prompt=True, help="output directory file")
@click.option('--logfile', type=click.File('w'), prompt=True, help="outputlogfile")
@click.option('--checkbarcodes/--no-checkbarcodes', default=True, help="check the barcodes file")
@click.option('--keepunassigned/--no-keepunassigned', default=True)
@click.option('--chunksize', default=100000)
def demultiplexfastq(forward_fastq, reverse_fastq, barcodefile, barcodelength, max_mismatches, outdirectory,
                     logfile, checkbarcodes, keepunassigned, chunksize):
    """
    Demultiplexing paired Fastq files with a barcode file.

    This is intended for use with MiSeq reads where the paired ends have
    barcodes inside of the Illlumina sequencing primers. The full barcode is
    a concatenation of barcodes on  the forward and reverse read. A two-column
    barcode file must be provided that gives the barcode-to-sample relationship.

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
    barcodes   = process_barcodefile(barcodefile, barcodelength, checkbarcodes=checkbarcodes)
    fastqs     = zip(SeqIO.parse(forward_fastq, 'fastq'), SeqIO.parse(reverse_fastq,'fastq'))

    #check to see if barcodes match and return a dictionary of values
    checked_fastqs = (check_barcode_fastq(fastqs=fastq,
                                          barcodedict=barcodes,
                                          barcodelength=barcodelength,
                                          maxdistance=max_mismatches) for fastq in fastqs)

    # chunk the fastqs to reduce the nubmer of files we write
    groups = grouper(checked_fastqs, chunksize)

    samples = 0
    mismatched_samples = 0
    barcode_mismatched = 0
    sampledict = defaultdict(int)

    if not os.path.exists(outdirectory):
        os.mkdir(outdirectory)

    for group in groups:
        # These dicts will acculate seqs which we will periodically write
        forwardseqlist = defaultdict(list)
        reverseseqlist = defaultdict(list)

        # grouper fills with `None` so we need to remove them
        fastqsnotnone = (fq for fq in group if fq)

        for fastqs in fastqsnotnone:
            sample = fastqs["sample"]
            fq = fastqs["forward_rec"]
            rq = fastqs["reverse_rec"]

            samples += 1
            if sample is None:
                mismatched_samples += 1
            if fastqs['barcode_distance'] > max_mismatches:
                barcode_mismatched += 1

            if shouldwrite(fastqs, barcodedistance=max_mismatches, keepunassigned=keepunassigned):

                # if keep unassigned is true and sample is None, write to Unassigned
                if sample is None:
                    sample = "Unassigned"

                sampledict[sample] += 1

                forwardseqlist[sample].append(fq)
                reverseseqlist[sample].append(rq)

        for sample, fqs in forwardseqlist.items():
            fqout = outdirectory + "/" + sample + "_F.fq"
            with open(fqout, "a") as handle:
                SeqIO.write(fqs, handle, "fastq")
        for sample, rqs in reverseseqlist.items():
            rqout = outdirectory + "/" + sample + "_R.fq"
            with open(rqout, "a") as handle:
                SeqIO.write(rqs, handle, "fastq")


    # write out the log here
    logfile.write("Demultiplex Log\n")
    logfile.write("Samples Processed: {}\n".format(samples))
    logfile.write("Samples With Mismatches less than {}: {}\n".format(max_mismatches, mismatched_samples))
    logfile.write("Counts by Sample:\n")
    for sample, count in sampledict.items():
        logfile.write("{}: {}\n".format(sample, count))

    print("Finished Demultiplexing")


def shouldwrite(fqdict, barcodedistance, keepunassigned, filterspacermismatch=False, printout=False):
    "do we want to write their record?"
    if printout:
        print(fqdict)
    if fqdict['sample'] is None:
        if keepunassigned is True:
            return True
        else:
            return False
    if filterspacermismatch and fqdict['spacermismatch'] is True:
        return False
    if fqdict['barcode_distance'] > barcodedistance:
        return False
    else:
        return True


def check_barcode_fastq(fastqs, barcodedict, barcodelength, maxdistance):
    "check for barcode matches and return a trimmed set of fastqs"
    samplematch      = None
    barcodedata      = None
    spacermismatch   = False
    barcode_distance = 0
    halfbarcode      = int(barcodelength/2)

    #get data
    fq, rq  = fastqs
    barcode = str(fq[:halfbarcode].seq) + str(rq[:halfbarcode].seq)

    #check for perfect match first:
    for	sample, samplebarcodedict in barcodedict.items():
        if samplebarcodedict['barcode'] == barcode:
            samplematch = sample
            barcodedata = samplebarcodedict

    #if not choose closest
    if not samplematch:
        for	sample, samplebarcodedict in barcodedict.items():
            hdist = hamdist(samplebarcodedict['barcode'], barcode)
            if hdist <= maxdistance:
                barcode_distance = hdist
                samplematch = sample
                barcodedata = samplebarcodedict

    # trim the sequences after checking the spacer sequence between the barcode and the primer
    fq = fq[halfbarcode:]
    rq = rq[halfbarcode:]

    # handle length spacer
    if barcodedata is not None:
        forward_spacer = barcodedata['forward_spacer']
        reverse_spacer = barcodedata['reverse_spacer']

        fseq = str(fq.seq)
        rseq = str(rq.seq)

        if fseq.startswith(forward_spacer):
            fq = fq[len(forward_spacer):]
        else:
            fq = fq[len(forward_spacer):]
            spacermismatch = True

        if rseq.startswith(reverse_spacer):
            rq = rq[len(reverse_spacer):]
        else:
            rq = rq[len(reverse_spacer):]
            spacermismatch = True

    samplewrite  = "Unassigned" if samplematch is None else samplematch

    fq.id =  "{}.{} barcode:{} barcodemismatches:{} spacermismatch: {}".format(
            samplewrite, fq.id, barcode, barcode_distance, str(spacermismatch))
    rq.id =  "{}.{} barcode:{} barcodemismatches:{} spacermismatch: {}".format(
            samplewrite, rq.id, barcode, barcode_distance, str(spacermismatch))


    # return updated values
    return {"sample": samplematch,
            "spacermismatch": spacermismatch,
            "forward_rec": fq,
            "reverse_rec": rq,
            "barcode_distance": barcode_distance}
