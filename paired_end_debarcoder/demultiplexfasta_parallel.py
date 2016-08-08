import click
import glob
import os
import shutil
import sys

from subprocess import call
from multiprocess import Pool

###################################################################
## Schemas

###################################################################
## Helper Functions

def which_split():
    "get split command using gpslit from homebrew on mac osx"

    if sys.platform == "darwin":
        #check for homebrew coreutils
        has_gsplit = shutil.which('gsplit')
        if has_gsplit is not None:
            return("gsplit")
        else:
            return("split")
    elif 'linux' in sys.platform:
        return("split")
    else:
        raise ValueError("Only Linux and MacOSX supported")


########################################################################################
## Public Functions

@click.command()
@click.option('--forward_fasta', prompt=True,help="name of the fasta forward file")
@click.option('--reverse_fasta', prompt=True, help="name of the fasta reverse file")
@click.option('--barcodefile', type=click.Path(exists=True), prompt=True,help="name of the barcode file")
@click.option('--barcodelength', type=click.INT, prompt=True, help="how long is the barcode")
@click.option('--outfile', prompt=True, help="output fasta file")
@click.option('--logfile', prompt=True, help="output log file")
@click.option('--max_mismatches', type=click.INT, default=1, help="maximum difference between sequence and barcode")
@click.option('--trimsize_forward',type=click.INT, default=1000)
@click.option('--trimsize_reverse',type=click.INT, default=1000)
@click.option('--includeshort/--no-includeshort', default=False)
@click.option('--spacersequence', default="N")
@click.option('--sampleindex', type=click.INT, default=1)
@click.option('--keepunassigned/--no-keepunassigned', default=True)
@click.option('--splitsize', type=click.INT, default=100000, help="size (in lines) to split fastq")
@click.option('--ncpus', type=click.INT, default=2, help="cpus to use")
def demultiplex_parallel(forward_fasta, reverse_fasta, barcodefile, barcodelength, outfile,logfile, max_mismatches,
                trimsize_forward, trimsize_reverse, includeshort, spacersequence, sampleindex, keepunassigned, splitsize, ncpus):
    """
    A wrapper for `demultiplex_fasta` that splitting the input fasta file and
    does the work on the pieces. Note that inputs are strings not Click.Files as is
    the case for the demultiplex_fasta single script.

    :param forward_fasta:
    :param reverse_fasta:
    :param barcodefile:
    :param barcodelength:
    :param outfile:
    :param logfile:
    :param max_mismatches:
    :param trimsize_forward:
    :param trimsize_reverse:
    :param includeshort:
    :param spacersequence:
    :param sampleindex:
    :return:
    """

    def makecallstring(forward_fasta,reverse_fasta, barcodefile, barcodelength,
                   outfile, logfile, max_mismatches, trimsize_forward,trimsize_reverse, spacersequence, keepunassigned):

        if keepunassigned is True:
            assignstring = "--keepunassigned"
        else:
            assignstring = "--no-keepunassigned"

        return """
        demultiplexfasta \
                --forward_fasta  {} \
                --reverse_fasta  {} \
                --barcodefile    {} \
                --barcodelength  {} \
                --outfile        {} \
                --logfile        {} \
                --max_mismatches {} \
                --trimsize_forward {} \
                --trimsize_reverse {} \
                --no-includeshort \
                {} \
                --spacersequence  {}
        """.format(forward_fasta,reverse_fasta, barcodefile,barcodelength,
                   outfile, logfile, max_mismatches, trimsize_forward,trimsize_reverse, assignstring, spacersequence)


    assert splitsize % 2 == 0


    splitcommand = which_split()
    print("Using the {} split commmand".format(splitcommand))

    #generate the split files nad check for length equivalency
    print("Splitting the Forward Fasta File, {}".format(forward_fasta))
    call("cat {} | {} -l {} - {}".format(forward_fasta, splitcommand, splitsize,"forward_"),shell=True)

    print("Splitting the Reverse Fastq File, {}".format(reverse_fasta))
    call("cat {} | {} -l {} - {}".format(reverse_fasta, splitcommand, splitsize,"reverse_"),shell=True)

    #get the names of the split files and zip them together
    split_files_forward = sorted(glob.glob("forward_*"))
    split_files_reverse = sorted(glob.glob("reverse_*"))
    split_outfiles      = sorted([x.replace("forward_","tempout_")    for x in split_files_forward])
    split_logfiles      = sorted([x.replace("forward_","logfileout_") for x in split_files_forward])

    assert len(split_files_forward) == len(split_files_reverse)

    callstrings = []
    for (f_fasta, r_fasta, out_file, log_file) in zip(split_files_forward,split_files_reverse,split_outfiles,split_logfiles):
        callstring = makecallstring(forward_fasta    = f_fasta,
                                    reverse_fasta    = r_fasta,
                                    barcodefile      = barcodefile,
                                    barcodelength    = barcodelength,
                                    outfile          = out_file,
                                    logfile          = log_file,
                                    max_mismatches   = max_mismatches,
                                    trimsize_forward = trimsize_forward,
                                    trimsize_reverse = trimsize_reverse,
                                    spacersequence   = spacersequence,
                                    keepunassigned   = keepunassigned)
        print(callstring)
        callstrings.append(callstring)

    #process the split files in parallel using multiprocessing
    print("Processing the Split Files in Parallel with {} cpus".format(ncpus))
    p = Pool(ncpus)
    results = p.imap_unordered(lambda x: call(x, shell=True), callstrings)
    for r in results:
        print("Processing split file.")

    #print(os.listdir('.'))
    print("cleaning up the split files....")
    p.map(os.remove, split_files_forward)
    p.map(os.remove, split_files_reverse)


    print("Concatenating the results to {}".format(outfile))
    call("cat tempout_* > {}".format(outfile), shell=True)

    print("Concatenating the log to results to {}".format(logfile))
    call("cat logfileout_*  > {}".format(logfile), shell=True)

    print("cleaning up the temporary files....")
    for ofile in split_outfiles + split_logfiles:
        try:
            os.remove(ofile)
        except:
            pass

    # check output filesize is not zero which will happen
    # if something went wrong with the splitting step due to,say,
    # an error with the mapping file
    if os.path.getsize(outfile) == 0:
        os.remove(outfile)
        print("Error with Your process.. aborting...")
        raise ValueError("Outputfile of size zero indicates an issues with your qiime setup")
