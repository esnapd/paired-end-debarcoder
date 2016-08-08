# content of test_demultiplexfasta.py.py
"""
This file tests the demultiplexfasta script. demutiplexfasta two fastqfiles and outputs a fastafile where samples
have been labelled by sample id.

"""
import pytest
import os
from Bio import SeqIO
from paired_end_debarcoder.demultiplexfastq import demultiplexfastq
from click.testing import CliRunner


def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename

#files
fq1 = fixname("data/smallfq1.fq")
fq2 = fixname("data/smallfq2.fq")
barcodesfile = fixname("data/smallbarcodes2.txt")

########################################################################################################################
########################################################################################################################
# Tests

def test_demultiplexfastq():
    "test basic demultiplexing"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outlogfile = "outlog.log"
        result = runner.invoke(demultiplexfastq,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 16,
                                '--max_mismatches', 0,
                                '--outdirectory', ".",
                                '--logfile', outlogfile,
                                '--no-keepunassigned'])

        #assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == {"outlog.log",
                                        "Sample1_F.fq", "Sample1_R.fq",
                                        "Sample2_F.fq", "Sample2_R.fq"}

        for (name, num) in [("Sample1_F.fq", 2), ("Sample1_R.fq", 2),
                            ("Sample2_F.fq", 1), ("Sample2_R.fq", 1)]:
            records = [r for r in SeqIO.parse(open(name,'r'),'fastq')]
            assert len(records) == num

    # this time all unaddigned
    with runner.isolated_filesystem():
        outlogfile = "outlog.log"
        result = runner.invoke(demultiplexfastq,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--barcodefile', barcodesfile,
                                '--barcodelength', 16,
                                '--max_mismatches', 0,
                                '--outdirectory', ".",
                                '--logfile', outlogfile])

        # assert the output is the correct length and sequence
        assert result.exit_code == 0
        assert set(os.listdir(".")) == {"outlog.log",
                                        "Sample1_F.fq", "Sample1_R.fq",
                                        "Sample2_F.fq", "Sample2_R.fq",
                                        "Unassigned_F.fq", "Unassigned_R.fq"}

