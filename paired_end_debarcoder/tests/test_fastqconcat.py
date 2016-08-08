# content of test_sample.py
"""
This file test the fastqconcat script. fastconcat takes two fastqfiles and outputs an outputfiel. Our tests use
example fastq files in form the data directory and write out an output file. We then use the output file to check
whether our commands are correct.

"""
import os
import random

from paired_end_debarcoder.fastq_concat import fastqconcat
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

from click.testing import CliRunner


def fixname(filename):
    "alter a filename to correspond to the test directory"
    currentdir = os.path.dirname(os.path.realpath(__file__))
    return currentdir + "/" + filename

#files
fq1 = fixname("data/fq1.fq")
fq2 = fixname("data/fq2.fq")

#sequences
fq1_first =  next(SeqIO.parse(open(fq1,'r'),'fastq'))
fq2_first =  next(SeqIO.parse(open(fq2,'r'),'fastq'))


# Tests
########################################################################################################################
########################################################################################################################
########################################################################################################################
def test_sequences():
    #basic fastq checking
    seq1 = str(fq1_first.seq)
    seq2 = str(fq2_first.seq)
    assert len(seq1) == 301
    assert len(seq2) == 301
    assert fq1_first.id == fq2_first.id
    assert seq1 == "ATAGCGACCTAGCGTACGCGATGTACACGTCCGGCTCCACAGGGACGCCCAAGGGGGTGATGAACGCGCACGGCGGCGTGGTCAACCGGCTGGCGTGGATGGAGGCCGAATACCGGCTCGGGGCCGGCGAGGCGGTGCTGCAGAAGACGCCGTACACCTTCGACGTGTCGGTGTGGGAGTTCTTCTGGCCGCTGATGACGGGCGCGCGGCTGGTGGTGGCGCGCCCCGGCGGCCACCGCGACCCCGGCTGCCTGGTCGAGACGATCGTCGCGGAGGGGATCACAACGCTTCACGTCGTCCC"
    assert seq2 == "CACTTCGAAGACGTCGCCGGTGCGGTACAGCCGCGCGCCCGGCTCGCCGAAGGGATCGGGGACGAAGCGCTCCGCCGTCTGCCCGGCGCGGTTCAGGTATCCGCGCGCCACCTGGATCCCGCCGATGTACAGCTCGCCCGCCACCCCCGGCGGCACGGGAGCCATCCTCCCATCCAGCAGGTAGGTCCGCACGTTCCCCATCGCCCGGCCGCGCGGGACGCCCCCGGATTCGCCCTCGGCGCCCCCCCCCACCGCAACCGCCACGGCGGCCTCCGCCGGGCGCGACCGGGTGGCGCGCCCC"

def test_fastqconcat_basic():
    "run fastqconcat and check that the sequences have been concatenated."
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fq"
        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 301,
                                '--keep_right', 301,
                                '--no-revcomp',
                                '--no-spacer',
                                '--spacercharacters', "NNNNNNNNNN"])

        #assert the output is the correct length and sequence
        assert not result.exception
        firstrecord = next(SeqIO.parse(open(outfq,'r'),'fastq')).seq
        assert len(firstrecord) == 602
        assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.seq)
        assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.seq)

def test_fastqconcat_revcomp():
    "Make sure that calling reverse complement yields the sequence that has been reverse complemented"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fq"
        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 301,
                                '--keep_right', 301,
                                '--revcomp',
                                '--no-spacer',
                                '--spacercharacters', "NNNNNNNNNN"])

        #assert the output is the correct length and sequence
        assert not result.exception
        firstrecord = next(SeqIO.parse(open(outfq,'r'),'fastq')).seq
        assert len(firstrecord) == 602
        assert str(firstrecord) == str(fq1_first.seq) + str(fq2_first.reverse_complement().seq)

def test_fastqconcat_trimming():
    "Keep only the first 250 characters of each"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fq"

        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 250,
                                '--keep_right', 250,
                                '--no-revcomp',
                                '--no-spacer',
                                '--spacercharacters', "NNNNNNNNNN"])

        #assert the output is the correct length and sequence
        assert not result.exception
        firstrecord = next(SeqIO.parse(open(outfq,'r'),'fastq')).seq
        assert len(firstrecord) == 500
        assert len(str(fq1_first.seq)[:250]) == 250
        assert len(str(fq2_first.seq)[:250]) == 250
        assert str(firstrecord) == str(fq1_first.seq)[:250] + str(fq2_first.seq)[:250]

def test_fastqconcat_discard():
    "When using a trim length greater than the sequence, discard the sequence"

    nucleotides = ['A','C','T','G']

    def randomfastq(name = "seq1",length=10):
        "make a fastq of length length with filler quality values"
        seq = "".join([random.choice(nucleotides) for x in range(length)])
        bioseq = Seq(seq, generic_dna)
        rec = SeqRecord(bioseq)
        rec.letter_annotations['phred_quality'] = [40] * len(bioseq)
        rec.name = name
        rec.id = name
        return rec

    runner=CliRunner()

    with runner.isolated_filesystem():
        "testing that single sequence is not included in the output"
        outfq = "out.fq"

        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 310,
                                '--keep_right', 310,
                                '--discard'])

        #output but no sequence
        assert not result.exception
        assert not os.path.exists(outfq)

    with runner.isolated_filesystem():
        "testing that single sequence is not included in the output"
        outfq = "out.fq"

        seq1f = randomfastq(name="seq1f", length=10)
        seq1r = randomfastq(name="seq1r", length=10)

        SeqIO.write(seq1f,"outf.fq","fastq")
        SeqIO.write(seq1r,"outr.fq","fastq")

        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', "outf.fq",
                                '--reverse_fastq', "outr.fq",
                                '--outfile', outfq,
                                '--keep_left', 11,
                                '--keep_right', 11,
                                '--discard'])

        #output but no sequence
        assert not result.exception
        assert not os.path.exists(outfq)

    with runner.isolated_filesystem():
        """
        Make Three sets of sequences:
        one where both are too short,
        one where both are long enough
        one there only one is long enough.

        Only where both are long enough shouuld we have output
        """
        outfq = "out.fq"

        seq1f = randomfastq(name="seq1f", length=10)
        seq2f = randomfastq(name="seq2f", length=12)
        seq3f = randomfastq(name="seq3f", length=12)

        seq1r = randomfastq(name="seq1r", length=10)
        seq2r = randomfastq(name="seq1r", length=12)
        seq3r = randomfastq(name="seq1r", length=10)

        SeqIO.write([seq1f,seq2f,seq3f],"outf.fq","fastq")
        SeqIO.write([seq1r,seq2r,seq3r],"outr.fq","fastq")

        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', "outf.fq",
                                '--reverse_fastq', "outr.fq",
                                '--outfile', outfq,
                                '--keep_left', 11,
                                '--keep_right', 11,
                                '--discard'])

        #output but no sequence
        assert not result.exception
        assert os.path.exists(outfq)
        records = [rec for rec in SeqIO.parse(open(outfq,'r'),'fastq')]
        assert len(records) == 1

def test_fastqconcat_spacers():
    "Test that the spacer is inserted"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fq"
        spacerchars = "NNNNNNNNNN"
        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 250,
                                '--keep_right', 250,
                                '--no-revcomp',
                                '--spacer',
                                '--spacercharacters', spacerchars])

        #assert the output is the correct length and sequence
        assert not result.exception
        firstrecord = next(SeqIO.parse(open(outfq,'r'),'fastq')).seq
        assert len(firstrecord) == 500 + len(spacerchars)
        assert str(firstrecord) == str(fq1_first.seq)[:250] + spacerchars +str(fq2_first.seq)[:250]

def test_fastqconcat_spacers_and_revcomp():
    "Check the spacer is inserted and the second sequence is reverse complemented"
    runner=CliRunner()
    with runner.isolated_filesystem():
        outfq = "out.fq"
        spacerchars = "NNNNNNNNNN"
        result = runner.invoke(fastqconcat,
                               ['--forward_fastq', fq1,
                                '--reverse_fastq', fq2,
                                '--outfile', outfq,
                                '--keep_left', 250,
                                '--keep_right', 250,
                                '--revcomp',
                                '--spacer',
                                '--spacercharacters', spacerchars])

        #assert the output is the correct length and sequence
        assert not result.exception
        firstrecord = next(SeqIO.parse(open(outfq,'r'),'fastq')).seq
        assert len(firstrecord) == 500 + len(spacerchars)
        assert str(firstrecord) == str(fq1_first.seq)[:250] + spacerchars +str(fq2_first.reverse_complement().seq)[-250:]
