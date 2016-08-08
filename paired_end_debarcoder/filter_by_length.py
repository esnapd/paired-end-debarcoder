import click
from Bio import SeqIO

@click.command()
@click.option('--fastafile', type = click.Path(exists=True), help="inputfastafile")
@click.option('--outputfasta',type = click.STRING, help="outputfile")
@click.option('--length', type=click.INT)
def filter_by_length(fastafile, outputfasta, length):
    """
    from an input Fastafile, filter reads not matching the length parameter
    """
    records = (rec for rec in SeqIO.parse(open(fastafile,'r'),'fasta') if len(rec.seq) == length)
    SeqIO.write(records, outputfasta,'fasta')
