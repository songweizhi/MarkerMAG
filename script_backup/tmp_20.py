
from sys import exit, stderr

try:
    from Bio import AlignIO, SeqIO
except ImportError:
    stderr.write("Error! BioPython is not detected!\n")
    exit(1)
