
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc

seq_in = 'ATGC'

seq_in_rc = get_rc(seq_in)

print(seq_in_rc)

print('ATGC'[::-1])