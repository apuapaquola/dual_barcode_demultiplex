#!/usr/bin/env python
'''Demultiplex fastq files with dual barcodes.

Uses fixed number of mismatches, allowing no indels.


--
Copyright (c) 2013 by Apua Paquola <apuapaquola@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''

import argparse
import sys
import gzip
import re
sys.path.append('/opt/readfq')
from readfq import readfq

def distance_k_variants(seq, distance, letters='ACGTN'):
    """Generator function that yields all up to 'distance'-nucleotide variants of seq, 
       no indels allowed.
    """
    yield(seq)
    if distance==0:
        return
    for i in range(len(seq)):
        for t in letters:
            if seq[i] != t:
                head=seq[:i]+t
                tail=seq[i+1:]
                for x in distance_k_variants(tail, distance-1, letters):
                    yield(head+x)


def read_barcode_dict(file_name):
    """Reads a tab-delimited file containing sample name and barcode into a
    dictionary indexed by barcode.
    Example file:
    A5-1    CGTGAT
    A5-2    ACATCG
    A5-3    GCCTAA
    A5-4    TGGTCA
    """
    f=open(file_name)
    d=dict()
    for line in f:
        a=line.split('\t')
        a[-1]=a[-1].rstrip()
        if a[1] in d:
            raise Exception("Barcode " + a[1] + " already present in dictionary")
        d[a[1]]=a[0]
    return(d)
    f.close()
        

def gapless_alignment(seq1, seq2):
    minlen=min(len(seq1),len(seq2))
    return(''.join([' |'[seq1[i]==seq2[i]] for i in range(minlen)]))

        
def expand_barcode_dict(bc, distance):
    """Given a barcode dictionary bc, returns a new dictionary containing all
    'distance'-nucleotide variants of each key referring to that key.
    Example:
    Input:
    {'AA': 'sample1'} 
    Output: expand_barcode_dict({'AA': 'sample1'}, distance=1)
    {'NA': 'AA', 'AN': 'AA', 'CA': 'AA', 'AG': 'AA', 'TA': 'AA', 'AA': 'AA', 'AC': 'AA', 'AT': 'AA', 'GA': 'AA'}
    """
    
    d=dict()
    for seq in bc.keys():  
        for seq_var in distance_k_variants(seq, distance):
            if seq_var in d:
                msg="Barcode conflict\n"
                msg+=bc[seq] + ": " +seq + " -> " + seq_var + "\n"
                msg+=seq+"\n"+gapless_alignment(seq, seq_var)+"\n"+seq_var+"\n"
                msg+="\n"
                msg+=bc[d[seq_var]] + ": " +d[seq_var] + " -> " + seq_var+"\n"
                msg+=d[seq_var]+"\n"+gapless_alignment(d[seq_var], seq_var)+"\n"+seq_var+"\n"
                raise Exception(msg)
            d[seq_var]=seq
    return(d)


def demultiplex_dual_barcodes():
    parser = argparse.ArgumentParser()
    parser.add_argument('--index1_file', required=True)
    parser.add_argument('--index1_max_distance', type=int, required=True)
    parser.add_argument('--index2_file', required=True)
    parser.add_argument('--index2_max_distance', type=int, required=True)
    parser.add_argument('--read1_fq', required=True)
    parser.add_argument('--read2_fq', required=True)
    parser.add_argument('--output_dir', required=True)
    #parser.add_argument('--output_file_prefix')
    args=parser.parse_args()

    bc=dict()
    bce=dict()
    for i in [1,2]:
        bc[i]=read_barcode_dict(args.__dict__['index'+str(i)+'_file'])
        bce[i]=expand_barcode_dict(bc[i], args.__dict__['index'+str(i)+'_max_distance'])



    outfiles=dict()
    def output_file(filename):
        if filename not in outfiles:
            outfiles[filename]=gzip.open(filename, 'wt')
        return outfiles[filename]


    f1=readfq(open(args.read1_fq))
    f2=readfq(open(args.read2_fq))

    for name1, seq1, qual1 in f1:
        m=re.search('#([ACGTN]+)_([ACGTN]+)/', name1)
        
        sample_id=dict()
        for i in [1,2]:
            index_seq=m.group(i)
            if index_seq in bce[i]:
                sample_id[i]=bc[i][bce[i][index_seq]]
            else:
                sample_id[i]=None
        
        if sample_id[1] is None or sample_id[2] is None:
            filename_suffix='Unassigned'
        else:
            filename_suffix=sample_id[1]+'_'+sample_id[2]
        
        base_output_filename=args.output_dir+'/'+filename_suffix
        
        f=output_file(base_output_filename+'_R1.fq.gz')
        f.write('\n'.join(['@'+name1,seq1,'+',qual1,'']))
        
        name2, seq2, qual2 = next(f2)
        f=output_file(base_output_filename+'_R2.fq.gz')
        f.write('\n'.join(['@'+name2,seq2,'+',qual2,'']))
        
        #print(output_filename)


    for f in outfiles.values():
        f.close()
        

if __name__ == "__main__":
    demultiplex_dual_barcodes()
