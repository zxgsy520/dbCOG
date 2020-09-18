#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import time
import random
import logging
import argparse

from Bio import SeqIO
from Bio import Entrez
from multiprocessing import Pool

LOG = logging.getLogger(__name__)

__version__ = "v1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


Entrez.email = "mingyan24@126.com"


def read_csv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(',')
        

def read_cog(file):

    data = {}

    for line in read_csv(file):
        if len(line)<=7:
            LOG.debug(line)
        else:
            data[line[2]] = line[7]

    return data


def download_seq(seqid, newid, sformat, dbase):

    handle = Entrez.efetch(db=dbase, id=seqid, rettype=sformat)

    record = SeqIO.read(handle, sformat)
    if newid != "":
        seqid = newid
    print('>%s|%s\n%s\n' % (seqid, record.id, record.seq))
    LOG.debug('%s download completed' % record.id)

    return 0


def download_cog(file, sformat, dbase, thread=10):

    data = read_cog(file)
   
    pool = Pool(processes=thread)

    for i in data:
        time.sleep(random.randint(0, 1))
        pool.apply_async(download_seq, (i, data[i], sformat, dbase))

    pool.close()
    pool.join()


def add_help_args(parser):
    parser.add_argument('input',
        help='Input the csv file of the cog database.')
    parser.add_argument('-f', '--format', metavar='STR', type=str, default='fasta',
        help='Set the downloaded file format, default=fasta')
    parser.add_argument('-db', '--dbase', metavar='STR', type=str, default='protein',
        help='Set database type(nucleotide or protein,default=nucleotide.')
    parser.add_argument('-t', '--thread', metavar='INT' ,type=int, default=10,
        help='Set the number of threads, default=10.')

    return parser


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
        download_cog_seqs.py -- Download the sequence of the cog database
attention:
        download_cog_seqs.py cog-20.cog.csv -f fasta -db protein >cog.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    download_cog(args.input, args.format, args.dbase, args.thread)


if __name__ == "__main__":
    main()
