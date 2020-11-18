#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import gzip
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "v1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta(file):

    LOG.info("reading message from %r" % file)

    if file.endswith("gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    seq = ''
    for line in fh:
        if type(line) == type(b''):
            line = line.decode('utf-8')

        line = line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith(">"):
            if seq!='':
                yield seq.split('\n')
            seq = "%s\n" % line.strip('>')
            continue
        seq += line

    if seq!='':
        yield seq.split('\n')
    fh.close()


def merge_cog_seq(files):

    for file in files:
        cog = file.split('/')[-1].split('.')[0]
        for seqid, seq in read_fasta(file):
            print(">%s|%s\n%s" % (cog, seqid, seq))

    return 0


def add_help_args(parser):

    parser.add_argument('input', nargs='+',
        help='Input COG sequence file.')

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
        merge_cog_seq.py -- Merge sequences from COG database
attention:
        merge_cog_seq.py COG*.fa.gz >COG.fasta
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    merge_cog_seq(args.input)


if __name__ == "__main__":
    main()
